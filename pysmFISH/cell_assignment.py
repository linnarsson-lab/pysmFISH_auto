from turtle import st
import numpy as np
import pandas as pd
import math
import gc
from skimage.segmentation import expand_labels
from scipy.spatial import KDTree
import numpy_groupies as npg
from numpy.lib.stride_tricks import sliding_window_view
import loompy
import zarr

class Cell_Assignment:
    
    def expand_labels_mask(self, mask, distance, out_file_name='expanded', 
                           chunk_size=5000, overlap_percentage=10, 
                           dtype = np.uint32):
        """Expand segmentation mask without overlapping labels.
        
        Expands the labels of a large segemented image in a memory efficient 
        way. The output (and ideally also the input) are Zarr arrays that
        are stored on disk and not loaded in memory at once.
        Uses the Scikit image expand_label() function under the hood.
        
        Args:
            mask (np.ndarray): Array with segmentation mask where the pixels of
                each individual cell is labeled with an unique integer. Ideally
                this should be a Zarr array.
            distance (int): Number of pixels to exapand the segmentation masks 
                with.
            out_file_name (str): Name of the output file. The output will be a 
                Zarr array. Name will be appended with '.zarr'.
                Defaults to 'expanded'
            chunk_size (int): Size of the chunks in the X axis. Size may differ
                depending on the shape of the dataset in the Y axis, pick it
                accordingly. Defaults to 5000.
            overlap_percentage (int): Percentage of overlap between the chunks. 
                Overlap is needed to correctly handle the expansion of labels
                wihtout overlapping with neighbouring cells. The overlap 
                percentage should be large enough so that the overlap contains
                at least 2 cell diameters. The overlap well be calculated as a
                percentage of the chunk_size. Defaults to 10. With the default
                settings the overlap will thus be 500 pixels. 
            dtype (type): Numpy dtype for output datasets. Defaults to 
                np.uint32.
        
        Returns:
            Out_file_name (str): Output is writtent to the file specified.
        
        """
        #Open memmap for results
        out_file_name = out_file_name + '.zarr'
        exp = zarr.open(out_file_name, mode='w', shape = mask.shape, chunks=(chunk_size, mask.shape[1]), dtype=dtype)
        #Change to uint64 when you get more than 4e9 cells.

        #Make chunks
        chunks = np.arange(0, ((math.ceil(mask.shape[0] / chunk_size)) + 1) * chunk_size, chunk_size)
        #Calculate overlap
        overlap_percentage = overlap_percentage / 100
        overlap = int(chunk_size * overlap_percentage)

        #Iterate over chunks and asign dots to cells
        for i, (view_min, view_max) in enumerate(sliding_window_view(chunks, 2)):
            print(i, view_min, view_max)
            if i == 0:
                view_min_extra = view_min    
                view_max_extra = view_max + overlap
            else:
                view_min_extra = view_min - overlap
                view_max_extra = view_max + overlap  
            
            #Expand labels
            img = mask[view_min_extra:view_max_extra, :]
            img = expand_labels(img, distance=distance)

            #Asign to large array
            if i == 0: #First chunk
                exp[view_min:view_max, :] = img[0:-overlap,:]

            elif i == chunks.shape[0] - 2: #Last chunck
                exp[view_min:, :] = img[overlap:]

            else: #Middle chunks
                exp[view_min:view_max, :] = img[overlap:-overlap,:]
            
        return out_file_name

    def assignment(self, mask, points, genes, offset = np.array([0,0]), 
                   chunk_offset=0, unique_genes=None, workers=-1):
        """Count molecules in cells.
        
        Count which points fall within a segmentation mask and return a gene by
        cell matrix.
        
        If running out of memory check the chunked version of the function.
        
        Args:
            mask (np.ndarray): Array with segmentation mask where the pixels of
                each individual cell are labeled with an unique positive 
                integer.
            points (np.ndarray): Arry with shape (X,2) with positions of the 
                detected molecules. These should match the mask coordinates, or 
                should match after applying the offset.
            genes (np.ndarray): Array with gene labels for the given points.
            offset (np.ndarray): This function asumes that the origin of both 
                the mask and the points co-localizes in the coordinates: (0,0).
                Using this "offset" parameter the origin of the mask can be 
                moved by passing an array with shape two for the X and Y 
                displacement. Defults to np.array([0,0]).
            chunk_offset (float): Use this only when using the chuncked version
                of this function. The chunk offset is used to move the mask 
                in X for every chunck processed.
            unique_genes (np.ndarray, optional): Array with unique genes in the 
                order that np.unique() would give. If given it will speed up 
                calculation. Defaults to None. 
            workers (int): Number of workers. If -1 is given it will use all
                processors. Defaults to -1.
        
        Returns:
            pandas.DataFrame: Dataframe with genes in rows, cells in columns 
                and counts per gene as values.
            np.ndarray: Array with unique genes.
            np.ndarray: Array with cell labels for each point. Points that fall
                outside cells are labeled -1.
        
        """
        
        #Get unique genes and cells if not given
        if type(unique_genes) == type(None):
            print('Calculating unique genes')
            unique_genes = np.unique(genes)
        
        #Get position of segmented cell pixels
        cell_idx = np.array(np.nonzero(mask)).T
        #Get ID or cell pixels
        cell_id = mask[cell_idx[:,0], cell_idx[:,1]]
        #Apply offset to mask coordinates
        cell_idx += offset + np.array([chunk_offset, 0])

        #Build tree
        #tree = KDTree(cell_idx) 
        tree = KDTree(cell_idx, leafsize=20, compact_nodes=False, balanced_tree=False) 
        #I would think it should be `cell_idx + np.array([0.5, 0.5])` so that it is the middle of the pixel, but this does not look correct in the plots.
        #Query the points
        #Minkowsky distance to infinite so that it is square to match the grid nature of the points
        r = tree.query_ball_point(points, p=np.inf, r=0.5, workers=workers)
        del tree
        gc.collect()

        #Find gene identity and cell identity of all points that fall inside a segmented cell
        gene_id = []
        in_cell_id = []
        all_labels = []
        for gi, ri in zip(genes, r):
            if ri != []:
                gene_id.append(gi)
                cell = cell_id[ri[0]]
                in_cell_id.append(cell)
                all_labels.append(cell)
            else:
                all_labels.append(-1)

        gene_id = np.asarray(gene_id)
        in_cell_id = np.asarray(in_cell_id)
        all_labels = np.asarray(all_labels)
        
        if in_cell_id.shape[0] > 0:
            #Group counts per cell
            #Numpy groupies can not handel it if there are missing lables:
            #For instance consider this array of cell IDs: [0,1,2,4]. In this case cell ID 3 is missing.
            #It needs to be filled. To prevent a "Future waring elementwise comparison failed" fill with empty arrays.
            agg = npg.aggregate_np(in_cell_id, gene_id, func='array', dtype=object, fill_value=np.array([], dtype=object))
            #Remove empty arrays (This are cell labels that do not exist).
            agg_shape = np.array([i.shape[0] for i in agg])
            agg_clean = agg[agg_shape > 0]

            # Count gene occurences for each cell ID
            cell_count = [np.unique(i, return_counts=True) for i in agg_clean]

            #Aggregate data in pandas dataframe
            unique_cell_id = np.unique(in_cell_id)
            df = pd.DataFrame(data=np.zeros((unique_genes.shape[0], unique_cell_id.shape[0])), index=unique_genes, columns=unique_cell_id, dtype=int)

            for (cc_g, cc_c), c_id in zip(cell_count, unique_cell_id):
                df.loc[cc_g, c_id] = cc_c

            return df, unique_genes, all_labels
        
        else:
            return None, unique_genes, all_labels
        
    def asignment_chunked(self, mask, points, genes, offset=np.array([0,0]), 
                          chunk_size_x=5000, unique_genes=None, workers=-1):
        """Count molecules in cells. Chunked version for memory efficiency.
        
        Count which points fall within a segmentation mask and return a gene by
        cell matrix.
        Chunking will be performed on the X axis only. So chunks will have a
        defined X dimention but not a defined Y dimention. Depending on the 
        shape of the dataset adjust the chunks_size.
        
        Args:
            mask (np.ndarray): Array with segmentation mask where the pixels of
                each individual cell are labeled with an unique positive 
                integer.
            points (np.ndarray): Arry with shape (X,2) with positions of the 
                detected molecules. These should match the mask coordinates, or 
                should match after applying the offset.
            genes (np.ndarray): Array with gene labels for the given points.
            offset (np.ndarray): This function asumes that the origin of both 
                the mask and the points co-localizes in the coordinates: (0,0). 
                Using this "offset" parameter the origin of the mask can be 
                moved by passing an array with shape two for the X and Y
                displacement. Defults to np.array([0,0]).
            chunk_size (int): Size of the chunks in the X axis. Size may differ
                depending on the shape of the dataset in the Y axis, pick it
                accordingly. Defaults to 5000.
            unique_genes (np.ndarray, optional): Array with unique genes in the 
                order that np.unique() would give. If given it will speed up 
                calculation. Defaults to None. 
            workers (int): Number of workers. If -1 is given it will use all
                processors. Do not expect a large performace increase.
                Defaults to -1.
        
        Returns:
            pandas.DataFrame: Dataframe with genes in rows, cells in columns 
                and counts per gene as values.
            np.ndarray: Array with unique genes.
            np.ndarray: Array with cell labels for each point. Points that fall
                outside cells are labeled -1.
        
        """
        #Get unique genes and cells if not given
        if type(unique_genes) == type(None):
            print('Calculating unique genes')
            unique_genes = np.unique(genes)
        
        #Make chunks
        chunks = np.arange(0, ((math.ceil(mask.shape[0] / chunk_size_x)) + 1) * chunk_size_x, chunk_size_x)
        
        #Iterate over chunks and asign dots to cells
        results = []
        point_cell_id = - np.ones(points.shape[0], dtype='int')
        for view_min, view_max in sliding_window_view(chunks, 2):
            filt = (points[:,0] > view_min ) & (points[:,0] < view_max)
            df, unique_genes, cell_id = self.assignment(mask[view_min : view_max+1, :], points[filt,:], genes[filt], 
                                            offset = np.array([0,0]), chunk_offset=view_min,
                                            unique_genes=unique_genes, workers=workers)
            results.append(df)
            point_cell_id[filt] = cell_id
            
        #Drop None results
        results = [r for r in results if type(r) != type(None)]
        #Merge all dataframes for all chunks
        total = pd.concat(results, axis=1)
        #Find all duplicated columns
        dup_filt = total.columns.duplicated(keep='first')
        duplicates = total.loc[:, dup_filt]
        #Remove duplicated columns but keep the first
        total = total.loc[:, ~dup_filt]
        #Sum the duplicated columns
        total.loc[:,duplicates.columns] = total.loc[:,duplicates.columns] + duplicates
        
        return total.sort_index(axis=1), unique_genes, point_cell_id

    def make_loom(self, df, ra={}, ca={}, fa={}, out_file_name = 'data.loom'):
        """Make loom file of dataframe.
        
        Args:
            df (pd.dataFrame): Pandas dataframes with genes in Index and
                cellIDs in columns.
            ra (dict): Row attributes. Dictionary with row attributes.
                Values should be in the same order and shape as the index.
                "Genes" will be automatically handled and do not need to be
                specified. Input should be numpy arrays or lists.
                Suggestion for standard names:
                    "Chromosome": Chromosome number of gene.
            ca (dict): Column attributes. Dictionary with colum attributes.
                Values should be in the same order and shape as the columns.
                "CellID" and "TotalMolecules" will be automatically handled
                and do not need to be specified. Input should be numpy 
                arrays or lists.
                Suggestion for standard names: 
                    "X" X coordinate of cell centoid.
                    "Y" Y coordinate of cell centroid.
                    "UMAP_1" UMAP component 1
                    "UMAP_2" UMAP component 2
                    "tSNE_1" tSNE component 1
                    "tSNE_2" tSNE component 2
                    "nucleus_size_um2" Size of nucleus in um2
                    "nucleus_size_px" Size of nucleus in pixels
            fa (dict): Global file atrributes.
            out_file_name (str): Name of the loom file to create. Should
                have the .loom extension. Defualts to "data.loom"
        
        """
        
        ra['Gene'] = df.index.to_numpy()
    
        ca['CellID'] = df.columns.to_numpy()
        ca['TotalMolecules'] = df.sum().to_numpy()
        
        loompy.create(out_file_name, df.values, row_attrs=ra, col_attrs=ca, file_attrs=fa)