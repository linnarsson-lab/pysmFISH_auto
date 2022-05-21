import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import KDTree
from skimage.draw import disk
import dask
from dask.diagnostics import ProgressBar
from typing import Union, Tuple
import pickle as pkl
try:
    import ripleyk
except ModuleNotFoundError:
    print('Please install module `ripleyk` using pip. Source: https://github.com/SamPIngram/RipleyK') 

class BeadAlignment:
    
    def __init__(self, initial_scale_factor:float =1, search_fraction:float =0.1,
                 initial_rotation: float=0, rotation_search_width: float=4,
                 search_radius:float=1000, samples:int = 25,               
                 focusing_bins:int =100, centering_mode: str = 'middle',
                 max_broad_sweeps:int =10, num_narrow_sweeps:int =2):
        """Initialize bead alignment class.
        
        Class to align two 2D point clouds that have (some) matching points.
        
        Initiate the class and then use the `.fit()` function to learn the
        transfrom parameters. Then any dataset can be transformed using the
        `.transfrom()` function.
        Args:
            initial_scale_factor (float, optional): Expected scale factor. 
                Number that when multiplied with the source data would bring 
                the source data approximately to the scale of the target. 
                Defaults to 1.
            search_fraction (float, optional): Fraction of initial_scale_factor
                that is used for the search range. Example; if the 
                initial_scale_factor is 1 and the search_fraction is 0.1, scale
                factors between 0.95 and 1.05 will be investigated.
                Defaults to 0.1.
            initial_rotation (float, optional): Initial rotation to center
                rotation search around. In degees. Defaults to 0.
            rotation_search_width (float, optional): Search range for rotation.
                Eample; if initial_rotation is 0 and rotation_search_width is
                4, rotations between -2 and 2 degree will be investigated.
                Defaults to 4.
            search_radius (float, optional): Radius within to look for matching
                points. Defaults to 1000.
            samples (int, optional): Number of samples within the search range.
                Will be made uneven because then the search range will include
                the initial_scale_factor. Defaults to 25.
            focusing_bins (int, optional): Number of bins used for focussing.
                Should be around 100. Defaults to 100.
            centering_mode (str, optional): Method to initially center the two 
                datasets. Options:
                `middle`: Will center the datasets based on their centers.
                `mean`: Will center the datasets based on their weighted 
                    centers. 
                Defaults to 'middle'.
            max_broad_sweeps (int, optional): maximum number of broad sweeps.
                When one sweep for scaling factors is not successfull it will
                increase the sampling density and range. Defaults to 10.
            num_narrow_sweeps (int, optional): Number of higher density sweeps
                once it has found a candidate scaling factor. Defaults to 2.
        """
        #Input
        self.initial_scale_factor = initial_scale_factor
        self.search_radius = search_radius
        self.search_fraction = search_fraction
        self.initial_rotation = initial_rotation
        self.rotation_search_width = rotation_search_width
        self.samples = samples
        self.focusing_bins = focusing_bins
        self.centering_mode = centering_mode
        #Change to uneven number of samples because then the initial_scale_factor is included
        if self.samples % 2 == 0:
            self.samples += 1
        self.max_broad_sweeps = max_broad_sweeps
        self.num_narrow_sweeps = num_narrow_sweeps

    def align(self, target: np.ndarray, source: np.ndarray, radius: float) -> Tuple[np.ndarray, np.ndarray]:
        """Radius neighbour search to find matching points in both datasets.
        Args:
            target (np.ndarray): Array with target point set.
            source (np.ndarray): Array with source point set.
            radius (float): Search radius. matching point should fall within 
                the radius from each other.
        Returns:
            np.ndarray: Jaggered array with distances for each target point to
                neighbouring source points. 
            np.ndarray: Jaggered array with angle in radians for each target 
                point to neighbouring source points. 
        """
        
        tree = KDTree(target, leaf_size=100)
        neighbours, dist = tree.query_radius(source, radius, return_distance=True)
        
        #Get coordinates of neighbours
        neigh_coord = [target[i] for i in neighbours]
        #Calculate X and Y delta between point and neighbours
        neigh_delta = [i-j for i,j in zip(neigh_coord, source)]
        #Calculate radians between point and neighbours
        neigh_rad = np.asarray([np.arctan2(i[:,0], i[:,1]) for i in neigh_delta], dtype='object')
        
        return dist, neigh_rad

    def flatten(self, dist: np.ndarray, radians: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Flatten jaggered arrays.
        
        Sepecifically used with the output of the `align()` function.
        Args:
            dist (np.ndarray): Jaggered array.
            radians (np.ndarray): Jaggered array
        Returns:
            Union[np.ndarray, np.ndarray]: Flattened arrays of the dist and
                radians input. 
        """
        dist_flat = np.hstack(dist.flatten())
        neigh_rad_flat = np.hstack(radians.flatten())
        return dist_flat, neigh_rad_flat

    def plot_align(self, dist: np.ndarray, radians: np.ndarray, s: float=0.5,
                   alpha: float=0.005):
        """Plot the alignment results
        Args:
            dist (np.ndarray): Jaggered array with distances for each target 
                point to neighbouring source points. 
            radians (np.ndarray): Jaggered array with angle in radians for each
                target point to neighbouring source points. 
            s (float, optional): Point size. Defaults to 0.5.
            alpha (float, optional): Alpha value of points. Defaults to 0.005.
        """
        
        #Flatten for plotting
        dist, radians = self.flatten(dist, radians)
        
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(projection='polar')
        ax.scatter(radians, dist, s=s, alpha=alpha, c='k')
        
    def pol2cart(self, rho: np.ndarray, phi: np.ndarray) -> np.ndarray:
        """Convert polar coordinates to cartesian.
        Args:
            rho (np.ndarray): Distances.
            phi (np.ndarray): Angles in radians. 
        Returns:
            np.ndarray: 2D array with cartesion coordinates.
        """
        x = rho * np.cos(phi)
        y = rho * np.sin(phi)
        return(np.array([x, y]).T)
    
    def make_mask(self, bins: int) -> np.ndarray:
        """Make boolean mask of a circle. 
        Args:
            bins (int): Number of bins used.
        Returns:
            np.ndarray: 2D array with boolean circular mask, where the center
                is True.
        """
        mask = np.zeros((bins, bins), dtype=np.uint8)
        r = int(bins/2)
        rr, cc = disk((r, r), r)
        mask[rr, cc] = 1
        return mask.astype('bool')
    
    def find_center(self, p: np.ndarray) -> np.ndarray:
        """Find the center of a dataset based on the extend of the data.
        
        Sensitive to outliers. 
        Args:
            p (np.ndarray): XY coordinates of points.
        Returns:
            np.ndarray: Center of dataset.
        """
        
        p_min = p.min(axis=0)
        p_max = p.max(axis=0)
        p_extend = p_max - p_min
        p_center = p_min + (0.5 * p_extend)
        return p_center

    def center_datasets(self, target: np.ndarray, source: np.ndarray, mode: Union[str, np.ndarray]='middle') -> Tuple[np.ndarray, np.ndarray]:
        """Center datasets as initial alignment.
        Args:
            target (np.ndarray): Array with target point set.
            source (np.ndarray): Array with source point set.
            mode (str, np.ndarray): Alignment mode:
                "middle": Calculates centers of both datasets and aligns these.
                "mean": Calculates weighted centers of both datasets and aligns
                    these
                If a numpy array is given with shape 2, the source dataset will
                be moved using these values.                 
        Returns:
            np.ndarray: moved transformed dataset.
            np.ndarray: Array with the XY distance the dataset has been moved.
        """

        if type(mode) == str:
    
            if mode.lower() == 'mean':
                #Center datasets (uses mean, which might not be optimal for weird datasets)
                target_mean = target.mean(axis=0)
                source_mean = source.mean(axis=0)

                delta = target_mean - source_mean
                source = source + delta
            
            elif mode.lower() == 'middle':
                t_center = self.find_center(target)
                s_center = self.find_center(source)
                
                delta = t_center - s_center
                source = source + delta

        elif type(mode) == np.ndarray:
            if mode.shape == (2,):
                source = source + mode
                delta = mode
            else:
                raise Exception(f'Transform array shoud have shape (2,), not: {mode.shape}')
        
        else:
            raise Exception(f'Invalid input: {mode.shape}')
        
        return source, delta

    def scale(self, source: np.ndarray, factor: float) -> np.ndarray:
        """Scale the source with a scale factor.
        Args:
            source (np.ndarray): Array with source point set.
            factor (float): Scaling factor.
        Returns:
            np.ndarray: Scaled dataset.
        """
        return source * factor
    
    def rotate(self, p: np.ndarray, origin: tuple=(0, 0), angle: float=0, degrees: bool=None) -> np.ndarray:
        """Rotate points around an origin with a certain angle.
        Args:
            p (np.ndarray): Array with xy coordinates
            origin (tuple, optional): Point to rotate around. Defaults to (0, 0).
            angle (float, optional): Angle in radians between -pi and pi.
                Defaults to 0.
            degrees (bool, optional): Angle to rotate in degree, if this is given 
                the function will ignore any "angle" input. Defaults to None.
        Returns:
            np.ndarray: Rotated points. 
        """
        if degrees != None:
            angle = np.deg2rad(degrees)
        R = np.array([[np.cos(angle), -np.sin(angle)],
                    [np.sin(angle),  np.cos(angle)]])
        o = np.atleast_2d(origin)
        p = np.atleast_2d(p)
        return np.squeeze((R @ (p.T-o.T) + o.T).T)
    
    def gini(self, x):
        sorted_x = np.sort(x)
        n = len(x)
        cumx = np.cumsum(sorted_x, dtype=float)
        return (n + 1 - 2 * np.sum(cumx) / cumx[-1]) / n
    
    def broad_evaluator(self, coord: np.ndarray, bins: int, mode: str='default') -> float:
        """Broad evaluator of the scale factor.
        Calculates how focused the matching points are after the scaling and
        neighbour search.
        For a more accurate evaluator see the `narrow_evaluator()` function.    
        
        Args:
            coord (np.ndarray): Cartesian coordinates of the neighbour search.
            bins (int): Number of bins for the histogram. Dependends on the 
                search radius. Should be around 100-500. 
        Returns:
            float: Focussing value, where higher values mean a more correct
                scale factor.
        """
        
        #Mask circle
        mask = self.make_mask(bins)
        
        #Histogram
        hist = np.histogram2d(coord[:,0], coord[:,1], bins=bins)
        
        # 
        data = hist[0][mask]
        mean = np.mean(data)
        std = np.std(data)
        peak = hist[0].max()
        
        if mode == 'scan':
            result = abs(peak - mean)
        elif mode == 'score':
            result = peak 
        elif mode == 'default':
            result = peak 

        return result
    
    def find_median(self, coord: np.ndarray, tolerance: float=1e-5) -> np.ndarray:
        """Iteratively find median of a point cloud.
        
        Each iteration it focuses the dataset around the median untill the 
        difference between iterations is below the tolerance.
        Args:
            coord (np.ndarray): Cartesian coordinates subset of the neighbour 
                search. Subset should be focused on the peak. 
            tolerance (float, optional): Tolerance. Defaults to 1e-5.
        Returns:
            np.ndarray: Coordinates of the peak. 
        """

        median = np.median(coord, axis=0)
        selection = coord
        for i in range(30):
            tree = KDTree(coord)
            std = np.std(selection, axis=0).mean()
            neigh = tree.query_radius(median.reshape(1, -1), 3*std)
            selection = coord[neigh[0]]
            new_median = np.median(selection, axis=0)
            if np.all(np.abs(median - new_median) < np.array([tolerance, tolerance])):
                return np.flip(new_median)
            else:
                median = new_median
        
        return np.flip(median)
      
    def narrow_evaluator(self, coord: np.ndarray, bins: int, find_offset: bool=False) -> Tuple[float, Union[np.ndarray, float]]:
        """Narrow evaluator of scale factor.
        Args:
            coord (np.ndarray): Cartesian coordinates of the neighbour search.
            bins (int): Number of bins for the histogram. Dependent on the 
                search radius. Should be around 100-500. 
            find_offset (bool, optional): If True, also find the XY offset. 
                Takes some computational time so it is advised to only perform
                in the last iteration. Defaults to False.
        Returns:
            float: Focussing value, where lower is better.
            Union[np.ndarray, float]: If `find_offset` is True retuns a numpy 
                array with the XY offset. If False, returns zero. 
        """

        #Histogram
        hist = np.histogram2d(coord[:,0], coord[:,1], bins=bins)
        
        x, y = np.where(hist[0] == np.amax(hist[0]))
        center_x = hist[1][x[0]]
        center_y = hist[2][y[0]]
        step = hist[1][1] - hist[1][0]
        
        #Select all point around focal point
        filt_x = np.logical_and(coord[:,0]>center_x - 3*step, coord[:,0]<center_x + 3*step)
        filt_y = np.logical_and(coord[:,1]>center_y - 3*step, coord[:,1]<center_y + 3*step)
        filt = np.logical_and(filt_x, filt_y)
        coord_filt = coord[filt]
        
        #Calculate clusterdness using Ripley's K
        x_extend = coord_filt[:,0].max() - coord_filt[:,0].min()
        y_extend = coord_filt[:,1].max() - coord_filt[:,1].min()
        r = step * 0.1
        result = ripleyk.calculate_ripley(r, [x_extend, y_extend], coord_filt[:,0], coord_filt[:,1], sample_shape='rectangle', CSR_Normalise=True)
        result = result * coord_filt.shape[0]
        
        #Measure nearest neighbour distance
        # Did not work because data is discreet due to the pixels and many points are on top of each other.
        # So that the NN-dist is often 0
        #tree = KDTree(coord_filt)
        #dist, _ = tree.query(coord_filt, return_distance=True, k=2)
        #result = dist[:,1].mean()
        
        #Get distance of lowest percentile
        # Also did not work because of the discreet nature of the data. 
        #perc = np.percentile(dist[:,1], 0.5)
        #dist_selec = dist[:,1][dist[:,1] <= perc]
        #print(dist_selec.shape, dist.shape, perc, (dist[:,1] < perc).sum())
        
        #Find offset    
        if find_offset:
            offset = self.find_median(coord_filt)
        
        #Find the offset of the data
        if find_offset:
            return result, offset
        else:
            return result, 0

    def eval_param_worker(self, target: np.ndarray, source: np.ndarray, factor: float, angle: float=0,
                          search_radius: float=2000, 
                          bins: int=100, centering_mode: Union[str, np.ndarray]='mean', find_offset: bool=False, 
                          evaluator: str='broad', plot: bool=False, mode: str='default'):
        """Evaluate a factor that transfroms the source to the target scale.
        Args:
            target (np.ndarray): Array with target point set.
            source (np.ndarray): Array with source point set.
            factor (float): Scale factor that when multiplied with the source,
                would bring the transformed source to the exact scale of the 
                target.
            angle (float): Angle to rotate the source data with. 
            search_radius (float, optional): Search radius. matching point 
                should fall within the radius from each other.
                Defaults to 2000.
            bins (int, optional): Number of bins for the histogram. Dependends 
                on the 
                search radius. Should be around 100-500. . Defaults to 100.
            mode (str, np.ndarray): Alignment mode:
                "middle": Calculates centers of both datasets and aligns these.
                "mean": Calculates weighted centers of both datasets and aligns
                    these
                If a numpy array is given with shape 2, the source dataset will
                be moved using these values.
             find_offset (bool, optional): If True, also find the XY offset. 
                Takes some computational time so it is advised to only perform
                in the last iteration. Defaults to False.
            evaluator (str, optional): Evaluator function to use. Either 
                'broad' or 'narrow'. Defaults to 'broad'.
            plot (bool, optional): If True plots output. Defaults to False.
        Returns:
            np.ndarray: Result of the evaluator function.
            If find_offset is True, also returns:
                np.ndarray: Found XY offset.
                np.ndarray: XY offset calculated by the centering function.
            The final XY offset is the above two offsets summed.
        """
        
        #Scale
        source = self.scale(source, factor)
        
        #Rotate
        source = self.rotate(source, origin=self.find_center(source), degrees=angle)
        
        #Center
        source, delta = self.center_datasets(target, source, mode=centering_mode)
        
        #Align
        dist, a = self.align(target, source, search_radius)
        
        #Plot
        if plot:
            self.plot_align(dist, a)
            
        #To cartesian
        dist_flat, a_flat = self.flatten(dist, a)
        coord = self.pol2cart(dist_flat, a_flat)
        
        #Evaluate
        if evaluator == 'broad':
            if mode == 'score':
                matched = self.broad_evaluator(coord, bins, mode='score')
                unmatched = len(source) - matched
                result = matched/unmatched # returns the overlap ratio, raw overlap between the two tested datasets 
            elif mode == 'scan':
                result = self.broad_evaluator(coord, bins, mode='scan')
                result = result
            else:
                result = self.broad_evaluator(coord, bins, mode='default')

        else:
            result, offset = self.narrow_evaluator(coord, bins, find_offset=find_offset)
        
        #Find the offset of the data
        if find_offset:
            return result, offset, delta
        else:
            return result

    def eval_param_worker_parallel(self, target: np.ndarray, source: np.ndarray, 
                                   scale_search_space: Union[float, np.ndarray], 
                                   angle_search_space:Union[float, np.ndarray], 
                                    search_radius: float=2000, bins: int=100, 
                                    centering_mode: Union[str, np.ndarray]='middle', find_offset: bool=False,
                                    evaluator: str='broad') -> dict:
        """Evaluate multiple scale factors in parallel.
        Args:
            target (np.ndarray): Array with target point set.
            source (np.ndarray): Array with source point set.
            search_space (np.ndarray): Array with scale factors to evaluate.
            search_radius (float, optional): Search radius. matching point 
                should fall within the radius from each other.
                Defaults to 2000.
            bins (int, optional): Number of bins for the histogram. Dependends 
                on the search radius. Should be around 100-500. 
                Defaults to 100.
            mode (str, np.ndarray): Alignment mode:
                "middle": Calculates centers of both datasets and aligns these.
                "mean": Calculates weighted centers of both datasets and aligns
                    these
                If a numpy array is given with shape 2, the source dataset will
                be moved using these values.
            find_offset (bool, optional): If True, also find the XY offset. 
                Takes some computational time so it is advised to only perform
                in the last iteration. Defaults to False.
            evaluator (str, optional): Evaluator function to use. Either 
                'broad' or 'narrow'. Defaults to 'broad'.
        Returns:
            dict: Dictionary with the follwing keys:
            'r': Evaluation results.
            If find_offset is True also contains:
            'c': XY offset found by evaluator.
            'd': XY offset found by centering function.
        """

        #Inputs and outputs
        # target_d = dask.delayed(target)           # MODIFIED (MARCOS) 
        # source_d = dask.delayed(source)
        
        if find_offset:
            results = {'r' : [],
                    'c' : [],
                    'd' : []}
        else:
            results = {'r' : []}
        
        #Evaluating scale factor or rotation angle    
        if type(scale_search_space) == np.ndarray:
            scale = True
            angle = False
            search_space = scale_search_space
        if type(angle_search_space) == np.ndarray:
            scale = False
            angle = True
            search_space = angle_search_space

        #Parallel iterate through options
        for i in search_space:
            if scale == True:
                r = dask.delayed(self.eval_param_worker)(target, source, i, angle_search_space, search_radius, bins, centering_mode, find_offset, evaluator)
            if angle == True:
                r = dask.delayed(self.eval_param_worker)(target, source, scale_search_space, i, search_radius, bins, centering_mode, find_offset, evaluator)
            if find_offset:
                results['r'].append(r[0])
                results['c'].append(r[1])
                results['d'].append(r[2])
            else:
                results['r'].append(r)

        #Compute          # ALTERED HERE REMOVED PROGRESS BAR (MARCOS)
        with ProgressBar():
            result = dask.compute(results)
            
        return result[0]

    def eval_results(self, result: dict, search_space: np.array, plot: bool=False, evaluator: str='broad', title: str=''):
        """Interpret results of evaluator.
        Args:
            result (dict): Dictionary with results from the 
                `eval_scale_worker_parallel()` function.
            search_space (np.ndarray): Array with scale factors that were
                evaluated.
            plot (bool, optional): If True plots results. Defaults to False.
            evaluator (str, optional): Evaluator function to use. Either 
                'broad' or 'narrow'. Defaults to 'broad'.
        Returns:
            bool: True if peak was found.
            float: New lower bound. To use when another iteration is perfomed.
                Will be narrower if peak is found, or wider if no peak is 
                found.
            float: New upper bound. To use when another iteration is perfomed.
                Will be narrower if peak is found, or wider if no peak is 
                found.
            int: New number of samples to take. Will be same as number of
                samples in provided search_space.
            float: Value of peak if found, otherwise 0.
            np.ndarray: XY offset found by evaluator.
            np.ndarray: XY offset found by centering function.
                The two XY offsets summed will be the final XY transformation.
        """

        if len(result) == 3:
            values = np.asarray(result['r'])
            centers = result['c']
            deltas = result['d']
        else:
            values = np.asarray(result['r'])
            centers = np.zeros_like(values)
            deltas = np.zeros_like(values)
                
        min_ss = search_space[0]
        max_ss = search_space[-1]
        range_ss = max_ss - min_ss
        half_range = 0.5 * range_ss
        n_samples = search_space.shape[0]
        step = search_space[1] - search_space[0]
        
        if plot:
            plt.figure(figsize=(7,2))
            plt.plot(values)
            plt.xticks(np.arange(n_samples), labels = search_space, rotation=-90)
            plt.ylabel('Score')
            plt.title(title)
        
        #Found a peak
        if values.sum() > 0:
            optimum = values.max() #Maximize intensity for broad sweep, or Ripley k for narrow sweep
            max_idx = np.where(values == optimum)[0][0]
            peak = search_space[max_idx]
            center = centers[max_idx]
            delta = deltas[max_idx]
            
            #Peak is located at the edge, so it might not be the best value yet
            if max_idx == 0 or max_idx == search_space.shape[0]:
                print('Peak located at edge')
                return False, peak - half_range, peak + half_range, n_samples, peak, center, delta
            
            #Peak is fully sampled, narrow search parameters
            else:
                return True, peak - 2*step, peak + 2*step, n_samples, peak, center, delta
        
        #Did not find a peak, widen search and make denser
        else:
            print('No peaks found. Widening the search to find the scaling factor. If this happends in the last cycles \
,double the width of the search and triple the number of samples.')
            return False, min_ss - half_range, max_ss + half_range, n_samples * 3, 0, 0, 0
        
    def find_transform(self, target: np.ndarray, source: np.ndarray, 
                          initial_scale_factor: float=1, scale_search_fraction: float=0.1, 
                          initial_rotation: float=0, rotation_search_width: float=4,
                          samples: int=25, max_broad_sweep: int=10, 
                        num_narrow_sweeps: int=2, search_radius: float=2000, bins: int=100,
                        centering_mode: Union[str, np.ndarray]='mean', plot: bool=False) -> Tuple[float, np.ndarray, np.ndarray]:
        """Find scale factor and XY offset.
        
        First performs one or multiple broad sweeps that approximately try to
        find the correct scale factor. The correct scale factor is the factor
        that when multiplied with the source data would bring the transformed
        source data to the same scale as the target data.
        Afterwards it performs a requested number of narrow sweeps to optimize
        the scale factor, and in the last iteration also find the XY offset.
        Togheter the scale factor and the XY offset form the function that can
        transfrom the source data to exactly match the target dataset.
        
        The search relies on an initial scaling factor that would be close to
        the actual value. Pick this carefully, toghether with the search range.
        However, the algorithm is quite robust and can have a broad range. Just
        increase the number of samples when doing a broader sweep.
        
        Args:
            target (np.ndarray): Array with target point set.
            source (np.ndarray): Array with source point set.
            initial_scale_factor (float, optional): Estimated scale factor.
                Ideally this is approximately correct. Defaults to 1.
            scale_search_fraction (float, optional): Used to calculate the 
                interval over which the initial search if performed. Example;
                if the initial_scale_factor is 1 and the search_fraction is 
                0.1, scale factors between 0.95 and 1.05 will be investigated.
                Defaults to 0.1.
            initial_rotation (float, optional): Initial rotation to center
                rotation search around. In degees. Defaults to 0.
            rotation_search_width (float, optional): Search range for rotation.
                Eample; if initial_rotation is 0 and rotation_search_width is
                4, rotations between -2 and 2 degree will be investigated.
                Defaults to 4.
            samples (int, optional): Number of samples to take between the 
                lower and upper bound, as defined by the `initial_scale_factor`
                and the `sarch_fraction`. When doing a broad search pick more 
                samples. Uneven numbers work better because then the 
                `initial_scale_factor` will the included in the searchspace.
                Defaults to 25.
            max_broad_sweep (int, optional): Maximum number of broad sweeps.
                If a broad sweep is not succesfull it will widen the search and
                make the search denser. Defaults to 10.
            num_narrow_sweeps (int, optional): Number of narrow sweeps to fine
                tune the scaling factor and XY offset. Defaults to 2.
            search_radius (float, optional): Search radius. matching point 
                should fall within the radius from each other.
                Defaults to 2000.
            bins (int, optional): Number of bins for the histogram. Dependends 
            on the search radius. Should be around 100-500. . Defaults to 100.
            centering_mode (Union[str, np.ndarray], optional):Alignment mode:
            "middle": Calculates centers of both datasets and aligns these.
            "mean": Calculates weighted centers of both datasets and aligns
                these.
                If a numpy array is given with shape 2, the source dataset will
                be moved using these values.
            plot (bool, optional): If true plots results. Defaults to False.
        Raises:
            Exception: When broad sweep does not find anything.
        Returns:
            float: Optimal scale factor.
            np.ndarray: XY offset found by evaluator.
            np.ndarray: XY offset found by centering algorithm. 
        """
        #Scale factor broad sweep
        half_range = initial_scale_factor * (0.5 * scale_search_fraction)
        min_search = initial_scale_factor - half_range
        max_search = initial_scale_factor + half_range
        
        count = 0
        while True:
            search_space = np.linspace(min_search, max_search, num=samples)
            result = self.eval_param_worker_parallel(target, 
                                                     source, 
                                                     search_space, 
                                                     0,
                                                     search_radius, 
                                                     bins, 
                                                     centering_mode, 
                                                     evaluator='broad')
            succes, min_search, max_search, samples, scale_factor, center, delta = self.eval_results(result,
                                                                                                     search_space, 
                                                                                                     plot, 
                                                                                                     evaluator='broad',
                                                                                                     title='Broad sweep scale factor')
            if succes == True:
                break
            count += 1
            if count >= max_broad_sweep:
                raise Exception('Bead alignment could not find the scale factor. Try again with different parameters.')
            
        #Rotation broad sweep
        a_half_range = rotation_search_width / 2
        a_min_search = initial_rotation - a_half_range
        a_max_search = initial_rotation + a_half_range
        count = 0
        while True:
            search_space = np.linspace(a_min_search, a_max_search, num=samples)
            result = self.eval_param_worker_parallel(target, 
                                                     source, 
                                                     scale_factor, 
                                                     search_space,
                                                     search_radius, 
                                                     bins, 
                                                     centering_mode, 
                                                     evaluator='broad')
            a_succes, a_min_search, a_max_search, a_samples, angle, a_center, a_delta = self.eval_results(result, 
                                                                                                          search_space, 
                                                                                                          plot, 
                                                                                                          evaluator='broad',
                                                                                                          title='Broad sweep rotation')
            if a_succes == True:
                break
            count += 1
            if count >= max_broad_sweep:
                raise Exception('Bead alignment could not find the rotation. Try again with different parameters.')


        #Narrow sweeps
        find_offset = False
        for i in range(num_narrow_sweeps):
            if i == num_narrow_sweeps-1:
                find_offset = True
                
            #Rotation
            search_space = np.linspace(a_min_search, a_max_search, num=a_samples)
            result = self.eval_param_worker_parallel(target, 
                                                     source, 
                                                     scale_factor, 
                                                     search_space,
                                                     search_radius, 
                                                     bins, 
                                                     centering_mode, 
                                                     evaluator='narrow')
            a_succes, a_min_search, a_max_search, a_samples, angle, a_center, a_delta = self.eval_results(result, 
                                                                                                          search_space, 
                                                                                                          plot, 
                                                                                                          evaluator='narrow',
                                                                                                          title='Narrow sweep rotation')
            
            #Scale factor
            search_space = np.linspace(min_search, max_search, num=samples)
            result = self.eval_param_worker_parallel(target, 
                                                     source, 
                                                     search_space, 
                                                     angle,
                                                     search_radius,
                                                     bins*(2+i), 
                                                     centering_mode, 
                                                     find_offset, 
                                                     evaluator='narrow')
            succes, min_search, max_search, samples, factor, center, delta = self.eval_results(result, 
                                                                                               search_space, 
                                                                                               plot, 
                                                                                               evaluator='narrow',
                                                                                               title='Narrow sweep scale factor')            
        
        return factor, center, delta, angle

    def fit(self, target: np.ndarray, source: np.ndarray, plot: bool=False):
        """Find the scale factor and XY offset.
        
        First performs one or multiple broad sweeps that approximately try to
        find the correct scale factor. The correct scale factor is the factor
        that when multiplied with the source data would bring the transformed
        source data to the same scale as the target data.
        Afterwards it performs a requested number of narrow sweeps to optimize
        the scale factor, and in the last iteration also find the XY offset.
        Togheter the scale factor and the XY offset form the function that can
        transfrom the source data to exactly match the target dataset.
        
        The search relies on an initial scaling factor that would be close to
        the actual value. Pick this carefully, toghether with the search range.
        However, the algorithm is quite robust and can have a broad range. Just
        increase the number of samples when doing a broader sweep.
        
        Output can be found as attributes:
        self.factor: Found scale factor.
        self.offset: Found XY offset
        Also available: (These summed are self.offset)
            self.center: XY offset found by centering function.
            self.delta: Xy offset found by evaluator.
        
        To change the parameters please change the parameters of the class:
        
            initial_scale_factor (float, optional): Expected scale factor. 
                Number that when multiplied with the source data would bring 
                the source data approximately to the scale of the target. 
                Defaults to 1.
            search_fraction (float, optional): Fraction of initial_scale_factor
                that is used for the search range. Example; if the 
                initial_scale_factor is 1 and the search_fraction is 0.1, scale
                factors between 0.95 and 1.05 will be investigated.
                Defaults to 0.1.
            initial_rotation (float, optional): Initial rotation to center
                rotation search around. In degees. Defaults to 0.
            rotation_search_width (float, optional): Search range for rotation.
                Eample; if initial_rotation is 0 and rotation_search_width is
                4, rotations between -2 and 2 degree will be investigated.
                Defaults to 4.
            search_radius (float, optional): Radius within to look for matching
                points. Defaults to 1000.
            samples (int, optional): Number of samples within the search range.
                Will be made uneven because then the search range will include
                the initial_scale_factor. Defaults to 25.
            focusing_bins (int, optional): Number of bins used for focussing.
                Should be around 100. Defaults to 100.
            centering_mode (str, optional): Method to initially center the two 
                datasets. Options:
                `middle`: Will center the datasets based on their centers.
                `mean`: Will center the datasets based on their weighted 
                    centers. 
                Defaults to 'middle'.
            max_broad_sweeps (int, optional): maximum number of broad sweeps.
                When one sweep for scaling factors is not successfull it will
                increase the sampling density and range. Defaults to 10.
            num_narrow_sweeps (int, optional): Number of higher density sweeps
                once it has found a candidate scaling factor. Defaults to 2.
        Args:
            target (np.ndarray): Array with target point set.
            source (np.ndarray): Array with source point set.
            plot (bool, optional): If True plots otimization results.
                Defaults to False.
        """
        factor, center, delta, angle = self.find_transform(target, source, self.initial_scale_factor, self.search_fraction,
                                                            self.initial_rotation, self.rotation_search_width,
                                                            self.samples, self.max_broad_sweeps, self.num_narrow_sweeps,
                                                            self.search_radius, self.focusing_bins, self.centering_mode,
                                                            plot=plot)
        
        self.factor = factor
        self.center = center
        self.delta = delta
        self.offset = center + delta
        self.angle = angle
        self.rotation_origin = self.find_center(source)
        
    def transform(self, X: np.ndarray) -> np.ndarray:
        """Transform 2D array with found transfrom function.
        Args:
            X (np.ndarray): 2D array with points to transfrom.
        Raises:
            Exception: If X is not a numpy array.
            Exception: If there are not 2 columns. Points should have X and Y
                coordinates in two columns. 
            Exception: If Transform parameters have not be found yet. Please 
                run `.fit()` method first.
        Returns:
            np.ndarray: Transformed coordinates.
        """
                
        if type(X) != np.ndarray:
            raise Exception(f'Input needs to by a numpy array, not: {type(X)}.')
        if X.shape[1] != 2:
            raise Exception(f'Input must have 2 columns, {X.shape[1]}.')
        
        if hasattr(self, 'factor') and hasattr(self, 'offset') and hasattr(self, 'angle') and hasattr(self, 'rotation_origin'):
            #Scale
            #X = (X * self.factor)
            #Rotate
            #X = self.rotate(X, origin=self.rotation_origin, degrees=self.angle)
            #Move
            #X = X + self.offset
            
            #Rotate
            X = self.rotate(X, origin=self.rotation_origin, degrees=self.angle)
            #Scale and move
            X = (X * self.factor) + self.offset
            
            
            
            return X
        else:
            raise Exception('Transform parameters not known. Please run `fit()` function first.')

        #---------------------------------------MARCOS-------------------------------------#

    def subsample(self, to_subsample, xlist, ylist):

        '''Subsample the given array into specified chunk.
        Args: 
            to_subsample: array to subsample.
        Returns:
            np.ndarray: subsampled array.
        '''

        height = max(to_subsample[:,0]) - min(to_subsample[:,0])
        width = max(to_subsample[:,1]) - min(to_subsample[:,1])

        subsampled = to_subsample[(to_subsample[:,0] > height*xlist[0]) & (to_subsample[:,0] < height*xlist[1]) & \
        (to_subsample[:,1] > width*ylist[0]) & (to_subsample[:,1] < width*ylist[1])]

        return subsampled

    def isnan(self, arr):

        '''Checks if the array contains NaN values.
        Args: 
            arr (np.ndarray): query array.
        Returns: 
            bool.
        '''

        asum = np.sum(arr)
        isnan = np.isnan(asum)
        return isnan

    def removenan(self, arr):

        '''Removes NaN values from array.
        Args: 
            arr (np.ndarray): Target array.
        Returns: 
            arr (np.ndarray): Array stripped from NaN values.
        '''

        arr = arr[~np.isnan(arr).any(axis=1)]
        return arr

    def is_subset(self, a, b):

        '''Checks if a is subset of b. 
        Args: 
            a (np.narray): array a.
            b (np.ndarray): array b.
        Returns:
            bool.
        '''

        b = np.unique(b)
        c = np.intersect1d(a,b)
        return c.size == b.size

    def plot_2d(self, arrays: list):

        '''Plotting function using matplotlib.
        Args:
            arrays (list): List of arrays (np.ndarray) to plot. 
        Returns:
            plt.scatter.
        '''

        for array in arrays:
            plt.scatter(array[:,0], array[:,1], s=1, label=str(array))

        plt.gca().set_aspect('equal')
        

    def plot_3d(self,a,b): 

        '''Interactive plotting function using open3d.
        Args:
            a, b (np.ndarray): arrays to plot. 
        Returns:
            Initializes open3d and plots the fed arrays.
        '''

        import open3d as o3

        zero = np.zeros((a.shape[0], 1))
        a = np.append(a, zero, axis=1)
        zero1 = np.zeros((b.shape[0], 1))
        b = np.append(b, zero1, axis=1)

        a[:,2]=0
        b[:,2]=0

        pcd_a = o3.geometry.PointCloud()
        pcd_a.points = o3.utility.Vector3dVector(a)

        pcd_b = o3.geometry.PointCloud()
        pcd_b.points = o3.utility.Vector3dVector(b)

        pcd_a.paint_uniform_color([1, 0, 0])
        pcd_b.paint_uniform_color([0, 1, 0])
        return(o3.visualization.draw_geometries([pcd_a, pcd_b]))

    def filter(self,source,target, max_dist):
            
        '''Two-way filter the target and source using a KDTree neighbour search. 
        Args: 
            source (np.ndarray): Source array.
            target (np.ndarray): Target array.
            max_dist (float): Neighbour search radius.
        Returns: 
            filteredS (np.ndarray): Filtered source array.
            filteredT (np.ndarray): Filtered target array.
        '''

        treeS = KDTree(target, leaf_size=100)
        distS, indexS = treeS.query(source, 1, return_distance=True)

        treeT = KDTree(source, leaf_size=100)
        distT, indexT = treeT.query(target, 1, return_distance=True)

        # Get filtered matrix
        filteredS = np.take(source, np.where(distS < max_dist)[0], axis=0)
        filteredT = np.take(target, np.where(distT < max_dist)[0], axis=0)

        return filteredS, filteredT

    def select_points(self,a,b):

        import open3d as o3

        print('press "shift" + click points to select')

        zero = np.zeros((a.shape[0], 1))
        a = np.append(a, zero, axis=1)
        zero1 = np.zeros((b.shape[0], 1))
        b = np.append(b, zero1, axis=1)

        a[:,2]=0
        b[:,2]=0

        pcd_a = o3.geometry.PointCloud()
        pcd_a.points = o3.utility.Vector3dVector(a)

        pcd_b = o3.geometry.PointCloud()
        pcd_b.points = o3.utility.Vector3dVector(b)

        pcd_a.paint_uniform_color([1, 0, 0])
        pcd_b.paint_uniform_color([0, 1, 0])

        vis = o3.visualization.draw_geometries_with_editing([pcd_a, pcd_b])
        vis.create_window()
        vis.add_geometry(pcd_a)
        vis.add_geometry(pcd_b)
        vis.run()
        vis.destroy_window()
        print('Selected points (index): ', vis.get_picked_points())

        return vis.get_picked_points()

    def bin(self, array: np.ndarray, chunk_size: list, overlap: list, area: list):

        '''Chunks array into smaller arrays.
        Args: 
            array (np.ndarray): array to chunk.
            chunk_size (list): list that contains the relative size of the chunk on X an Y coordinates: "[X, Y]". 
                               Example: "[0.2, 0.2]" for chunks 0.2 by 0.2 of the lenght and width, respectively.
            overlap (list): list with relative overlap allowed for each chunk, for X and Y: "[X overlap, Y overlap]".
            area (list): Area chunked: [[minX, maxX], [minY, maxY]]".
            Returns: 
                chunks (list): list of chunked arrays (np.ndarray).
        '''

        scanning_array = self.subsample(array, area[0], area[1])
        chunks = []
        x_divisions = [[0,chunk_size[0]]]
        y_divisions = [[0,chunk_size[1]]]
        xmin, xmax = 0, chunk_size[0]
        ymin, ymax = 0, chunk_size[1]
        xsteps = int(1/(chunk_size[0]*(1-overlap[0])))-1
        ysteps = int(1/(chunk_size[1]*(1-overlap[1])))-1

        for x_division in range(xsteps):
            if xmax < 1:
                xmin += chunk_size[0]*(1-overlap[0])
                xmax += chunk_size[0]*(1-overlap[0])
                x_divisions.append([xmin, xmax])

        for y_division in range(ysteps):
            if ymax < 1:
                ymin += chunk_size[1]*(1-overlap[1])
                ymax += chunk_size[1]*(1-overlap[1])
                y_divisions.append([ymin, ymax])
        
        for xrange in x_divisions:
            for yrange in y_divisions:
                chunk = dask.delayed(self.subsample)(array, xrange, yrange)
                chunks.append(chunk)            

        with ProgressBar():
            chunks = dask.compute(chunks)[0]

        new = []
        for index, chunk in enumerate(chunks):
            if self.is_subset(scanning_array, chunk) == True:
                new.append(chunks[index]) 
        chunks = new

        return chunks

    def scan(self, target: np.ndarray, source: np.ndarray, chunk_size: list, plot: bool=False, ref_chunk: list=[[0,1],[0.4,0.5]], overlap: list=[0, 0], area: list=[[0, 1], [0, 1]]):

        '''Scans all source chunks against a specified target reference chunk and 
        returns the best match using the broad evaluator.
        Score provided is the matched/unmatched ratio on the search radius provided.
        Args: 
            target (np.ndarray): target array.
            source (np.ndarray): source array.
            chunk_size (list): list with the X and Y relative sizes of 
                               each chunk: "[Xsize, Ysize]."
            plot (bool): If true, plots each iteration for each chunk 
                         and its position in the source sample.
            ref_chunk (list): list of lists with the relative dimensions 
                              of the reference chunk: "[[minX, maxX], [minY, maxY]]".
            overlap (list): list with relative overlap allowed for 
                            each chunk, for X and Y: "[X overlap, Y overlap]".
            area (list): Area chunked: [[minX, maxX], [minY, maxY]]".
        Returns:
            best_match (np.ndarray): Highest-scoring array. 
            '''

        print("Binning...")
        source_chunks = self.bin(source, chunk_size, overlap=overlap, area=area)
        results = []
        ref = self.subsample(target, ref_chunk[0], ref_chunk[1])


        print("Scanning...")
        for index, chunk in enumerate(source_chunks):
            if len(chunk)/len(ref) >= 0.3: # threshold of 0.1 means that chunk needs to be at least 10% of the length of the reference 
                result = dask.delayed(self.eval_param_worker)(ref, chunk, factor=2/3, evaluator="broad", centering_mode='middle', mode='score')
                results.append(result)
                if plot == True:
                    result = dask.compute(result)[0]
                    print("Scanning chunk: {}".format(index))
                    self.plot_2d([source, chunk])
                    plt.show()
                    print('Score for chunk {}: {}'.format(index, result))
            else: 
                results.append(0) # empty tuple for dask compatibility
                if plot == True:
                    print("Chunk {} is invalid".format(index))

        with ProgressBar():
            results = dask.compute(results)[0]
        index_max = max(range(len(results)), key=results.__getitem__)
        best_match = source_chunks[index_max]
        print('Best match: chunk {}'.format(index_max))

        return best_match
            

class AlignmentPipeline(BeadAlignment):

    def __init__(self, scanning_search_radius: float=2000, filtered_search_radius: float=250, 
                 filtering_radius_1: float=70, filtering_radius_2: float=20, scanning_chunk_size: list=[1, 0.1], 
                 manual_subsampling: bool=False, sub_source: list=[[0,1],[0.4,0.5]],
                 ref_chunk: list=[[0,1],[0.4,0.5]], 
                 overlap: list=[0, 0], scanning_area: list=[[0,1],[0,1]]):

                 '''Initialize AlignmentPipeline class. 
                 
                 Child-class of BeadAlignment that uses its functions to 
                 heuristically employ the distance-angle algorithm in a 
                 sequential fashion.
                 It begins by first centering the datasets by automatically 
                 scanning through the source in a selected area to find 
                 homologous regions to
                 a specified reference in the target dataset. Alternatively, 
                 the user can manually select potentially homologous regions 
                 by setting manual_subsampling to False. 
                 After centering the dataset, the "filter" function will 
                 eliminate noise that lies outside the specified search radius,
                 execute the main alignment algorithm again, then filter a 
                 second time, followed by a final alignment.

                 The alignment can be executed by using the ".run_alignment_pipeline()"
                 function.

                 Args: 
                    scanning_search_radius (float): Search radius for the scanning model to 
                                                    center the datasets. Defaults to 2000.
                    filtered_search_radius (float): Search radius for the (post) filtered 
                                                    model to do the "main" and "adjustment" 
                                                    alignments. Defaults to 250.
                    filtering_radius_1 (float): Radius of coordinates to keep after centering. 
                                                Defaults to 70.
                    filtering_radius_2 (float): Radius of coordinates to keep after main 
                                                alignment round. Defaults to 20.
                    scanning_chunk_size (list): Relative size of each chunk: [Xsize, Ysize].
                    manual_subsampling (bool): True for performing manual subsampling 
                                               reflecting homologous regions instead of 
                                               scanning to find them. Defaults to False.
                    sub_source (list): List of lists with X and Y relative dimensions to 
                                       subsample the source
                                       if manual_subsampling = True.
                    ref_chunk (list): List of lists with X and Y dimensions to be used as the
                                      reference chunk in the target sample for both the scanning 
                                      and manual subsampling mode.
                    overlap (list): List with relative overlap allowed for each chunk, 
                                    for X and Y: "[X overlap, Y overlap]".
                    scanning_area (list): Scanned area. List of lists with X and Y dimensions.

                 '''

                 #Input
                 self.scanning_search_radius = scanning_search_radius
                 self.filtered_search_radius = filtered_search_radius
                 self.filtering_radius_1 = filtering_radius_1
                 self.filtering_radius_2 = filtering_radius_2
                 self.scanning_chunk_size = scanning_chunk_size
                 self.manual_subsampling = manual_subsampling
                 self.sub_source = sub_source
                 self.ref_chunk = ref_chunk
                 self.overlap = overlap
                 self.scanning_area = scanning_area

    def center_datasets_pipeline(self, target: np.ndarray, source: np.ndarray, plot: bool=False, score_centering: bool=False):

        '''
        Workhorse function to initially center the datasets. If manual_subsampling = False, it 
        will scan through the scanning area to find the region homologous to the 
        reference area in the target sample. If manual_subsampling = True, it will
        assume the sub_source provided corresponds to the homologous region to 
        the reference chunk.

        Args:
            target (np.ndarray): Target array.
            source (np.ndarray): Source array.
            plot (bool): True for spatially plotting each iteration. 
            score_centering (bool): Set to True to score the centering (match/unmatched ratio).
        Returns:
            centered (np.ndarray): Centered source array.
            model_centered (class): Returns the centering model to be able to use ".transform()".
        '''

        model_centering = BeadAlignment(initial_scale_factor=(2/3), search_radius=self.scanning_search_radius, search_fraction=0.1, centering_mode='middle')

        # Scan or perform manual subsampling
        if self.manual_subsampling == True:
            print("Subsampling manually...")
            ref = self.subsample(target, self.ref_chunk[0], self.ref_chunk[1])
            best_match = model_centering.subsample(source, self.sub_source[0], self.sub_source[1])
            if plot == True:
                self.plot_2d([target, source+50000, ref, best_match+50000])
                plt.show()
            print("Centering...")
            model_centering.fit(ref, best_match)

        else:
            ref = self.subsample(target, self.ref_chunk[0], self.ref_chunk[1])
            if plot == True:
                print('!!WARNING!!: Plotting hinders parallelization')
                print("Reference target:")
                plt.scatter(target[:,0], target[:,1], label='target', c='g', s=1)
                plt.scatter(ref[:,0], ref[:,1], label='ref_area_target', c='r', s=1)
                plt.gca().set_aspect('equal')
                plt.show()
            best_match = model_centering.scan(target, source, self.scanning_chunk_size, plot, ref_chunk=self.ref_chunk, overlap=self.overlap, area=self.scanning_area)
            print("Centering...")
            model_centering.fit(ref, best_match)

        # Transform based on subsampling to obtain initial centered dataset
        centered = model_centering.transform(source)
        print("Source centered!")
        if score_centering == True:
            centering_score = self.eval_param_worker(target, centered, factor=1, evaluator="broad", centering_mode=np.array([0,0]), mode='score')
            print("Centering score: {}".format(centering_score))
            if centering_score < 0.01:
                print('The matching ratio is low after centering, the alignment will likely fail. Try increasing the scanning chunk size.')
        if plot == True:
            self.plot_2d([target, centered])
            plt.show()

        return centered, model_centering

    def run_pipeline(self, target: np.ndarray, source: np.ndarray, plot: bool=False, score_centering: bool=False):

        '''
        Runs the alignment pipeline. First it removes possible NaN values, then it centers
        the source dataset according to target, followed by a first filtering round. After, 
        it performs the main alignment, filters it again, and concludes by adjusting the
        points by using the second filtering round to bring the points closer together.

        Args: 
            target (np.ndarray): Target array.
            source (np.ndarray): Source array.
            plot (bool): True for spatially plotting each iteration. 
            score_centering (bool): Set to True to score the centering (match/unmatched ratio).
        Returns: 
            transformed (np.ndarray): Transformed dataset. 
            global_transform (function): Transform 2D array with found transfrom function.
        '''

        # 0. Remove possible NaN values
        target = self.removenan(target)
        source = self.removenan(source)

        # 1. Center datasets (pipeline mode)
        centered, model_centering = self.center_datasets_pipeline(target, source, plot, score_centering)

        # 2. filtering based on centered alignment
        centered_filtered, target_filtered = self.filter(centered, target, self.filtering_radius_1)
        print("Filtered based on centered alignment")

        # 3. Whole dataset (filtered)
        print("Main alignment round...")
        model = BeadAlignment(initial_scale_factor=1, search_radius=self.filtered_search_radius, search_fraction=0.1, centering_mode=np.array([0,0]))
        # 3.1. 
        model.fit(target_filtered, centered_filtered)
        # 3.2. 
        transformed = model.transform(centered)


        # 4. 2nd filtering
        transformed_filtered2, target_filtered2 = model.filter(transformed,target,self.filtering_radius_2)
        print("Filtered based on second alignment")
        # 4.1.
        model_filtered = BeadAlignment(initial_scale_factor=1, search_radius=self.filtered_search_radius, search_fraction=0.1, centering_mode=np.array([0,0]))
        # 4.2.
        print("Adjustment alignment round...")
        model_filtered.fit(target_filtered2, transformed_filtered2)
        # 4.3.
        transformed = model_filtered.transform(transformed)

        print("Done!")

        def global_transform(array: np.ndarray):

            # return model_filtered.transform(model.transform(model_centering.transform(array)))   
            return model.transform(model_centering.transform(array))   
            
        return transformed, global_transform

    def plot_2d(self, arrays):
        return super().plot_2d(arrays)

    def plot_3d(self, a, b):
        return super().plot_3d(a, b)

    def plot_initial(self, target, source):

        '''
        Plots target reference area, source scanning area, target, source, subsampled
        source area (if manual_subsampling = True).

        Args:
            target (np.ndarray): Target array.
            source (np.ndarray): Source array.
        Returns:
            plt.scatter
        '''

        ref_chunk = self.subsample(target, self.ref_chunk[0], self.ref_chunk[1])
        if self.manual_subsampling == True:
            sub_source = self.subsample(source, self.sub_source[0], self.sub_source[1])
        scanning_area = self.subsample(source, self.scanning_area[0], self.scanning_area[1])
        if self.manual_subsampling == True:    
            self.plot_2d([target, source + 50000, ref_chunk, sub_source+50000])
        else: 
            self.plot_2d([target, source + 50000, ref_chunk, scanning_area+50000])

#----------------------------------------------END---------------------------------------#