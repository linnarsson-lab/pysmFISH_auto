from ast import excepthandler
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import KDTree, KernelDensity
from typing import Tuple, Union
import pandas as pd
import os
from copy import deepcopy


def align(
    target: np.ndarray, source: np.ndarray, radius: float
) -> Tuple[np.ndarray, np.ndarray]:
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

    # Get coordinates of neighbours
    neigh_coord = [target[i] for i in neighbours]
    # Calculate X and Y delta between point and neighbours
    neigh_delta = [i - j for i, j in zip(neigh_coord, source)]
    # Calculate radians between point and neighbours
    neigh_rad = np.asarray(
        [np.arctan2(i[:, 0], i[:, 1]) for i in neigh_delta], dtype="object"
    )

    return dist, neigh_rad


def flatten(dist: np.ndarray, radians: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
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


def plot_align(
    dist: np.ndarray, radians: np.ndarray, s: float = 5, alpha: float = 0.001
):
    """Plot the alignment results

    Args:
        dist (np.ndarray): Jaggered array with distances for each target
            point to neighbouring source points.
        radians (np.ndarray): Jaggered array with angle in radians for each
            target point to neighbouring source points.
        s (float, optional): Point size. Defaults to 0.5.
        alpha (float, optional): Alpha value of points. Defaults to 0.005.
    """

    # Flatten for plotting
    dist, radians = flatten(dist, radians)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(projection="polar")
    ax.scatter(radians, dist, s=s, alpha=alpha, c="k")


def pol2cart(rho: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """Convert polar coordinates to cartesian.

    Args:
        rho (np.ndarray): Distances.
        phi (np.ndarray): Angles in radians.

    Returns:
        np.ndarray: 2D array with cartesion coordinates.
    """
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return np.array([x, y]).T


def search_points(center: np.ndarray, r: float) -> np.ndarray:
    """Generate 8 points around center, including center.

    Includes center position. will generate 8 points around center at 45 and 90
    degree.

    Args:
        center (np.ndarray): Array with X,Y coordinates of center.
        r (float): Radius.

    Returns:
        np.ndarray: Array with shape (9,2) with X, Y coordinates of points.
    """

    # Offset for points at right angles
    xy0 = np.array([r, 0])
    xy2 = np.array([0, r])
    xy4 = np.array([-r, 0])
    xy6 = np.array([0, -r])

    # Offset for points at 45degree
    a = np.deg2rad(45)
    x1 = np.cos(a) * r
    y1 = np.sin(a) * r
    xy1 = np.array([x1, y1])
    xy3 = np.array([-x1, y1])
    xy5 = np.array([-x1, -y1])
    xy7 = np.array([x1, -y1])

    points = np.array(
        [
            center,
            center + xy0,
            center + xy1,
            center + xy2,
            center + xy3,
            center + xy4,
            center + xy5,
            center + xy6,
            center + xy7,
        ]
    )
    return points


def gradient_acent(
    xy: np.ndarray,
    bandwidth: float = 5,
    search_radius: float = 20,
    max_iterations: int = 500,
    resolution: float = 0.05,
    rtol: float = 1e-3,
    plot: bool = False,
    alpha: float = 0.1,
) -> np.ndarray:
    """Gradient acent on 2D point set.

    Finds highest density area in a set of points.

    Args:
        xy (np.ndarray): Array with shape (n,2) with points.
        bandwidth (float, optional): Bandwidth for KDE with gaussian kernel.
            Should be large enough to completely smooth the points. Meaning
            that there are no local maxima.
            Defaults to 5.
        search_radius (float, optional): Search radius used for the align
            function. Defaults to 20.
        max_iterations (int, optional): Maximum number of iterations. If the
            resolution is small, more iterations will be needed. Defaults to 500.
        resolution (float, optional): Increments per iteration. If resolution
            is small, it will take more iterations. Defaults to 0.05.
        rtol (float, optional): Relative tolerance for sklearn KDE function.
            Defaults to 1e-3.
        plot (bool, optional): If True, plots output and optimization history.
            Defaults to False.
        alpha

    Returns:
        np.ndarray: Array with coordinates of the highest density area of the
            point set.
    """

    # Make KDE
    kde = KernelDensity(
        bandwidth=bandwidth, kernel="gaussian", rtol=rtol, algorithm="ball_tree"
    )
    kde.fit(xy)

    # Initial scoring to find starting point
    step = search_radius / 5
    xiyi = np.meshgrid(
        np.arange(-search_radius, search_radius, step),
        np.arange(-search_radius, search_radius, step),
    )
    xi = xiyi[0].ravel()
    yi = xiyi[1].ravel()
    xy_i = np.stack((xi, yi), axis=1)
    # Data is a circle so remove points outside circle
    r = np.sqrt((xy_i[:, 0]) ** 2 + (xy_i[:, 1]) ** 2)
    outside = r < search_radius
    xy_i = xy_i[outside]
    # Score samples
    score = kde.score_samples(xy_i)
    score_max = score.max()

    # Gradient acent
    start_index = np.where(score == score_max)[0][0]
    center = xy_i[start_index]

    centers = [center]
    scores = [score_max]

    for i in range(max_iterations):
        test_points = search_points(center, r=resolution)
        new_score = kde.score_samples(test_points)
        new_score_max = new_score.max()
        if new_score_max <= score_max:
            break
        else:
            score_max = new_score_max
            center = test_points[np.where(new_score == new_score_max)[0][0]]

            centers.append(center)
            scores.append(score_max)

    cent_final = centers[-1]

    if plot:
        fig, axes = plt.subplots(figsize=(20, 10), ncols=2)
        ax0, ax1 = axes

        ax0.scatter(xy[:, 0], xy[:, 1], s=50, alpha=alpha, label="Input data")
        c = np.array(centers)
        ax0.plot(c[:, 0], c[:, 1], c="k", label="Optimization path")

        ax0.scatter(cent_final[0], cent_final[1], label="Result")
        ax0.set_aspect("equal")
        ax0.set_title(f"Data and result, Final translation: {cent_final}")
        ax0.legend(loc="upper right")

        ax1.plot(scores)
        ax1.set_title(f"Scores, {i} iterations")
        ax1.set_xlabel("Iterations")
        ax1.set_ylabel("Score")

    return cent_final


def find_offset(
    A: np.ndarray,
    B: np.ndarray,
    radius: float = 20,
    max_iterations: int = 500,
    bandwidth: float = 5,
    resolution: float = 0.05,
    initial_offset: np.ndarray = np.array([0, 0]),
    rtol: float = 1e-3,
    plot: bool = False,
    alpha: float = 0.1,
    integer: bool = False,
) -> np.array:
    """Find the offset between two point sets with (some) matching dots.

    Args:
        A (np.ndarray): Point set A.
        B (np.ndarray): Point set B.
        radius (float, optional): Search radius to find matching points between
            A and B. Defaults to 20.
        max_iterations (int, optional): Maximum number of iterations. If the
            resolution is small, more iterations will be needed. Defaults to
            500.
        bandwidth (float, optional): Bandwidth for KDE with gaussian kernel.
            Should be large enough to completely smooth the points.
            Defaults to 5.
        resolution (float, optional): Increments per iteration. If resolution
            is small, it will take more iterations. Defaults to 0.05.
        initial_offset (np.ndarray, optional): If B has already been moved,
            pass that offset so that the offsets can be combined to get the
            total offset value. Defaults to np.array([0,0]).
        rtol (float, optional): Relative tolerance for sklearn KDE function.
            Defaults to 1e-3.
        plot (bool, optional): If True, plots output and optimization history.
            Defaults to False.
        alpha
        integer (bool, optional): If True, transforms the output to integers.
            Defaults to False.

    Returns:
        np.array: Array with X, Y offset values.
    """

    # Find neighbours
    dist, neigh_rad = align(A, B, radius)
    # Flatten output
    dist_flat, rad_flat = flatten(dist, neigh_rad)
    # Translate to Cartesian coordinates
    xy = pol2cart(dist_flat, rad_flat)
    # Find center
    offset = gradient_acent(
        xy,
        bandwidth=bandwidth,
        search_radius=radius,
        max_iterations=max_iterations,
        resolution=resolution,
        rtol=rtol,
        plot=plot,
        alpha=alpha,
    )
    # Flip offset
    offset = np.flip(offset)
    # Change to integer
    if integer:
        offset = np.round(offset).astype("int")
    # Add any previously applied offset.
    offset += initial_offset
    return offset


def find_matches(A: np.ndarray, B: np.ndarray, id_A: np.ndarray, id_B: np.ndarray, offset: np.ndarray, 
                 max_radius: int=10, remove_distinct_genes:bool=False) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Find matching points between two point sets.
    If there are matching points the mean position is returned. 
    
    Args:
        A (np.ndarray): Point set A.
        B_transformed (np.ndarray): Transfromed piont set B so that it aligns
            with A. 
            Use the `find_offset()` function to find the values to transform
            B so that it aligns with A. 
        id_A    
        id_B
        max_radius (int, optional): Max radius within which pionts are 
            considered the same. Defaults to 10.
    Returns:
        Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]: 
            final_positions: Array with all unique points of A, all unique 
                points of B, and the mean position of the matching points of
                A and B.
            A_merged: Array with all new positions of A.
            B_merged: Array with all new positions of B.
            A_filt: Boolean array with all points of A that are also in B.
            B_filt: Boolean array with all points of B that are also in A.
    """
    
    B_transformed = B + offset
    #Find nearest neighbours
    tree = KDTree(A)
    dist, neigh = tree.query(B_transformed, return_distance=True)
    neigh_f, dist_f = neigh.ravel(), dist.ravel()
    
    #Find pairs of matching points that should be merged
    pairs = np.stack((np.arange(B_transformed.shape[0]), neigh[:,0]), axis=1)

    #Exclude points that have multiple matches, keep only those with the shortest distance
    to_del = []
    u, c = np.unique(neigh, return_counts=True)
    dup = u[c > 1]
    for d in dup:
        where = np.where(neigh_f==d)[0]
        filt = dist_f[where] > dist_f[where].min()
        to_del.append(where[filt])

    #Exclude points outside radius
    to_del.append(np.where(dist_f > max_radius)[0])

    #Exclude duplicates and invalid matches
    pairs = np.delete(pairs, np.unique(np.concatenate(to_del)), axis=0)

    #Exclude pairs that do not have the same identity
    #print('pairs',pairs)
    if remove_distinct_genes:
        filt = id_A[pairs[:,1]] == id_A[pairs[:,1]]
        pairs = pairs[filt]
    else:
        filt = id_A[pairs[:,1]] == id_B[pairs[:,0]]
        pairs = pairs[filt]
    
    to_del = []
    u, c = np.unique(pairs[:,1], return_counts=True)
    dup = u[c > 1]
    #print('duplicates',dup)
    for d in dup:
        where = np.where(pairs[:,1]==d)[0]
        to_del.append(where[1:])
        
    #Exclude duplicates and invalid matches
    if len(to_del) > 0:
        pairs = np.delete(pairs, np.unique(np.concatenate(to_del)), axis=0)
    #Filter input datasets
    filt_A = ~np.isin(np.arange(A.shape[0]), pairs[:,1])
    unique_A = A[filt_A]
    filt_B = ~np.isin(np.arange(B.shape[0]), pairs[:,0])
    unique_B = B[filt_B]
    merged_positions = np.mean(np.array([B[pairs[:,0]], A[pairs[:,1]]]), axis=0)

    #Final output of points unique to A, points unique to B and the merged positions of points that are found in both
    final_positions = np.vstack((unique_A, unique_B, merged_positions))
    A_merged = np.vstack((unique_A, merged_positions))
    B_merged = np.vstack((unique_B, merged_positions))
    
    return final_positions, A_merged, B_merged, ~filt_A, ~filt_B

###############################################################################
##  pysmFISH specific functions ##
###############################################################################


def combine_boolean(l: list) -> np.ndarray:
    """Combine multiple boolean arrays with logical and.

    Args:
        l (list): List with arrays. Shapes should match.

    Returns:
        np.ndarray: Combined boolean array.
    """
    filt = np.asarray(l[0])
    if len(l) > 1:
        for ll in l[1:]:
            filt = np.logical_and(filt, np.asarray(ll))
    return filt


def load_parquet(fname: str, columns: list):
    """Load specific columns of a parquet file

    Args:
        fname (str): Filename.
        columns (list): List of column names to fetch.

    Returns:
        Pandas dataframe: Pandas dataframe with requested columns.
    """
    return pd.read_parquet(fname, columns=columns)


def select_overlap(
    rc: np.ndarray, corner_coordinates: Union[list, np.ndarray]
) -> np.ndarray:
    """Return boolean array of points that fall between 4 extrema.

    Args:
        rc (np.ndarray): Array with shape (n,2) with row, column coordinates.
        corner_coordinates (Union[list, np.ndarray]): List or array with
            extrema in the order: [row_min, row_max, col_min, col_max].

    Returns:
        np.ndarray: Boolean array with True for points that fall within
            extrema.
    """
    filt0 = rc[:, 0] >= corner_coordinates[0]
    filt1 = rc[:, 0] <= corner_coordinates[1]
    filt2 = rc[:, 1] >= corner_coordinates[2]
    filt3 = rc[:, 1] <= corner_coordinates[3]
    filt = combine_boolean([filt0, filt1, filt2, filt3])
    return filt

def clip_borders(dataframe:pd.array, bounds:int=0,
        intensity_quantile:float = 0.5) -> pd.array:
    """Clip borders of dataframe to avoid edge effects.
    Args:
        dataframe (pd.array): Dataframe with shape (n,m,p).
        bounds (pd.array): Length to clip the outer frame in pixels.
            Defaults to 0.
        intensity_quantile (float): Quantile to clip the outer frame dots.
    Returns:
        pd.array: Dataframe with shape (n-2,m-2,p).
    """
    if dataframe.shape[0] > 0:

        center_dots = dataframe[(dataframe['r_px_original'] >bounds) & 
                                        (dataframe['r_px_original'] < 2048 -bounds) & 
                                        (dataframe['c_px_original'] > bounds) & 
                                        (dataframe['c_px_original'] < 2048-bounds) ]

        center_dots = center_dots[center_dots.hamming_distance <= 2/16]
        
        return center_dots


def clip_overlap(
    corner_coordinates: Union[list, np.ndarray], fov0: np.ndarray, fov1: np.ndarray
) -> Tuple[str, float, str]:
    """Calculate position to clip two adjacent FOVs.

    Args:
        corner_coordinates (Union[list, np.ndarray]): List or array with
            extrema in the order: [row_min, row_max, col_min, col_max].
        fov0 (np.ndarray): (x,y) coordinates of FOV 0.
        fov1 (np.ndarray): (x,y) coordinates of FOV 1.

    Returns:
        Tuple[str, float, str]:
        index - String  name of column to perform the clipping on.
        cutoff - Float widht the cutoff value.
        orrientation - String describing how FOV 0 is positioned in respect to
            FOV 1.
    """
    r_min, r_max, c_min, c_max = corner_coordinates
    dr = r_max - r_min
    dc = c_max - c_min

    # FOVs overlap left-right
    if dc > dr:
        column = "r_px_microscope_stitched"  # Select on Y (row) coordinate
        cutoff = r_min + (0.05 * dr)
        if fov0[0] > fov1[0]:
            orrientation = "fov0_right_of_fov1"
        else:
            orrientation = "fov0_left_of_fov1"

    # FOVs overlap top-bottom
    else:
        column = "c_px_microscope_stitched"  # Select on X (column) coordinate
        cutoff = c_min + (0.25 * dc)
        if fov0[1] > fov1[1]:
            orrientation = "fov0_above_fov1"
        else:
            orrientation = "fov0_under_fov1"

    return column, cutoff, orrientation


def clean_microscope_stitched(
    fovs: Union[list, np.ndarray],
    fov_df: dict,
    overlapping_regions: dict,
    fov_coords: Union[np.ndarray, None] = None,
    mode: str = "merge",
    bead_channel: str = "Europium",
    max_hamming_dist: float = 3 / 16,
    matching_dot_radius: float = 5,
    out_folder: str = "",
    exp_name: str = "",
    remove_distinct_genes:bool=False,
    clip_size:int = 0, 
    verbose: bool = False,
) -> Tuple[str, str]:
    """Remove overlapping dots in microscope stitched data.

    Caveat: If the same gene name is present in different channels the overlap
    removal likely merges these points. Make sure they have unique names.

    Args:
        fovs (Union[list, np.ndarray]): List or array with FOV numbers.
        fov_df (dict): Dictionary with FOV numbers as keys and dataframe with
            decoded point data as values. Dataframes should contain the
            following columns: `channel`, `hamming_distance`,
            `r_px_microscope_stitched`, `c_px_microscope_stitched`,
            `decoded_genes` and `round_num`.
            It is advised to also include `dot_id` so that output can be
            matched to input.
        overlapping_regions (dict): Dictionary with FOV numbers as keys and as
            value a dictionary with a tuple of neighbouring FOV IDs as keys,
            like (1, 30). And the 4 extrema that define the overlapping region
            between the two FOVs. These extrema should be in the format:
            [row_min, row_max, col_min, col_max]
        fov_coords (Union[np.ndarray, None], optional): When mode is `merge`
            coordinates of fovs need to be given. Coordinates should be a
            (n,2) array in the same order as `fovs`. Defaults to None.
        mode (str, optional): Mode to handle overlap. Options: `merge` or
            `clip`. If `merge`, will discard points that are found in both FOVs
            and keep unique points. If `clip`, will clip a FOV up to the middle
            of the overlap. Defaults to 'merge'.
        bead_channel (str, optional): Name of bead channel to exclude.
            If nothing needs to be excluded set to a nonsense name.
            Defaults to 'Europium'.
        max_hammin_dist (float, optional): Maximum hamming distance for a point
            to be considered valid. Defaults to 3/16.
        matching_dot_radius (float, optional): Radius to look for matching
            points for removal. Search will be conducted after alignment.
            Defaults to 5.
        out_folder (str, optional): Path to output folder. Needs to have
            trailing slash. Defaults to '/'.
        exp_name (str, optional): Name of experiment to name the output files.
            Defaults to ''.
        remove_distinct_genes (bool, optional): If True overlapping dots from
            different genes will be removed. Defaults to False.
        clip_size (int, optional): If more than 0 the FOVs outmost pixels will
            be clipped to avoid edge effects. Defaults to 0.
        verbose (bool, optional): If True prints progress. Defaults to False.

    Returns:
        Tuple[str, str]
        filename (str): Path and name of merged RNA output file.
        filename (str): Path and name of merge large beads output file.

        Output is saved:
        - For each FOV a reduced dataframe is saved in the output folder with
          the extension: _microscope_merged_fov_X.parquet where X is the FOV
          number. This dataframe contains all RNA points with new locations
          if they had a match with an adjacent FOV. In the other FOV these
          points will have been deleted.
        - Merged dataframe with all the points after merging. File extension
          will be: _RNA_microscope_stitched.parquet. and located in the 
          out_folder. The file will contain the following columns: 
          `r_px_microscope_stitched`, `c_px_microscope_stitched`,
          `hamming_distance` and `decoded_genes`.
        - Merged dataframe with all the large bead coordinates after merging.
          File extention will be: _large_beads_microscope_stitched.parquet. The
          file will contain the the following columns: 
          `r_px_microscope_stitched` and `c_px_microscope_stitched`.

    """

    valid_modes = ['clip', 'merge']
    if mode not in valid_modes:
        raise Exception(f'Dot removal mode not recognized. Input mode: {mode} is not part of implemented modes: {valid_modes}')
    if mode == 'clip':
        clip_size=0
        
    # Fetch all large beads (Does not run bead removal in overlapping regions)
    large_beads_list = []
    for i in fovs:
        bead_data = fov_df[i]
        large_bead_filt = combine_boolean([bead_data.channel == bead_channel, 
                                           bead_data.mapped_beads_type == "large",
                                           bead_data.round_num == 1])
        bead_data = bead_data.loc[large_bead_filt]
        large_beads_list.append(bead_data)
    large_beads_df = pd.concat(large_beads_list)
    large_beads_df = large_beads_df.dropna(subset=['r_px_microscope_stitched','c_px_microscope_stitched','mapped_beads_type'])

    # Merge overlapping points.
    for i in fovs:
        # Get neighbours of FOV
        fov_neigbours = overlapping_regions[i]

        # Check if FOV has neighbours
        if fov_neigbours:
            key, values = zip(*fov_neigbours.items())
            # Iterate through neighbours
            for k, v in zip(key, values):
                pair, corner_coordinates = k, v
                fov0, fov1 = pair
                if verbose:
                    print(
                        f"FOV :{i}/{fovs[-1]}. Pair: {fov0} - {fov1}           ",
                        end="\r",
                    )

                # Get FOV 0 RNA
                d0 = fov_df[fov0]
                # Load all valid RNA signal
                d0_filt_rna = np.logical_and(
                    d0.channel != bead_channel, d0.hamming_distance < max_hamming_dist
                )
                d0_filt = deepcopy(d0_filt_rna)
                d0_rc = d0.loc[
                    d0_filt, ["r_px_microscope_stitched", "c_px_microscope_stitched"]
                ].to_numpy()
                # Select RNA in overlap
                filt0_overlap = select_overlap(d0_rc, corner_coordinates)
                d0_rc_overlap = d0_rc[filt0_overlap]
                d0_filt[d0_filt] = filt0_overlap
                # Get gene ID of each point
                d0_id = d0.loc[d0_filt, "decoded_genes"].to_numpy()

                # Get FOV 1 RNA
                d1 = fov_df[fov1]
                # Load all valid RNA signal
                d1_filt_rna = np.logical_and(
                    d1.channel != bead_channel, d1.hamming_distance < max_hamming_dist
                )
                d1_filt = deepcopy(d1_filt_rna)
                d1_rc = d1.loc[
                    d1_filt, ["r_px_microscope_stitched", "c_px_microscope_stitched"]
                ].to_numpy()
                # Select RNA in overlap
                filt1_overlap = select_overlap(d1_rc, corner_coordinates)
                d1_rc_overlap = d1_rc[filt1_overlap]
                d1_filt[d1_filt] = filt1_overlap
                # Get gene ID of each point
                d1_id = d1.loc[d1_filt, "decoded_genes"].to_numpy()

                if mode == "merge":
                    # Check if there are overlapping points
                    if d0_rc_overlap.shape[0] > 10 and d1_rc_overlap.shape[0] > 10:

                        # Find offset between FOVs
                        offset = find_offset(
                            d0_rc_overlap, d1_rc_overlap, resolution=0.05
                        )
                        # Find matching dots between FOVs
                        (
                            final_positions,
                            d0_merged,
                            d1_merged,
                            d0_match_filt,
                            d1_match_filt,
                        ) = find_matches(
                            d0_rc_overlap,
                            d1_rc_overlap,
                            d0_id,
                            d1_id,
                            offset,
                            max_radius=matching_dot_radius,
                            remove_distinct_genes=remove_distinct_genes
                        )

                        # Merged points will be put only in dataset 0, and deleted from dataset 1
                        # Correct new positions of merged points in dataset 0
                        df0_cleaned = d0[d0_filt_rna]
                        df0_cleaned.loc[
                            filt0_overlap,
                            ["r_px_microscope_stitched", "c_px_microscope_stitched"],
                        ] = d0_merged
                        fov_df[fov0] = df0_cleaned

                        # Delet positions in dataset 1 that are merged and should only be in dataset 0
                        df1_cleaned = d1[d1_filt_rna]
                        index_to_drop = df1_cleaned[filt1_overlap].index[d1_match_filt]
                        df1_cleaned = df1_cleaned.drop(index_to_drop)
                        fov_df[fov1] = df1_cleaned

                    # Overlap could not be determined, just delete non-(valid-)RNA points
                    else:
                        # Dataset 0
                        df0_cleaned = d0[d0_filt_rna]
                        fov_df[fov0] = df0_cleaned

                        # Dataset 1
                        df1_cleaned = d1[d1_filt_rna]
                        fov_df[fov1] = df1_cleaned

                elif mode == "clip":
                    column, cutoff, orrientation = clip_overlap(
                        corner_coordinates, fov_coords[fov0], fov_coords[fov1]
                    )
                    # Clean RNA
                    df0_cleaned = d0[d0_filt_rna]
                    df1_cleaned = d1[d1_filt_rna]
                    # Get clips
                    if orrientation == "fov0_under_fov1":
                        df0_clip_filt = df0_cleaned.loc[:, column] < cutoff
                        df1_clip_filt = df1_cleaned.loc[:, column] > cutoff
                    elif orrientation == "fov0_above_fov1":
                        df0_clip_filt = df0_cleaned.loc[:, column] > cutoff
                        df1_clip_filt = df1_cleaned.loc[:, column] < cutoff
                    elif orrientation == "fov0_right_of_fov1":
                        df0_clip_filt = df0_cleaned.loc[:, column] > cutoff
                        df1_clip_filt = df1_cleaned.loc[:, column] < cutoff
                    elif orrientation == "fov0_left_of_fov1":
                        df0_clip_filt = df0_cleaned.loc[:, column] < cutoff
                        df1_clip_filt = df1_cleaned.loc[:, column] > cutoff
                    # Clip datasets
                    fov_df[fov0] = df0_cleaned.loc[df0_clip_filt]
                    fov_df[fov1] = df1_cleaned.loc[df1_clip_filt]

        # FOV has no neighbours. Just delete beads and non-valid-RNA points
        else:
            d = fov_df[i]
            d_filt_rna = np.logical_and(
                d.channel != bead_channel, d.hamming_distance < max_hamming_dist
            )
            fov_df[i] = d[d_filt_rna]

    # Save individual FOV RNA results
    for i in fovs:
        fname = os.path.join(
            out_folder, (exp_name + f"_microscope_merged_fov_{i}.parquet")
        )
        fov_df[i].to_parquet(fname)

    # Merge RNA results
    #merged_df = pd.concat([fov_df[i] for i in fovs])
    merged_df = pd.concat([clip_borders(fov_df[i], clip_size) for i in fovs])
    merged_df = merged_df.dropna(subset=['r_px_microscope_stitched','c_px_microscope_stitched','decoded_genes'])

    # Save RNA results
    fname_rna_merge = os.path.join(out_folder, (exp_name + "_RNA_microscope_stitched.parquet"),)

    try:
        merged_df.loc[
            :,
            [
                "r_px_microscope_stitched",
                "c_px_microscope_stitched",
                "hamming_distance",
                "decoded_genes",
                "channel",
                "round_num",
                "raw_barcodes",
                "mapped_beads_type",
                "bit_1_intensity",	
                "bit_2_intensity",	
                "bit_3_intensity",
                "bit_4_intensity",
                "bit_5_intensity",
                "bit_6_intensity",
                "bit_7_intensity",
                "bit_8_intensity",
                "bit_9_intensity",
                "bit_10_intensity",	
                "bit_11_intensity",	
                "bit_12_intensity",	
                "bit_13_intensity",
                "bit_14_intensity",
                "bit_15_intensity",
                "bit_16_intensity",
            ],
        ].to_parquet(fname_rna_merge)
    except:
        merged_df.loc[
            :,
            [
                "r_px_microscope_stitched",
                "c_px_microscope_stitched",
                "hamming_distance",
                "decoded_genes",
            ],
        ].to_parquet(fname_rna_merge)
    
    # Save Large bead results
    fname_large_beads = os.path.join(out_folder, (exp_name + "_large_beads_microscope_stitched.parquet"))
    large_beads_df.loc[:,["r_px_microscope_stitched", "c_px_microscope_stitched"]].to_parquet(fname_large_beads)
    
    if verbose:
        print(f"Merged RNA output saved in: {fname_rna_merge}")
        print(f"Merged Large bead output saved in: {fname_large_beads}")

    return fname_rna_merge, fname_large_beads


def clean_microscope_stitched2(
    fovs: Union[list, np.ndarray],
    fov_df: dict,
    overlapping_regions: dict,
    fov_coords: Union[np.ndarray, None] = None,
    mode: str = "merge",
    bead_channel: str = "Europium",
    max_hamming_dist: float = 3 / 16,
    matching_dot_radius: float = 5,
    out_folder: str = "",
    exp_name: str = "",
    remove_distinct_genes:bool=False,
    clip_size:int = 0, 
    verbose: bool = False,
) -> Tuple[str, str]:
    """Remove overlapping dots in microscope stitched data.

    Caveat: If the same gene name is present in different channels the overlap
    removal likely merges these points. Make sure they have unique names.

    Args:
        fovs (Union[list, np.ndarray]): List or array with FOV numbers.
        fov_df (dict): Dictionary with FOV numbers as keys and dataframe with
            decoded point data as values. Dataframes should contain the
            following columns: `channel`, `hamming_distance`,
            `r_px_microscope_stitched`, `c_px_microscope_stitched`,
            `decoded_genes` and `round_num`.
            It is advised to also include `dot_id` so that output can be
            matched to input.
        overlapping_regions (dict): Dictionary with FOV numbers as keys and as
            value a dictionary with a tuple of neighbouring FOV IDs as keys,
            like (1, 30). And the 4 extrema that define the overlapping region
            between the two FOVs. These extrema should be in the format:
            [row_min, row_max, col_min, col_max]
        fov_coords (Union[np.ndarray, None], optional): When mode is `merge`
            coordinates of fovs need to be given. Coordinates should be a
            (n,2) array in the same order as `fovs`. Defaults to None.
        mode (str, optional): Mode to handle overlap. Options: `merge` or
            `clip`. If `merge`, will discard points that are found in both FOVs
            and keep unique points. If `clip`, will clip a FOV up to the middle
            of the overlap. Defaults to 'merge'.
        bead_channel (str, optional): Name of bead channel to exclude.
            If nothing needs to be excluded set to a nonsense name.
            Defaults to 'Europium'.
        max_hammin_dist (float, optional): Maximum hamming distance for a point
            to be considered valid. Defaults to 3/16.
        matching_dot_radius (float, optional): Radius to look for matching
            points for removal. Search will be conducted after alignment.
            Defaults to 5.
        out_folder (str, optional): Path to output folder. Needs to have
            trailing slash. Defaults to '/'.
        exp_name (str, optional): Name of experiment to name the output files.
            Defaults to ''.
        remove_distinct_genes (bool, optional): If True overlapping dots from
            different genes will be removed. Defaults to False.
        clip_size (int, optional): If more than 0 the FOVs outmost pixels will
            be clipped to avoid edge effects. Defaults to 0.
        verbose (bool, optional): If True prints progress. Defaults to False.

    Returns:
        Tuple[str, str]
        filename (str): Path and name of merged RNA output file.
        filename (str): Path and name of merge large beads output file.

        Output is saved:
        - For each FOV a reduced dataframe is saved in the output folder with
          the extension: _microscope_merged_fov_X.parquet where X is the FOV
          number. This dataframe contains all RNA points with new locations
          if they had a match with an adjacent FOV. In the other FOV these
          points will have been deleted.
        - Merged dataframe with all the points after merging. File extension
          will be: _RNA_microscope_stitched.parquet. and located in the 
          out_folder. The file will contain the following columns: 
          `r_px_microscope_stitched`, `c_px_microscope_stitched`,
          `hamming_distance` and `decoded_genes`.
        - Merged dataframe with all the large bead coordinates after merging.
          File extention will be: _large_beads_microscope_stitched.parquet. The
          file will contain the the following columns: 
          `r_px_microscope_stitched` and `c_px_microscope_stitched`.

    """

    valid_modes = ['clip', 'merge']
    if mode not in valid_modes:
        raise Exception(f'Dot removal mode not recognized. Input mode: {mode} is not part of implemented modes: {valid_modes}')
    if mode == 'clip':
        clip_size=0
        
    # Fetch all large beads (Does not run bead removal in overlapping regions)
    large_beads_list = []
    for i in fovs:
        bead_data = fov_df[i]
        large_bead_filt = combine_boolean([bead_data.channel == bead_channel, 
                                           bead_data.mapped_beads_type == "large",
                                           bead_data.round_num == 1])
        bead_data = bead_data.loc[large_bead_filt]
        large_beads_list.append(bead_data)
    large_beads_df = pd.concat(large_beads_list)
    large_beads_df = large_beads_df.dropna(subset=['r_px_microscope_stitched','c_px_microscope_stitched','mapped_beads_type'])

    # Merge overlapping points.
    for reremove in [1,2]:
        for i in fovs:
            # Get neighbours of FOV
            fov_neigbours = overlapping_regions[i]

            # Check if FOV has neighbours
            if fov_neigbours:
                key, values = zip(*fov_neigbours.items())
                # Iterate through neighbours
                for k, v in zip(key, values):
                    pair, corner_coordinates = k, v
                    fov0, fov1 = pair
                    if verbose:
                        print(
                            f"FOV :{i}/{fovs[-1]}. Pair: {fov0} - {fov1}           ",
                            end="\r",
                        )

                    # Get FOV 0 RNA
                    d0 = fov_df[fov0]
                    # Load all valid RNA signal
                    d0_filt_rna = np.logical_and(
                        d0.channel != bead_channel, d0.hamming_distance < max_hamming_dist
                    )
                    d0_filt = deepcopy(d0_filt_rna)
                    d0_rc = d0.loc[
                        d0_filt, ["r_px_microscope_stitched", "c_px_microscope_stitched"]
                    ].to_numpy()
                    # Select RNA in overlap
                    filt0_overlap = select_overlap(d0_rc, corner_coordinates)
                    d0_rc_overlap = d0_rc[filt0_overlap]
                    d0_filt[d0_filt] = filt0_overlap
                    # Get gene ID of each point
                    d0_id = d0.loc[d0_filt, "decoded_genes"].to_numpy()

                    # Get FOV 1 RNA
                    d1 = fov_df[fov1]
                    # Load all valid RNA signal
                    d1_filt_rna = np.logical_and(
                        d1.channel != bead_channel, d1.hamming_distance < max_hamming_dist
                    )
                    d1_filt = deepcopy(d1_filt_rna)
                    d1_rc = d1.loc[
                        d1_filt, ["r_px_microscope_stitched", "c_px_microscope_stitched"]
                    ].to_numpy()
                    # Select RNA in overlap
                    filt1_overlap = select_overlap(d1_rc, corner_coordinates)
                    d1_rc_overlap = d1_rc[filt1_overlap]
                    d1_filt[d1_filt] = filt1_overlap
                    # Get gene ID of each point
                    d1_id = d1.loc[d1_filt, "decoded_genes"].to_numpy()

                    if mode == "merge":
                        # Check if there are overlapping points
                        if d0_rc_overlap.shape[0] > 10 and d1_rc_overlap.shape[0] > 10:

                            # Find offset between FOVs
                            offset = find_offset(
                                d0_rc_overlap, d1_rc_overlap, resolution=0.05
                            )
                            # Find matching dots between FOVs
                            (
                                final_positions,
                                d0_merged,
                                d1_merged,
                                d0_match_filt,
                                d1_match_filt,
                            ) = find_matches(
                                d0_rc_overlap,
                                d1_rc_overlap,
                                d0_id,
                                d1_id,
                                offset,
                                max_radius=matching_dot_radius*reremove,
                                remove_distinct_genes=remove_distinct_genes
                            )

                            # Merged points will be put only in dataset 0, and deleted from dataset 1
                            # Correct new positions of merged points in dataset 0
                            df0_cleaned = d0[d0_filt_rna]
                            df0_cleaned.loc[
                                filt0_overlap,
                                ["r_px_microscope_stitched", "c_px_microscope_stitched"],
                            ] = d0_merged
                            fov_df[fov0] = df0_cleaned

                            # Delet positions in dataset 1 that are merged and should only be in dataset 0
                            df1_cleaned = d1[d1_filt_rna]
                            index_to_drop = df1_cleaned[filt1_overlap].index[d1_match_filt]
                            df1_cleaned = df1_cleaned.drop(index_to_drop)
                            fov_df[fov1] = df1_cleaned

                        # Overlap could not be determined, just delete non-(valid-)RNA points
                        else:
                            # Dataset 0
                            df0_cleaned = d0[d0_filt_rna]
                            fov_df[fov0] = df0_cleaned

                            # Dataset 1
                            df1_cleaned = d1[d1_filt_rna]
                            fov_df[fov1] = df1_cleaned

                    elif mode == "clip":
                        column, cutoff, orrientation = clip_overlap(
                            corner_coordinates, fov_coords[fov0], fov_coords[fov1]
                        )
                        # Clean RNA
                        df0_cleaned = d0[d0_filt_rna]
                        df1_cleaned = d1[d1_filt_rna]
                        # Get clips
                        if orrientation == "fov0_under_fov1":
                            df0_clip_filt = df0_cleaned.loc[:, column] < cutoff
                            df1_clip_filt = df1_cleaned.loc[:, column] > cutoff
                        elif orrientation == "fov0_above_fov1":
                            df0_clip_filt = df0_cleaned.loc[:, column] > cutoff
                            df1_clip_filt = df1_cleaned.loc[:, column] < cutoff
                        elif orrientation == "fov0_right_of_fov1":
                            df0_clip_filt = df0_cleaned.loc[:, column] > cutoff
                            df1_clip_filt = df1_cleaned.loc[:, column] < cutoff
                        elif orrientation == "fov0_left_of_fov1":
                            df0_clip_filt = df0_cleaned.loc[:, column] < cutoff
                            df1_clip_filt = df1_cleaned.loc[:, column] > cutoff
                        # Clip datasets
                        fov_df[fov0] = df0_cleaned.loc[df0_clip_filt]
                        fov_df[fov1] = df1_cleaned.loc[df1_clip_filt]

            # FOV has no neighbours. Just delete beads and non-valid-RNA points
            else:
                d = fov_df[i]
                d_filt_rna = np.logical_and(
                    d.channel != bead_channel, d.hamming_distance < max_hamming_dist
                )
                fov_df[i] = d[d_filt_rna]

    # Save individual FOV RNA results
    for i in fovs:
        fname = os.path.join(
            out_folder, (exp_name + f"_microscope_merged_fov_{i}.parquet")
        )
        fov_df[i].to_parquet(fname)

    # Merge RNA results
    #merged_df = pd.concat([fov_df[i] for i in fovs])
    merged_df = pd.concat([clip_borders(fov_df[i], clip_size) for i in fovs])
    merged_df = merged_df.dropna(subset=['r_px_microscope_stitched','c_px_microscope_stitched','decoded_genes'])

    # Save RNA results
    fname_rna_merge = os.path.join(out_folder, (exp_name + "_RNA_microscope_stitched.parquet"),)

    try:
        merged_df.loc[
            :,
            [
                "r_px_microscope_stitched",
                "c_px_microscope_stitched",
                "hamming_distance",
                "decoded_genes",
                #"channel",
                #"round_num",
                #"raw_barcodes",
                #"mapped_beads_type",
                #"bit_1_intensity",	
                #"bit_2_intensity",	
                #"bit_3_intensity",
                #"bit_4_intensity",
                #"bit_5_intensity",
                #"bit_6_intensity",
                #"bit_7_intensity",
                #"bit_8_intensity",
                #"bit_9_intensity",
                #"bit_10_intensity",	
                #"bit_11_intensity",	
                #"bit_12_intensity",	
                #"bit_13_intensity",
                #"bit_14_intensity",
                #"bit_15_intensity",
                #"bit_16_intensity",
            ],
        ].to_parquet(fname_rna_merge)
    except:
        merged_df.loc[
            :,
            [
                "r_px_microscope_stitched",
                "c_px_microscope_stitched",
                "hamming_distance",
                "decoded_genes",
            ],
        ].to_parquet(fname_rna_merge)
    
    # Save Large bead results
    fname_large_beads = os.path.join(out_folder, (exp_name + "_large_beads_microscope_stitched.parquet"))
    large_beads_df.loc[:,["r_px_microscope_stitched", "c_px_microscope_stitched"]].to_parquet(fname_large_beads)
    
    if verbose:
        print(f"Merged RNA output saved in: {fname_rna_merge}")
        print(f"Merged Large bead output saved in: {fname_large_beads}")

    return fname_rna_merge, fname_large_beads