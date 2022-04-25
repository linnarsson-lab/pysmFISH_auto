import pickle
import itertools
import numpy as np
import zarr

from pathlib import Path
from skimage import measure


def load_segmented_data(fov_subdataset, experiment_path):
    experiment_path = Path(experiment_path)
    experiment_name = fov_subdataset.experiment_name
    segmented_file_path = (
        experiment_path / "fresh_tissue" / "fresh_tissue_segmented_nuclei_data.zarr"
    )
    st = zarr.DirectoryStore(segmented_file_path)
    root = zarr.group(store=st, overwrite=False)
    grp_name = (
        experiment_name + "_segmetented_fresh_tissue_fov_" + str(fov_subdataset.fov_num)
    )

    dataset_name = "segmentation_mask_fov_" + str(fov_subdataset.fov_num)

    mask = root[grp_name][dataset_name][...]
    return mask


def load_stitched_segmented_data(fov_subdataset, segmentation_output_path):
    segmented_region_img = pickle.load(
        open(
            segmentation_output_path
            / ("registered_objs_dict_fov_" + str(fov_subdataset.fov_num) + ".pkl"),
            "rb",
        )
    )

    obj_segmented_img = {
        el: coords_dict["stitched_coords"]
        for (el, coords_dict) in segmented_region_img.items()
    }
    return obj_segmented_img


def load_stitched_segmented_data_full_dict(fov_subdataset, segmentation_output_path):
    segmented_region_dict = pickle.load(
        open(
            segmentation_output_path
            / ("registered_objs_dict_fov_" + str(fov_subdataset.fov_num) + ".pkl"),
            "rb",
        )
    )
    return segmented_region_dict


def stitch_using_coords_general_segmented_objects(
    fov, obj_dict, tile_corners_coords_pxl, reference_corner_fov_position, metadata
):
    """
    Tiles are placed directly on the position indicated by the microscope
    coords
    """

    r_coords = tile_corners_coords_pxl[fov, 0]
    c_coords = tile_corners_coords_pxl[fov, 1]

    if obj_dict:

        if reference_corner_fov_position == "top-left":
            for el, coords_dict in obj_dict.items():
                coords_dict["stitched_coords"] = np.vstack(
                    [
                        r_coords + coords_dict["original_coords"][:, 0],
                        c_coords + coords_dict["original_coords"][:, 1],
                    ]
                ).T
                coords_dict["stitched_centroid"] = np.vstack(
                    [
                        r_coords + coords_dict["original_centroid"][0],
                        c_coords + coords_dict["original_centroid"][1],
                    ]
                ).T

        elif reference_corner_fov_position == "top-right":
            for el, coords_dict in obj_dict.items():
                coords_dict["stitched_coords"] = np.vstack(
                    [
                        r_coords + coords_dict["original_coords"][:, 0],
                        c_coords
                        - (
                            metadata["img_width"] - coords_dict["original_coords"][:, 1]
                        ),
                    ]
                ).T
                coords_dict["stitched_centroid"] = np.vstack(
                    [
                        r_coords + coords_dict["original_centroid"][0],
                        c_coords
                        - (metadata["img_width"] - coords_dict["original_centroid"][1]),
                    ]
                ).T

        elif reference_corner_fov_position == "bottom_left":
            for el, coords_dict in obj_dict.items():
                coords_dict["stitched_coords"] = np.vstack(
                    [
                        r_coords
                        + (
                            metadata["img_height"]
                            - coords_dict["original_coords"][:, 0]
                        ),
                        c_coords + coords_dict["original_coords"][:, 1],
                    ]
                ).T
                coords_dict["stitched_centroid"] = np.vstack(
                    [
                        r_coords
                        + (
                            metadata["img_height"] - coords_dict["original_centroid"][0]
                        ),
                        c_coords + coords_dict["original_centroid"][1],
                    ]
                ).T

    return obj_dict


def register_coords_obj(
    fov_subdataset,
    experiment_fpath,
    segmentation_output_path,
    tile_corners_coords_global_pxl,
    metadata,
    reference_corner_fov_position="top-right",
):

    mask = load_segmented_data(fov_subdataset, experiment_fpath)
    segmentation_output_path = Path(experiment_fpath) / "fresh_tissue" / "segmentation"

    segmented_regions = measure.regionprops(mask)
    segmented_regions_dict = {}
    for prop in segmented_regions:
        segmented_regions_dict[str(fov_subdataset.fov_num) + "-" + str(prop.label)] = {}
        segmented_regions_dict[str(fov_subdataset.fov_num) + "-" + str(prop.label)][
            "original_coords"
        ] = prop.coords
        segmented_regions_dict[str(fov_subdataset.fov_num) + "-" + str(prop.label)][
            "stitched_coords"
        ] = np.nan
        segmented_regions_dict[str(fov_subdataset.fov_num) + "-" + str(prop.label)][
            "original_centroid"
        ] = prop.centroid
        segmented_regions_dict[str(fov_subdataset.fov_num) + "-" + str(prop.label)][
            "stitched_centroid"
        ] = prop.centroid
        segmented_regions_dict[str(fov_subdataset.fov_num) + "-" + str(prop.label)][
            "area"
        ] = prop.area
        segmented_regions_dict = stitch_using_coords_general_segmented_objects(
            fov_subdataset.fov_num,
            segmented_regions_dict,
            tile_corners_coords_global_pxl,
            reference_corner_fov_position,
            metadata,
        )
    pickle.dump(
        segmented_regions_dict,
        open(
            segmentation_output_path
            / ("registered_objs_dict_fov_" + str(fov_subdataset.fov_num) + ".pkl"),
            "wb",
        ),
    )


#    return segmented_regions_dict


def objects_overlapping_region(obj_segmented_img, chunk_coords):
    objs_list = list()
    for lab in obj_segmented_img.keys():
        if (
            (obj_segmented_img[lab][:, 0] >= chunk_coords[0]).any()
            and (obj_segmented_img[lab][:, 0] <= chunk_coords[1]).any()
            and (obj_segmented_img[lab][:, 1] >= chunk_coords[2]).any()
            and (obj_segmented_img[lab][:, 1] <= chunk_coords[3]).any()
        ):
            objs_list.append(lab)
    return objs_list


def overlapping_objs(
    obj_segmented_img1,
    obj_segmented_img2,
    chunk_coords,
    min_overlapping_pixels_segmentation=20,
):
    objs_list_img1 = objects_overlapping_region(obj_segmented_img1, chunk_coords)
    objs_list_img2 = objects_overlapping_region(obj_segmented_img2, chunk_coords)

    all_combinations = list(itertools.product(objs_list_img1, objs_list_img2))

    intersection = []
    overlapping_sizes = {}
    for couple in all_combinations:
        obj1_set = set(map(tuple, obj_segmented_img1[couple[0]].astype(int)))
        obj2_set = set(map(tuple, obj_segmented_img2[couple[1]].astype(int)))
        result = not obj1_set.isdisjoint(obj2_set)
        if result:
            # Check minimum number of pixels to correct for potential error caused by the rounding
            number_overlapping_pixels = len(obj1_set.intersection(obj2_set))
            if number_overlapping_pixels >= min_overlapping_pixels_segmentation:
                intersection.append(result)
                overlapping_sizes[couple] = number_overlapping_pixels
            else:
                intersection.append(not result)
        else:
            intersection.append(result)
    if (np.array(intersection)).any():
        all_combinations = np.array(all_combinations)
        intersecting_cpls = all_combinations[intersection]
    else:
        intersecting_cpls = []
    return intersecting_cpls, overlapping_sizes


def correct_multiple_overlapping_objects(intersecting_cpls, overlapping_sizes):

    intersecting_cpls_tpl_trimmed = []
    intersecting_cpls_tpl = [tuple(x) for x in intersecting_cpls]

    intersecting_cpls_removed_overlapping = intersecting_cpls_tpl.copy()
    tracker = intersecting_cpls_tpl.copy()

    combined_intersecting_cpls_tpl = []
    tracker = []
    reducer = intersecting_cpls_tpl.copy()
    for tpl in intersecting_cpls_tpl:
        if tpl not in tracker:
            combined = [tpl]
            tracker.append(tpl)
            reducer.remove(tpl)
            runner = reducer
            for cpl in runner:
                if set(tpl).intersection(cpl):
                    tracker.append(cpl)
                    combined.append(cpl)
                    reducer.remove(cpl)
            combined_intersecting_cpls_tpl.append(combined)

    # combined_intersecting_cpls_tpl = set(combined_intersecting_cpls_tpl)
    for grp in combined_intersecting_cpls_tpl:
        if len(grp) > 1:
            # Select the object with largest overlap (maybe not best criteria)
            selected = {your_key: overlapping_sizes[your_key] for your_key in grp}
            counts = list(selected.values())
            keys = list(selected.keys())
            selected_cpl = keys[counts.index(max(counts))]
            intersecting_cpls_tpl_trimmed.append(selected_cpl)
        else:
            intersecting_cpls_tpl_trimmed.append(grp)

    return intersecting_cpls_removed_overlapping, overlapping_sizes


def merge_objects(obj_segmented_img1, obj_segmented_img2, intersecting_cpls):
    merged_objs = {}
    for cpl in intersecting_cpls:
        obj1 = obj_segmented_img1[cpl[0]].astype(int)
        obj2 = obj_segmented_img2[cpl[1]].astype(int)
        merged_obj = np.vstack((obj1, obj2))
        merged_obj = np.unique(merged_obj, axis=0)
        merged_objs[cpl[0]] = merged_obj
    return merged_objs


def adjust_overlapping_objects(
    obj_segmented_img1,
    obj_segmented_img2,
    chunk_coords,
    segmentation_output_path,
    min_overlapping_pixels_segmentation=20,
):

    intersecting_cpls, overlapping_sizes = overlapping_objs(
        obj_segmented_img1,
        obj_segmented_img2,
        chunk_coords,
        min_overlapping_pixels_segmentation,
    )

    # Correct for multiple overlapping large objects
    intersecting_cpls, overlapping_sizes = correct_multiple_overlapping_objects(
        intersecting_cpls, overlapping_sizes
    )

    objs = merge_objects(obj_segmented_img1, obj_segmented_img2, intersecting_cpls)

    # for cpl in intersecting_cpls:
    #   del obj_segmented_img1[cpl[0]]
    #   del obj_segmented_img2[cpl[1]]
    #   obj_segmented_img1[cpl[0]] = objs[cpl[0]]

    return intersecting_cpls, objs


def remove_overlapping_obj(
    fov_subdataset,
    nuclei_org_tiles,
    segmentation_output_path,
    min_overlapping_pixels_segmentation=20,
):

    fov = fov_subdataset.fov_num
    overlapping_dict = nuclei_org_tiles.overlapping_regions[fov]
    all_obj_segmented = {}
    output_cpl = {}
    corner_cpl = []

    # Correct overlapping objects at the intersection between the images
    for cpl, chunk_coords in overlapping_dict.items():
        obj_segmented_img1 = load_stitched_segmented_data(
            cpl[0], segmentation_output_path
        )
        obj_segmented_img2 = load_stitched_segmented_data(
            cpl[1], segmentation_output_path
        )
        intersecting_cpls, objs = adjust_overlapping_objects(
            obj_segmented_img1,
            obj_segmented_img2,
            chunk_coords,
            segmentation_output_path,
            min_overlapping_pixels_segmentation,
        )
        all_obj_segmented[cpl[0]] = obj_segmented_img1
        all_obj_segmented[cpl[1]] = obj_segmented_img2
        output_cpl[cpl] = (intersecting_cpls, objs)
        corner_cpl.append(cpl[1])

    # Check if processing fov has one single neighbour
    if len(overlapping_dict.keys()) == 2:
        # Correct overlapping objects at the corner
        intersecting_cpls, objs = adjust_overlapping_objects(
            all_obj_segmented[corner_cpl[0]],
            all_obj_segmented[corner_cpl[1]],
            chunk_coords,
            segmentation_output_path,
            min_overlapping_pixels_segmentation,
        )
        all_obj_segmented[corner_cpl[0]] = obj_segmented_img1
        all_obj_segmented[corner_cpl[1]] = obj_segmented_img2
        output_cpl[tuple(corner_cpl)] = (intersecting_cpls, objs)

    # Format the output
    all_cpls = []
    all_intersecting_step = []
    all_obj_step = {}
    for cpl, (intersecting_cpls, objs) in output_cpl.items():
        all_intersecting_step.append(intersecting_cpls)
        all_obj_step.update(objs)
        all_cpls.append(cpl)

    unraveled_intersecting_step = [
        val for big_list in all_intersecting_step for el in big_list for val in el
    ]

    # To use in case to sort by fov for parallelisation
    # all_couples = set([el for cpl in all_cpls for el in cpl])
    # for el in all_couples:
    #   intersecting_step_dict[el] = [val for val in unraveled_intersecting_step if str(el) in val]

    return unraveled_intersecting_step, all_obj_step


def create_label_image(
    experiment_fpath,
    segmentation_output_path,
    dataset_nuclei,
    nuclei_org_tiles,
    nuclei_adjusted_coords,
    metadata_nuclei,
    client,
    min_overlapping_pixels_segmentation=20,
):

    # Get objects properties and register to the data
    nuclei_org_tiles.tile_corners_coords_pxl = nuclei_adjusted_coords
    reference_corner_fov_position = nuclei_org_tiles.reference_corner_fov_position

    all_futures = []
    for idx, fov_subdataset in dataset_nuclei.dataset.iterrows():
        future = client.submit(
            register_coords_obj,
            fov_subdataset,
            experiment_fpath,
            segmentation_output_path,
            nuclei_adjusted_coords,
            metadata_nuclei,
            reference_corner_fov_position="top-right",
        )
        all_futures.append(future)
    _ = client.gather(all_futures)

    # Remove overlapping objects paralle
    all_futures = []
    for idx, fov_subdataset in dataset_nuclei.dataset.iterrows():
        future = client.submit(
            remove_overlapping_obj,
            fov_subdataset,
            nuclei_org_tiles,
            segmentation_output_path,
            min_overlapping_pixels_segmentation,
        )
        all_futures.append(future)
    all_outputs = client.gather(all_futures)

    # tmp_saving
    pickle.dump(
        all_outputs, open(segmentation_output_path / ("output_removal_tmp.pkl"), "wb")
    )

    all_obj_original = {}
    for idx, fov_subdataset in dataset_nuclei.dataset.iterrows():
        segmented_region_img = load_stitched_segmented_data_full_dict(
            fov_subdataset, segmentation_output_path
        )
        all_obj_original.update(segmented_region_img)

    all_remove_objs_list = []
    all_add_objs = {}
    for rem, obj in all_outputs:
        all_remove_objs_list.append(rem)
        all_add_objs.update(obj)
    all_remove_objs_list = [val for el in all_remove_objs_list for val in el]

    for cell_id, coords in all_add_objs.items():
        all_obj_original[cell_id]["stitched_coords"] = coords

    all_remove_objs_list_tmp = all_remove_objs_list.copy()
    for obj in all_obj_original:
        if obj in all_remove_objs_list_tmp:
            all_remove_objs_list_tmp.remove(obj)

    for obj in all_remove_objs_list_tmp:
        all_obj_original.pop(obj, None)

    pickle.dump(
        all_obj_original,
        open(segmentation_output_path / ("segmented_objects_dict.pkl"), "wb"),
    )

    # Create image
    # Calculate image size

    row_min = nuclei_adjusted_coords[:, 0].min()
    col_min = nuclei_adjusted_coords[:, 1].min()
    row_max = nuclei_adjusted_coords[:, 0].max()
    col_max = nuclei_adjusted_coords[:, 1].max()

    nrows = np.int(row_max - row_min + metadata_nuclei["img_height"])
    ncols = np.int(col_max - col_min + metadata_nuclei["img_width"])

    img = np.zeros([nrows, ncols])

    for idx, data_dict in enumerate(all_obj_original.items()):
        coords = data_dict[1]["stitched_coords"]
        coords = coords - [row_min, col_min]
        coords = coords.astype(int)
        img[coords[:, 0], coords[:, 1]] = idx

    zarr_fpath = segmentation_output_path / "image_segmented_labels.zarr"
    store = zarr.DirectoryStore(zarr_fpath, "w")
    grp = zarr.group(store=store, overwrite=True)
    grp.create_dataset(name="segmented_labels_image", data=img)
