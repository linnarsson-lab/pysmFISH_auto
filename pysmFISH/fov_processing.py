"""
Utility functions and computational graphs used to run processing for a single fov
"""

import pickle
import gc
import dask
import sys
import pandas as pd
import numpy as np
import zarr
from typing import *
from pathlib import Path
from dask import delayed
from dask import dataframe as dd
from dask.base import tokenize
from skimage import img_as_uint


import pysmFISH
from pysmFISH import dots_calling
from pysmFISH import io
from pysmFISH import fovs_registration
from pysmFISH import barcodes_analysis
from pysmFISH import stitching
from pysmFISH import preprocessing
from pysmFISH import configuration_files
from pysmFISH import utils
from pysmFISH import microscopy_file_parsers
from pysmFISH import data_models
from pysmFISH import segmentation_NN
from pysmFISH.logger_utils import selected_logger


def combine_steps(*args):
    pass


def single_fov_round_processing_eel(
    fov_subdataset: pd.Series,
    analysis_parameters: dict,
    running_functions: dict,
    dark_img: np.ndarray,
    experiment_fpath: str,
    preprocessed_zarr_fpath: str,
    save_steps_output: bool = False,
    start_from_preprocessed_imgs: bool = False,
) -> Tuple[pd.DataFrame, Tuple]:

    """Function to run eel processing and counting on a single fov

    Args:
        fov_subdataset (pd.Series): Contains all the metadata relative to a fov.
        analysis_parameters (dict): Processing parameters.
        running_functions (dict): Preprocessing and counting function to run on
            the selected fov.
        dark_img (np.ndarray): Image for the correction of dark signal of the
            camera.
        experiment_fpath (str): Path to the experiment to process
        preprocessed_zarr_fpath (str): Path to the zarr container where the preprocessed
            images will be saved
        save_steps_output (bool, optional): Determine if to save the intermediate
            processing steps. Defaults to False.

    Returns:
        Tuple[pd.DataFrame,Tuple]: counts, filt_out
    """

    logger = selected_logger()
    experiment_fpath = Path(experiment_fpath)

    experiment_name = fov_subdataset.experiment_name
    pipeline = fov_subdataset.pipeline
    processing_type = fov_subdataset.processing_type
    zarr_grp_name = fov_subdataset.grp_name

    raw_data_location = Path(fov_subdataset.raw_data_location)
    parsed_raw_data_fpath = raw_data_location.parent

    if processing_type == "fish":
        processing_parameters = analysis_parameters["fish"]
        filtering_fun = running_functions["fish_channels_preprocessing"]
        counting_fun = running_functions["fish_channels_dots_calling"]

    elif "beads" in processing_type:
        processing_parameters = analysis_parameters[processing_type]
        filtering_fun = running_functions["reference_channels_preprocessing"]
        counting_fun = running_functions["reference_channels_dots_calling"]

    elif "nuclei" in processing_type:
        processing_parameters = analysis_parameters[processing_type]
        filtering_fun = running_functions["reference_channels_preprocessing"]

    elif processing_type != "staining":
        processing_parameters = analysis_parameters[processing_type]
        filtering_fun = running_functions["reference_channels_preprocessing"]

    if start_from_preprocessed_imgs:
        # Load already filtered data
        filt_out = io.load_general_zarr(
            fov_subdataset, preprocessed_zarr_fpath, tag="preprocessed_data"
        )
        filt_out = ((filt_out[0],), filt_out[1])

    else:

        filt_out = getattr(pysmFISH.preprocessing, filtering_fun)(
            zarr_grp_name, parsed_raw_data_fpath, processing_parameters, dark_img
        )

    counts = getattr(pysmFISH.dots_calling, counting_fun)(
        filt_out[0][0], fov_subdataset, processing_parameters
    )

    if save_steps_output:

        # Save the file as zarr
        store = zarr.DirectoryStore(preprocessed_zarr_fpath)
        root = zarr.group(store=store, overwrite=False)
        tag_name = (
            experiment_name
            + "_"
            + fov_subdataset.channel
            + "_round_"
            + str(fov_subdataset.round_num)
            + "_fov_"
            + str(fov_subdataset.fov_num)
        )
        dgrp = root.create_group(tag_name, overwrite=True)
        for k, v in filt_out[1].items():
            dgrp.attrs[k] = v
        fov_name = "preprocessed_data_fov_" + str(fov_subdataset.fov_num)
        dgrp.attrs["fov_name"] = fov_name
        img = utils.convert_to_uint16(
            filt_out[0][-1]
        )  # Must change to zero to save final processed image
        dset = dgrp.create_dataset(
            fov_name, data=img, shape=img.shape, chunks=None, overwrite=True
        )

        # counts.to_parquet(raw_counts_path / (fname + '.parquet'),index=False)

    # return counts, (fov_subdataset.channel,fov_subdataset.round_num,img)
    return counts, filt_out


def single_fov_round_processing_serial_nuclei(
    fov_subdataset: pd.Series,
    analysis_parameters: dict,
    running_functions: dict,
    dark_img: np.ndarray,
    experiment_fpath: str,
    preprocessed_zarr_fpath: str,
    save_steps_output=False,
    start_from_preprocessed_imgs=False,
) -> Tuple[np.ndarray, pd.Series]:

    """Function to run serial processing of nuclei
    Some of the input variable are not used but I wanted to keep the same type of input
    for all the single fov processing functions.

    Args:
        fov_subdataset (pd.Series): Contains all the metadata relative to a fov.
        analysis_parameters (dict): Processing parameters.
        running_functions (dict): Preprocessing and counting function to run on
            the selected fov.
        dark_img (np.ndarray): Image for the correction of dark signal of the
            camera.
        experiment_fpath (str): Path to the experiment to process
        preprocessed_zarr_fpath (str): Path to the zarr container where the preprocessed
            images will be saved
        save_steps_output (bool, optional): Determine if to save the intermediate
            processing steps. Defaults to False.

    Returns:
        Tuple[np.ndarray,pd.Series]: (img,fov_subdataset)
    """
    logger = selected_logger()
    experiment_fpath = Path(experiment_fpath)
    experiment_name = fov_subdataset.experiment_name

    if start_from_preprocessed_imgs:

        store = zarr.DirectoryStore(preprocessed_zarr_fpath)
        root = zarr.group(store=store, overwrite=False)
        tag_name = (
            experiment_name
            + "_"
            + fov_subdataset.channel
            + "_round_"
            + str(fov_subdataset.round_num)
            + "_fov_"
            + str(fov_subdataset.fov_num)
        )

        # Load already filtered data
        filt_out = io.load_general_zarr(
            fov_subdataset, preprocessed_zarr_fpath, tag="preprocessed_data"
        )
        filt_out = ((filt_out[0],), filt_out[1])

        return (filt_out[0][0], fov_subdataset)

    else:

        # Path of directory where to save the intermediate results
        filtered_img_path = experiment_fpath / "results"
        raw_counts_path = experiment_fpath / "results"

        pipeline = fov_subdataset.pipeline
        processing_type = fov_subdataset.processing_type
        zarr_grp_name = fov_subdataset.grp_name

        raw_data_location = Path(fov_subdataset.raw_data_location)
        parsed_raw_data_fpath = raw_data_location.parent

        processing_parameters = analysis_parameters[processing_type]

        filt_out = getattr(
            pysmFISH.preprocessing,
            running_functions["reference_channels_preprocessing"],
        )(zarr_grp_name, parsed_raw_data_fpath, processing_parameters)

        if save_steps_output:

            # Save the file as zarr
            store = zarr.DirectoryStore(preprocessed_zarr_fpath)
            root = zarr.group(store=store, overwrite=False)
            tag_name = (
                experiment_name
                + "_"
                + fov_subdataset.channel
                + "_round_"
                + str(fov_subdataset.round_num)
                + "_fov_"
                + str(fov_subdataset.fov_num)
            )
            dgrp = root.create_group(tag_name, overwrite=True)
            for k, v in filt_out[1].items():
                dgrp.attrs[k] = v
            fov_name = "preprocessed_data_fov_" + str(fov_subdataset.fov_num)
            dgrp.attrs["fov_name"] = fov_name
            img = utils.convert_to_uint16(filt_out[0][-1])
            dset = dgrp.create_dataset(
                fov_name, data=img, shape=img.shape, chunks=None, overwrite=True
            )

        return (img, fov_subdataset)


def processing_barcoded_eel_fov_graph(
    experiment_fpath: str,
    analysis_parameters: dict,
    running_functions: dict,
    tiles_org,
    metadata: dict,
    grpd_fovs: pd.DataFrame,
    save_intermediate_steps: bool,
    preprocessed_image_tag: str,
    client,
    chunks_size: int,
    save_bits_int: int,
    start_from_preprocessed_imgs: False,
):

    """Processing graph for eel type of experiments. Run all the
    steps that can be applied to a single FOV.
    1) Filtering and counting
    2) Register all imaging rounds
    3) Identification of the barcodes (decoding)
    4) Stitching using the stage coords
    5) Generate an output file with the counts for visualisation

    IMPORTANT:
    Because some of the processing steps take quite a bit of time it is necessary
    to process the FOV in chunks to avoid that the processes will fail (workers get
    lost and not comunicate with the scheduler).

    Args:
        experiment_fpath (str): Path to the experiment to process
        analysis_parameters (dict): Parameters to use for the processing
        running_functions (dict): Function to run for preprocessing and counting
        tiles_org (tile_org): Organization of the tiles in the large image
        metadata (dict): Metadata describing the experiment
        grpd_fovs (pd.DataFrame): Database of the experiment grouped by fov
        save_intermediate_steps (bool): Save preprocessed images and raw counts
        preprocessed_image_tag (str): Tag to label the preprocessed images zarr container
        client (distributed.Client): Dask Client that take care of the graph processing
        chunks_size (int): Number of FOVs to process in one go
        save_bits_int (int): Save the intensity of the barcodes (also negative barcodes)
                        and the position of the bits that are flipped
        start_from_preprocessed_imgs (bool): Run the processing starting from the counting
                using preprocessed images. default: False
    """

    experiment_fpath = Path(experiment_fpath)
    io.create_empty_zarr_file(experiment_fpath.as_posix(), preprocessed_image_tag)
    preprocessed_zarr_fpath = experiment_fpath / (
        experiment_fpath.stem + "_" + preprocessed_image_tag + ".zarr"
    )

    microscopy_file_parsers.create_dark_img(experiment_fpath, metadata)

    dark_img = preprocessing.load_dark_image(experiment_fpath)

    # did this conversion to avoid to pass self to dask
    # analysis_parameters = analysis_parameters
    # running_functions = running_functions
    # tile_corners_coords_pxl = tile_corners_coords_pxl

    dark_img = delayed(dark_img)
    analysis_parameters = delayed(analysis_parameters)
    running_functions = delayed(running_functions)
    tile_corners_coords_pxl = delayed(tiles_org.tile_corners_coords_pxl)

    list_all_channels = metadata["list_all_channels"]
    stitching_channel = metadata["stitching_channel"]
    fish_channels = set(list_all_channels) ^ set([stitching_channel])

    codebook_dict = configuration_files.load_codebook(experiment_fpath, metadata)
    codebook_dict = delayed(codebook_dict)

    all_processing = []
    all_filtered_images = []
    all_fovs = list(grpd_fovs.groups.keys())
    chunks = [
        all_fovs[x : x + chunks_size] for x in range(0, len(all_fovs), chunks_size)
    ]
    for chunk in chunks:
        all_processing = []
        all_counts_fov = {}
        for fov_num in chunk:
            all_filtered_images = {}
            all_counts_fov = {}
            all_counts_fov_concat = {}

            fov_group = grpd_fovs.get_group(fov_num)
            channel_grpd = fov_group.groupby("channel")

            all_filtered_imges = {}
            for channel_proc in list_all_channels:
                all_filtered_images[channel_proc] = {}
                all_counts_fov[channel_proc] = []
                group = channel_grpd.get_group(channel_proc)
                for index_value, fov_subdataset in group.iterrows():
                    round_num = fov_subdataset.round_num
                    channel = fov_subdataset.channel

                    fov = fov_subdataset.fov_num
                    experiment_name = fov_subdataset.experiment_name
                    dask_delayed_name = (
                        "filt_count_"
                        + experiment_name
                        + "_"
                        + channel
                        + "_round_"
                        + str(round_num)
                        + "_fov_"
                        + str(fov)
                        + "-"
                        + tokenize()
                    )
                    fov_out = delayed(
                        single_fov_round_processing_eel, name=dask_delayed_name, nout=2
                    )(
                        fov_subdataset,
                        analysis_parameters,
                        running_functions,
                        dark_img,
                        experiment_fpath,
                        preprocessed_zarr_fpath,
                        save_steps_output=save_intermediate_steps,
                        start_from_preprocessed_imgs=start_from_preprocessed_imgs,
                        dask_key_name=dask_delayed_name,
                    )
                    counts, filt_out = fov_out[0], fov_out[1]

                    all_counts_fov[channel_proc].append(counts)

                    if save_bits_int:
                        if channel_proc != fov_subdataset.stitching_channel:
                            all_filtered_images[channel_proc][
                                round_num
                            ] = filt_out  # store it if it gets too big

                name = (
                    "concat_"
                    + experiment_name
                    + "_"
                    + channel_proc
                    + "_"
                    + "_fov_"
                    + str(fov)
                    + "-"
                    + tokenize()
                )
                all_counts_fov_concat[channel_proc] = delayed(pd.concat, name=name)(
                    all_counts_fov[channel_proc], axis=0, ignore_index=True
                )

            if save_intermediate_steps:

                for channel in list_all_channels:
                    name = (
                        "save_raw_counts_"
                        + experiment_name
                        + "_"
                        + channel
                        + "_"
                        + "_fov_"
                        + str(fov)
                        + "-"
                        + tokenize()
                    )

                    saved_raw_counts = delayed(
                        all_counts_fov_concat[channel].to_parquet, name=name
                    )(
                        Path(experiment_fpath)
                        / "results"
                        / (
                            experiment_name
                            + "_raw_counts_channel_"
                            + channel
                            + "_fov_"
                            + str(fov)
                            + ".parquet"
                        ),
                        index=False,
                    )

                    all_processing.append(saved_raw_counts)

            # name = 'register_' +experiment_name + '_' + stitching_channel + '_' \
            #                     + '_fov_' +str(fov) + '-' + tokenize()
            # registered_counts = delayed(fovs_registration.beads_based_registration,name=name)(all_counts_fov_concat[stitching_channel],
            #                                     analysis_parameters)

            registration_stitching_channel_output = delayed(
                fovs_registration.beads_based_registration_stitching_channel, name=name
            )(all_counts_fov_concat[stitching_channel], analysis_parameters, metadata)

            stitching_channel_df, all_rounds_shifts, all_rounds_matching_dots = (
                registration_stitching_channel_output[0],
                registration_stitching_channel_output[1],
                registration_stitching_channel_output[2],
            )

            stitched_coords_reference_df = delayed(
                stitching.stitch_using_coords_general, name=name
            )(
                stitching_channel_df,
                tile_corners_coords_pxl,
                tiles_org.reference_corner_fov_position,
                metadata,
                tag="microscope_stitched",
            )
            all_stitched_coords = []
            all_stitched_coords.append(stitched_coords_reference_df)

            for processing_channel in fish_channels:

                # Register fish
                name = (
                    "register_fish_channels_"
                    + experiment_name
                    + "_"
                    + processing_channel
                    + "_"
                    + "_fov_"
                    + str(fov)
                    + "-"
                    + tokenize()
                )

                registered_counts = delayed(
                    fovs_registration.beads_based_registration_fish, name=name
                )(
                    all_counts_fov_concat[processing_channel],
                    all_rounds_shifts,
                    all_rounds_matching_dots,
                    analysis_parameters,
                )

                # Decoded fish
                name = (
                    "decode_"
                    + experiment_name
                    + "_"
                    + processing_channel
                    + "_"
                    + "_fov_"
                    + str(fov)
                    + "-"
                    + tokenize()
                )
                decoded = delayed(
                    barcodes_analysis.extract_barcodes_NN_fast_multicolor, name=name
                )(
                    registered_counts,
                    analysis_parameters,
                    codebook_dict[processing_channel],
                    metadata,
                )

                if save_bits_int:

                    # all_filtered_images = {}
                    # group_bits = channel_grpd.get_group(processing_channel)
                    # for index_value, fov_subdataset in group_bits.iterrows():
                    #     round_num = fov_subdataset.round_num
                    #     name = 'load_filtered_image_' +experiment_name + '_' \
                    #             + '_fov_' +str(fov_num) + '-' + tokenize()
                    #     filt_out = delayed(io.load_general_zarr,name=name)(fov_subdataset,preprocessed_zarr_fpath,tag='preprocessed_data')
                    # all_filtered_images[round_num] = filt_out

                    name = (
                        "combine_shifted_images_"
                        + experiment_name
                        + "_"
                        + "_fov_"
                        + str(fov_num)
                        + "-"
                        + tokenize()
                    )

                    combined_shift_images = delayed(
                        fovs_registration.combine_register_filtered_image_single_channel,
                        name=name,
                    )(
                        all_filtered_images[processing_channel],
                        metadata,
                        all_rounds_shifts,
                    )

                    name = (
                        "extract_dots_intensities_"
                        + experiment_name
                        + "_"
                        + "_fov_"
                        + str(fov_num)
                        + "-"
                        + tokenize()
                    )
                    extracted_intensities = delayed(
                        barcodes_analysis.extract_dots_images, name=name
                    )(decoded[1], combined_shift_images, experiment_fpath, metadata)

                    # Stitch to the microscope reference coords
                    name = (
                        "stitch_to_mic_coords_"
                        + experiment_name
                        + "_"
                        + processing_channel
                        + "_"
                        + "_fov_"
                        + str(fov)
                        + "-"
                        + tokenize()
                    )
                    stitched_coords = delayed(
                        stitching.stitch_using_coords_general, name=name
                    )(
                        extracted_intensities,
                        tile_corners_coords_pxl,
                        tiles_org.reference_corner_fov_position,
                        metadata,
                        tag="microscope_stitched",
                    )

                else:
                    stitched_coords = delayed(
                        stitching.stitch_using_coords_general, name=name
                    )(
                        decoded[1],
                        tile_corners_coords_pxl,
                        tiles_org.reference_corner_fov_position,
                        metadata,
                        tag="microscope_stitched",
                    )

                all_stitched_coords.append(stitched_coords)

            name = "concat_" + experiment_name + "_fov_" + str(fov) + "-" + tokenize()
            all_stitched_coords = delayed(pd.concat, name=name)(
                all_stitched_coords, axis=0, ignore_index=True
            )

            name = (
                "save_df_"
                + experiment_name
                + "_"
                + channel
                + "_"
                + "_fov_"
                + str(fov)
                + "-"
                + tokenize()
            )
            saved_file = delayed(all_stitched_coords.to_parquet, name=name)(
                Path(experiment_fpath)
                / "results"
                / (experiment_name + "_decoded_fov_" + str(fov) + ".parquet"),
                index=False,
            )

            all_processing.append(saved_file)

        _ = dask.compute(*all_processing)
        client.run(gc.collect)

    io.consolidate_zarr_metadata(preprocessed_zarr_fpath)


# TODO: Needs to be adjusted
def processing_barcoded_eel_fov_graph_from_decoding(
    experiment_fpath: str,
    analysis_parameters: dict,
    tiles_org,
    metadata: dict,
    grpd_fovs: pd.DataFrame,
    client,
    chunks_size: int,
):

    """Processing graph for eel type of experiments. Runs the decoding of
    a preprocessed experiment

    3) Identification of the barcodes (decoding)
    4) Stitching using the stage coords
    5) Generate an output file with the counts for visualisation

    IMPORTANT:
    Because some of the processing steps take quite a bit of time it is necessary
    to process the FOV in chunks to avoid that the processes will fail (workers get
    lost and not comunicate with the scheduler).

    Args:
        experiment_fpath (str): Path to the experiment to process
        analysis_parameters (dict): Parameters to use for the processing
        tiles_org (tile_org): Organization of the tiles in the large image
        metadata (dict): Metadata describing the experiment
        grpd_fovs (pd.DataFrame): Database of the experiment grouped by fov
        client (distributed.Client): Dask Client that take care of the graph processing
        chunks_size (int): Number of FOVs to process in one go

    """

    experiment_fpath = Path(experiment_fpath)

    analysis_parameters = delayed(analysis_parameters)
    tile_corners_coords_pxl = delayed(tiles_org.tile_corners_coords_pxl)

    list_all_channels = metadata["list_all_channels"]
    stitching_channel = metadata["stitching_channel"]
    fish_channels = set(list_all_channels) ^ set([stitching_channel])

    codebook_dict = configuration_files.load_codebook(experiment_fpath, metadata)
    codebook_dict = delayed(codebook_dict)

    all_processing = []
    all_fovs = list(grpd_fovs.groups.keys())
    chunks = [
        all_fovs[x : x + chunks_size] for x in range(0, len(all_fovs), chunks_size)
    ]
    for chunk in chunks:
        all_processing = []
        all_counts_fov = {}
        for fov_num in chunk:

            counts_fpath = list(
                (experiment_fpath / "results").glob(
                    "*_decoded_fov_" + str(fov_num) + ".parquet"
                )
            )[0]
            experiment_name = experiment_fpath.stem

            name = (
                "load_counts_"
                + experiment_name
                + "_"
                + "_fov_"
                + str(fov_num)
                + "-"
                + tokenize()
            )
            counts_fov = delayed(pd.read_parquet, name=name)(counts_fpath)

            # all_stitched_coords = []
            # for processing_channel in fish_channels:

            #     # Decoded fish
            #     name = 'decode_' +experiment_name + '_' + processing_channel + '_' \
            #                         + '_fov_' +str(fov_num) + '-' + tokenize()
            #     decoded = delayed(barcodes_analysis.extract_barcodes_NN_fast_multicolor,name=name)(counts_fov,
            #                                                             analysis_parameters,codebook_dict[processing_channel],
            #                                                             metadata)

            #     # Stitch to the microscope reference coords
            #     name = 'stitch_to_mic_coords_' +experiment_name + '_' + processing_channel + '_' \
            #                         + '_fov_' +str(fov_num) + '-' + tokenize()
            #     stitched_coords = delayed(stitching.stitch_using_coords_general,name=name)(decoded[1],
            #                                                     tile_corners_coords_pxl,tiles_org.reference_corner_fov_position,
            #                                                     metadata,tag='microscope_stitched')

            #     all_stitched_coords.append(stitched_coords)

            # name = 'concat_' +experiment_name + \
            #                         '_fov_' + str(fov_num) + '-' + tokenize()
            # all_stitched_coords = delayed(pd.concat,name=name)(all_stitched_coords,axis=0,ignore_index=True)

            # name = 'save_df_' +experiment_name + '_' \
            #                     + '_fov_' +str(fov_num) + '-' + tokenize()
            # saved_file = delayed(all_stitched_coords.to_parquet,name=name)(Path(experiment_fpath) / 'results'/ (experiment_name + \
            #                 '_decoded_fov_' + str(fov_num) + '.parquet'),index=False)

            # all_processing.append(saved_file)
            all_processing.append(counts_fov)

        _ = dask.compute(*all_processing)


def processing_barcoded_eel_fov_starting_from_registration_graph(
    experiment_fpath: str,
    analysis_parameters: dict,
    running_functions: dict,
    tiles_org,
    metadata: dict,
    grpd_fovs: pd.DataFrame,
    preprocessed_image_tag: str,
    client,
    chunk_size: int,
    save_bits_int: int,
):
    """Processing graph for runnning analysis of eel type experiments
    skipping the preprocessing and counting. It is useful when there
    are issue with the registration of the different rounds. The
    processing restart with the registration of the differen rounds.
    1) Register all imaging rounds
    2) Identification of the barcodes (decoding)
    3) Stitching using the stage coords
    4) Generate an output file with the counts for visualisation

    IMPORTANT:
    Because some of the processing steps take quite a bit of time it is necessary
    to process the FOV in chunks to avoid that the processes will fail (workers get
    lost and not comunicate with the scheduler).

    Args:
        experiment_fpath (str): Path to the experiment to process
        analysis_parameters (dict): Parameters to use for the processing
        running_functions (dict): Function to run for preprocessing and counting
        tiles_org (tile_org): Organization of the tiles in the large image
        metadata (dict): Metadata describing the experiment
        grpd_fovs (pd.DataFrame): Database of the experiment grouped by fov
        preprocessed_image_tag (str): Tag to label the preprocessed images zarr container
        client (distributed.Client): Dask Client that take care of the graph processing
        chunks_size (int): Number of FOVs to process in one go
        save_bits_int (int): Save the intensity of the barcodes (also negative barcodes)
                        and the position of the bits that are flipped
    """

    experiment_fpath = Path(experiment_fpath)
    experiment_name = experiment_fpath.stem
    preprocessed_zarr_fpath = experiment_fpath / (
        experiment_fpath.stem + "_" + preprocessed_image_tag + ".zarr"
    )

    list_all_channels = metadata["list_all_channels"]
    stitching_channel = metadata["stitching_channel"]
    fish_channels = set(list_all_channels) ^ set([stitching_channel])
    total_rounds = metadata["total_rounds"]
    all_processing = []

    analysis_parameters = delayed(analysis_parameters)
    tile_corners_coords_pxl = delayed(tiles_org.tile_corners_coords_pxl)
    codebook_dict = configuration_files.load_codebook(experiment_fpath, metadata)
    codebook_dict = delayed(codebook_dict)

    all_fovs = list(grpd_fovs.groups.keys())
    chunks = [all_fovs[x : x + chunk_size] for x in range(0, len(all_fovs), chunk_size)]

    for chunk in chunks:
        all_processing = []
        for fov_num in chunk:

            fov_group = grpd_fovs.get_group(fov_num)
            channel_grpd = fov_group.groupby("channel")

            stitching_channel_fpath = list(
                (experiment_fpath / "results").glob(
                    "*_raw_counts_channel_"
                    + stitching_channel
                    + "_fov_"
                    + str(fov_num)
                    + ".parquet"
                )
            )[0]

            name = (
                "load_counts_"
                + experiment_name
                + "_"
                + "_fov_"
                + str(fov_num)
                + "-"
                + tokenize()
            )
            all_counts_fov = delayed(pd.read_parquet, name=name)(
                stitching_channel_fpath
            )

            name = (
                "regitration_stitching_"
                + experiment_name
                + "_"
                + "_fov_"
                + str(fov_num)
                + "-"
                + tokenize()
            )
            registration_stitching_channel_output = delayed(
                fovs_registration.beads_based_registration_stitching_channel, name=name
            )(all_counts_fov, analysis_parameters, metadata)

            stitching_channel_df, all_rounds_shifts, all_rounds_matching_dots = (
                registration_stitching_channel_output[0],
                registration_stitching_channel_output[1],
                registration_stitching_channel_output[2],
            )

            name = (
                "stitching_to_microscope_"
                + experiment_name
                + "_"
                + "_fov_"
                + str(fov_num)
                + "-"
                + tokenize()
            )
            stitched_coords_reference_df = delayed(
                stitching.stitch_using_coords_general, name=name
            )(
                stitching_channel_df,
                tile_corners_coords_pxl,
                tiles_org.reference_corner_fov_position,
                metadata,
                tag="microscope_stitched",
            )
            all_stitched_coords = []
            all_stitched_coords.append(stitched_coords_reference_df)

            for processing_channel in fish_channels:

                # Register fish
                name = (
                    "load_fish_channels_"
                    + experiment_name
                    + "_"
                    + processing_channel
                    + "_"
                    + "_fov_"
                    + str(fov_num)
                    + "-"
                    + tokenize()
                )

                channel_fpath = list(
                    (experiment_fpath / "results").glob(
                        "*_raw_counts_channel_"
                        + processing_channel
                        + "_fov_"
                        + str(fov_num)
                        + ".parquet"
                    )
                )[0]

                name = (
                    "load_counts_"
                    + experiment_name
                    + "_"
                    + "_fov_"
                    + str(fov_num)
                    + "-"
                    + tokenize()
                )
                fish_counts_fov = delayed(pd.read_parquet, name=name)(channel_fpath)

                registered_counts = delayed(
                    fovs_registration.beads_based_registration_fish, name=name
                )(
                    fish_counts_fov,
                    all_rounds_shifts,
                    all_rounds_matching_dots,
                    analysis_parameters,
                )

                # Decoded fish
                name = (
                    "decode_"
                    + experiment_name
                    + "_"
                    + processing_channel
                    + "_"
                    + "_fov_"
                    + str(fov_num)
                    + "-"
                    + tokenize()
                )
                decoded = delayed(
                    barcodes_analysis.extract_barcodes_NN_fast_multicolor, name=name
                )(
                    registered_counts,
                    analysis_parameters,
                    codebook_dict[processing_channel],
                    metadata,
                )

                if save_bits_int:

                    # all_filtered_images = {}
                    # group_bits = channel_grpd.get_group(processing_channel)
                    # for index_value, fov_subdataset in group_bits.iterrows():
                    #     round_num = fov_subdataset.round_num
                    #     name = 'load_filtered_image_' +experiment_name + '_' \
                    #             + '_fov_' +str(fov_num) + '-' + tokenize()
                    #     filt_out = delayed(io.load_general_zarr,name=name)(fov_subdataset,preprocessed_zarr_fpath,tag='preprocessed_data')
                    # all_filtered_images[round_num] = filt_out

                    name = (
                        "combine_shifted_images_"
                        + experiment_name
                        + "_"
                        + "_fov_"
                        + str(fov_num)
                        + "-"
                        + tokenize()
                    )

                    combined_shift_images = delayed(
                        fovs_registration.combine_register_filtered_image_single_channel,
                        name=name,
                    )(
                        all_filtered_images[processing_channel],
                        metadata,
                        all_rounds_shifts,
                    )

                    name = (
                        "extract_dots_intensities_"
                        + experiment_name
                        + "_"
                        + "_fov_"
                        + str(fov_num)
                        + "-"
                        + tokenize()
                    )
                    extracted_intensities = delayed(
                        barcodes_analysis.extract_dots_images, name=name
                    )(decoded[1], combined_shift_images, experiment_fpath, metadata)

                    # Stitch to the microscope reference coords
                    name = (
                        "stitch_to_mic_coords_"
                        + experiment_name
                        + "_"
                        + processing_channel
                        + "_"
                        + "_fov_"
                        + str(fov_num)
                        + "-"
                        + tokenize()
                    )
                    stitched_coords = delayed(
                        stitching.stitch_using_coords_general, name=name
                    )(
                        extracted_intensities,
                        tile_corners_coords_pxl,
                        tiles_org.reference_corner_fov_position,
                        metadata,
                        tag="microscope_stitched",
                    )

                else:
                    stitched_coords = delayed(
                        stitching.stitch_using_coords_general, name=name
                    )(
                        decoded[1],
                        tile_corners_coords_pxl,
                        tiles_org.reference_corner_fov_position,
                        metadata,
                        tag="microscope_stitched",
                    )

                all_stitched_coords.append(stitched_coords)

            name = (
                "concat_" + experiment_name + "_fov_" + str(fov_num) + "-" + tokenize()
            )
            all_stitched_coords = delayed(pd.concat, name=name)(
                all_stitched_coords, axis=0, ignore_index=True
            )

            name = (
                "save_df_"
                + experiment_name
                + "_"
                + "_fov_"
                + str(fov_num)
                + "-"
                + tokenize()
            )
            saved_file = delayed(all_stitched_coords.to_parquet, name=name)(
                Path(experiment_fpath)
                / "results"
                / (experiment_name + "_decoded_fov_" + str(fov_num) + ".parquet"),
                index=False,
            )

            all_processing.append(saved_file)

        _ = dask.compute(*all_processing)


def combine_filtered_images(
    output_list: list, experiment_fpath: str, metadata: pd.DataFrame, save: bool = False
):
    """Function used to combine all the filtered images for a fov/channel in a single
        image stack
    Args:
        output_list (list): list containing the output of preprocessing
        experiment_fpath (str): path to the experiment to process
        metadata (pd.DataFrame): dataframe containing the metadata
        save (bool, optional): Determine if the filtered images should be stored Defaults to False.
    Returns:
        img_stack (np.ndarray): image stack of all the images for a fov. The position in the
                stack correspond to round_num-1
    """
    experiment_fpath = Path(experiment_fpath)

    img_stack = np.zeros(
        [metadata["total_rounds"], metadata["img_width"], metadata["img_height"]]
    )

    for img, img_meta in output_list:
        round_num = img_meta.round_num
        img_stack[round_num - 1, :, :] = img

    if save:
        # Add conversion to more compress ftype
        img_meta = output_list[0][1]
        channel = img_meta.channel
        fov = img_meta.fov_num
        fpath = (
            experiment_fpath
            / "results"
            / (
                experiment_fpath.stem
                + "_"
                + channel
                + "_combined_img_fov_"
                + fov
                + ".npy"
            )
        )
        np.save(fpath, img_stack)

    return img_stack


def processing_serial_fish_fov_graph(
    experiment_fpath: str,
    analysis_parameters: dict,
    running_functions: dict,
    tiles_org,
    metadata: dict,
    grpd_fovs: pd.DataFrame,
    save_intermediate_steps: bool,
    preprocessed_image_tag: str,
    client,
    chunks_size: int,
    start_from_preprocessed_imgs=False,
):
    """Processing graph for serial type of experiments.

    Args:

        experiment_fpath (str): Path to the experiment to process
        analysis_parameters (dict): Parameters to use for the processing
        running_functions (dict): Function to run for preprocessing and counting
        tiles_org (tile_org): Organization of the tiles in the large image
        metadata (dict): Metadata describing the experiment
        grpd_fovs (pd.DataFrame): Database of the experiment grouped by fov
        save_intermediate_steps (bool): Save preprocessed images and raw counts
        preprocessed_image_tag (str): Tag to label the preprocessed images zarr container
        client (distributed.Client): Dask Client that take care of the graph processing
        chunks_size (int): Number of FOVs to process in one go

    """

    experiment_fpath = Path(experiment_fpath)
    io.create_empty_zarr_file(experiment_fpath, preprocessed_image_tag)
    preprocessed_zarr_fpath = experiment_fpath / (
        experiment_fpath.stem + "_" + preprocessed_image_tag + ".zarr"
    )

    microscopy_file_parsers.create_dark_img(experiment_fpath, metadata)

    dark_img = preprocessing.load_dark_image(experiment_fpath)

    # did this conversion to avoid to pass self to dask
    analysis_parameters = analysis_parameters
    running_functions = running_functions
    tile_corners_coords_pxl = tiles_org.tile_corners_coords_pxl

    dark_img = delayed(dark_img)
    analysis_parameters = delayed(analysis_parameters)
    running_functions = delayed(running_functions)
    tile_corners_coords_pxl = delayed(tile_corners_coords_pxl)

    all_processing = []
    all_filtered_images = []
    all_fovs = list(grpd_fovs.groups.keys())
    chunks = [
        all_fovs[x : x + chunks_size] for x in range(0, len(all_fovs), chunks_size)
    ]
    for chunk in chunks:
        all_processing = []
        all_filtered_images = []
        for fov_num in chunk:
            group = grpd_fovs.get_group(fov_num)
            # for fov_num, group in grpd_fovs:
            all_counts_fov = []
            all_nuclei_fov = []
            for index_value, fov_subdataset in group.iterrows():
                round_num = fov_subdataset.round_num
                channel = fov_subdataset.channel
                fov = fov_subdataset.fov_num
                stitching_type = fov_subdataset.stitching_type
                experiment_name = fov_subdataset.experiment_name
                processing_type = fov_subdataset.processing_type

                if processing_type == "nuclei":
                    dask_delayed_name = (
                        "filt_"
                        + experiment_name
                        + "_"
                        + channel
                        + "_round_"
                        + str(round_num)
                        + "_fov_"
                        + str(fov)
                        + "-"
                        + tokenize()
                    )

                    out_nuclei = delayed(
                        single_fov_round_processing_serial_nuclei,
                        name=dask_delayed_name,
                    )(
                        fov_subdataset,
                        analysis_parameters,
                        running_functions,
                        dark_img,
                        experiment_fpath,
                        preprocessed_zarr_fpath,
                        save_steps_output=save_intermediate_steps,
                        start_from_preprocessed_imgs=start_from_preprocessed_imgs,
                    )
                    all_nuclei_fov.append(out_nuclei)

                else:
                    dask_delayed_name = (
                        "filt_count_"
                        + experiment_name
                        + "_"
                        + channel
                        + "_round_"
                        + str(round_num)
                        + "_fov_"
                        + str(fov)
                        + "-"
                        + tokenize()
                    )
                    fov_out = delayed(
                        single_fov_round_processing_eel, name=dask_delayed_name
                    )(
                        fov_subdataset,
                        analysis_parameters,
                        running_functions,
                        dark_img,
                        experiment_fpath,
                        preprocessed_zarr_fpath,
                        save_steps_output=save_intermediate_steps,
                        start_from_preprocessed_imgs=start_from_preprocessed_imgs,
                    )

                    counts, filt_out = fov_out[0], fov_out[1]
                    all_counts_fov.append(counts)

                    # if channel != fov_subdataset.stitching_channel:
                    #     all_filtered_images.append(filt_out)

            name = (
                "concat_fish_"
                + experiment_name
                + "_"
                + channel
                + "_"
                + "_fov_"
                + str(fov)
                + "-"
                + tokenize()
            )
            all_counts_fov = delayed(pd.concat, name=name)(
                all_counts_fov, axis=0, ignore_index=True
            )

            if stitching_type == "nuclei":

                name = (
                    "create_nuclei_stack"
                    + experiment_name
                    + "_"
                    + channel
                    + "_"
                    + "_fov_"
                    + str(fov)
                    + "-"
                    + tokenize()
                )
                filtered_nuclei_stack = delayed(combine_filtered_images, name=name)(
                    all_nuclei_fov, experiment_fpath, metadata
                )

                name = (
                    "register_"
                    + experiment_name
                    + "_"
                    + channel
                    + "_"
                    + "_fov_"
                    + str(fov)
                    + "-"
                    + tokenize()
                )
                registered_counts = delayed(
                    fovs_registration.nuclei_based_registration, name=name
                )(all_counts_fov, filtered_nuclei_stack, analysis_parameters)

            else:

                name = (
                    "register_"
                    + experiment_name
                    + "_"
                    + channel
                    + "_"
                    + "_fov_"
                    + str(fov)
                    + "-"
                    + tokenize()
                )
                registered_counts = delayed(
                    fovs_registration.beads_based_registration, name=name
                )(all_counts_fov, analysis_parameters)

            name = (
                "stitch_to_mic_coords_"
                + experiment_name
                + "_"
                + channel
                + "_"
                + "_fov_"
                + str(fov)
                + "-"
                + tokenize()
            )

            stitched_coords = delayed(stitching.stitch_using_coords_general, name=name)(
                registered_counts,
                tile_corners_coords_pxl,
                tiles_org.reference_corner_fov_position,
                metadata,
                tag="microscope_stitched",
            )

            name = (
                "register_and_combine_filt_imgs"
                + experiment_name
                + "_"
                + channel
                + "_"
                + "_fov_"
                + str(fov)
                + "-"
                + tokenize()
            )

            # combined_images = delayed(fovs_registration.combine_register_filtered_images,name=name)(all_filtered_images,stitched_coords,
            #                                                                                 fov_subdataset.stitching_channel)

            name = (
                "save_df_"
                + experiment_name
                + "_"
                + channel
                + "_"
                + "_fov_"
                + str(fov)
                + "-"
                + tokenize()
            )
            saved_file = delayed(stitched_coords.to_parquet, name=name)(
                Path(experiment_fpath)
                / "results"
                / (experiment_name + "_decoded_fov_" + str(fov) + ".parquet"),
                index=False,
            )

            all_processing.append(saved_file)

        _ = dask.compute(*all_processing)

    io.consolidate_zarr_metadata(preprocessed_zarr_fpath)


def single_fov_fresh_tissue_beads(
    processing_tag: str,
    fov_subdataset: pd.Series,
    analysis_parameters: dict,
    running_functions: dict,
    dark_img: np.ndarray,
    experiment_fpath: str,
    preprocessed_zarr_fpath: str,
    save_steps_output: bool = True,
):
    """[summary]

    Args:
        processing_tag (str): name describing the processing
        fov_subdataset (pd.Series): metadata corresponding to the specific fov
        analysis_parameters (dict): paramters used for processing
        running_functions (dict): functions to run preprocessing and detection of the beads
        dark_img (np.ndarray): background of the camera
        experiment_fpath (str): path to the experiment to process.
        preprocessed_zarr_fpath (str): path to the zarr container with the preprocessed images
        save_steps_output (bool, optional): Determine if to save preprocessing data. Defaults to True.

    Returns:
        [Tuple[pd.DataFrame, Tuple]]: counts, filt_out if processing_tag == 'beads'
        [Tuple]: filt_out if processing_tag == 'nuclei'
    """

    logger = selected_logger()
    experiment_fpath = Path(experiment_fpath)

    experiment_name = fov_subdataset.experiment_name
    zarr_grp_name = fov_subdataset.grp_name

    parsed_raw_data_fpath = (
        experiment_fpath / "fresh_tissue" / (experiment_name + "_img_data.zarr")
    )

    processing_parameters = analysis_parameters["fresh-tissue"][processing_tag]

    if processing_tag == "beads":
        filtering_fun = running_functions["fresh_sample_reference_preprocessing"]
        counting_fun = running_functions["fresh_sample_reference_dots_calling"]

        filt_out = getattr(pysmFISH.preprocessing, filtering_fun)(
            zarr_grp_name, parsed_raw_data_fpath, processing_parameters, dark_img
        )

        counts = getattr(pysmFISH.dots_calling, counting_fun)(
            filt_out[0][0], fov_subdataset, processing_parameters
        )

        if save_steps_output:

            # Save the file as zarr
            store = zarr.DirectoryStore(preprocessed_zarr_fpath)
            root = zarr.group(store=store, overwrite=False)
            tag_name = (
                experiment_name
                + "_fresh_tissue_"
                + processing_tag
                + "_fov_"
                + str(fov_subdataset.fov_num)
            )
            dgrp = root.create_group(tag_name, overwrite=True)
            for k, v in filt_out[1].items():
                dgrp.attrs[k] = v
            fov_name = "preprocessed_data_fov_" + str(fov_subdataset.fov_num)
            dgrp.attrs["fov_name"] = fov_name
            img = utils.convert_to_uint16(filt_out[0][-1])
            dset = dgrp.create_dataset(
                fov_name, data=img, shape=img.shape, chunks=None, overwrite=True
            )

    elif processing_tag == "nuclei":
        filtering_fun = running_functions["fresh_sample_nuclei_preprocessing"]
        filt_out = getattr(pysmFISH.preprocessing, filtering_fun)(
            zarr_grp_name, parsed_raw_data_fpath, processing_parameters
        )

        if save_steps_output:

            # Save the file as zarr
            store = zarr.DirectoryStore(preprocessed_zarr_fpath)
            root = zarr.group(store=store, overwrite=False)
            tag_name = (
                experiment_name
                + "_fresh_tissue_"
                + processing_tag
                + "_fov_"
                + str(fov_subdataset.fov_num)
            )
            dgrp = root.create_group(tag_name, overwrite=True)
            for k, v in filt_out[1].items():
                dgrp.attrs[k] = v
            fov_name = "preprocessed_data_fov_" + str(fov_subdataset.fov_num)
            dgrp.attrs["fov_name"] = fov_name
            img = img_as_uint(filt_out[0][-1])
            dset = dgrp.create_dataset(
                fov_name, data=img, shape=img.shape, chunks=None, overwrite=True
            )

    if processing_tag == "beads":
        return counts, filt_out

    elif processing_tag == "nuclei":
        return filt_out


def make_fresh_beads_count_like_eel(data, eel_metadata):
    data["r_px_registered"] = data["r_px_original"]
    data["c_px_registered"] = data["c_px_original"]
    data["hamming_distance"] = 0
    data["decoded_genes"] = "beads"
    data["machine"] = eel_metadata["machine"]
    return data


def segmentation_NN_fov(
    img: np.ndarray,
    fov_subdataset: pd.Series,
    segmented_file_path: str,
    fresh_tissue_segmentation_engine: str,
    diameter_size: int,
    model,
):

    experiment_name = fov_subdataset.experiment_name
    segmented_file_path = Path(segmented_file_path)

    nuclei_segmentation = segmentation_NN.Segmenation_NN(
        fresh_tissue_segmentation_engine, diameter_size, model
    )
    mask = nuclei_segmentation.segment(img)

    if isinstance(mask, np.ndarray):
        # Save the file as zarr
        store = zarr.DirectoryStore(segmented_file_path)
        root = zarr.group(store=store, overwrite=False)
        tag_name = (
            experiment_name
            + "_segmetented_fresh_tissue_fov_"
            + str(fov_subdataset.fov_num)
        )
        dgrp = root.create_group(tag_name, overwrite=True)
        fov_name = "segmentation_mask_fov_" + str(fov_subdataset.fov_num)

        dset = dgrp.create_dataset(
            fov_name,
            data=mask,
            shape=mask.shape,
            chunks=None,
            overwrite=True,
        )


def process_fresh_sample_graph(
    experiment_fpath: str,
    running_functions: dict,
    analysis_parameters: dict,
    client,
    chunks_size: int,
    tag_ref_beads: str,
    tag_nuclei: str,
    eel_metadata: dict,
    fresh_tissue_segmentation_engine: str,
    diameter_size: int,
    parsing: bool = True,
    overlapping_percentage: int = 5,
    save_steps_output: bool = True,
):
    """Processing graph for the low magnification images of the
    tissues nuclei acquired before eel and used for the
    segmentation and identification of the cells

    1) Parsing of the raw images (if required)
    2) Create the fresh nuclei images dataset
    3) Preprocessing and counting (beads) or
       preprocessing (nuclei)

    Args:
        experiment_fpath (str): path to the experiment to process
        running_functions (dict): function used to run preprocessing and counting
        analysis_parameters (dict): parameters used to run the analysis
        client (distributed.Client): dask client coordinating the processing
        chunks_size (int): processing in chunks
        tag_ref_beads (str): str in the files name used to identify the images of the beads
        tag_nuclei (str): str in the files name used to identify the images of the nuclei
        eel_metadata (dict): overall experiment info
        parsing (bool, optional): Determine if the images need to be parsed. Defaults to True.
        save_steps_output (bool): save the processed data
    """

    logger = selected_logger()
    all_parsing = []
    if fresh_tissue_segmentation_engine == "stardist":
        from stardist.models import StarDist2D

        model = StarDist2D.from_pretrained("2D_versatile_fluo")
        # model_path = (
        #     Path(experiment_fpath) / "fresh_tissue" / "stardist_2D_versatile_fluo.zip"
        # )
        # model.export_TF(fname=model_path)
    else:
        # TODO fix the same thing for cellpose
        pass

    if parsing:
        presence_nuclei = 0
        presence_beads = 0
        all_fresh_tissue_fpath = list(
            (Path(experiment_fpath) / "fresh_tissue").glob("*.nd2")
        )
        if all_fresh_tissue_fpath:
            for fpath in all_fresh_tissue_fpath:
                if tag_ref_beads in fpath.stem:
                    parsed_beads_fpath = (
                        experiment_fpath
                        / "fresh_tissue"
                        / (fpath.stem + "_img_data.zarr")
                    )
                    parsing_future = client.submit(
                        microscopy_file_parsers.nikon_nd2_parser_simple_mfov, fpath
                    )
                    all_parsing.append(parsing_future)
                    presence_beads = 1
                elif tag_nuclei in fpath.stem:
                    parsed_nuclei_fpath = (
                        experiment_fpath
                        / "fresh_tissue"
                        / (fpath.stem + "_img_data.zarr")
                    )
                    parsing_future = client.submit(
                        microscopy_file_parsers.nikon_nd2_parser_simple_mfov, fpath
                    )
                    all_parsing.append(parsing_future)
                    presence_nuclei = 1

            if presence_beads and presence_nuclei:
                _ = client.gather(all_parsing)
                io.consolidate_zarr_metadata(parsed_beads_fpath)
                io.consolidate_zarr_metadata(parsed_nuclei_fpath)
            else:
                if not presence_beads:
                    logger.error(f"missing fresh-tissue beads file")
                    sys.exit(f"missing fresh-tissue beads file")
                elif not presence_nuclei:
                    logger.error(f"missing fresh-tissue nuclei file")
                    sys.exit(f"missing fresh-tissue nuclei file")
                else:
                    logger.error(f"missing fresh-tissue beads and nuclei files")
                    sys.exit(f"missing fresh-tissue beads and nuclei files")

    else:
        all_fresh_tissue_fpath = list(
            (Path(experiment_fpath) / "fresh_tissue").glob("*.zarr")
        )
        if all_fresh_tissue_fpath:
            for fpath in all_fresh_tissue_fpath:
                if tag_ref_beads in fpath.stem:
                    parsed_beads_fpath = fpath
                    presence_beads = 1
                elif tag_nuclei in fpath.stem:
                    parsed_nuclei_fpath = fpath
                    presence_nuclei = 1

        if not presence_beads:
            logger.error(f"missing fresh-tissue beads parsed file")
            sys.exit(f"missing fresh-tissue beads parsed file")
        elif not presence_nuclei:
            logger.error(f"missing fresh-tissue nuclei parsed file")
            sys.exit(f"missing fresh-tissue nuclei parsed file")
        elif (not presence_nuclei) & (not presence_beads):
            logger.error(f"missing fresh-tissue beads and nuclei parsed files")
            sys.exit(f"missing fresh-tissue beads and nuclei parsed files")

    # Create dataset

    ds_beads = data_models.Dataset()
    ds_nuclei = data_models.Dataset()
    ds_beads.create_full_dataset_from_zmetadata(parsed_beads_fpath)
    ds_nuclei.create_full_dataset_from_zmetadata(parsed_nuclei_fpath)

    ds_beads.dataset["processing_type"] = "undefined"
    ds_beads.dataset["overlapping_percentage"] = overlapping_percentage
    ds_beads.dataset["machine"] = eel_metadata["machine"]
    ds_beads.dataset["round_num"] = 1
    ds_beads.save_dataset(ds_beads.dataset, ds_beads.dataset_fpath)

    ds_nuclei.dataset["processing_type"] = "undefined"
    ds_nuclei.dataset["overlapping_percentage"] = overlapping_percentage
    ds_nuclei.dataset["machine"] = eel_metadata["machine"]
    ds_nuclei.dataset["round_num"] = 1
    ds_nuclei.save_dataset(ds_nuclei.dataset, ds_nuclei.dataset_fpath)

    beads_grpd_fovs = ds_beads.dataset.groupby("fov_num")
    nuclei_grpd_fovs = ds_nuclei.dataset.groupby("fov_num")

    # In this case I fake a dark image. It must be collected from the
    # robofish system
    img_width = ds_nuclei.dataset.iloc[0].img_width
    img_height = ds_nuclei.dataset.iloc[0].img_height
    dark_img = np.zeros([img_width, img_height])

    base_path = experiment_fpath / "fresh_tissue"
    utils.create_dir(base_path / "results")
    # nuclei_base = parsed_nuclei_fpath.stem.split('_img_data.zarr')[0]
    # beads_base = parsed_beads_fpath.stem.split('_img_data.zarr')[0]
    nuclei_filtered_fpath = base_path / (
        base_path.stem + "_nuclei_preprocessed_img_data.zarr"
    )
    io.create_empty_zarr_file(base_path.as_posix(), tag="nuclei_preprocessed_img_data")
    beads_filtered_fpath = base_path / (
        base_path.stem + "_beads_preprocessed_img_data.zarr"
    )
    io.create_empty_zarr_file(base_path.as_posix(), tag="beads_preprocessed_img_data")

    io.create_empty_zarr_file(base_path.as_posix(), tag="segmented_nuclei_data")
    segmented_file_path = base_path / (base_path.stem + "_segmented_nuclei_data.zarr")

    client.run(gc.collect)

    processing_tag = "beads"

    metadata = ds_beads.collect_metadata(ds_beads.dataset)

    round_num = 1
    stitching_channel = "Europium"
    tiles_org = stitching.organize_square_tiles(
        base_path, ds_beads.dataset, metadata, round_num
    )
    tiles_org.run_tiles_organization()

    all_fovs = list(beads_grpd_fovs.groups.keys())
    chunks = [
        all_fovs[x : x + chunks_size] for x in range(0, len(all_fovs), chunks_size)
    ]
    for chunk in chunks:
        all_processing = []
        all_counts_beads = []
        for fov_num in chunk:
            fov_subdataset = beads_grpd_fovs.get_group(fov_num).iloc[
                0
            ]  # Olny one round of imaging

            round_num = fov_subdataset.round_num
            channel = fov_subdataset.channel
            fov = fov_subdataset.fov_num
            experiment_name = fov_subdataset.experiment_name
            dask_delayed_name = "filt_count_beads_fov" + str(fov) + "_" + tokenize()
            fov_out = delayed(single_fov_fresh_tissue_beads, name=dask_delayed_name)(
                processing_tag,
                fov_subdataset,
                analysis_parameters,
                running_functions,
                dark_img,
                experiment_fpath,
                preprocessed_zarr_fpath=beads_filtered_fpath,
                save_steps_output=save_steps_output,
                dask_key_name=dask_delayed_name,
            )
            counts, filt_out = fov_out[0], fov_out[1]

            # name = 'concat_all_counts_beads_fresh_tissue'+ '-' + tokenize()
            # all_counts_fov = delayed(pd.concat,name=name)(all_counts_beads,axis=0,ignore_index=True)

            name = "add missing fields" + "-" + tokenize()
            counts_adj = delayed(make_fresh_beads_count_like_eel, name=name)(
                counts, eel_metadata
            )

            # Stitch to the microscope reference coords
            name = (
                "stitch_to_mic_coords_"
                + experiment_name
                + "_"
                + "_fov_"
                + str(fov)
                + "-"
                + tokenize()
            )
            stitched_coords = delayed(stitching.stitch_using_coords_general, name=name)(
                counts_adj,
                tiles_org.tile_corners_coords_pxl,
                tiles_org.reference_corner_fov_position,
                metadata,
                tag="microscope_stitched",
            )

            # Add registration and recalculation of all the coords
            name = (
                "save_df_beads_fresh_tissue"
                + experiment_name
                + "_"
                + channel
                + "_"
                + "_fov_"
                + str(fov)
                + "-"
                + tokenize()
            )
            saved_file = delayed(stitched_coords.to_parquet, name=name)(
                base_path
                / "results"
                / (
                    experiment_name
                    + "_counts_beads_fresh_tissue_decoded_fov_"
                    + str(fov_num)
                    + ".parquet"
                ),
                index=False,
            )

            all_processing.append(saved_file)
        _ = dask.compute(all_processing)
        client.run(gc.collect)

    processing_tag = "nuclei"
    metadata = ds_nuclei.collect_metadata(ds_beads.dataset)

    all_fovs = list(nuclei_grpd_fovs.groups.keys())
    chunks = [
        all_fovs[x : x + chunks_size] for x in range(0, len(all_fovs), chunks_size)
    ]
    for chunk in chunks:
        scattered_mdel = client.scatter(model)
        all_processing_nuclei = []
        for fov_num in chunk:
            fov_subdataset = nuclei_grpd_fovs.get_group(fov_num).iloc[
                0
            ]  # Olny one round of imaging
            round_num = fov_subdataset.round_num
            channel = fov_subdataset.channel
            fov = fov_subdataset.fov_num
            experiment_name = fov_subdataset.experiment_name
            dask_delayed_name = "filt_nuclei_fov" + str(fov) + "_" + tokenize()
            fov_out = delayed(single_fov_fresh_tissue_beads, name=dask_delayed_name)(
                processing_tag,
                fov_subdataset,
                analysis_parameters,
                running_functions,
                dark_img,
                experiment_fpath,
                preprocessed_zarr_fpath=nuclei_filtered_fpath,
                save_steps_output=save_steps_output,
                dask_key_name=dask_delayed_name,
            )
            dask_delayed_name = "segment_nuclei_fov" + str(fov) + "_" + tokenize()
            mask_out = delayed(segmentation_NN_fov, name=dask_delayed_name)(
                fov_out[0][-1],
                fov_subdataset,
                segmented_file_path,
                fresh_tissue_segmentation_engine,
                diameter_size,
                model=scattered_mdel,
            )

            all_processing_nuclei.append(mask_out)
            # all_processing_nuclei.append(fov_out)
        # end = delayed(combine_steps)(saved_file,all_processing_nuclei)

        _ = dask.compute(all_processing_nuclei)
        client.run(gc.collect)

    return ds_beads, ds_nuclei, metadata


# TODO Remove functions


# def collect_bits_intensity_graph(dataset:pd.Dataframe, experiment_fpath:str,
#                             grpd_fovs,metadata:Dict,chunks_size,codebooks,client):

#     # Write to add stitching reference to the dataset to make in easy to run the
#     # bits analysis


#     # Need to run for fov

#     bit_channels = list(set(metadata['list_all_channels']).difference(metadata['stitching_channel'])))
#     experiment_fpath = Path(experiment_fpath)
#     experiment_name = metadata['experiment_name']
#     filtered_images_path =  Path(experiment_fpath) / (metadata['experiment_name'] + 'preprocessed_img_data.zarr')

#     all_fovs = list(grpd_fovs.groups.keys())
#     chunks = [all_fovs[x:x+chunks_size] for x in range(0, len(all_fovs), chunks_size)]
#     for chunk in chunks:
#         all_processing = []
#         all_filtered_images = []
#         for fov_num in chunk:
#             group = grpd_fovs.get_group(fov_num)
#             grpd_channel = group.groupby('channel')

#             # Load counts
#             name = 'Load_counts' +experiment_name + '_' + \
#                                         + '_fov_' +str(fov_num) + '-' + tokenize()

#             counts_fpath = experiment_fpath / 'results' / (experiment_name + '_decoded_fov_' +str(fov_num) + '.parquet')
#             counts_df = delayed(pd.read_parquet,name=name)(counts_fpath)


#             for channel, channel_grp in grpd_channel:
#                 all_filtered_images = []
#                 for index_value, fov_subdataset in channel_grp.iterrows():
#                     # Create ((img,), metadata) list to match the one used in the eel graph in order
#                     # to used the same set of functions

#                     name =  'load_processed_images_' + fov_subdataset.grp_name+ '-' + tokenize()

#                     img = delayed(io.load_general_zarr,name=name)(fov_subdataset,filtered_images_path)
#                     filt_out = ((img,), metadata)
#                     all_filtered_images.append(filt_out)


#                 name = 'register_and_combine_filt_imgs' +experiment_name + '_' + channel + '_' \
#                                         + '_fov_' +str(fov_num) + '-' + tokenize()

#                 combined_images = delayed(fovs_registration.combine_register_filtered_images,name=name)(all_filtered_images,counts_df,
#                                                                                                 fov_subdataset.stitching_channel)


#                 name = 'extract_barcodes_int' +experiment_name + '_' + channel + '_' \
#                                     + '_fov_' +str(fov_num) + '-' + tokenize()
#                 barcodes_max = delayed(barcodes_analysis.extract_dots_images,name=name)(counts_df,combined_images,
#                                                                                         experiment_fpath)


#                 name = 'extract_bit_flip_direction' +experiment_name + '_' + channel + '_' \
#                                     + '_fov_' +str(fov) + '-' + tokenize()
#                 flip_direction = delayed(barcodes_analysis.define_flip_direction,name=name)(codebook,
#                                                                                         experiment_fpath,
#                                                                                         stitched_coords)


# """Collected the intensity of the 1 and 0 bits and the flipping direction.
#     This is an extra function to run after the data are collected to avoid
#     to slow down the initial dots counting.

# Args:
#     dataset ([type]): [description]
# """
