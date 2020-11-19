""" 
Set of tasks to run prefect jobs 
"""


from pathlib import Path

# prefect related imports
import prefect
from prefect import task
from prefect.engine import signals

# pysmFISH imports
from pysmFISH.utilities_tasks import load_raw_images
from pysmFISH.io import save_images_metadata, save_dots_data
from pysmFISH.dots_calling import osmFISH_peak_based_detection
from pysmFISH.preprocessing_tasks import preprocessing_dot_raw_image, load_dark_image
from pysmFISH.logger_utils import setup_logger_prefect_UI


# testing import
from pathlib import Path










@task(task_run_name=lambda **kwargs: f"fish-preprocessing-{kwargs['zarr_grp_name']}")
def single_fish_filter_count(zarr_grp_name,
                    parsed_raw_data_fpath,
                    experiment_fpath,
                    FlatFieldKernel,FilteringSmallKernel, 
                    LaplacianKernel,
                    min_distance,
                    min_obj_size,
                    max_obj_size,
                    num_peaks_per_label):

    logger = prefect.context.get("logger")
    # logger = setup_logger_prefect_UI()
    try:
        raw_fish_images_meta = load_raw_images(zarr_grp_name,
                                    parsed_raw_data_fpath)
    except:
        logger.error(f'cannot load {zarr_grp_name} raw fish image')
        signals.FAIL(f'cannot load {zarr_grp_name} raw fish image')
    else:
        logger.info(f'loaded {zarr_grp_name} raw fish image')
        try:
            dark_img = load_dark_image(experiment_fpath)
        except:
            logger.error(f'cannot load dark reference fish image')
            signals.FAIL(f'cannot load dark reference fish image')
        else:
            logger.info('loaded dark reference image')
            try:
                filtered_fish_images_metadata = preprocessing_dot_raw_image(raw_fish_images_meta,dark_img,
                                        FlatFieldKernel,FilteringSmallKernel, 
                                        LaplacianKernel)
            except:
                  logger.error(f'cannot filter fish image {zarr_grp_name}')
                  signals.FAIL(f'cannot filter fish image {zarr_grp_name}')

            else:
                logger.info(f'filtered fish image {zarr_grp_name}')
                try:
                    save_images_metadata(filtered_fish_images_metadata)
                except:
                    logger.error(f'cannot save fish image {zarr_grp_name}')
                    signals.FAIL(f'cannot save fish image {zarr_grp_name}')
            
                else:
                    logger.info(f'saved filtered fish image {zarr_grp_name}')
                    try:
                        fish_counts = osmFISH_peak_based_detection(filtered_fish_images_metadata,
                                                    min_distance,
                                                    min_obj_size,
                                                    max_obj_size,
                                                    num_peaks_per_label)
                    except:
                        logger.error(f'cannot count dots in fish image {zarr_grp_name}')
                        signals.FAIL(f'cannot count dots in fish image {zarr_grp_name}')

                    else:
                        logger.info(f'counted dots in fish image {zarr_grp_name}')
                        try:
                            save_dots_data(fish_counts)
                        except:
                            logger.error(f'cannot save the counts of fish image {zarr_grp_name}')
                            signals.FAIL(f'cannot save the counts of fish image {zarr_grp_name}')
                        else:
                            logger.info(f'completed preprocessing and counting fish image {zarr_grp_name}')


@task(task_run_name=lambda **kwargs: f"beads-preprocessing-{kwargs['zarr_grp_name']}")
def single_beads_filter_count(zarr_grp_name,
                    parsed_raw_data_fpath,
                    experiment_fpath,
                    FlatFieldKernel,FilteringSmallKernel, 
                    LaplacianKernel,
                    min_distance,
                    min_obj_size,
                    max_obj_size,
                    num_peaks_per_label):

    logger = setup_logger_prefect_UI()
    try:
        raw_beads_images_meta = load_raw_images(zarr_grp_name,
                                    parsed_raw_data_fpath)
    except:
        logger.error(f'cannot load {zarr_grp_name} raw beads image')
        signals.FAIL(f'cannot load {zarr_grp_name} raw beads image')
    else:
        logger.info(f'loaded {zarr_grp_name} raw beads image')
        try:
            dark_img = load_dark_image(experiment_fpath)
        except:
            logger.error(f'cannot load dark reference beads image')
            signals.FAIL(f'cannot load dark reference beads image')
        else:
            logger.info('loaded dark reference image')
            try:
                filtered_beads_images_metadata = preprocessing_dot_raw_image(raw_beads_images_meta,dark_img,
                                        FlatFieldKernel,FilteringSmallKernel, 
                                        LaplacianKernel)
            except:
                  logger.error(f'cannot filter beads image {zarr_grp_name}')
                  signals.FAIL(f'cannot filter beads image {zarr_grp_name}')

            else:
                logger.info(f'filtered beads image {zarr_grp_name}')
                try:
                    save_images_metadata(filtered_beads_images_metadata)
                except:
                    logger.error(f'cannot save beads image {zarr_grp_name}')
                    signals.FAIL(f'cannot save beads image {zarr_grp_name}')
            
                else:
                    logger.info(f'saved filtered beads image {zarr_grp_name}')
                    try:
                        beads_counts = osmFISH_peak_based_detection(filtered_beads_images_metadata,
                                                    min_distance,
                                                    min_obj_size,
                                                    max_obj_size,
                                                    num_peaks_per_label)
                    except:
                        logger.error(f'cannot count dots in beads image {zarr_grp_name}')
                        signals.FAIL(f'cannot count dots in beads image {zarr_grp_name}')

                    else:
                        logger.info(f'counted dots in beads image {zarr_grp_name}')
                        try:
                            save_dots_data(beads_counts)
                        except:
                            logger.error(f'cannot save the counts of beads image {zarr_grp_name}')
                            signals.FAIL(f'cannot save the counts of beads image {zarr_grp_name}')
                        else:
                            logger.info(f'completed preprocessing and counting beads image {zarr_grp_name}')