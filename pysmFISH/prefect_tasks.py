""" 
Set of tasks to run prefect jobs 
"""


from pathlib import Path

# prefect related imports
import prefect
from prefect import task
from prefect import Task
from prefect.engine import signals

# pysmFISH imports
from pysmFISH.utilities_tasks import load_raw_images
from pysmFISH.io import save_images_metadata, save_dots_data
from pysmFISH.dots_calling import osmFISH_peak_based_detection
from pysmFISH.preprocessing_tasks import preprocessing_dot_raw_image, load_dark_image
from pysmFISH.logger_utils import setup_logger_prefect_UI


# testing import
from pathlib import Path



class single_fish_filter_count(Task):
    """
    Task to:
    - preprocess the fish images
    - count the dots and save the data
    
    Args:
    -----
        zarr_grp_name: str
            group representing the image to process
        parsed_raw_data_fpath: str
            path to the zarr file containing the parsed images
        FlatFieldKernel: np.ndarray
            size of the kernel use to remove the backgroun (ex. [10,10])
        FilteringSmallKernel: np.ndarray
            size of the kernel use to smooth the dots (ex. [8,8])
        LaplacianKernel: np.ndarray
            size of the kernel use enhance dots signal (ex. [1,1])
        min_distance: int
            minimum distance between two dots
        min_obj_size: int
            minimum size of a dot
        max_obj_size: int
            max size of a dot or a cluster of dots (depending on num_peaks_per_label)
        num_peaks_per_label
            max number of peaks called in a masked object
    
    """

    def run(self,
            zarr_grp_name,
            parsed_raw_data_fpath,
            FlatFieldKernel,
            FilteringSmallKernel, 
            LaplacianKernel,
            min_distance,
            min_obj_size,
            max_obj_size,
            num_peaks_per_label):


        """
        Function to:
        - preprocess the fish images
        - count the dots and save the data
        
        Args:
        -----
            zarr_grp_name: str
                group representing the image to process
            parsed_raw_data_fpath: str
                path to the zarr file containing the parsed images
            FlatFieldKernel: np.ndarray
                size of the kernel use to remove the backgroun (ex. [10,10])
            FilteringSmallKernel: np.ndarray
                size of the kernel use to smooth the dots (ex. [8,8])
            LaplacianKernel: np.ndarray
                size of the kernel use enhance dots signal (ex. [1,1])
            min_distance: int
                minimum distance between two dots
            min_obj_size: int
                minimum size of a dot
            max_obj_size: int
                max size of a dot or a cluster of dots (depending on num_peaks_per_label)
            num_peaks_per_label
                max number of peaks called in a masked object
        
        """

        parsed_raw_data_fpath = Path(parsed_raw_data_fpath)
        experiment_fpath = parsed_raw_data_fpath.parent
        
        try:
            raw_fish_images_meta = load_raw_images(zarr_grp_name,
                                        parsed_raw_data_fpath)
        except:
            self.logger.error(f'cannot load {zarr_grp_name} raw fish image')
            signals.FAIL(f'cannot load {zarr_grp_name} raw fish image')
        else:
            self.logger.info(f'loaded {zarr_grp_name} raw fish image')
            try:
                # This may chnaged if the image will be store in shoji
                dark_img = load_dark_image(experiment_fpath)
            except:
                self.logger.error(f'cannot load dark reference fish image')
                signals.FAIL(f'cannot load dark reference fish image')
            else:
                self.logger.info('loaded dark reference image')
    
                filtered_fish_images_metadata = preprocessing_dot_raw_image(raw_fish_images_meta,dark_img,
                                        FlatFieldKernel,FilteringSmallKernel, 
                                        LaplacianKernel)

                save_images_metadata(filtered_fish_images_metadata)

                fish_counts = osmFISH_peak_based_detection(filtered_fish_images_metadata,
                                                        min_distance,
                                                        min_obj_size,
                                                        max_obj_size,
                                                        num_peaks_per_label)
                save_dots_data(fish_counts)

                # except:
                #     self.logger.error(f'cannot filter fish image {zarr_grp_name}')
                #     signals.FAIL(f'cannot filter fish image {zarr_grp_name}')

                # else:
                #     self.logger.info(f'filtered fish image {zarr_grp_name}')
                #     try:
                #         save_images_metadata(filtered_fish_images_metadata)
                #     except:
                #         self.logger.error(f'cannot save fish image {zarr_grp_name}')
                #         signals.FAIL(f'cannot save fish image {zarr_grp_name}')
                
                #     else:
                #         self.logger.info(f'saved filtered fish image {zarr_grp_name}')
                #         try:
                #             fish_counts = osmFISH_peak_based_detection(filtered_fish_images_metadata,
                #                                         min_distance,
                #                                         min_obj_size,
                #                                         max_obj_size,
                #                                         num_peaks_per_label)
                #         except:
                #             self.logger.error(f'cannot count dots in fish image {zarr_grp_name}')
                #             signals.FAIL(f'cannot count dots in fish image {zarr_grp_name}')

                #         else:
                #             self.logger.info(f'counted dots in fish image {zarr_grp_name}')
                #             try:
                #                 save_dots_data(fish_counts)
                #             except:
                #                 self.logger.error(f'cannot save the counts of fish image {zarr_grp_name}')
                #                 signals.FAIL(f'cannot save the counts of fish image {zarr_grp_name}')
                #             else:
                #                 self.logger.info(f'completed preprocessing and counting fish image {zarr_grp_name}')



class single_beads_filter_count(Task):
    """
    Task to:
    - preprocess the beads images
    - count the dots and save the data
    
    Args:
    -----
        zarr_grp_name: str
            group representing the image to process
        parsed_raw_data_fpath: str
            path to the zarr file containing the parsed images
        FlatFieldKernel: np.ndarray
            size of the kernel use to remove the backgroun (ex. [10,10])
        FilteringSmallKernel: np.ndarray
            size of the kernel use to smooth the dots (ex. [8,8])
        LaplacianKernel: np.ndarray
            size of the kernel use enhance dots signal (ex. [1,1])
        min_distance: int
            minimum distance between two dots
        min_obj_size: int
            minimum size of a dot
        max_obj_size: int
            max size of a dot or a cluster of dots (depending on num_peaks_per_label)
        num_peaks_per_label
            max number of peaks called in a masked object
    
    """
    def run(self, zarr_grp_name,
                        parsed_raw_data_fpath,
                        FlatFieldKernel,
                        FilteringSmallKernel, 
                        LaplacianKernel,
                        min_distance,
                        min_obj_size,
                        max_obj_size,
                        num_peaks_per_label):

        """
        Function to:
        - preprocess the fish images
        - count the dots and save the data
        
        Args:
        -----
            zarr_grp_name: str
                group representing the image to process
            parsed_raw_data_fpath: str
                path to the zarr file containing the parsed images
            FlatFieldKernel: np.ndarray
                size of the kernel use to remove the backgroun (ex. [10,10])
            FilteringSmallKernel: np.ndarray
                size of the kernel use to smooth the dots (ex. [8,8])
            LaplacianKernel: np.ndarray
                size of the kernel use enhance dots signal (ex. [1,1])
            min_distance: int
                minimum distance between two dots
            min_obj_size: int
                minimum size of a dot
            max_obj_size: int
                max size of a dot or a cluster of dots (depending on num_peaks_per_label)
            num_peaks_per_label
                max number of peaks called in a masked object
        
        """

        parsed_raw_data_fpath = Path(parsed_raw_data_fpath)
        experiment_fpath = parsed_raw_data_fpath.parent

        try:
            raw_beads_images_meta = load_raw_images(zarr_grp_name,
                                        parsed_raw_data_fpath)
        except:
            self.logger.error(f'cannot load {zarr_grp_name} raw beads image')
            signals.FAIL(f'cannot load {zarr_grp_name} raw beads image')
        else:
            self.logger.info(f'loaded {zarr_grp_name} raw beads image')
            try:
                dark_img = load_dark_image(experiment_fpath)
            except:
                self.logger.error(f'cannot load dark reference beads image')
                signals.FAIL(f'cannot load dark reference beads image')
            else:
                self.logger.info('loaded dark reference image')
                
                filtered_beads_images_metadata = preprocessing_dot_raw_image(raw_beads_images_meta,dark_img,
                                        FlatFieldKernel,FilteringSmallKernel, 
                                        LaplacianKernel)

                save_images_metadata(filtered_beads_images_metadata)

                beads_counts = osmFISH_peak_based_detection(filtered_beads_images_metadata,
                                                    min_distance,
                                                    min_obj_size,
                                                    max_obj_size,
                                                    num_peaks_per_label)
                save_dots_data(beads_counts)

                # except:
                #       logger.error(f'cannot filter beads image {zarr_grp_name}')
                #       signals.FAIL(f'cannot filter beads image {zarr_grp_name}')

                # else:
                #     logger.info(f'filtered beads image {zarr_grp_name}')
                #     try:
                #         save_images_metadata(filtered_beads_images_metadata)
                #     except:
                #         logger.error(f'cannot save beads image {zarr_grp_name}')
                #         signals.FAIL(f'cannot save beads image {zarr_grp_name}')
                
                #     else:
                #         logger.info(f'saved filtered beads image {zarr_grp_name}')
                #         try:
                #             beads_counts = osmFISH_peak_based_detection(filtered_beads_images_metadata,
                #                                         min_distance,
                #                                         min_obj_size,
                #                                         max_obj_size,
                #                                         num_peaks_per_label)
                #         except:
                #             logger.error(f'cannot count dots in beads image {zarr_grp_name}')
                #             signals.FAIL(f'cannot count dots in beads image {zarr_grp_name}')

                #         else:
                #             logger.info(f'counted dots in beads image {zarr_grp_name}')
                #             try:
                #                 save_dots_data(beads_counts)
                #             except:
                #                 logger.error(f'cannot save the counts of beads image {zarr_grp_name}')
                #                 signals.FAIL(f'cannot save the counts of beads image {zarr_grp_name}')
                #             else:
                #                 logger.info(f'completed preprocessing and counting beads image {zarr_grp_name}')