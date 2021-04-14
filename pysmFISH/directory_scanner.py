from pathlib import Path
import os

from flows.general_flow_single_fov import general_flow
from pysmFISH.utils import free_space
from pysmFISH.logger_utils import selected_logger


def chuparana():
    print(f'executing function in runs')
    time.sleep(10)


def end_processing_file(path_destination, completion_pattern='processing_completed.txt'):
    
    logger = selected_logger()
    
    fname = Path(path_destination) / completion_pattern
    try:
        open(fname, 'a').close()
    except OSError:
        logger.error('Failed to create the processing termination file')
    else:
        logger.info('processing termination files created')



def experiment_processing_runner(path_source:str,
                                         flag_file_key:str,
                                         completion_pattern:str,
                                         min_free_space:int=1000,
                                         monitor_frequency:int=43200):

    logger = selected_logger()
    path_source = Path(path_source)

    try:

        while True:
            files = list(path_source.glob('*' + flag_file_key))
            files.sort(key=lambda x: os.path.getmtime(x))
            logger.error(f'processing files {files}')
            logger.error(f'experiments to process {len(files)}')
            # Check that the previous file has been processed
            completion_file = list(path_source.glob(completion_pattern))

            if len(completion_file):
                logger.error(f'completion file present')
                # Check if there is enough space to run processing
                allowed_space = free_space(path_source, min_free_space)
                for flagged_file_path in files:
                    if allowed_space:
                        try:
                            experiment_name = (flagged_file_path.name).split(flag_file_key)[0]
                            experiment_fpath = flagged_file_path.parent / experiment_name
                            os.stat(experiment_fpath)
                        except:
                            logger.error(f'No experiment folder matching the flagged file')
                            logger.error(f'The flagged file is removed')
                            os.remove(flagged_file_path.as_posix())
                        else:
                            experiment_fpath = flagged_file_path.parent / experiment_name
                            os.remove(completion_file[0].as_posix())
                            try:
                                print(f'RUN PROCESSING {experiment_fpath.stem}')
                                chuparana()
                                status = 'SUCCESS'
                            except:
                                logger.error(f'Broken processing')
                                os.remove(flagged_file_path.as_posix())
                                experiment_fpath.rename(experiment_fpath.parent / (experiment_fpath.stem + '_PROCESSING_FAILED'))
                            else:
                                
                                if status == 'SUCCESS':
                                    allowed_space = free_space(path_source, min_free_space)
                                    end_processing_file(path_source,
                                                completion_pattern='processing_completed.txt')
                                    os.remove(flagged_file_path.as_posix())
                                
            else:
                logger.error(f'completion file is missing')

            time.sleep(monitor_frequency)

    except KeyboardInterrupt:
        logger.info(f'scanning interrupted')
        pass
        # log output

