from prefect.tasks.notifications import github
from prefect import task

@task(name='input_files_error')
def report_input_files_errors(git_processing_issues_repo:str, experiment_name:str):
    title = 'Processing of ' + experiment_name + 'failed because error in the input file'
    body = 'Before re run the processing of experiment ' + experiment_name + 'you need \
            to make sure that all the input files are correct'
