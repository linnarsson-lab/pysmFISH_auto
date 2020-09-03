from prefect.tasks.github.issues import OpenGitHubIssue
from prefect import task
from prefect.triggers import all_successful, all_failed

from path import Path


@task(name='input_files_error', trigger=all_failed)
def report_input_files_errors(git_processing_issues_repo:str, experiment_fpath:str,token:str):
    experiment_fpath = Path(experiment_fpath)
    experiment_name = experiment_fpath.stem
    title = 'Processing of ' + experiment_name + ' failed because error in the input file'
    body = 'Before re run the processing of experiment ' + experiment_name + ' you need \
            to make sure that all the input files are correct'
    issue = OpenGitHubIssue(repo=git_processing_issues_repo,title=title,body=body)
    issue.run(repo=git_processing_issues_repo,title=title,body=body,token=token)


if __name__ == "__main__":
    from prefect import Flow
    from prefect.engine.executors.dask import LocalDaskExecutor
    with Flow("test_running") as flow:
        git_processing_issues_repo = 'linnarsson-lab/smFISH_processing_issues'
        experiment_name = 'test_issues'
        token='1e3d1e426d46007f876cc5550bd336729b0a0686'
        report_input_files_errors(git_processing_issues_repo,experiment_name,token)

    executor = LocalDaskExecutor()
    flow_state = flow.run(executor=executor)
