import click

from cli_prep.commands import setup_processing_env
from cli_folder_scanner.commands import folder_monitoring_utilities
# from cli_run_flows.commands import flows_runner


@click.group()
def pysmFISH_cli():
    pass

pysmFISH_cli.add_command(setup_processing_env)
pysmFISH_cli.add_command(folder_monitoring_utilities)
# pysmFISH_cli.add_command(flows_runner)


if __name__ == '__main__':
    pysmFISH_cli()