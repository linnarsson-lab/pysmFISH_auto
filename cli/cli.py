import click

from cli_viz.commands import viz_utilities
from cli_data_org.commands import data_organization_utilities
from cli_prep.commands import setup_processing_env
from cli_run_flows.commands import flows_runner
from cli_folder_scanner.commands import folder_monitoring_utilities


@click.group()
def pysmFISH_cli():
    pass

pysmFISH_cli.add_command(viz_utilities)
pysmFISH_cli.add_command(data_organization_utilities)
pysmFISH_cli.add_command(setup_processing_env)
pysmFISH_cli.add_command(flows_runner)
pysmFISH_cli.add_command(folder_monitoring_utilities)



if __name__ == '__main__':
    pysmFISH_cli()