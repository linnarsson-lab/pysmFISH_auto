import click

from cli_tranfer_files.commands import data_transfer
from cli_viz import vizcounts

@click.group()
def pysmFISH_cli():
    pass

pysmFISH_cli.add_command(data_transfer)
pysmFISH_cli.add_command(vizcounts)

if __name__ == '__main__':
    pysmFISH_cli()