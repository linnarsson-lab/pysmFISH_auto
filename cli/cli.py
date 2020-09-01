import click

from cli_tranfer_files.commands import data_transfer

@click.group()
def pysmFISH_cli():
    pass

pysmFISH_cli.add_command(data_transfer)

if __name__ == '__main__':
    pysmFISH_cli()