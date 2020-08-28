import click

from cli_tranfer_files.commands import transfer_data_to_processing_hd

@click.group()
def pysmFISH_cli():
    pass

pysmFISH_cli.add_command(transfer_data_to_processing_hd)

if __name__ == '__main__':
    pysmFISH_cli()