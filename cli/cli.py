import click

from cli_viz.commands import viz_utilities

@click.group()
def pysmFISH_cli():
    pass

pysmFISH_cli.add_command(viz_utilities)

if __name__ == '__main__':
    pysmFISH_cli()