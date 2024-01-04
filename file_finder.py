import os
import click

@click.command()
@ckick.argument('file_path')

def main():
    files = os.listdir(file_path)
    for file in files
