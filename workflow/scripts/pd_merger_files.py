from os import path
from functools import reduce, partial
import pandas as pd
from os import path
import click


@click.command()
@click.option(
    "-f",
    "--file_to_merge",
    multiple=True,
    help="each file to merge, files should be identically formatted.",
)
@click.option(
    "-o",
    "--output_file",
    help="Output for the merged file",
)
@click.option(
    "-n",
    "--name",
    multiple=True,
    help="Each name for each column after merge.",
)
@click.option("--index", default=0, help="pd.read_csv index argumet", type=int)
@click.option(
    "--merge_column", default=1, help="Column from each file to merg", type=int
)
@click.option("--sep", default="\t", help="CSV seperator")
@click.option("--header/--no-header", default=True, help="is there a header?")
def pa_merge_file(
    file_to_merge=(), name=(), index=0, merge_column=None, sep="\t", header=True, output_file='merge.txt'
):
    if header:
        files = [pd.read_csv(i, sep=sep, index_col=index) for i in file_to_merge]
    else:
        files = [
            pd.read_csv(i, sep=sep, index_col=index, header=header)
            for i in file_to_merge
        ]
    files = [i.iloc[:, merge_column] for i in files]
    if len(name) > 1:
        for i, j in zip(files, name):
            i.name = j
    else:
        for i, j in zip(files, file_to_merge):
            i.name = path.basename(j)

    merged_out=reduce(partial(pd.merge, left_index=True, right_index=True, how="outer"), files)
    merged_out.to_csv(output_file, sep=sep, index=True)




if __name__ == "__main__":
    pa_merge_file()

