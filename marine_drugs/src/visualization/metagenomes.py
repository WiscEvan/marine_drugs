#!/usr/bin/env python

import argparse
import logging
import os

import pandas as pd
import plotly.graph_objects as go

from Bio.SeqIO import parse
from typing import Union


logger = logging.getLogger(__name__)

logging.basicConfig(level=logging.INFO)


def is_valid_suffix(out: str) -> bool:
    """Check whether `out` has valid extension.
    Choices: [png, pdf, jpeg]

    Parameters
    ----------
    out : str
        Filename

    Returns
    -------
    bool
        Whether `out` contains extension in choices or not
    """
    suffixes = {".png", ".pdf", ".jpeg"}
    __, ext = os.path.splitext(out)
    if ext in suffixes:
        return True
    return False


def get_base_stats(assembly) -> pd.DataFrame:
    """Retrieve assembly GC content, length and coverage for each contig.

    Returns
    -------
    pd.DataFrame
        [description]
    """
    content = [
        {
            "contig": record.id,
            "GC": GC(record.seq),
            "length": len(record.seq),
            "coverage": float(record.id.rsplit("_cov_", 1)[-1]),
        }
        for record in parse(assembly, "fasta")
    ]
    return pd.DataFrame(content).set_index("contig")


def plot_gc_vs_coverage(
    assembly_or_dataframe: Union[str, pd.DataFrame], out: str = None, title: str = None
) -> go.Figure:
    """Returns plotly 2D scatter figure object after retrieving `assembly`
    GC contents, coverages and lengths. Y-axis is GC content and coverage is X-axis.
    Markers are sized by contig length

    Note
    ----
    To be used later in notebook with either::

    >>>fig = plot_gc_vs_coverage(...)
    >>>fig.show()

    or::

    >>>from IPython.display import Image
    >>>fig = plot_gc_vs_coverage(...)
    >>>img_bytes = fig.to_image("jpeg", width=1000)
    >>>Image(img_bytes)

    Parameters
    ----------
    assembly : [type]
        [description]
    out : str, optional
        [description], by default None
    title : str, optional
        [description], by default None

    Returns
    -------
    go.Figure
        Plotly Figure graph object of 2D scatterplot (GC% vs Coverage).
    """
    if isinstance(assembly_or_dataframe, pd.DataFrame):
        df = assembly_or_dataframe
        assert "coverage" in df.columns
        assert "length" in df.columns
        assert "GC" in df.columns
    elif isinstance(assembly_or_dataframe, str):
        exts = [".tsv", ".tab", ".csv", ".fasta", ".fna", ".fn"]
        seps = ["\t", "\t", ",", None, None, None]
        for ext, sep in zip(exts, seps):
            if ext in assembly_or_dataframe:
                if sep:
                    df = pd.read_csv(assembly_or_dataframe, sep=sep, index_col="contig")
                else:
                    df = get_base_stats(assembly)
                break
    fig = go.Figure()
    trace = go.Scattergl(
        x=df.coverage,
        y=df.GC,
        name=os.path.basename(assembly),
        mode="markers",
        marker_size=df.length,
    )
    fig.add_trace(trace)
    fig.update_xaxes(title_text="Coverage")
    # Set y-axes titles
    fig.update_yaxes(title_text="GC Content (%)")
    title = os.path.basename(assembly) if not title else title
    fig.update_layout(
        title_text=f"{title}<br>Num. Sequences: {df.shape[0]:,}",
        template="simple_white",
    )
    if out:
        if not is_valid_suffix(out):
            logger.warning(f"{out} file format not supported.. Skipping writing")
        else:
            fig.write_image(out, width=1000)
    return fig


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate plot given input data")
    parser.add_argument(
        "--metagenome-stats",
        help="Path to metagenome contigs stats containing contig, coverage, GC and length columns. (must contain .tsv, .tab, .csv extension in filename)",
    )
    parser.add_argument(
        "--assembly",
        help="Path to metagenome assembly (must contain .fasta, .fna or .fn extension in filename)",
    )
    parser.add_argument("--title", help="Title text to place in figure header")
    parser.add_argument("out", help="</path/to/write/figure.png>")
    args = parser.parse_args()
    if not args.metagenome_stats and not args.assembly:
        raise ValueError("One of --assembly or --metagenome-stats is required")

    if args.metagenome_stats:
        plot_gc_vs_coverage(
            assembly_or_dataframe=args.metagenome_stats, out=args.out, title=args.title
        )
    else:
        plot_gc_vs_coverage(
            assembly_or_dataframe=args.assembly, out=args.out, title=args.title
        )
