#!/usr/bin/env python

import argparse
import logging
import os

import multiprocessing as mp
import numpy as np
import pandas as pd

from Bio.SeqIO import parse
from Bio.SeqUtils import GC

from pathlib import Path
from typing import List
from functools import lru_cache


logger = logging.getLogger(__name__)

logging.basicConfig(level=logging.DEBUG)


class Metagenome:
    def __init__(self, assembly: str):
        self.assembly = assembly

    @property
    @lru_cache(maxsize=None)
    def filename(self) -> str:
        # Example location: data/assemblies/FL2015-34.fasta
        return os.path.basename(self.assembly)

    @property
    @lru_cache(maxsize=None)
    def size(self) -> int:
        df = self.get_base_stats()
        return df.length.sum()

    @property
    @lru_cache(maxsize=None)
    def records(self) -> list:
        return [record for record in parse(self.assembly, "fasta")]

    @lru_cache(maxsize=None)
    def get_base_stats(self) -> pd.DataFrame:
        """Retrieve metagenome GC content, length and coverage for each contig.

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
            for record in self.records
        ]
        return pd.DataFrame(content).set_index("contig")

    def fragmentation_metric(self, quality_measure=0.50) -> float:
        """Describes the quality of assembled genomes that are fragmented in
        contigs of different length.

        For more information see:
            http://www.metagenomics.wiki/pdf/definition/assembly/n50

        Parameters
        ----------
        quality_measure : 0 < float < 1
            Description of parameter `quality_measure` (the default is .50).
            I.e. default measure is N50, but could use .1 for N10 or .9 for N90

        Returns
        -------
        int
            Minimum contig length to cover `quality_measure` of genome (i.e. length
            weighted median)

        """
        if not 0 < quality_measure < 1:
            raise ValueError(
                f"quality measure must be b/w 0 and 1! given: {quality_measure}"
            )
        target_size = self.size * quality_measure
        lengths = []
        for length in sorted([len(rec.seq) for rec in self.records], reverse=True):
            lengths.append(length)
            if sum(lengths) > target_size:
                return length

    def describe(self) -> pd.DataFrame:
        df = self.get_base_stats()
        length_weighted_gc = np.average(a=df.GC, weights=df.length / df.length.sum())
        length_weighted_coverage = np.average(
            a=df.coverage, weights=df.length / df.length.sum()
        )
        return {
            "Assembly": self.filename,
            "Num. Sequences": df.shape[0],
            "Size (bp)": self.size,
            "N50": self.fragmentation_metric(0.5),
            "N10": self.fragmentation_metric(0.1),
            "N90": self.fragmentation_metric(0.9),
            "Largest sequence": df.length.idxmax(),
            "Length Weighted Avg. GC content": length_weighted_gc,
            "Min GC content": df.GC.min(),
            "Max GC content": df.GC.max(),
            "Std. dev. GC content": df.GC.std(),
            "Length Weighted Avg. Coverage": length_weighted_coverage,
            "Min Coverage": df.coverage.min(),
            "Max Coverage": df.coverage.max(),
            "Std. dev. Coverage": df.coverage.std(),
        }


def assembly_statter(assembly):
    logger.info(f"describing {assembly}")
    return Metagenome(assembly).describe()


def pool_stats(assemblies: List[str], cpus: int = mp.cpu_count()) -> pd.DataFrame:
    """Multiprocessing metagenome statter"""

    cpus = len(assemblies) if cpus > len(assemblies) else cpus
    pool = mp.Pool(cpus)
    logger.debug(
        f" Describing stats for {len(assemblies):,} metagenomes using {cpus} cpus."
    )
    stats = pool.map(assembly_statter, assemblies)
    pool.close()
    pool.join()
    return pd.DataFrame(stats)


def base_statter(args):
    assembly, outdir = args
    logger.info(f" Describing {assembly} base stats.")
    metagenome = Metagenome(assembly)
    df = metagenome.get_base_stats()
    suffix = Path(metagenome.filename).suffix
    filename = metagenome.filename.replace(suffix, ".stats.tsv")
    out = os.path.join(outdir, filename)
    df.to_csv(out, sep="\t", index=True, header=True)
    logger.info(f" Finished describing {assembly} base stats.")
    return out


def pool_base_statter(
    assemblies: List[str], outdir: str, cpus: int = mp.cpu_count()
) -> List[str]:

    cpus = len(assemblies) if cpus > len(assemblies) else cpus
    pool = mp.Pool(cpus)
    logger.debug(
        f" Retrieving base stats for {len(assemblies):,} metagenomes using {cpus} cpus."
    )
    args = [(assembly, outdir) for assembly in assemblies]
    stats = pool.map(base_statter, args)
    pool.close()
    pool.join()
    return stats


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Stat sponge metagenome(s)")
    parser.add_argument("assembly", nargs="+")
    parser.add_argument(
        "out", help="path to write sponge metagenomes table of statistics"
    )
    parser.add_argument(
        "--write-base-stats",
        help="Write base statistics for each contig in metagenome to provided directory. (File naming will follow assembly.stats.tsv)",
    )

    args = parser.parse_args()
    if not os.path.exists(args.out):
        if len(args.assembly) > 1:
            df = pool_stats(args.assembly)
        else:
            df = pd.DataFrame([assembly_statter(args.assembly)])
        df.to_csv(args.out, sep="\t", index=False, header=True)
        logger.info(f"Wrote stats of {df.shape[0]} sponges to {args.out}")
    else:
        logger.info(f"{args.out} already exists. Skipping...")

    if args.write_base_stats:
        if not os.path.exists(args.write_base_stats):
            logger.info(f"Created {args.write_base_stats}")
            os.makedirs(args.write_base_stats)
        if len(args.assembly) > 1:
            filepaths = pool_base_statter(
                assemblies=args.assembly, outdir=args.write_base_stats
            )
        else:
            arg = (args.assembly, args.write_base_stats)
            filepaths = [base_statter(arg)]
