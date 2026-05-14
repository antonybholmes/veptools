# -*- coding: utf-8 -*-
"""
Generate a tss distribution for a region file

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import gzip
import pandas as pd
from .utils import NA

DBSNP_COL = "dbSNP_RSID"


class GnomadAnnotator:
    def __init__(self, vcf: str):
        """
        Annotate a MAF file with gnomAD annotations from a VCF file.
        The VCF file should have been preprocessed to have an ID field in the INFO column
        that matches the format "chr_pos_ref/alt" for exact matching.
        Arguments:
            vcf:    path to the matching VCF that had gnomAD annotations added to it
                    (can be gzipped). This is not the original gnomAD VCF.
        """
        self._vcf = vcf
        self._rsmap = {}

    def _load(self):
        """
        Lazy load the VCF file and build a mapping from the ID field to the gnomAD annotations.
        The ID field should be in the format "chr_pos_ref/alt" for exact matching.
        """
        if len(self._rsmap) > 0:
            return

        if self._vcf.endswith(".gz"):
            open_func = gzip.open(self._vcf, "rt")
        else:
            open_func = open(self._vcf, "r")
        with open_func as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                # chrom = fields[0]
                # pos = int(fields[1])
                rsid = fields[2]
                # ref = fields[3]
                # alt = fields[4]
                info = fields[7]

                # will contain VEP_ID=, plus gnomAD annotations like gnomAD_AF=, gnomAD_AFR_AF=, etc.
                info_tokens = info.split(";")

                id = info_tokens[0].split("=")[1] if "=" in info_tokens[0] else ""

                if id == "" or not rsid.startswith("rs"):
                    continue

                self._rsmap[id] = rsid

    def annotate(self, fin: str, fout: str, chunksize: int = 200000):

        self._load()

        pc = 0
        chunk = 1
        first = True

        for df in pd.read_csv(
            fin,
            sep="\t",
            header=0,
            keep_default_na=False,
            chunksize=chunksize,
        ):
            print(f"Processing chunk {chunk}...")
            chunk += 1

            df[DBSNP_COL] = NA

            for index, row in df.iterrows():
                pc += 1

                chr = row["Chromosome"]
                start = row["Start_Position"]
                ref = row["Reference_Allele"]
                alt = row["Tumor_Seq_Allele2"]

                if chr == "chrMT":
                    chr = "chrM"

                # exact
                id = f"{chr}_{start}_{ref}/{alt}"

                if id in self._rsmap:
                    df.at[index, DBSNP_COL] = self._rsmap[id]

            print(f"Processed {pc} splice site variants")

            df.to_csv(
                fout,
                sep="\t",
                header=first,
                mode="w" if first else "a",
                index=False,
            )

            first = False
