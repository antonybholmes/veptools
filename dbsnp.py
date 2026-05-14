# -*- coding: utf-8 -*-
"""
Generate a tss distribution for a region file

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import gzip
import numpy as np
import pandas as pd
from .utils import SEP, NA

DBSNP_COL = "dbSNP_RSID"


class DBSNPAnnotator:
    def __init__(self, vcf: str):
        self._vcf = vcf
        self._rsmap = {}

    def _load(self):
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
                chrom = fields[0]
                pos = int(fields[1])
                rsid = fields[2]
                ref = fields[3]
                alt = fields[4]
                info = fields[7]
                id = info.split("=")[1] if "=" in info else ""

                if id == "" or not rsid.startswith("rs"):
                    continue

                self._rsmap[id] = rsid

    def annotate(self, fin: str, fout: str, chunksize: int = 200000):

        self._load()

        pc = 0
        chunk = 1
        first = True

        # out = "bcca2024-16se_73primary_29cl_20icg_hg19.vep_annotated.splice.maf.v4.vcf.txt"  # re.sub(r"\.[^.]+$", "", f) + ".splice_site_annotated.xlsx"

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

                # if chr == "chrMT":
                #    chr = "chrM"

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
