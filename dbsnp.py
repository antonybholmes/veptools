# -*- coding: utf-8 -*-
"""
Generate a tss distribution for a region file

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import gzip
import pandas as pd

from .vcf import vcf_record_info_fields
from .utils import NA

DBSNP_COL = "dbSNP_RSID"

DBSNP_COLS = [
    {"id": "RS", "header": "dbSNP_RSID"},  # dbSNP RSID
]


class DBSNPAnnotator:
    def __init__(self, vcf: str):
        self._vcf = vcf
        self._dbsnp_map = {}

    def _load(self):
        if len(self._dbsnp_map) > 0:
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

                info = fields[7]

                data = vcf_record_info_fields(info)

                id = data.get("VEP_ID", "")

                if id == "":
                    continue

                self._dbsnp_map[id] = data

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

                for c in DBSNP_COLS:
                    df.at[index, c["header"]] = self._dbsnp_map.get(id, {}).get(
                        c["id"], NA
                    )

            print(f"Processed {pc} splice site variants")

            df.to_csv(
                fout,
                sep="\t",
                header=first,
                mode="w" if first else "a",
                index=False,
            )

            first = False
