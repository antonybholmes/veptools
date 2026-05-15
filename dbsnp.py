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

                # print(id, data)

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

            for c in DBSNP_COLS:
                df[c["header"]] = NA

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
                    v = self._dbsnp_map.get(id, {}).get(c["id"], NA)

                    # prefix with "rs" if it's an RSID and doesn't already start with "rs"
                    if c["id"] == "RS" and v != NA and not str(v).startswith("rs"):
                        v = f"rs{v}"

                    df.at[index, c["header"]] = v

            print(f"Processed {pc} splice site variants")

            df.to_csv(
                fout,
                sep="\t",
                header=first,
                mode="w" if first else "a",
                index=False,
            )

            first = False
