import collections
import os
import re
from logging import info
from urllib.parse import unquote

import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Font

SEP = "|"
NA = "."

COLUMNS = [
    "VEP_Gene_Symbol",
    "VEP_Biotype",
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Location",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "Variant_Type",
    "VEP_HGVSc",
    "VEP_HGVSp",
    "VEP_Exon",
    "VEP_Exons",
    "VEP_Total_Exons",
    "VEP_Variant_Classification",
    "Sample",
    "Tumor_Sample_Barcode",
    "Dataset",
    "Annotation_Database",
    "db",
    "Notes",
    "VEP_Gene_ID",
    "VEP_Is_Hugo_Gene",
    "VEP_Transcript",
    "VEP_Canonical",
    "VEP_Is_Canonical",
    "MANE_RefSeq",
    "MANE_status",
    "VEP_Variant_Severity",
    "Hugo_Symbol (Rahul - not really Hugo)",
    "Variant_Classification",
    "Protein_Change",
    "amino_acid_change",
    "AAChange",
    "AA",
    "cDNAChange",
    "codon_change",
    "t_alt_count",
    "t_depth",
    "VAF",
    "CCDS",
    "CCDS_AA_Length",
    "VEP_Secondary_Gene_Symbol",
    "VEP_Secondary_HGVSp",
    "VEP_Secondary_HGVSc",
    "VEP_Secondary_Variant_Classification",
    "VEP_Secondary_Variant_Severity",
    "VEP_Secondary_Gene_ID",
    "VEP_Secondary_Transcript",
    "VEP_Secondary_Exon",
    "VEP_Secondary_Exons",
    "VEP_Secondary_Total_Exons",
    "VEP_Secondary_Canonical",
    "VEP_Secondary_Biotype",
    "Secondary_CCDS",
    "Secondary_CCDS_AA_Length",
    "Splice_Protein_Change (AH)",
    "Splice_Transcript_ID",
    "Splice_Transcript_Type (1=canonical, 2=longest)",
    "Splice_Transcript_Type (1=exact, 2=canonical, 3=longest)",
    "Splice_Nearest_Exon_Number",
    "Splice_Nearest_Exon_Dist_bp",
    "Is_Strict_Splice_+/-_2bp",
    "Splice_Annotation_Database",
]


def reorder_columns(fin: str, fout: str, chunk_size: int = 200000):

    first = True
    chunk = 1

    for df in pd.read_csv(
        fin,
        sep="\t",
        header=0,
        keep_default_na=False,
        chunksize=chunk_size,
    ):
        print(f"Processing chunk {chunk}...")
        chunk += 1

        df.rename(
            columns={
                "DNAChange": "cDNAChange",
                "AAChange": "Protein_Change",
                "Sample": "Tumor_Sample_Barcode",
                "VEP_Annotation_Database": "db",
            },
            inplace=True,
        )

        current_columns = df.columns.tolist()

        # move columns that are in desired_order to the front, in the order specified by desired_order
        new_order = [col for col in COLUMNS if col in current_columns]
        new_order += [col for col in current_columns if col not in new_order]
        df = df[new_order]

        df.to_csv(
            fout,
            sep="\t",
            header=first,
            mode="w" if first else "a",
            index=False,
        )
        first = False
