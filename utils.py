from typing import Callable

import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Font

# optional: enforce specific column types
dtype_map = {
    "VEP_Gene_Symbol": "string",
    "Hugo_Symbol (Rahul - not really Hugo)": "string",
}


def maf_to_excel(df, fout):
    df = df.copy()

    for col, typ in dtype_map.items():
        if col in df.columns:
            df[col] = df[col].astype(typ)

    df.to_excel(fout, index=False)

    arial = Font(name="Arial")
    # define font
    header_font = Font(name="Arial", bold=True)

    # reopen workbook
    wb = load_workbook(fout)
    ws = wb.active

    for col in ws.columns:
        for cell in col:
            cell.font = arial

    # apply to header row
    for cell in ws[1]:
        cell.font = header_font

    wb.save(fout)


def chr_sizes(assembly: str = "hg19"):
    df_sizes = pd.read_csv(
        f"/ifs/archive/cancer/Lab_RDF/scratch_Lab_RDF/ngs/references/ucsc/{assembly}_chromosome_sizes.txt",
        sep="\t",
        header=0,
    )

    sizes = [
        {"chr": row["chrom"], "size": row["size"]} for i, row in df_sizes.iterrows()
    ]

    return sizes


def chunk_table(
    fin: str, fout: str, fn: Callable[[pd.DataFrame], None], chunksize=200000
):
    first = True
    chunk = 1

    for df in pd.read_csv(
        fin,
        sep="\t",
        header=0,
        keep_default_na=False,
        chunksize=chunksize,
    ):
        print(f"Processing chunk {chunk}...")
        chunk += 1

        fn(df)

        if fout:
            df.to_csv(
                fout,
                sep="\t",
                header=first,
                mode="w" if first else "a",
                index=False,
            )

        first = False
