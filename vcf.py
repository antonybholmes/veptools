import os

import gal
import libdna

from veptools import utils

ASSEMBLY = "hg19"
ASSEMBLY_DIR = "/ifs/archive/cancer/Lab_RDF/scratch_Lab_RDF/ngs/dna/hg19"


class VCFMaker:
    def __init__(self, assembly: str = "hg19", id_mode: str = "vep_id"):
        self.dna = libdna.DNA4Bit(ASSEMBLY_DIR)
        self.sizes = utils.chr_sizes(assembly)
        self.id_mode = id_mode

    def write_vcf(self, df, fout):
        df = df.sort_values(by=["Chromosome", "Start_Position"])

        # write VCF

        used = set()
        skips = 0

        with open(fout, "w") as f:
            # header
            print("##fileformat=VCFv4.2", file=f)
            print(
                '##INFO=<ID=VEP_ID,Number=1,Type=String,Description="VEP Identifier">',
                file=f,
            )

            for r in self.sizes:
                print(f"##contig=<ID={r['chr']},length={r['size']}>", file=f)

            print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", file=f)

            for i, row in df.iterrows():
                chrom = row["Chromosome"]
                pos = row["Start_Position"]
                # default to same as original
                id = "."
                ref = row["Reference_Allele"]
                alt = row["Tumor_Seq_Allele2"]
                sample = "."

                if "Tumor_Sample_Barcode" in row:
                    sample = row["Tumor_Sample_Barcode"]
                elif "Sample" in row:
                    sample = row["Sample"]
                else:
                    sample = "."

                new_ref = ref
                new_alt = alt
                new_pos = pos

                # need to fix for insertions and deletions

                if ref == "-":
                    # insertion we need to add base before
                    loc = gal.genomic.Location(chrom, pos, pos)
                    base = self.dna.dna(loc).upper()

                    # print(loc, base)

                    new_ref = base
                    new_alt = base + alt

                if alt == "-":
                    # deletion we need to report base before
                    new_pos = pos - 1

                    loc = gal.genomic.Location(chrom, new_pos, new_pos)
                    base = self.dna.dna(loc).upper()

                    new_alt = base
                    new_ref = base + ref

                # make a VEP id to track
                if self.id_mode == "vep_id":
                    id = f"{chrom}_{pos+1 if ref == '-' else pos}_{ref}/{alt}"

                    # in vep_id mode we skip
                    if id in used:
                        # print(f"Duplicate variant {id} at row {i}, skipping")
                        skips += 1
                        continue
                else:
                    id = sample

                print(
                    f"{chrom}\t{new_pos}\t{id}\t{new_ref}\t{new_alt}\t.\tPASS\tVEP_ID={id}",
                    file=f,
                )

                used.add(id)

        print(
            f"Finished writing VCF with {len(used)} variants, skipped {skips} duplicates"
        )

    def split_by_sample(self, df, dir: str = "output", write_header: bool = True):
        df = df.sort_values(by=["Chromosome", "Start_Position"])

        # write VCF

        used = set()
        skips = 0

        os.makedirs(dir, exist_ok=True)

        samples = df["Tumor_Sample_Barcode"].unique()

        print("total", df.shape[0])

        total = 0

        for sample in samples:
            df_sample = df[df["Tumor_Sample_Barcode"] == sample]

            total += df_sample.shape[0]

            fout = os.path.join(dir, f"{sample}.vcf")

            with open(fout, "w") as f:
                # header
                if write_header:
                    print("##fileformat=VCFv4.2", file=f)
                    print(
                        '##INFO=<ID=VEP_ID,Number=1,Type=String,Description="VEP Identifier">',
                        file=f,
                    )

                    for r in self.sizes:
                        print(f"##contig=<ID={r['chr']},length={r['size']}>", file=f)

                    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", file=f)

                for i, row in df_sample.iterrows():
                    chrom = row["Chromosome"]
                    pos = row["Start_Position"]
                    # default to same as original
                    id = "."
                    ref = row["Reference_Allele"]
                    alt = row["Tumor_Seq_Allele2"]
                    sample = "."

                    if "Tumor_Sample_Barcode" in row:
                        sample = row["Tumor_Sample_Barcode"]
                    elif "Sample" in row:
                        sample = row["Sample"]
                    else:
                        sample = "."

                    new_ref = ref
                    new_alt = alt
                    new_pos = pos

                    # need to fix for insertions and deletions

                    if ref == "-":
                        # insertion we need to add base before
                        loc = gal.genomic.Location(chrom, pos, pos)
                        base = self.dna.dna(loc).upper()

                        # print(loc, base)

                        new_ref = base
                        new_alt = base + alt

                    if alt == "-":
                        # deletion we need to report base before
                        new_pos = pos - 1

                        loc = gal.genomic.Location(chrom, new_pos, new_pos)
                        base = self.dna.dna(loc).upper()

                        new_alt = base
                        new_ref = base + ref

                    # make a VEP id to track
                    if self.id_mode == "vep_id":
                        id = f"{chrom}_{pos+1 if ref == '-' else pos}_{ref}/{alt}"

                        # in vep_id mode we skip
                        if id in used:
                            # print(f"Duplicate variant {id} at row {i}, skipping")
                            skips += 1
                            continue
                    else:
                        id = sample

                    print(
                        f"{chrom}\t{new_pos}\t{id}\t{new_ref}\t{new_alt}\t.\tPASS\tID={id}",
                        file=f,
                    )

                    used.add(id)

        print(
            f"Finished writing VCF with {len(used)} variants, skipped {skips} duplicates"
        )

        print("running", total)
