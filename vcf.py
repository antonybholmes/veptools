import gzip
import os

import gal
import libdna

from reorder import SEP
from veptools import utils

ASSEMBLY = "hg19"
ASSEMBLY_DIR = "/ifs/archive/cancer/Lab_RDF/scratch_Lab_RDF/ngs/dna/hg19"


def chr_to_int(ch):
    ch = ch.replace("chr", "")
    if ch.isdigit():
        return int(ch)
    return {"X": 23, "Y": 24, "MT": 25}.get(ch, 100)


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
        df["chr_order"] = df["Chromosome"].map(chr_to_int)

        df = df.sort_values(by=["Tumor_Sample_Barcode", "chr_order", "Start_Position"])

        # write VCF

        os.makedirs(dir, exist_ok=True)

        print("total", df.shape[0])

        total = 0
        current_sample = None
        f = None

        dna_lookup = {}

        # make header
        header = "##fileformat=VCFv4.2"
        header += '\n##INFO=<ID=ID,Number=1,Type=String,Description="Sample Id">'

        for r in self.sizes:
            header += f"\n##contig=<ID={r['chr']},length={r['size']}>"
        header += "\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"

        for i, (index, row) in enumerate(df.iterrows()):
            sample = row["Tumor_Sample_Barcode"]

            if sample != current_sample:
                print(f"Processing sample {sample} at row {i}, total samples {total}")

                if f is not None:
                    f.close()

                fout = os.path.join(dir, f"{sample}.vcf")
                f = open(fout, "w")

                current_sample = sample
                total += 1

                # header
                if write_header:
                    print(header, file=f)

            chrom = row["Chromosome"]
            pos = row["Start_Position"]
            ref = row["Reference_Allele"]
            alt = row["Tumor_Seq_Allele2"]

            new_ref = ref
            new_alt = alt
            new_pos = pos

            # need to fix for insertions and deletions

            if ref == "-":
                # insertion we need to add base before

                key = f"{chrom}:{pos}"

                if key not in dna_lookup:
                    loc = gal.genomic.Location(chrom, pos, pos)
                    base = self.dna.dna(loc).upper()
                    dna_lookup[key] = base

                base = dna_lookup[key]  # gal.genomic.Location(chrom, pos, pos)
                # base = self.dna.dna(loc).upper()

                # print(loc, base)

                new_ref = base
                new_alt = base + alt
            elif alt == "-":
                # deletion we need to report base before
                new_pos = pos - 1

                key = f"{chrom}:{new_pos}"

                if key not in dna_lookup:
                    loc = gal.genomic.Location(chrom, new_pos, new_pos)
                    base = self.dna.dna(loc).upper()
                    dna_lookup[key] = base

                base = dna_lookup[key]

                new_alt = base
                new_ref = base + ref
            else:
                # snv
                pass

            print(
                f"{chrom}\t{new_pos}\t{sample}\t{new_ref}\t{new_alt}\t.\tPASS\tID={sample}",
                file=f,
            )

            if i % 100000 == 0:
                print(
                    f"Processed {i} rows, current sample {sample}, total samples {total}"
                )

        print(
            f"Finished writing VCF with {total} samples, total variants {df.shape[0]}"
        )

        print("running", total)


def extract_vcf_info_fields(
    vcf_file: str,
) -> tuple[list[str], dict[str, dict[str, str]]]:
    """
    Extract INFO field definitions from a VCF file and return a list of
    INFO field ids and a mapping from id to its definition
    (number, type, description).
    """

    if not os.path.exists(vcf_file):
        raise FileNotFoundError(f"VCF file {vcf_file} does not exist")

    if vcf_file.endswith(".gz"):
        open_func = gzip.open(vcf_file, "rt")
    else:
        open_func = open(vcf_file, "r")

    ret = []

    with open_func as f:
        for line in f:
            if not line.startswith("##"):
                break

            if line.startswith("##INFO="):
                # extract id, number, type, description
                id = None
                number = None
                type = None
                description = None
                tokens = line.strip().split("=", 1)[1].strip("<>").split(",")
                for token in tokens:
                    key, value = token.split("=", 1)
                    key = key.strip()
                    value = value.strip().strip('"')
                    if key == "ID":
                        id = value
                    elif key == "Number":
                        number = value
                    elif key == "Type":
                        type = value
                    elif key == "Description":
                        description = value

                if id is not None:
                    ret.append(
                        {
                            "id": id,
                            "number": number,
                            "type": type,
                            "description": description,
                        }
                    )

    field_map = {field["id"]: field for field in ret}

    return ret, field_map
