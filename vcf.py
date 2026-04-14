import gal
import libdna

from veptools import utils

ASSEMBLY = "hg19"
ASSEMBLY_DIR = "/ifs/archive/cancer/Lab_RDF/scratch_Lab_RDF/ngs/dna/hg19"


class VCFMaker:
    def __init__(self, assembly: str = "hg19"):
        self.dna = libdna.DNA4Bit(ASSEMBLY_DIR)
        self.sizes = utils.chr_sizes(assembly)

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
                id = f"{chrom}_{pos+1 if ref == '-' else pos}_{ref}/{alt}"

                if id in used:
                    # print(f"Duplicate variant {id} at row {i}, skipping")
                    skips += 1
                    continue

                print(
                    f"{chrom}\t{new_pos}\t{id}\t{new_ref}\t{new_alt}\t.\tPASS\tVEP_ID={id}",
                    file=f,
                )

                used.add(id)

        print(
            f"Finished writing VCF with {len(used)} variants, skipped {skips} duplicates"
        )
