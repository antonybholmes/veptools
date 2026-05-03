import collections

from logging import info
from pathlib import Path
from typing import Union
from urllib.parse import unquote
from .utils import (
    load_hugo,
    get_is_hugo_gene,
    load_transcripts,
    load_ccds_lengths,
    NA,
    SEP,
)

import pandas as pd

VEP_VERSION = "Ensembl_VEP_115"
HUGO_PATH = Path(__file__).parent / "res" / "hugo_approved_20260409.tsv"
MANE_PATH = Path(__file__).parent / "res" / "grch38" / "MANE.GRCh38.v1.5.summary.txt"
SEVERITY_PATH = Path(__file__).parent / "res" / "severity.tsv"


def load_severity() -> dict[str, dict]:
    """
    Load severity map for given assembly which is just a map of severity ids to their descriptions.
    https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html?redirect=no
    """

    df = pd.read_csv(str(SEVERITY_PATH), sep="\t", header=0, keep_default_na=False)

    severity_map = {}

    for _, row in df.iterrows():
        term = row["term"]
        severity = row["severity"]
        impact = row["impact"]

        severity_map[term] = {
            "severity": severity,
            "impact": impact,
        }

    return severity_map


CONSEQUENCE_SEVERITY_MAP = load_severity()

# https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html?redirect=no
# CONSEQUENCE_SEVERITY = {
#     "transcript_ablation": 1,
#     "splice_acceptor_variant": 2,
#     "splice_donor_variant": 3,
#     "stop_gained": 4,
#     "frameshift_variant": 5,
#     "stop_lost": 6,
#     "start_lost": 7,
#     "transcript_amplification": 8,
#     "feature_elongation": 9,
#     "feature_truncation": 10,
#     "inframe_insertion": 11,
#     "inframe_deletion": 11,
#     "missense_variant": 12,
#     "protein_altering_variant": 13,
#     "splice_donor_5th_base_variant": 14,
#     "splice_region_variant": 15,
#     "splice_donor_region_variant": 16,
#     "splice_polypyrimidine_tract_variant": 17,
#     "incomplete_termiNAl_codon_variant": 18,
#     "start_retained_variant": 19,
#     "stop_retained_variant": 20,
#     "synonymous_variant": 21,
#     "coding_sequence_variant": 22,
#     "mature_miRNA_variant": 23,
#     "5_prime_UTR_variant": 24,
#     "3_prime_UTR_variant": 25,
#     "non_coding_transcript_exon_variant": 26,
#     "intron_variant": 27,
#     "NMD_transcript_variant": 28,
#     "non_coding_transcript_variant": 29,
#     "coding_transcript_variant": 30,
#     "upstream_gene_variant": 31,
#     "downstream_gene_variant": 32,
#     "TFBS_ablation": 33,
#     "TFBS_amplification": 34,
#     "TF_binding_site_variant": 35,
#     "regulatory_region_ablation": 36,
#     "regulatory_region_amplification": 37,
#     "regulatory_region_variant": 38,
#     "intergenic_variant": 39,
#     "sequence_variant": 40,
# }

# change 3 letter amino acid NAmes to 1 letter
AA_THREE_TO_ONE_MAP = {
    "ala": "A",
    "cys": "C",
    "asp": "D",
    "glu": "E",
    "phe": "F",
    "gly": "G",
    "his": "H",
    "ile": "I",
    "lys": "K",
    "leu": "L",
    "met": "M",
    "asn": "N",
    "pro": "P",
    "gln": "Q",
    "arg": "R",
    "ser": "S",
    "thr": "T",
    "val": "V",
    "trp": "W",
    "tyr": "Y",
    "ter": "*",
}


# def lookup_gene(
#     id: str,
#     gene_lookup_map: dict[str, str],
#     previous_gene_lookup_map: dict[str, str],
#     alias_gene_lookup_map: dict[str, str],
# ) -> str:
#     id_lower = id.lower()
#     if id_lower in gene_lookup_map:
#         return gene_lookup_map[id_lower]
#     elif id_lower in previous_gene_lookup_map:
#         return previous_gene_lookup_map[id_lower]
#     elif id_lower in alias_gene_lookup_map:
#         return alias_gene_lookup_map[id_lower]
#     else:
#         return id


def extract_csq_header(vcf_file: str) -> tuple[list[str], dict[str, int]]:
    fields = []
    with open(vcf_file) as f:
        for line in f:
            if line.startswith("##INFO=<ID=CSQ"):
                # Extract the Format: part
                start = line.find("Format:") + len("Format:")
                end = line.find('">', start)
                fields = [x.strip() for x in line[start:end].split(SEP)]
                break

    return fields, {field: i for i, field in enumerate(fields)}


def find_vcf_header_line(vcf_file: str) -> int:
    with open(vcf_file) as f:
        for i, line in enumerate(f):
            if line.startswith("#CHROM"):
                return i
    return -1


def format_hgvs(hgvs: str) -> str:
    # VEP likes to url encode certain characters in the HGVS strings, so we need to decode them
    hgvs = unquote(hgvs)

    if ":" in hgvs:
        # remove protein change prefix if present, as we only care about
        # the actual change for annotation purposes
        hgvs = hgvs.split(":")[1]

    return hgvs


def format_hgvsp(hgvsp: str) -> str:
    hgvsp = format_hgvs(hgvsp)

    # replace three letter amino acid codes with one letter codes for consistency and easier parsing later
    hgvsp = hgvsp.lower()

    for three, one in AA_THREE_TO_ONE_MAP.items():
        hgvsp = hgvsp.replace(three, one)

    return hgvsp


def get_highest_severity(consequences: list[str]) -> int:
    # default non existent low priority consequence to 1000 so that it gets sorted to the end
    consequence = {"severity": 1000, "impact": NA}

    # try to assign the most severe consequence based on the predefined severity ranking,
    # if no consequences match, will stay at 1000
    for c in consequences:
        c = c.strip()
        if c in CONSEQUENCE_SEVERITY_MAP:
            s = CONSEQUENCE_SEVERITY_MAP[c]["severity"]
            if s < consequence["severity"]:
                consequence = CONSEQUENCE_SEVERITY_MAP[c]

    return consequence


def parse_csq_with_severity(
    csq: str,
    header_fields: list[str],
) -> list[dict]:
    """
    Parse a VEP CSQ field and return canonical transcripts sorted by
    if canonical and by severity etc.

    Parameters:
        csq_string (str): The CSQ field from VEP VCF (comma-separated transcripts)
        header_fields (list of str): Column names from CSQ header

    Returns:
        tuple of lists of dicts: Each dict contains transcript info + 'severity_rank', sorted by severity
    """
    transcripts = csq.split(",")
    results = []

    # Map field names to their index
    field_index = {name: i for i, name in enumerate(header_fields)}

    for t in transcripts:
        fields = t.split(SEP)
        # Only consider canonical transcripts

        # is_canonical = False
        # is_protein_coding = False
        # has_protein_change = False

        canonical_flag = fields[field_index.get("CANONICAL", -1)]

        is_canonical = canonical_flag == "YES"

        biotype = fields[field_index.get("BIOTYPE", -1)]
        is_protein_coding = "protein_coding" in biotype
        hgvsp = fields[field_index.get("HGVSp", -1)]
        hgvsc = fields[field_index.get("HGVSc", -1)]

        hgvsc = format_hgvs(hgvsc)
        hgvsp = format_hgvsp(hgvsp)

        has_protein_change = hgvsp.startswith("p.")  # hgvsp != ""
        is_nonsense = has_protein_change and ("Ter" in hgvsp or "*" in hgvsp)

        consequences = fields[field_index.get("Consequence", 1)].split("&")

        severity = get_highest_severity(consequences)

        transcript_info = {name: fields[i] for i, name in enumerate(header_fields)}

        # if hgvsp != "":
        #    print(hgvsp)

        transcript_info["biotype"] = biotype
        transcript_info["hgvsc"] = hgvsc
        transcript_info["hgvsp"] = hgvsp
        transcript_info["is_canonical"] = is_canonical
        transcript_info["is_protein_coding"] = is_protein_coding
        transcript_info["has_protein_change"] = has_protein_change
        transcript_info["is_nonsense"] = is_nonsense
        transcript_info["severity_rank"] = severity["severity"]
        transcript_info["severity_impact"] = severity["impact"]

        transcript_info["gene_id"] = transcript_info.get("Gene", NA)
        transcript_info["gene_symbol"] = transcript_info.get("SYMBOL", NA)

        results.append(transcript_info)

    # Sort transcripts by should have hugo symbol, then protein coding, then canonical, then severity
    # results.sort(
    #     key=lambda x: (
    #         not x["is_hugo_gene"],
    #         not x["is_protein_coding"],
    #         not x["is_canonical"],
    #         not x["has_protein_change"],
    #         x["severity_rank"],
    #     )
    # )

    # is_primary = len(results) > 0 and results[0]["is_hugo_gene"]
    # primary = results[0] if is_primary else []
    # secondary = results[1:] if is_primary else results

    return results


def extract_csq(info_field: str) -> str:
    """
    Extract the CSQ field from a VCF INFO column.
    Returns empty string if not present.
    """
    if "CSQ=" not in info_field:
        return ""

    csq_part = info_field.split("CSQ=")[1]
    # Stop at the next semicolon if present
    csq_part = csq_part.split(";")[0]
    return csq_part


def make_vep_id(chrom: str, start: int, ref: str, alt: str) -> str:
    vep_id = f"{chrom}_{start + 1 if ref == '-' else start}_{ref}/{alt}"
    return vep_id


def blank_val(v: Union[str, bool, int]) -> Union[str, int]:
    # check if str is int, if so return int, otherwise return NA if blank, otherwise return original string

    if isinstance(v, bool):
        return int(v)

    if isinstance(v, int):
        # treat -1 as missing
        if v == -1:
            return NA

        return v

    v = str(v).strip()

    if v.isdigit():
        return int(v)

    if v == "" or v == NA or v == "-1":
        return NA

    return v


def load_gene_lookup_maps() -> tuple[dict[str, str], dict[str, str], dict[str, str]]:
    gene_lookup_map, previous_gene_lookup_map, alias_gene_lookup_map = load_hugo(
        str(HUGO_PATH)
    )

    return gene_lookup_map, previous_gene_lookup_map, alias_gene_lookup_map


def load_transcript_map(assembly: str = "hg19") -> dict[str, dict]:
    """
    Load CCDS map for given assembly which is just a map of symbols to CCDS ids.
    """
    if assembly == "hg38":
        path = (
            Path(__file__).parent
            / "res"
            / assembly
            / "gencode_v48_basic_transcripts.tsv"
        )
    else:
        path = (
            Path(__file__).parent
            / "res"
            / assembly
            / "gencode_v48lift37_basic_transcripts.tsv"
        )

    print(f"Loading CCDS map from {path}...")

    transcript_map = {}

    load_transcripts(str(path), transcript_map)

    print("x", transcript_map.get("ENST00000445750"))

    # augment with v19 transcripts for hg19, which have some ccds annotations missing from v48lift37
    if assembly == "hg19":
        path = (
            Path(__file__).parent
            / "res"
            / assembly
            / "gencode_transcripts_v19_grch37.tsv"
        )

        print(f"Loading v19 CCDS map from {path}...")

        load_transcripts(str(path), transcript_map)

        path = (
            Path(__file__).parent
            / "res"
            / assembly
            / "gencode_transcripts_v10_grch37.tsv"
        )

        print(f"Loading v10 CCDS map from {path}...")

        # v10_transcript_map = {}
        load_transcripts(str(path), transcript_map)

        # transcript_map = transcript_map | v19_transcript_map | v10_transcript_map

    return transcript_map


# def load_v19_transcript_map(assembly: str = "hg19") -> dict[str, dict]:
#     """
#     Load CCDS map for given assembly which is just a map of symbols to CCDS ids.
#     """

#     path = (
#         Path(__file__).parent / "res" / assembly / "gencode_transcripts_v19_grch37.tsv"
#     )

#     print(f"Loading CCDS map from {path}...")

#     transcript_map = load_transcripts(str(path))

#     return transcript_map


def load_ccds_length_map(assembly: str = "hg19") -> dict[str, dict]:
    """
    Load CCDS length map for given assembly which is just a map of CCDS ids to their amino acid lengths.
    """
    path = Path(__file__).parent / "res" / "CCDS_aa_lengths.tsv"

    print(f"Loading CCDS length map from {path}...")

    ccds_length_map = load_ccds_lengths(str(path))

    if assembly == "hg19":
        # add multiple versions ccds lengths as well, which are missing some ccds that are in v48lift37
        # in the hope we can match old and withdrawn ccds ids that are missing from the main map.

        path = (
            Path(__file__).parent / "res" / assembly / "CCDS_aa_lengths.v15.grch37.tsv"
        )

        print(f"Loading hg19 CCDS length map from {path}...")

        v19_ccds_length_map = load_ccds_lengths(str(path))

        path = (
            Path(__file__).parent
            / "res"
            / assembly
            / "CCDS_aa_lengths.vhs37-1.grch37.tsv"
        )

        print(f"Loading hg19 CCDS length map from {path}...")

        vhs371_ccds_length_map = load_ccds_lengths(str(path))

        path = (
            Path(__file__).parent
            / "res"
            / assembly
            / "CCDS_aa_lengths.vhs37-3.grch37.tsv"
        )

        print(f"Loading hg19 CCDS length map from {path}...")

        vhs373_ccds_length_map = load_ccds_lengths(str(path))

        # add this info to the main ccds length map
        ccds_length_map = (
            ccds_length_map
            | v19_ccds_length_map
            | vhs373_ccds_length_map
            | vhs371_ccds_length_map
        )

    return ccds_length_map


def load_mane_map():
    df = pd.read_csv(MANE_PATH, sep="\t", header=0, keep_default_na=False)

    mane_map = {}

    for _, row in df.iterrows():
        ensembl_id = row["Ensembl_Gene"].split(".")[0]
        symbol = row["symbol"]
        refseq = row["RefSeq_nuc"].split(".")[0]
        transcript = row["Ensembl_nuc"].split(".")[0]
        status = row["MANE_status"]

        data = {
            "ensembl": ensembl_id,
            "gene_symbol": symbol,
            "refseq": refseq,
            "transcript": transcript,
            "status": status,
        }

        mane_map[ensembl_id] = data
        mane_map[refseq] = data
        mane_map[transcript] = data
        mane_map[symbol] = data

    return mane_map


class VEPAnnotation:
    def __init__(
        self,
        vcf_file: str,
        assembly: str = "hg19",
    ):
        self.vcf_file = vcf_file
        self.assembly = assembly

        self.header_line = find_vcf_header_line(vcf_file)
        self.csq_fields, self.csq_field_map = extract_csq_header(vcf_file)
        self.hugo_info = load_hugo(str(HUGO_PATH))

        # self.transcript_v19_map = load_v19_transcript_map()
        self.transcript_map = load_transcript_map(assembly)
        self.ccds_length_map = load_ccds_length_map()
        self.mane_map = load_mane_map()

        self.annotation_map = None

    def parse_vcf(self) -> dict[str, list[dict]]:
        if self.annotation_map is not None:
            return self.annotation_map

        self.annotation_map = collections.defaultdict(
            lambda: collections.defaultdict(list)
        )

        chunk = 1
        # chunk vcf and build annotation map
        for df_vep in pd.read_csv(
            self.vcf_file,
            sep="\t",
            header=self.header_line,
            keep_default_na=False,
            chunksize=200000,
        ):
            print(f"Processing VEP chunk {chunk}...")
            chunk += 1
            for _, row in df_vep.iterrows():
                vep_id = row["ID"]

                csq = extract_csq(row["INFO"])

                if csq == "":
                    continue

                transcripts = parse_csq_with_severity(
                    csq,
                    self.csq_fields,
                )

                # add extra annotation for sorting
                for transcript in transcripts:
                    transcript_id = transcript.get("Feature", NA)

                    # add hugo gene symbol and is_hugo_gene flag for sorting
                    hugo_info = get_is_hugo_gene(
                        transcript["gene_id"],
                        transcript["gene_symbol"],
                        self.hugo_info,
                    )

                    # symbol will be either the original symbol or the hugo symbol if found,
                    # but we keep the original symbol in the gene_symbol field for reference
                    transcript["gene_symbol"] = hugo_info["symbol"]
                    transcript["hgnc_id"] = hugo_info["hgnc_id"]
                    transcript["is_hugo_gene"] = hugo_info["is_hugo_gene"]

                    gencode_info = self.transcript_map.get(transcript_id, {})

                    transcript["exons"] = gencode_info.get("exons", NA)

                    gencode_is_canonical = gencode_info.get("is_canonical", 0)

                    # use VEP or GENCODE decide if canonical
                    if gencode_is_canonical:
                        transcript["is_canonical"] = 1

                    ccds = gencode_info.get("ccds", NA)

                    # if self.assembly == "hg19":
                    #     v19_info = self.transcript_v19_map.get(transcript_id, {})
                    #     if ccds == NA:
                    #         ccds = v19_info.get("ccds", ccds)

                    transcript["ccds"] = ccds

                    ccds_info = self.ccds_length_map.get(ccds, {})

                    # ccds must always be numeric for sorting
                    ccds_aa_length = ccds_info.get("aa_length", -1)

                    # we don't use na here as we need to sort by this so
                    # cannot intermix strings and ints
                    transcript["ccds_aa_length"] = ccds_aa_length

                    # if transcript["ccds"] != NA:
                    #     print(
                    #         f"Transcript {transcript_id} has CCDS {transcript['ccds']} with AA length {transcript['ccds_aa_length']}"
                    #     )
                    #     exit(0)

                    transcript["has_ccds"] = int(transcript["ccds"] != NA)
                    transcript["has_mane"] = 0
                    transcript["mane_refseq"] = NA
                    transcript["mane_status"] = NA
                    mane_info = self.mane_map.get(transcript_id, None)
                    if mane_info:
                        transcript["has_mane"] = 1
                        transcript["mane_refseq"] = mane_info["refseq"]
                        transcript["mane_status"] = mane_info["status"]

                # sort by which we think might be priority, we want
                # mutations in a ccds to be more likely to be considered primary,
                # and if multiple ccds, then the longest one, as a proxy for most complete annotation.
                # Also want to prioritize transcripts with protein changes and with hugo symbols,
                # as these are more likely to be relevant for our analysis
                transcripts.sort(
                    key=lambda t: (
                        # if has hugo symbol, more likely to be primary
                        not t["is_hugo_gene"],
                        # if protein coding, more likely to be primary
                        not t["is_protein_coding"],
                        # if has ccds, more likely to be primary
                        not t["has_ccds"],
                        # if has mane, more likely to be primary
                        not t["has_mane"],
                        # if canonical, more likely to be primary
                        not t["is_canonical"],
                        # if nothing else try longest CCDS length, as a proxy for most complete transcript
                        -t["ccds_aa_length"],
                        # protein changes should come first
                        not t["has_protein_change"],
                        # nonsense mutations should come first
                        not t["is_nonsense"],
                        # pick the one with most severe consequence
                        t["severity_rank"],
                    )
                )

                for ti, transcript in enumerate(transcripts):
                    transcript_id = transcript.get("Feature", NA)

                    vep_exon_info = transcript.get("EXON", NA)
                    exon_num = NA
                    exon_total = transcript["exons"]
                    consequences = SEP.join(
                        transcript.get("Consequence", NA).split("&")
                    )

                    if vep_exon_info != NA and "/" in vep_exon_info:
                        # become strings but will be
                        # turned into nums by blank_val
                        exon_num, exon_total = vep_exon_info.split("/")

                    annotation = {
                        "gene_id": blank_val(transcript["gene_id"]),
                        "gene_symbol": blank_val(transcript["gene_symbol"]),
                        "hgnc_id": blank_val(transcript["hgnc_id"]),
                        "is_hugo_gene": blank_val(transcript["is_hugo_gene"]),
                        # "hugo_gene_symbol": blank_val(hugo_gene_symbol),
                        "is_canonical": blank_val(transcript["is_canonical"]),
                        "is_protein_coding": blank_val(transcript["is_protein_coding"]),
                        "has_protein_change": blank_val(
                            transcript["has_protein_change"]
                        ),
                        "is_nonsense": blank_val(transcript["is_nonsense"]),
                        "transcript_id": blank_val(transcript_id),
                        "exon": blank_val(exon_num),
                        "exons": blank_val(exon_total),
                        "biotype": blank_val(transcript["biotype"]),
                        "consequence": blank_val(consequences),
                        "severity": blank_val(transcript["severity_rank"]),
                        "hgvsp": blank_val(transcript["hgvsp"]),
                        "hgvsc": blank_val(transcript["hgvsc"]),
                        "ccds": blank_val(transcript["ccds"]),
                        "ccds_aa_length": blank_val(transcript["ccds_aa_length"]),
                        "mane_refseq": blank_val(transcript["mane_refseq"]),
                        "mane_status": blank_val(transcript["mane_status"]),
                    }

                    # if transcript_id == "ENST00000332831":
                    #     print(
                    #         f"Found transcript ENST00000332831 with annotation: {annotation}"
                    #     )

                    # if is_protein_coding:
                    #     # matcher = re.search(r":(c\..+)", info[10])
                    #     hgvscs = transcript.get("HGVSc", NA)
                    #     hgvsps = transcript.get("HGVSp", NA)

                    #     hgvscs = decode_hgvs(hgvscs)
                    #     hgvsps = decode_hgvs(hgvsps)

                    #     if ":" in hgvscs:
                    #         hgvscs = hgvscs.split(":")[1]

                    #     # matcher = re.search(r":(p\..+)", info[11])

                    #     # print(f"Raw HGVSp: {hgvsp}")

                    #     if ":" in hgvsps:
                    #         hgvsps = hgvsps.split(":")[1]

                    #     protein_lc = hgvsps.lower()

                    #     for three, one in AA_THREE_TO_ONE_MAP.items():
                    #         protein_lc = protein_lc.replace(three, one)

                    #     annotation["hgvsp"] = blank_val(protein_lc)
                    #     annotation["hgvsc"] = blank_val(hgvscs)

                    # mode = 0 if ti == 0 and is_hugo_gene and is_protein_coding else 1
                    # first item wins primary annotation
                    primary_mode = 0 if ti == 0 else 1

                    self.annotation_map[vep_id][primary_mode].append(annotation)
            # break

        print("Finished parsing VCF and building annotation map.")
        return self.annotation_map

    def annotate_maf(self, fin: str, fout: str, chunksize: int = 200000):
        # annotate maf
        self.parse_vcf()

        # add new columns for VEP annotation

        print("Adding VEP annotations to MAF...")

        first = True
        chunk = 1
        for df in pd.read_csv(
            fin,
            sep="\t",
            header=0,
            keep_default_na=False,
            chunksize=chunksize,
        ):
            print(f"Processing VEP chunk {chunk}...")
            chunk += 1

            # delete col
            # df = df.drop(columns=["VEP_HGVSp"])
            if first:
                df.rename(
                    columns={"Hugo_Symbol": "Hugo_Symbol (Rahul - not really Hugo)"},
                    inplace=True,
                )

            df["VEP_Gene_ID"] = NA
            # df["VEP_Hugo_Gene_Symbol"] = NA
            df["VEP_Gene_Symbol"] = NA
            df["VEP_Is_Hugo_Gene"] = NA
            df["VEP_Biotype"] = NA
            df["VEP_HGVSp"] = NA
            df["VEP_HGVSc"] = NA
            df["VEP_Variant_Classification"] = NA
            df["VEP_Variant_Severity"] = NA
            df["VEP_Transcript"] = NA
            df["VEP_Exon"] = NA
            df["VEP_Total_Exons"] = NA
            df["VEP_Is_Canonical"] = NA
            df["CCDS"] = NA
            df["CCDS_AA_Length"] = NA
            df["MANE_RefSeq"] = NA
            df["MANE_status"] = NA

            df["VEP_Secondary_Gene_ID"] = NA
            df["VEP_Secondary_Gene_Symbol"] = NA
            # df["VEP_Secondary_Hugo_Gene_Symbol"] = NA
            df["VEP_Secondary_Biotype"] = NA
            df["VEP_Secondary_HGVSp"] = NA
            df["VEP_Secondary_HGVSc"] = NA
            df["VEP_Secondary_Variant_Classification"] = NA
            df["VEP_Secondary_Variant_Severity"] = NA
            df["VEP_Secondary_Transcript"] = NA
            df["VEP_Secondary_Exon"] = NA
            df["VEP_Secondary_Total_Exons"] = NA
            df["VEP_Secondary_Canonical"] = NA
            df["Secondary_CCDS"] = NA
            df["Secondary_CCDS_AA_Length"] = NA
            df["VEP_Annotation_Database"] = VEP_VERSION

            for i, row in df.iterrows():
                chrom = row["Chromosome"]
                start = row["Start_Position"]
                ref = row["Reference_Allele"]
                alt = row["Tumor_Seq_Allele2"]

                vep_id = make_vep_id(chrom, start, ref, alt)

                annotations = self.annotation_map[vep_id][0]

                # if primary annotation exists, add it, otherwise add secondary annotations concateNAted with SEP
                if len(annotations) > 0:
                    annotation = annotations[0]
                    transcript_id = annotation["transcript_id"]

                    df.at[i, "VEP_Gene_ID"] = annotation["gene_id"]
                    df.at[i, "VEP_Gene_Symbol"] = annotation["gene_symbol"]
                    # df.at[i, "VEP_Hugo_Gene_Symbol"] = annotation["hugo_gene_symbol"]

                    df.at[i, "VEP_Is_Hugo_Gene"] = annotation["is_hugo_gene"]

                    df.at[i, "VEP_Biotype"] = annotation["biotype"]
                    df.at[i, "VEP_HGVSp"] = annotation["hgvsp"]
                    df.at[i, "VEP_HGVSc"] = annotation["hgvsc"]
                    df.at[i, "VEP_Variant_Classification"] = annotation["consequence"]
                    df.at[i, "VEP_Variant_Severity"] = annotation["severity"]
                    df.at[i, "VEP_Transcript"] = transcript_id
                    df.at[i, "VEP_Exon"] = annotation["exon"]
                    df.at[i, "VEP_Total_Exons"] = annotation["exons"]
                    df.at[i, "VEP_Is_Canonical"] = int(annotation["is_canonical"])
                    df.at[i, "CCDS"] = annotation["ccds"]
                    df.at[i, "CCDS_AA_Length"] = annotation["ccds_aa_length"]
                    mane_info = self.mane_map.get(transcript_id, None)
                    if mane_info:
                        df.at[i, "MANE_RefSeq"] = mane_info["refseq"]
                        df.at[i, "MANE_status"] = mane_info["status"]

                #
                # secondary annotations
                #

                annotations = self.annotation_map[vep_id][1]

                if len(annotations) > 0:
                    hgvsps = SEP.join([a["hgvsp"] for a in annotations])
                    hgvscs = SEP.join([a["hgvsc"] for a in annotations])

                    consequences = SEP.join([a["consequence"] for a in annotations])

                    severities = SEP.join([str(a["severity"]) for a in annotations])

                    gene_ids = SEP.join([a["gene_id"] for a in annotations])
                    gene_symbols = SEP.join([a["gene_symbol"] for a in annotations])
                    # hugo_gene_symbols = SEP.join(
                    #    [a["hugo_gene_symbol"] for a in annotations]
                    # )
                    transcript_ids = SEP.join([a["transcript_id"] for a in annotations])
                    ccdss = [
                        self.transcript_map.get(a["transcript_id"], {}).get("ccds", NA)
                        for a in annotations
                    ]

                    exons = SEP.join([str(a["exon"]) for a in annotations])
                    exons_total = SEP.join([str(a["exons"]) for a in annotations])
                    is_canonical = SEP.join(
                        [str(a["is_canonical"]) for a in annotations]
                    )

                    biotypes = SEP.join([a["biotype"] for a in annotations])

                    df.at[i, "VEP_Secondary_Gene_ID"] = gene_ids
                    df.at[i, "VEP_Secondary_Gene_Symbol"] = gene_symbols
                    # df.at[i, "VEP_Secondary_Hugo_Gene_Symbol"] = hugo_gene_symbols
                    df.at[i, "VEP_Secondary_Biotype"] = biotypes
                    df.at[i, "VEP_Secondary_HGVSp"] = hgvsps
                    df.at[i, "VEP_Secondary_HGVSc"] = hgvscs
                    df.at[i, "VEP_Secondary_Variant_Classification"] = consequences
                    df.at[i, "VEP_Secondary_Variant_Severity"] = severities
                    df.at[i, "VEP_Secondary_Transcript"] = transcript_ids
                    df.at[i, "VEP_Secondary_Exon"] = exons
                    df.at[i, "VEP_Secondary_Total_Exons"] = exons_total
                    df.at[i, "VEP_Secondary_Canonical"] = is_canonical
                    df.at[i, "Secondary_CCDS"] = SEP.join(ccdss)
                    df.at[i, "Secondary_CCDS_AA_Length"] = SEP.join(
                        [
                            str(
                                blank_val(
                                    self.ccds_length_map.get(ccds, {}).get(
                                        "aa_length", NA
                                    )
                                )
                            )
                            for ccds in ccdss
                        ]
                    )

            df.to_csv(
                fout,
                sep="\t",
                mode="w" if first else "a",
                header=first,
                index=False,
            )

            first = False

        print("Finished adding VEP annotations to MAF.")


# map ensembl to gene symbol
# df_genes = pd.read_csv(
#     "/ifs/archive/cancer/Lab_RDF/scratch_Lab_RDF/ngs/references/gencode/grch37/gencode_v48lift37_basic_genes.tsv",
#     sep="\t",
#     header=0,
#     keep_default_na=False,
# )


# gene_lookup_map = {row["gene_id"]: row["gene_symbol"] for _, row in df_genes.iterrows()}
# gene_lookup_map = gene_lookup_map | (
#     {row["gene_symbol"]: row["gene_symbol"] for _, row in df_genes.iterrows()}
# )

# df_hugo = df_genes[df_genes["hgnc_id"] != ""]
# hugo_genes = set(df_hugo["Approved symbol"]) | set(df_hugo["HGNC ID"])


# VEP likes to url encode certain characters in the HGVS strings
# def decode_hgvs(encoded_str: str) -> str:
#    return unquote(encoded_str)
