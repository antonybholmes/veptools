# -*- coding: utf-8 -*-
"""
Generate a tss distribution for a region file

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sqlite3

import numpy as np
import pandas as pd

BASES_FROM_CDS_STILL_CONSIDERED_SPLICE = 2
SEP = "|"

SPLICE_PROTEIN_CHANGE_COL = "Splice_Protein_Change (AH)"
SPLICE_TRANSCRIPT_COL = "Splice_Transcript_Type (1=exact, 2=canonical, 3=longest)"

STRICT_SPLICE_COL = f"Is_Strict_Splice_+/-_{BASES_FROM_CDS_STILL_CONSIDERED_SPLICE}bp"


FIND_MATCHING_TRANSCRIPT = """SELECT DISTINCT
    g.id as gid,
    g.gene_id,
    g.gene_symbol,
    g.chr,
    g.strand,
    t.id as tid,
    t.transcript_id
    FROM genes as g
    JOIN transcripts AS t ON t.gene_id = g.id
    WHERE
        t.transcript_id = :transcript
    """

FIND_CANONICAL_TRANSCRIPTS = """SELECT DISTINCT
    g.id as gid,
    g.gene_id,
    g.gene_symbol,
    g.chr,
    g.strand,
    t.id as tid,
    t.transcript_id
    FROM genes as g
    JOIN transcripts AS t ON t.gene_id = g.id
    WHERE
        g.chr = :chr AND
        t.start <= :position AND
        t.end >= :position AND
        t.is_canonical = 1;
    """


FIND_LONGEST_TRANSCRIPTS = """SELECT DISTINCT
    g.id as gid,
    g.gene_id,
    g.gene_symbol,
    g.chr,
    g.strand,
    t.id as tid,
    t.transcript_id,
    t.start,
    t.end,
    t.end - t.start AS length
    FROM genes as g
    JOIN transcripts AS t ON t.gene_id = g.id
    WHERE
        g.chr = :chr AND   
        t.start <= :position AND
        t.end >= :position
    ORDER BY length DESC
    LIMIT 1;
    """

# distance to cloeset exon
FIND_CLOSEST_EXON_SQL = """SELECT DISTINCT
    e.id, 
    e.start,
    e.end,
    e.exon_number,
    :start <= e.end AND :end >= e.start as within_exon,
    MIN(ABS(e.start - :start), ABS(e.end - :start), ABS(e.start - :end), ABS(e.end - :end)) AS min_abs_dist,
    (cds.id IS NOT NULL) AS has_cds
    FROM exons AS e
    JOIN transcripts AS t ON t.id = e.transcript_id
    LEFT JOIN cds AS cds ON cds.exon_id = e.id
    WHERE t.id = :transcript_id
    ORDER BY min_abs_dist ASC
    LIMIT 1;    
"""

FIND_CLOSEST_CDS_SQL = """SELECT DISTINCT
    cds.id,
    cds.start,
    cds.end,
    e.id as exon_id,
    e.exon_number,
    cds.offset,
    :start >= cds.start AND :end <= cds.end as within_cds,
    MIN(ABS(cds.start - :start), ABS(cds.end - :start), ABS(cds.start - :end), ABS(cds.end - :end)) AS min_abs_dist
    FROM cds
    JOIN exons AS e ON e.id = cds.exon_id
    JOIN transcripts AS t ON t.id = e.transcript_id
    WHERE t.id = :transcript_id
    ORDER BY min_abs_dist ASC
    LIMIT 1;
"""


# FIND_CLOSEST_CDS_POS_SQL = """SELECT DISTINCT
#     cds.id,
#     cds.start,
#     cds.end,
#     cds.exon_number,
#     cds.offset,
#     :position >= cds.start AND :position <= cds.end as within_cds
#     FROM cds
#     JOIN transcripts AS t ON t.id = cds.transcript_id
#     WHERE
#         (cds.end <= :position OR (:position >= cds.start AND :position <= cds.end)) AND
#         t.id = :transcript_id
#     ORDER BY cds.exon_number DESC
#     LIMIT 1;
# """

# FIND_CLOSEST_CDS_NEG_SQL = """SELECT DISTINCT
#     cds.id,
#     cds.start,
#     cds.end,
#     cds.exon_number,
#     cds.offset,
#     :position >= cds.start AND :position <= cds.end as within_cds
#     FROM cds
#     JOIN transcripts AS t ON t.id = cds.transcript_id
#     WHERE
#         (cds.start >= :position OR (:position >= cds.start AND :position <= cds.end)) AND
#         t.id = :transcript_id
#     ORDER BY cds.exon_number DESC
#     LIMIT 1;
# """

# CDS_SQL = """SELECT DISTINCT
#     cds.id,
#     cds.exon_id,
#     cds.start,
#     cds.end,
#     cds.exon_number,
#     cds.offset
#     FROM cds as cds
#     JOIN transcripts AS t ON t.id = cds.transcript_id
#     JOIN genes AS g ON g.id = t.gene_id
#     WHERE
#         cds.exon_number = :exon_number AND
#         t.id = :transcript_id
#     ORDER BY cds.exon_number;
# """

# CDS_BEFORE_SQL = """SELECT DISTINCT
#     cds.id,
#     cds.exon_id,
#     cds.exon_number,
#     cds.offset
#     FROM cds as cds
#     JOIN transcripts AS t ON t.id = cds.transcript_id
#     JOIN genes AS g ON g.id = t.gene_id
#     WHERE
#         cds.exon_number <= :exon_number AND
#         t.id = :transcript_id
#     ORDER BY cds.exon_number;
# """


DB = "/home/antony/development/ngs/references/gencode/grch37/gencode.v49lift37.basic.grch37.db"


class SpliceAnnotator:
    def __init__(self, db_path=DB):
        self.conn = sqlite3.connect(db_path)
        self.conn.row_factory = sqlite3.Row  # access columns by name
        self.cursor = self.conn.cursor()

        self.cursor.execute("SELECT id, name FROM chromosomes;")
        self.chr_map = {row[1]: row[0] for row in self.cursor.fetchall()}

        self.cursor.execute("SELECT id, species, genome, version FROM metadata;")
        self.metadata = self.cursor.fetchone()

    def close(self):
        self.conn.close()

    def _process_splice_sites(
        self, df, index, start, end, transcripts, transcript_type
    ):

        annotations = []

        for transcript in transcripts:
            # print("transcript:", transcript, start)

            # sum lengths
            cds_pos = 0
            nearest_aa = -1
            protein_change = "5'UTR"

            # use the gene symbol we found in the db for
            # reporting rather than the input one which can be mangled
            # because of Excel etc
            gene_symbol = transcript["gene_symbol"]

            strand = transcript["strand"]
            transcript_db_id = transcript["tid"]
            transcript_id = transcript["transcript_id"]

            # sql = FIND_CLOSEST_CDS_POS_SQL if strand == 1 else FIND_CLOSEST_CDS_NEG_SQL

            self.cursor.execute(
                FIND_CLOSEST_CDS_SQL,
                {
                    "start": start,
                    "end": end,
                    "transcript_id": transcript_db_id,
                },
            )

            cds = self.cursor.fetchone()

            if cds is not None:
                cds_start = cds["start"]
                cds_end = cds["end"]
                # cds_exon_number = cds["exon_number"]
                within_cds = cds["within_cds"]

                # number of bases from cds 1 to this cds start
                # so we can calculate cds length without summing all previous cds
                # for cds 1 this can be -2, -1, or 0 depending on splice site position
                # Offset has phase of CDS 1 built in
                offset = cds["offset"]
                pos_in_cds = 0

                if within_cds == 1:
                    # we are within the CDS so we can just
                    # calculate the amino acid position directly
                    # this is zero based
                    pos_in_cds = (start - cds_start) if strand == 1 else (cds_end - end)
                else:
                    # we are between CDS regions so we add the full length of this CDS
                    if (strand == 1 and start > cds_end) or (
                        strand == -1 and end < cds_start
                    ):
                        # we only add the CDS length if our test site is after (positive strand)
                        # or before (negative strand) the CDS region
                        pos_in_cds = cds_end - cds_start  # + 1

                # total cds length must include offset  in CDS 1, which can be -2, -1, or 0 depending on splice site position
                # since a phase of 1 or 2 will shorten the total CDS length by 1 or 2 bases respectively,
                # which will affect the nearest amino acid calculation, we need to add the offset to the total CDS length
                cds_pos = pos_in_cds + offset

                # for cds 1 this can be -2, -1, or 0 depending on splice site position
                # so the length must be fixed to be at least 0 for this special case
                cds_pos = max(0, cds_pos)

                nearest_aa = cds_pos // 3 + 1

                # (total_cds_length + 2) // 3

                protein_change = f"p.X{nearest_aa}_splice"
            else:
                print(
                    f"No CDS found for transcript {transcript_id} for gene {gene_symbol} at position {start}-{end}."
                )

            #
            # Lets determine closest exon to splice site and test if we
            # are within a strict splice site distance to it

            # print(
            #     "Finding closest exon for splice site at position:",
            #     start,
            #     transcript_db_id,
            #     FIND_CLOSEST_EXON_SQL,
            # )

            self.cursor.execute(
                FIND_CLOSEST_EXON_SQL,
                {
                    "start": start,
                    "end": end,
                    "transcript_id": transcript_db_id,
                },
            )

            closest_exon = self.cursor.fetchone()

            if closest_exon is None:
                print(
                    f"No closest exon found for gene {gene_symbol} at position {start}-{end} splice"
                )
                continue

            min_abs_dist = closest_exon["min_abs_dist"]
            is_strict_splice = min_abs_dist <= BASES_FROM_CDS_STILL_CONSIDERED_SPLICE

            # if the closest exon does not have a CDS then we need to
            # think it is a 5'UTR mutation that has no effect on protein.
            # If there is no closest cds then we cannot annotate the protein change
            # and if the closest exon has no cds then it is likely a UTR exon
            # and we should not annotate a protein change either
            if cds is None or closest_exon["id"] != cds["exon_id"]:
                # closest exon has no CDS so likely 5' UTR
                protein_change = "5'UTR"
                nearest_aa = -1

            annotations.append(
                {
                    "gene_symbol": gene_symbol,
                    "protein_change": protein_change,
                    "aa_position": nearest_aa,
                    "nearest_exon": closest_exon["exon_number"],
                    "nearest_exon_dist": min_abs_dist,
                    "transcript_id": transcript_id,
                    "transcript_type": transcript_type,
                    "is_strict_splice": 1 if is_strict_splice else 0,
                }
            )

        if len(annotations) == 0:
            # cds starts after splice site so likely 5' UTR
            annotations = [
                {
                    "gene_symbol": gene_symbol,
                    "protein_change": "5'UTR",
                    "aa_position": "",
                    "nearest_exon": -1,
                    "nearest_exon_dist": -1,
                    "transcript_id": "",
                    "transcript_type": "",
                    "is_strict_splice": 0,
                }
            ]

        # df.at[index, "Gene_Symbol"] = SEP.join(
        #     sorted(set([a["gene_symbol"] for a in annotations]))
        # )

        # remove blank protein changes
        protein_changes = [
            a["protein_change"] for a in annotations if a["protein_change"] != ""
        ]

        unique_protein_changes = set(protein_changes)

        # print("changes:", protein_changes)

        # special case where we only have 1 annotation
        # so we will not repeat it in the splice change column
        if len(unique_protein_changes) < 2:
            protein_changes = list(sorted(unique_protein_changes))

        df.at[index, SPLICE_PROTEIN_CHANGE_COL] = SEP.join(protein_changes)

        df.at[index, "Splice_Transcript_ID"] = SEP.join(
            [a["transcript_id"] for a in annotations]
        )
        df.at[index, SPLICE_TRANSCRIPT_COL] = SEP.join(
            [str(a["transcript_type"]) for a in annotations]
        )
        df.at[index, "Splice_Nearest_Exon_Number"] = SEP.join(
            [str(a["nearest_exon"]) for a in annotations]
        )
        df.at[index, "Splice_Nearest_Exon_Dist_bp"] = SEP.join(
            [str(a["nearest_exon_dist"]) for a in annotations]
        )

        df.at[index, STRICT_SPLICE_COL] = SEP.join(
            [str(a["is_strict_splice"]) for a in annotations]
        )

    def _process_other_sites(self, df, index, start, end, transcripts, transcript_type):

        annotations = []

        for transcript in transcripts:
            # print("transcript:", transcript, start)

            # use the gene symbol we found in the db for
            # reporting rather than the input one which can be mangled
            # because of Excel etc
            gene_symbol = transcript["gene_symbol"]

            transcript_db_id = transcript["tid"]
            transcript_id = transcript["transcript_id"]

            #
            # Lets determine closest exon to splice site and test if we
            # are within a strict splice site distance to it

            self.cursor.execute(
                FIND_CLOSEST_EXON_SQL,
                {
                    "start": start,
                    "end": end,
                    "transcript_id": transcript_db_id,
                },
            )

            closest_exon = self.cursor.fetchone()

            if closest_exon is None:
                print(
                    f"No closest exon found for gene {gene_symbol} at position {start}-{end}"
                )
                continue

            min_abs_dist = closest_exon["min_abs_dist"]
            exon_number = closest_exon["exon_number"]

            annotations.append(
                {
                    "gene_symbol": gene_symbol,
                    "nearest_exon": exon_number,
                    "nearest_exon_dist": min_abs_dist,
                    "transcript_id": transcript_id,
                    "transcript_type": transcript_type,
                }
            )

        if len(annotations) == 0:
            # cds starts after splice site so likely 5' UTR
            annotations = [
                {
                    "gene_symbol": gene_symbol,
                    "nearest_exon": -1,
                    "nearest_exon_dist": -1,
                    "transcript_id": "",
                    "transcript_type": "",
                }
            ]

        # special case where we only have 1 annotation
        # so we will not repeat it in the splice change column

        # df.at[index, "Gene_Symbol"] = SEP.join(
        #     sorted(set([a["gene_symbol"] for a in annotations]))
        # )

        df.at[index, "Splice_Transcript_ID"] = SEP.join(
            [a["transcript_id"] for a in annotations]
        )

        df.at[index, SPLICE_TRANSCRIPT_COL] = SEP.join(
            [str(a["transcript_type"]) for a in annotations]
        )

        df.at[index, "Splice_Nearest_Exon_Number"] = SEP.join(
            [str(a["nearest_exon"]) for a in annotations]
        )
        df.at[index, "Splice_Nearest_Exon_Dist_bp"] = SEP.join(
            [str(a["nearest_exon_dist"]) for a in annotations]
        )

    def annotate_splice_sites(self, fin: str, fout: str, chunksize: int = 200000):

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

            df[SPLICE_PROTEIN_CHANGE_COL] = ""
            df["Splice_Transcript_ID"] = ""

            idx = np.where(df.columns == "Splice_Transcript_ID")[0][0]
            df.insert(idx + 1, SPLICE_TRANSCRIPT_COL, "")

            # df["Splice_Nearest_AA_Position"] = ""
            df["Splice_Nearest_Exon_Number"] = ""

            # find index of Splice_Nearest_Exon column
            if "Splice_Nearest_Exon_Dist_bp" not in df.columns:
                idx = np.where(df.columns == "Splice_Nearest_Exon_Number")[0][0]
                df.insert(idx + 1, "Splice_Nearest_Exon_Dist_bp", "")

            df[STRICT_SPLICE_COL] = ""

            df["Splice_Annotation_Database"] = self.metadata["version"]

            for index, row in df.iterrows():
                pc += 1

                variant_classification = row["VEP_Variant_Classification"]

                if "splice" not in variant_classification.lower():
                    continue

                # if row["Hugo_Symbol"] != "AC005019.3":
                # continue

                start = row["Start_Position"]
                end = row["End_Position"]
                gene_symbol = row["VEP_Gene_Symbol"]
                chr = row["Chromosome"]
                transcript = row["VEP_Transcript"]

                if chr == "chrMT":
                    chr = "chrM"

                chr_id = self.chr_map.get(chr, "")

                # exact
                transcript_type = "1"

                self.cursor.execute(
                    FIND_MATCHING_TRANSCRIPT,
                    {"chr": chr_id, "position": start, "transcript": transcript},
                )
                transcripts = self.cursor.fetchall()

                # print("Annotating splice site for row:", index, pc, gene_symbol, len(transcripts))

                # c for canonical

                if not transcripts:
                    print(
                        f"No exact transcripts found for gene {gene_symbol} {transcript} at position {start}-{end}, trying longest transcript..."
                    )

                    self.cursor.execute(
                        FIND_CANONICAL_TRANSCRIPTS,
                        {"chr": chr_id, "position": start},
                    )
                    transcripts = self.cursor.fetchall()

                    # for canonical since we can't find the exact transcript we will just report the canonical one
                    transcript_type = "2"

                    if not transcripts:
                        print(
                            f"No canonical transcripts found for gene {gene_symbol} {transcript} at position {start}-{end}, trying longest transcript..."
                        )

                        # try longest transcript if no canonical found
                        self.cursor.execute(
                            FIND_LONGEST_TRANSCRIPTS, {"chr": chr_id, "position": start}
                        )
                        transcripts = self.cursor.fetchall()

                        # l for longest
                        transcript_type = "3"

                        if not transcripts:
                            print(
                                f"No transcripts found for gene {gene_symbol} at position {start}-{end}. This will be skipped."
                            )
                            continue

                # if "splice" in row["Variant_Classification"].lower():
                self._process_splice_sites(
                    df, index, start, end, transcripts, transcript_type
                )
                # else:
                #    process_other_sites(
                #        index, row, chr, start, end, transcripts, transcript_type
                #    )

                if pc % 1000 == 0:
                    print(f"Processed {pc} rows of {df.shape[0]}")

            print(f"Processed {pc} splice site variants")

            df.to_csv(
                fout,
                sep="\t",
                header=first,
                mode="w" if first else "a",
                index=False,
            )

            first = False
