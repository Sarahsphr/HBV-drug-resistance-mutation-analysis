# count_RT_mutations.py
# Counts known RT mutations per genotype from aligned FASTA files
# and saves results as a CSV table.
# No reference sequence is used.

import re
import csv
from pathlib import Path
from Bio import SeqIO

#PATH TO YOUR RT ALIGNMENTS 

ALIGN_DIR = Path("data/Alig_RT")

#OUTPUT CSV 
OUTPUT_CSV = ALIGN_DIR / "RT_mutation_counts.csv"

#MUTATIONS TO COUNT 
KNOWN_MUTATIONS = [
    "L180M+M204V",
    "L180M+M204V+T184L",
    "L180M+M204V+S202G",
    "L180M+M204V+M250V",
    "L180M+M204V+T184G+S202I",
    "L180M+A181C+M204V",
    "L180M+A181C+M204V+M250V",
    "L180M+A181C+S202G+M204V",
    "L180M+A181V+M204V",
    "L180M+M204V+T184A",
    "L180M+M204V+M250V+I169T",
    "A181V",
    "A181T",
    "A181V+N236T",
    "A181T+N236T",
    "M204I",
    "N236T",
    "A194T",
    "A194T+L180M+M204V",
    "S106C+H126Y+D134E",
    "S106C+H126Y+D134E+L269I",
]

#REQUIRED GENOTYPE ORDER 
GENOTYPE_ORDER = [
    "A", "B", "C", "D", "E", "F", "G", "H", "I", "J",
    "AB", "AC", "AD", "AE", "AF", "AG",
    "BC", "BD", "CD", "CG", "DE", "DF", "EF", "FG", "FH", "GH"
]

MUT_RE = re.compile(r"^([A-Z\*])(\d+)([A-Z\*])$")


def parse_mut(mut: str):
    """
    Convert mutation string into parts.

    Example:
      L180M -> [('L', 180, 'M')]
      L180M+M204V -> [('L', 180, 'M'), ('M', 204, 'V')]

    The first amino acid is kept only for reading the mutation name.
    It is NOT used for reference checking.
    """
    mut = mut.strip().upper().replace(" ", "")
    parts = mut.split("+")
    parsed = []

    for part in parts:
        m = MUT_RE.match(part)
        if not m:
            raise ValueError(
                f"Bad mutation format: {mut} "
                f"(expected like L180M or L180M+M204V)"
            )

        ref_aa, pos, alt_aa = m.group(1), int(m.group(2)), m.group(3)
        parsed.append((ref_aa, pos, alt_aa))

    return parsed


def has_mutation_direct(seq_aln: str, mut_parts) -> bool:
    """
    Direct mutation counting.

    For L180M:
      - check sequence position 180 directly
      - count only if amino acid is M

    For combined mutations:
      - all listed mutant amino acids must be present

    Gaps and unknown residues are treated as no mutation.
    """
    seq_aln = seq_aln.upper()

    for ref_aa, pos, alt_aa in mut_parts:
        col = pos - 1

        if col >= len(seq_aln):
            return False

        aa = seq_aln[col]

        if aa in {"-", "X", "?"}:
            return False

        if aa != alt_aa:
            return False

    return True


def genotype_from_filename(name: str) -> str:
    """
    Examples:
      GTA_alig_RT.fas -> A
      GTAB_alig_RT.fas -> AB
      GTJ_alig_RT.fas -> J
    """
    prefix = name.split("_")[0].upper()
    gt = prefix.replace("GT", "", 1)

    if not gt:
        raise ValueError(f"Could not parse genotype from filename: {name}")

    return gt


def analyze_alignment(file_path: Path):
    records = list(SeqIO.parse(str(file_path), "fasta"))

    if not records:
        return 0, {m: 0 for m in KNOWN_MUTATIONS}

    parsed = {m: parse_mut(m) for m in KNOWN_MUTATIONS}
    counts = {m: 0 for m in KNOWN_MUTATIONS}

    # Count every sequence in the file.
    # No reference sequence is excluded.
    for rec in records:
        seq_aln = str(rec.seq).upper()

        for mut, mut_parts in parsed.items():
            if has_mutation_direct(seq_aln, mut_parts):
                counts[mut] += 1

    return len(records), counts


def main():
    if not ALIGN_DIR.exists():
        raise FileNotFoundError(f"Folder not found: {ALIGN_DIR}")

    files = sorted(ALIGN_DIR.glob("GT*_alig_RT*.fas*"))

    if not files:
        raise FileNotFoundError(
            f"No RT alignment files found in {ALIGN_DIR} "
            f"matching GT*_alig_RT*.fas*"
        )

    results = {}
    skipped_files = []

    for f in files:
        try:
            gt = genotype_from_filename(f.name)
            total, counts = analyze_alignment(f)

            results[gt] = {
                "Total Sequence": total,
                **counts
            }

        except Exception as e:
            skipped_files.append((f.name, str(e)))

    print("\n========== RT mutation counts ==========\n")

    for gt in GENOTYPE_ORDER:
        if gt in results:
            row = results[gt]
            print(f"Genotype {gt}")
            print(f"  Total Sequence: {row['Total Sequence']}")

            for mut in KNOWN_MUTATIONS:
                print(f"  {mut}: {row[mut]}")

            print()

    extra_gts = sorted(set(results.keys()) - set(GENOTYPE_ORDER))

    for gt in extra_gts:
        row = results[gt]
        print(f"Genotype {gt}")
        print(f"  Total Sequence: {row['Total Sequence']}")

        for mut in KNOWN_MUTATIONS:
            print(f"  {mut}: {row[mut]}")

        print()

    with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)

        writer.writerow(["Genotype", "Total Sequence"] + KNOWN_MUTATIONS)

        for gt in GENOTYPE_ORDER:
            if gt in results:
                row = results[gt]
                writer.writerow(
                    [gt, row["Total Sequence"]] +
                    [row[m] for m in KNOWN_MUTATIONS]
                )

        for gt in extra_gts:
            row = results[gt]
            writer.writerow(
                [gt, row["Total Sequence"]] +
                [row[m] for m in KNOWN_MUTATIONS]
            )

    print(f"CSV saved to:\n{OUTPUT_CSV}\n")

    if skipped_files:
        print("Skipped files:")

        for fname, err in skipped_files:
            print(f"  {fname}: {err}")


if __name__ == "__main__":
    main()
