#!/usr/bin/env python3
import sys, re, glob
from pathlib import Path

# --- Helpers ---------------------------------------------------------------

def norm_txt(s: str) -> str:
    """Normalize text to UTF-8 friendly and collapse odd spaces."""
    if s is None:
        return ""
    # Replace common non-breaking/figure spaces with normal spaces
    s = s.replace("\u00A0"," ").replace("\u2007"," ").replace("\u202F"," ")
    # Collapse whitespace and strip
    s = re.sub(r"\s+", " ", s).strip()
    return s

def read_tsv_safely(p):
    """Read tab file without pandas, tolerant to encodings."""
    import csv
    # try utf-8 first, fallback to latin-1
    for enc in ("utf-8", "utf-8-sig", "latin-1"):
        try:
            with open(p, "r", encoding=enc, errors="replace", newline="") as fh:
                rows = list(csv.reader(fh, delimiter="\t"))
            return rows
        except Exception:
            continue
    raise RuntimeError(f"Could not read {p} with UTF-8 or latin-1.")

def find_count_col(header):
    """Return index of the total-reads column (case/variant tolerant)."""
    candidates = [r"^total_reads$", r"^Total_reads$", r"^TOTAL_READS$",
                  r"^reads$", r"^Total$", r"^TOTAL$"]
    for i, h in enumerate(header):
        n = norm_txt(h)
        for pat in candidates:
            if re.match(pat, n, flags=re.IGNORECASE):
                return i
    # If nothing obvious, last column is often counts
    return len(header) - 1

def find_org_col(header):
    for i, h in enumerate(header):
        if norm_txt(h).lower() in ("organism_name", "organism", "taxon", "taxon_name"):
            return i
    # fallback: second column
    return 1

def find_name_col(header):
    for i, h in enumerate(header):
        if norm_txt(h).lower() in ("cenote_name", "filename", "file", "sample"):
            return i
    # fallback: first column
    return 0

# --- Grouping rules (ORDER MATTERS: first match wins) ---------------------

# Each entry: ("OutputGroupName", [list of substrings to match (case-insensitive)])
RULES = [
    ("Polyomaviridae",
        ["betapolyomavirus", "polyomaviridae"]),

    ("Anelloviridae",
        ["anelloviridae"]),

    ("Papillomaviridae",
        ["betapapillomavirus", "papillomaviridae", "gammapapillomavirus"]),

    ("Mastadenovirus",
        ["mastadenovirus"]),

    ("Parvoviridae",
        ["parvoviridae", "dependoparvovirus"]),

    ("Herpesviridae",
        ["orthoherpesviridae", "cytomegalovirus", "rhadinovirus"]),

    ("Hepadnaviridae",
        ["orthohepadnavirus"]),

    ("Bacteriophage",
        ["phage", "aliceevansviridae", "petitvirales", "microviridae",
         "autographiviridae", "autolykiviridae", "bacillales",
         "staphylococcaceae", "rountreeviridae",
         "orlajensenviridaem", "orlajensenviridae",  # include both spellings
         "kyanoviridae", "inoviridae", "herelleviridae", "glaedevirus",
         "emdodecavirus", "dubowvirus", "demerecviridae", "crassvirales",
         "casjensviridae", "salasmaviridae", "schitoviridae", "steigviridae",
         "mycoplasmataceae", "mycoplasmatales", "tubulavirales",
         "winoviridae", "peduoviridae"]),

    ("Cressdnaviricota",
        ["cirlivirales", "circoviridae", "vilyaviridae", "rohanvirales",
         "sepolyvirales", "cremevirales", "geminiviridae", "begomovirus",
         "geplafuvirales", "genomoviridae", "gemykrogvirus",
         "nenyaviridae", "smacoviridae", "porprismacovirus"]),

    ("Caudoviricetes",
        ["bievrevirus"]),

    ("Varidnaviria",
        ["imitervirales", "pimascovirales", "algavirales", "adintovirus"]),

    # Keep "Other" last so it doesn't capture specific families first
    ("Other",
        ["virus"]),
]

def apply_rules(name: str):
    """Return group string or None (to drop unassigned)."""
    q = norm_txt(name).lower()
    for group, needles in RULES:
        if any(n in q for n in needles):
            return group
    return None  # drop if no rule matches

# --- Main -----------------------------------------------------------------

def process_one(tsv_path: Path):
    rows = read_tsv_safely(tsv_path)
    if not rows:
        return

    # Normalize header
    header = [norm_txt(x) for x in rows[0]]
    i_name = find_name_col(header)
    i_org  = find_org_col(header)
    i_cnt  = find_count_col(header)

    out_rows = {}  # group -> count sum
    sample_name = Path(tsv_path.name).name  # keep input filename as CENOTE_NAME

    for r in rows[1:]:
        if not r or len(r) <= max(i_name, i_org, i_cnt):
            continue
        org = norm_txt(r[i_org])
        # coerce count
        try:
            cnt = float(norm_txt(r[i_cnt]).replace(",", ""))
        except Exception:
            cnt = 0.0
        grp = apply_rules(org)
        if grp is None:
            # drop anything unassigned
            continue
        out_rows[grp] = out_rows.get(grp, 0.0) + cnt

    # Write output
    out_path = Path(str(tsv_path).replace("_CENOTE2_contig_abundance.tsv",
                                          "_CENOTE2_GROUP.abundance.tsv"))
    with open(out_path, "w", encoding="utf-8") as out:
        out.write("CENOTE_NAME\tORGANISM_NAME\tTotal_reads\n")
        # sort safely by group name
        for grp in sorted(out_rows.keys(), key=lambda s: s.lower()):
            cnt = out_rows[grp]
            out.write(f"{sample_name}\t{grp}\t{int(round(cnt))}\n")
    print(f"Wrote: {out_path}")

def main():
    files = sorted(glob.glob("*_CENOTE2_contig_abundance.tsv"))
    if not files:
        print("No *_CENOTE2_contig_abundance.tsv files found in this folder.", file=sys.stderr)
        sys.exit(1)
    for f in files:
        process_one(Path(f))

if __name__ == "__main__":
    main()

