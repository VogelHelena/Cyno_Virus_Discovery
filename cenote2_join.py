#!/usr/bin/env python3
import pandas as pd
import glob
import os
import re
from typing import Dict

def extract_sample_prefix_from_pileup(basename: str) -> str:
    if basename.endswith("_pileup3_cov.txt"):
        return basename[: -len("_pileup3_cov.txt")]
    m = re.match(r"(.+?)_pileup3_cov", basename)
    if m:
        return m.group(1)
    return os.path.splitext(basename)[0]

def extract_sample_prefix_from_contig_summary(basename: str) -> str:
    for pat in ["_c_nh2_CONTIG_SUMMARY.tsv", "_CONTIG_SUMMARY.tsv"]:
        if basename.endswith(pat):
            basename = basename[: -len(pat)]
            break
    basename = re.sub(r"_c_nh2$", "", basename)
    return basename

def _read_table(path: str) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep='\t', engine='python')
    except Exception:
        return pd.read_csv(path, sep=r'\s+', engine='python')

def load_pileup_files() -> Dict[str, pd.DataFrame]:
    pileup_map = {}
    for path in glob.glob("*_pileup3_cov.txt"):
        base = os.path.basename(path)
        sample = extract_sample_prefix_from_pileup(base)
        df = _read_table(path)
        df.columns = [c.strip() for c in df.columns]
        if '#ID' not in df.columns:
            for cand in ['ID', 'Contig', 'Contig_ID', 'ORIGINAL_NAME']:
                if cand in df.columns:
                    df['#ID'] = df[cand]
                    break
        if '#ID' in df.columns:
            df['#ID'] = df['#ID'].astype(str).str.strip().str.replace(r'^"|"$', '', regex=True)
        df['__sample_prefix'] = sample
        df['__source_file'] = base
        pileup_map[sample] = df
    return pileup_map

def load_contig_summaries() -> Dict[str, pd.DataFrame]:
    contig_map = {}
    candidates = set(glob.glob("*_c_nh2/*_CONTIG_SUMMARY.tsv"))
    candidates.update(glob.glob("**/*_CONTIG_SUMMARY.tsv", recursive=True))
    for path in sorted(candidates):
        base = os.path.basename(path)
        parent = os.path.basename(os.path.dirname(path))
        sample = extract_sample_prefix_from_contig_summary(base)
        if sample == "" or sample.endswith("_CONTIG_SUMMARY"):
            sample = extract_sample_prefix_from_contig_summary(parent)
        df = _read_table(path)
        df.columns = [c.strip() for c in df.columns]
        # normalize string columns
        for col in ['ORIGINAL_NAME', 'CENOTE_NAME', 'ORGANISM_NAME', 'END_FEATURE', 'ORF_CALLER', 'HALLMARK_NAMES']:
            if col in df.columns:
                df[col] = df[col].astype(str).str.strip().str.replace(r'^"|"$', '', regex=True)
        # avoid column collision with pileup LENGTH by renaming to CENOTE_LENGTH
        if 'LENGTH' in df.columns:
            df = df.rename(columns={'LENGTH': 'CENOTE_LENGTH'})
        df['__sample_prefix'] = sample
        df['__source_file'] = os.path.join(os.path.basename(os.path.dirname(path)), base)
        contig_map[sample] = df
    return contig_map

def main(output_dir: str = "."):
    pileups = load_pileup_files()
    contigs = load_contig_summaries()

    common_samples = sorted(set(pileups.keys()) & set(contigs.keys()))
    if not common_samples:
        print("No matching sample prefixes between pileup and contig summary files.")
        print("Pileup samples:", sorted(pileups.keys()))
        print("Contig summary samples:", sorted(contigs.keys()))
        return

    for sample in common_samples:
        df_p = pileups[sample].copy()
        df_c = contigs[sample].copy()

        if '#ID' not in df_p.columns:
            raise ValueError(f"Pileup for sample {sample} is missing '#ID' column. Columns: {df_p.columns.tolist()}")

        # Columns to pull from Cenote summary
        cenote_cols = ['ORIGINAL_NAME', 'CENOTE_NAME', 'ORGANISM_NAME', 'END_FEATURE',
                       'CENOTE_LENGTH', 'ORF_CALLER', 'NUM_HALLMARKS', 'HALLMARK_NAMES', '__source_file']

        cenote_cols = [c for c in cenote_cols if c in df_c.columns]
        df_c_sub = df_c[cenote_cols].rename(columns={'__source_file': 'CENOTE_SOURCE'})

        merged = df_p.merge(
            df_c_sub,
            left_on='#ID',
            right_on='ORIGINAL_NAME',
            how='left',
            validate='many_to_one'
        )

        merged['PILEUP_CONTIG_ID'] = merged['#ID'].astype(str)
        merged['CONTIG_ORIGINAL_NAME'] = merged['ORIGINAL_NAME'].where(merged['ORIGINAL_NAME'].notna(), merged['PILEUP_CONTIG_ID'])

        # order pileup columns and exclude '#ID'
        pileup_cols = [c for c in df_p.columns if not c.startswith('__') and c != '#ID']

        # final columns
        front_cols = ['__sample_prefix', 'CENOTE_NAME', 'CONTIG_ORIGINAL_NAME',
                      'ORGANISM_NAME', 'END_FEATURE', 'CENOTE_LENGTH', 'ORF_CALLER',
                      'NUM_HALLMARKS', 'HALLMARK_NAMES']
        tail_cols = ['CENOTE_SOURCE', 'PILEUP_CONTIG_ID']

        out_cols = [c for c in front_cols if c in merged.columns] + pileup_cols + [c for c in tail_cols if c in merged.columns]
        out_all = merged[out_cols].copy().rename(columns={'__sample_prefix': 'SAMPLE'})

        out_name = f"{sample}_CENOTE2_contig_abundance.tsv"
        out_all.to_csv(os.path.join(output_dir, out_name), sep='\t', index=False)

        # matched-only
        matched = merged[merged['CENOTE_NAME'].notna()].copy()
        out_matched = matched[out_cols].copy().rename(columns={'__sample_prefix': 'SAMPLE'})
        matched_name = f"{sample}_CENOTE2_matched_contigs.tsv"
        out_matched.to_csv(os.path.join(output_dir, matched_name), sep='\t', index=False)

        total = len(merged)
        mcount = len(matched)
        print(f"[{sample}] total={total} matched={mcount} ({(mcount/total*100 if total else 0):.1f}%)")
        print(f" -> wrote {out_name} (all contigs)")
        print(f" -> wrote {matched_name} (matched contigs only)")

if __name__ == "__main__":
    main()

