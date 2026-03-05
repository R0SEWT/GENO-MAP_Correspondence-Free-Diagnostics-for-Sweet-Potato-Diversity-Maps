#!/usr/bin/env python3
"""Centralised loader for DArT / DArTSeq genotype and metrics files.

Implements:
  - ADR-0003: missingness filtering (sample/marker thresholds, max-impute cap)
  - ADR-0004: sentinel normalisation (Chrom/ChromPos/AlnEvalue), chromosome
    name standardisation, ``is_mapped`` flag, genotype NaN handling
  - ADR-0002: duplicate guard (Wild SNP == LowDensity SNP)

Usage as a library::

    from load_dart import load_genotypes, load_metrics, filter_missingness

Usage as a quick self-test::

    python scripts/load_dart.py --test
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Column classification rules (shared with EDA notebook)
# ---------------------------------------------------------------------------

CAT_RULES: Dict[str, List[str]] = {
    "ID": ["alleleid", "cloneid", "markerid", "snpid", "uid"],
    "Secuencia": [
        "allelesequence",
        "trimmedsequence",
        "snp",
        "trimmedsequenceref",
    ],
    "Genómica": [
        "chrom",
        "chrompos",
        "alncount",
        "alncnt",
        "alnevalue",
        "strand",
        "genomeposition",
    ],
    "Calidad": [
        "callrate",
        "oneratioref",
        "oneratiosnp",
        "freqhomref",
        "freqhomsnp",
        "freqhet",
        "picref",
        "picsnp",
        "avgpic",
        "avgcountref",
        "avgcountsnp",
        "ratioavgcountrefavgcountsnp",
        "reproducibility",
        "repavg",
        "rowsum",
        "polymorphicinformationcontent",
        "readcountref",
        "readcountsnp",
    ],
}


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------


def md5_file(path: Path, chunk_size: int = 65536) -> str:
    """Compute MD5 hex digest of a file."""
    h = hashlib.md5()
    with open(path, "rb") as f:
        while chunk := f.read(chunk_size):
            h.update(chunk)
    return h.hexdigest()


def detect_separator(path: Path, n_bytes: int = 8192) -> str:
    """Auto-detect field separator (``,``, ``;``, or ``\\t``)."""
    with open(path, "r", newline="", encoding="utf-8-sig") as f:
        sample = f.read(n_bytes)
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",;\t")
        return dialect.delimiter
    except csv.Error:
        return "," if sample.count(",") >= sample.count(";") else ";"


def classify_column(col_name: str) -> str:
    """Classify a column into ID / Secuencia / Genómica / Calidad / Muestra / Otro."""
    low = col_name.lower().replace("_", "").replace(" ", "")
    for cat, keywords in CAT_RULES.items():
        for kw in keywords:
            if kw in low:
                return cat
    if col_name and col_name[0].isdigit():
        return "Muestra"
    return "Otro"


def _strip_bom(cols: List[str]) -> List[str]:
    """Strip possible BOM characters from column names."""
    return [c.strip().strip("\ufeff") for c in cols]


def _normalise_chrom(s: pd.Series) -> pd.Series:
    """Normalise chromosome names: ``Chr01`` → ``chr1``, keep ``*`` as-is."""
    def _norm(v: object) -> object:
        if pd.isna(v):
            return v
        v_str = str(v).strip()
        if v_str == "*":
            return v_str
        m = re.match(r"[Cc]hr0*(\d+)", v_str)
        if m:
            return f"chr{m.group(1)}"
        return v_str
    return s.map(_norm)


# ---------------------------------------------------------------------------
# Genotype loading
# ---------------------------------------------------------------------------


def detect_format(header: List[str]) -> str:
    """Detect CSV layout: ``sample_columns`` or ``marker_metrics``."""
    if header and header[0].strip() == "Sample_code":
        return "sample_columns"
    if any("AlleleID" in h for h in header):
        return "marker_metrics"
    return "sample_columns"


def load_genotypes(
    path: Path,
    *,
    max_rows: int | None = None,
    max_sample_cols: int | None = None,
    seed: int = 42,
    fmt: str = "auto",
) -> Tuple[pd.DataFrame, List[str], Dict[str, Dict[str, str]]]:
    """Load a DArT genotype file and return a *samples × markers* matrix.

    Parameters
    ----------
    path : Path
        CSV/TSV genotype file.
    max_rows : int | None
        Limit marker rows read (``nrows``).  ``None`` reads all.
    max_sample_cols : int | None
        Subsample sample columns (for EDA).  ``None`` keeps all.
    seed : int
        RNG seed for column subsampling.
    fmt : str
        ``"auto"`` (default), ``"sample_columns"``, or ``"marker_metrics"``.

    Returns
    -------
    X : pd.DataFrame
        Shape *(n_samples, n_markers)*.  Genotype values are numeric
        (0/1/2 for SNP, 0/1 for SilicoDArT) with ``NaN`` for missing.
    sample_ids : list[str]
        Sample identifiers (column headers in original).
    sample_meta : dict
        Per-sample metadata extracted from non-marker rows (``sample_columns``
        format only).
    """
    sep = detect_separator(path)

    # Read header
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.reader(f, delimiter=sep)
        header = _strip_bom(next(reader))

    if fmt == "auto":
        fmt = detect_format(header)

    if fmt == "sample_columns":
        return _load_sample_columns(path, sep, header, max_rows, max_sample_cols, seed)
    else:
        return _load_marker_metrics(path, sep, header, max_rows, max_sample_cols, seed)


def _load_sample_columns(
    path: Path,
    sep: str,
    header: List[str],
    max_rows: int | None,
    max_sample_cols: int | None,
    seed: int,
) -> Tuple[pd.DataFrame, List[str], Dict[str, Dict[str, str]]]:
    """Layout where col-0 = row label, cols 1..n = sample IDs."""
    col_cats = {col: classify_column(col) for col in header}
    sample_cols = [col for col, cat in col_cats.items() if cat == "Muestra"]
    meta_cols = [col for col, cat in col_cats.items() if cat != "Muestra"]

    if max_sample_cols and len(sample_cols) > max_sample_cols:
        rng = np.random.default_rng(seed)
        sample_cols = list(rng.choice(sample_cols, max_sample_cols, replace=False))

    usecols = meta_cols[:1] + sample_cols
    df = pd.read_csv(
        path,
        sep=sep,
        nrows=max_rows,
        usecols=usecols,
        low_memory=False,
        dtype=str,
        encoding="utf-8-sig",
    )
    df.columns = _strip_bom(list(df.columns))

    first_col = df.columns[0]
    labels = df[first_col].astype(str)
    marker_mask = labels.str.contains(r"[0-9]", na=False)
    meta_rows = df.loc[~marker_mask]
    marker_rows = df.loc[marker_mask].copy()

    # Extract per-sample metadata from non-marker rows
    sample_meta: Dict[str, Dict[str, str]] = {sid: {} for sid in sample_cols}
    for _, row in meta_rows.iterrows():
        key = str(row.iloc[0]).strip()
        for sid in sample_cols:
            if sid in row.index:
                val = row[sid]
                if pd.notna(val) and str(val).strip():
                    sample_meta[sid][key] = str(val).strip()

    # Build genotype matrix
    geno = marker_rows[[c for c in sample_cols if c in marker_rows.columns]].copy()
    geno = geno.replace(["-", "", " "], np.nan)
    geno = geno.apply(pd.to_numeric, errors="coerce")
    geno.index = marker_rows[first_col].values
    X = geno.T  # samples × markers
    sample_ids = list(X.index.astype(str))

    return X, sample_ids, sample_meta


def _load_marker_metrics(
    path: Path,
    sep: str,
    header: List[str],
    max_rows: int | None,
    max_sample_cols: int | None,
    seed: int,
) -> Tuple[pd.DataFrame, List[str], Dict[str, Dict[str, str]]]:
    """Layout with leading metadata columns, sample columns start where the
    header begins with a digit."""
    sample_start = next(
        (i for i, name in enumerate(header) if name and name[0].isdigit()), None
    )
    if sample_start is None:
        raise ValueError(f"Cannot detect sample columns in {path}")

    sample_cols = header[sample_start:]
    if max_sample_cols and len(sample_cols) > max_sample_cols:
        rng = np.random.default_rng(seed)
        sample_cols = list(rng.choice(sample_cols, max_sample_cols, replace=False))

    marker_id_col = header[0]
    usecols = [marker_id_col] + sample_cols
    df = pd.read_csv(
        path,
        sep=sep,
        nrows=max_rows,
        usecols=usecols,
        low_memory=False,
        dtype=str,
        encoding="utf-8-sig",
    )
    df.columns = _strip_bom(list(df.columns))

    geno = df[[c for c in sample_cols if c in df.columns]].copy()
    geno = geno.replace(["-", "", " "], np.nan)
    geno = geno.apply(pd.to_numeric, errors="coerce")
    geno.index = df[marker_id_col].astype(str).values
    X = geno.T  # samples × markers
    sample_ids = list(X.index.astype(str))
    sample_meta: Dict[str, Dict[str, str]] = {s: {} for s in sample_ids}

    return X, sample_ids, sample_meta


# ---------------------------------------------------------------------------
# Metrics loading  (ADR-0004)
# ---------------------------------------------------------------------------


def load_metrics(path: Path) -> pd.DataFrame:
    """Load a DArT metrics/report file with sentinel normalisation.

    Applies ADR-0004 rules:
      - ``ChromPos = 0`` → ``NaN``
      - ``AlnEvalue = 999`` → ``NaN``
      - ``AlnCnt = 0`` (when Chrom is ``*``) → ``NaN``
      - ``Chrom`` names normalised (``Chr01`` → ``chr1``)
      - Boolean column ``is_mapped`` added
    """
    sep = detect_separator(path)
    df = pd.read_csv(path, sep=sep, low_memory=False, encoding="utf-8-sig")
    df.columns = _strip_bom(list(df.columns))

    # --- Sentinel normalisation ---
    chrom_cols = [c for c in df.columns if c.lower().startswith("chrom") and "pos" not in c.lower()]
    chrompos_cols = [c for c in df.columns if re.match(r"(?i)chrom.?pos", c)]
    alnevalue_cols = [c for c in df.columns if re.match(r"(?i)aln.?evalue", c)]
    alncnt_cols = [c for c in df.columns if re.match(r"(?i)aln.?cnt", c)]

    # is_mapped: True when *none* of the Chrom columns equals "*"
    if chrom_cols:
        is_mapped = ~df[chrom_cols[0]].astype(str).str.strip().eq("*")
        df.insert(len(chrom_cols), "is_mapped", is_mapped)
        for c in chrom_cols:
            df[c] = _normalise_chrom(df[c].astype(str))
            df.loc[df[c] == "*", c] = np.nan

    for c in chrompos_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
        df.loc[df[c] == 0, c] = np.nan

    for c in alnevalue_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
        df.loc[df[c] == 999, c] = np.nan

    for c in alncnt_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
        # Zero alignment count is sentinel only for unmapped markers
        if "is_mapped" in df.columns:
            df.loc[~df["is_mapped"] & (df[c] == 0), c] = np.nan

    return df


# ---------------------------------------------------------------------------
# Missingness filtering  (ADR-0003)
# ---------------------------------------------------------------------------


def filter_missingness(
    X: pd.DataFrame,
    *,
    sample_thresh: float = 0.50,
    marker_thresh: float = 0.50,
    max_impute_frac: float = 0.30,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Apply ADR-0003 missingness thresholds.

    Parameters
    ----------
    X : pd.DataFrame
        Genotype matrix *samples × markers*.
    sample_thresh : float
        Drop samples with more than this fraction missing.
    marker_thresh : float
        Drop markers with more than this fraction missing.
    max_impute_frac : float
        Flag markers above this threshold — they should NOT be imputed.

    Returns
    -------
    X_clean : pd.DataFrame
        Filtered matrix.
    report : dict
        ``samples_dropped``, ``markers_dropped``, ``markers_no_impute``,
        ``shape_before``, ``shape_after``, etc.
    """
    shape_before = X.shape

    # --- Drop samples exceeding threshold ---
    sample_miss = X.isna().mean(axis=1)
    bad_samples = sample_miss[sample_miss > sample_thresh].index.tolist()
    X = X.drop(index=bad_samples, errors="ignore")

    # --- Drop markers exceeding threshold ---
    marker_miss = X.isna().mean(axis=0)
    bad_markers = marker_miss[marker_miss > marker_thresh].index.tolist()
    X = X.drop(columns=bad_markers, errors="ignore")

    # --- Flag markers that should NOT be imputed (still above max_impute_frac) ---
    remaining_miss = X.isna().mean(axis=0)
    no_impute_markers = remaining_miss[remaining_miss > max_impute_frac].index.tolist()

    report: Dict[str, Any] = {
        "shape_before": list(shape_before),
        "shape_after": list(X.shape),
        "samples_dropped": len(bad_samples),
        "samples_dropped_ids": bad_samples,
        "markers_dropped": len(bad_markers),
        "markers_no_impute": len(no_impute_markers),
        "markers_no_impute_ids": no_impute_markers,
        "sample_thresh": sample_thresh,
        "marker_thresh": marker_thresh,
        "max_impute_frac": max_impute_frac,
        "missing_rate_before": float(
            pd.isna(pd.DataFrame(index=range(shape_before[0]), columns=range(shape_before[1]))).sum().sum()
            if shape_before[0] == 0 else 0
        ),
        "missing_rate_after": float(X.isna().sum().sum() / max(X.shape[0] * X.shape[1], 1)),
    }

    return X, report


# ---------------------------------------------------------------------------
# Duplicate guard  (ADR-0002)
# ---------------------------------------------------------------------------

# The LowDensity SNP and Wild SNP files are identical (MD5: ae014cb2a8a17cc215f73b545c571f8b).
KNOWN_DUPLICATE_MD5 = "ae014cb2a8a17cc215f73b545c571f8b"
WILD_SNP_FILENAME = "01_Report_DSp25-515_SNPs_Filtered_by _reads.csv"


def check_duplicate_guard(path: Path) -> bool:
    """Return True if *path* is the known Wild/LowDensity SNP duplicate.

    Does NOT block loading — just warns so callers can decide.
    """
    if path.name == WILD_SNP_FILENAME and "Wild" in str(path):
        return True
    return False


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------


def _self_test() -> None:
    """Quick smoke test using Global SNP (if available)."""
    data_dir = Path(__file__).resolve().parent.parent / "data"
    global_snp = data_dir / "10.21223P30BVZYY_Genetic_diversity" / "SNP_Genotypes.csv"
    global_silico_metrics = data_dir / "10.21223P30BVZYY_Genetic_diversity" / "SNP_metrics.csv"

    print("=" * 60)
    print("load_dart.py — self-test")
    print("=" * 60)

    # --- Genotype loading ---
    if global_snp.exists():
        print(f"\n[1] Loading genotypes (subsample): {global_snp.name}")
        X, sids, smeta = load_genotypes(global_snp, max_rows=500, max_sample_cols=200)
        print(f"    Shape: {X.shape}  (samples={len(sids)})")
        assert X.shape[0] == len(sids), "row count mismatch"
        assert not (X == "-").any().any(), "sentinel '-' still present"
        vals = set(X.values.flatten()) - {np.nan}
        numeric_vals = {v for v in vals if not np.isnan(v)} if vals else set()
        print(f"    Unique geno values (excl NaN): {sorted(numeric_vals)[:10]}")

        # --- Missingness filtering ---
        print("\n[2] Filtering missingness (thresh=0.50)")
        X_clean, report = filter_missingness(X, sample_thresh=0.50)
        print(f"    Before: {report['shape_before']}  After: {report['shape_after']}")
        print(f"    Samples dropped: {report['samples_dropped']}")
        print(f"    Markers dropped: {report['markers_dropped']}")
        print(f"    Markers no-impute (>30%): {report['markers_no_impute']}")
    else:
        print(f"[skip] Genotype file not found: {global_snp}")

    # --- Metrics loading ---
    if global_silico_metrics.exists():
        print(f"\n[3] Loading metrics: {global_silico_metrics.name}")
        mdf = load_metrics(global_silico_metrics)
        print(f"    Shape: {mdf.shape}")
        if "is_mapped" in mdf.columns:
            mapped_pct = mdf["is_mapped"].mean() * 100
            print(f"    is_mapped: {mapped_pct:.1f}%")
        chrom_cols_found = [c for c in mdf.columns if "chrom" in c.lower() and "pos" not in c.lower()]
        if chrom_cols_found:
            sample_vals = mdf[chrom_cols_found[0]].dropna().unique()[:10]
            print(f"    Chrom values (sample): {list(sample_vals)}")
    else:
        print(f"[skip] Metrics file not found: {global_silico_metrics}")

    # --- Duplicate guard ---
    wild_snp = data_dir / "10.21223P33VYY8C_Wild" / WILD_SNP_FILENAME
    if wild_snp.exists():
        print(f"\n[4] Duplicate guard: {wild_snp.name}")
        is_dup = check_duplicate_guard(wild_snp)
        print(f"    Is known duplicate path: {is_dup}")

    print("\n" + "=" * 60)
    print("Self-test complete.")
    print("=" * 60)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="DArT loader — self-test or show info")
    parser.add_argument("--test", action="store_true", help="Run self-test with local data")
    parser.add_argument("--info", type=Path, help="Show file info (separator, format, shape)")
    args = parser.parse_args()

    if args.test:
        _self_test()
    elif args.info:
        sep = detect_separator(args.info)
        with open(args.info, "r", encoding="utf-8-sig") as f:
            header = _strip_bom(next(csv.reader(f, delimiter=sep)))
        fmt = detect_format(header)
        print(f"File:       {args.info}")
        print(f"Separator:  {repr(sep)}")
        print(f"Format:     {fmt}")
        print(f"Columns:    {len(header)}")
        print(f"Header[:5]: {header[:5]}")
    else:
        parser.print_help()
