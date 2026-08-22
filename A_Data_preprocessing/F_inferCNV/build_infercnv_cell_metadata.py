#!/usr/bin/env python3
"""Build per-cell InferCNV metadata for Seurat AddMetaData (CMETS inferCNV subclusters notebook).

Env:
  INFERCNV_ROOT — sample directories with out/<cell_groupings> (default: 20241001 results).
  OUT_CSV — output path (default: Revesion/infercnv_cell_metadata.csv).
  INFERCNV_PRED_CNV_GENES — basename under each sample out/; if unset, first match wins:
      HMM_CNV_predictions...Pnorm_0.5.pred_cnv_genes.dat, then
      17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_genes.dat
  INFERCNV_CELL_GROUPINGS — basename of cell_groupings file (default: 17_HMM_...cell_groupings).
  INFERCNV_PER_CELL_CNV_TSV — optional tab file (cell_id, infercnv_cnv_call) from
      build_infercnv_per_cell_cnv_metadata_table.R → merged as column infercnv_cnv_call.
"""

from __future__ import annotations

import csv
import glob
import os

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INFERCNV_ROOT = os.environ.get(
    "INFERCNV_ROOT",
    "/diskmnt/Projects/MetNet_analysis_2/Colorectal/InferCNV/20241001_inferCNV_results",
)
OUT_CSV = os.environ.get(
    "OUT_CSV",
    os.path.join(_SCRIPT_DIR, "infercnv_cell_metadata.csv"),
)
PRED_CANDIDATES = [
    (
        "HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters."
        "Pnorm_0.5.pred_cnv_genes.dat"
    ),
    "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_genes.dat",
]
CG_NAME = os.environ.get(
    "INFERCNV_CELL_GROUPINGS",
    "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings",
)
# Discrete HMM state used when a tumor subcluster has no rows in pred (or gene missing); matches inferCNV neutral.
NEUTRAL_HMM_STATE = int(os.environ.get("INFERCNV_NEUTRAL_HMM_STATE", "3"))
# All genes: infercnv_*_gain_hmm TRUE iff discrete HMM state >= this (matches infercnv_kras_gain_tumor_table.R).
GAIN_MIN_HMM_STATE = int(os.environ.get("INFERCNV_GAIN_MIN_HMM_STATE", "4"))

# Same default panel as infercnv_kras_gain_tumor_table.R (no TH1L). Override: INFERCNV_GENES=KRAS,EGFR
_GENES_ENV = os.environ.get("INFERCNV_GENES", "").strip()
if _GENES_ENV:
    GENES_OF_INTEREST = [g.strip() for g in _GENES_ENV.split(",") if g.strip()]
else:
    GENES_OF_INTEREST = [
        "KRAS",
        "AURKA",
        "EGFR",
        "MYC",
        "CDK8",
        "CDX2",
        "ERBB2",
        "PIK3CA",
        "ROS1",
    ]

GENE_SET = set(GENES_OF_INTEREST)


def infercnv_sample_from_path(cg_path: str) -> str:
    parts = cg_path.replace("\\", "/").split("/")
    for i, x in enumerate(parts):
        if x == "out" and i > 0:
            return parts[i - 1]
    raise ValueError(cg_path)


def gene_column_triplet(gene: str) -> tuple[str, str, str]:
    """HMM state, gain flag, region — KRAS keeps legacy *_chr12_region name."""
    if gene == "KRAS":
        return (
            "infercnv_kras_hmm_state",
            "infercnv_kras_gain_hmm",
            "infercnv_kras_chr12_region",
        )
    g = gene.lower()
    return (
        f"infercnv_{g}_hmm_state",
        f"infercnv_{g}_gain_hmm",
        f"infercnv_{g}_region",
    )


def resolve_pred_path(out_dir: str) -> str | None:
    explicit = os.environ.get("INFERCNV_PRED_CNV_GENES", "").strip()
    if explicit:
        p = os.path.join(out_dir, explicit)
        return p if os.path.isfile(p) else None
    for name in PRED_CANDIDATES:
        p = os.path.join(out_dir, name)
        if os.path.isfile(p):
            return p
    return None


def load_per_cell_cnv_calls(tsv_path: str) -> dict[str, str]:
    """cell_id -> infercnv_cnv_call (last wins on duplicates)."""
    if not tsv_path or not os.path.isfile(tsv_path):
        return {}
    out: dict[str, str] = {}
    with open(tsv_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            cid = (row.get("cell_id") or "").strip()
            if not cid:
                continue
            call = row.get("infercnv_cnv_call")
            out[cid] = "" if call is None else str(call).strip()
    return out


def build_fieldnames(include_cnv_call: bool) -> list[str]:
    head = [
        "cell",
        "infercnv_sample_id",
        "infercnv_cell_group",
        "infercnv_lineage",
        "infercnv_subcluster",
        "infercnv_is_tumor",
        "infercnv_tumor_clone",
    ]
    if include_cnv_call:
        head.append("infercnv_cnv_call")
    for gene in GENES_OF_INTEREST:
        head.extend(gene_column_triplet(gene))
    head.append("infercnv_hmm_leiden_subclusters")
    return head


def load_gene_states_by_group(pred_path: str) -> tuple[dict[str, dict[str, dict]], set[str]]:
    """(cell_group_name -> gene -> {state, gene_region_name}, all cell_group_names in pred file)."""
    out: dict[str, dict[str, dict]] = {}
    groups_in_pred: set[str] = set()
    if not os.path.isfile(pred_path):
        return out, groups_in_pred
    with open(pred_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            cg = row.get("cell_group_name") or ""
            if cg:
                groups_in_pred.add(cg)
            gsym = row.get("gene") or ""
            if gsym not in GENE_SET:
                continue
            try:
                st = int(row["state"])
            except (TypeError, ValueError):
                continue
            if cg not in out:
                out[cg] = {}
            out[cg][gsym] = {
                "state": st,
                "gene_region_name": row.get("gene_region_name", "") or "",
            }
    return out, groups_in_pred


def gene_values_for_cell_group(
    genes_by_cg: dict[str, dict[str, dict]],
    cg: str,
    is_tumor: bool,
    groups_in_pred: set[str],
) -> dict[str, str]:
    """Flat dict of column_name -> value for all genes in GENES_OF_INTEREST."""
    flat: dict[str, str] = {}
    if not is_tumor:
        for gene in GENES_OF_INTEREST:
            c_state, c_gain, c_region = gene_column_triplet(gene)
            flat[c_state] = ""
            flat[c_gain] = ""
            flat[c_region] = ""
        return flat

    # Subcluster absent from pred file → neutral HMM state at all genes of interest.
    subcluster_in_pred = cg in groups_in_pred
    gmap = genes_by_cg.get(cg, {})
    for gene in GENES_OF_INTEREST:
        c_state, c_gain, c_region = gene_column_triplet(gene)
        gd = gmap.get(gene)
        if not subcluster_in_pred or not gd:
            st = NEUTRAL_HMM_STATE
            flat[c_state] = str(st)
            flat[c_gain] = "TRUE" if st >= GAIN_MIN_HMM_STATE else "FALSE"
            flat[c_region] = ""
            continue
        st = gd["state"]
        flat[c_state] = str(st)
        flat[c_gain] = "TRUE" if st >= GAIN_MIN_HMM_STATE else "FALSE"
        flat[c_region] = gd.get("gene_region_name", "") or ""
    return flat


def main() -> None:
    pattern = os.path.join(INFERCNV_ROOT, "*", "out", CG_NAME)
    cg_files = sorted(glob.glob(pattern))
    if not cg_files:
        raise SystemExit(f"No cell_groupings files under {INFERCNV_ROOT}")

    cnv_tsv = os.environ.get("INFERCNV_PER_CELL_CNV_TSV", "").strip()
    cnv_by_cell = load_per_cell_cnv_calls(cnv_tsv)
    include_cnv_call = bool(cnv_by_cell) or bool(cnv_tsv)
    fieldnames = build_fieldnames(include_cnv_call)

    n_written = 0
    n_missing_pred = 0
    with open(OUT_CSV, "w", newline="") as outf:
        w = csv.DictWriter(outf, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()

        for cg_path in cg_files:
            sample_id = infercnv_sample_from_path(cg_path)
            out_dir = os.path.dirname(cg_path)
            pred_path = resolve_pred_path(out_dir)
            if not pred_path:
                n_missing_pred += 1
            genes_by_cg, groups_in_pred = load_gene_states_by_group(pred_path or "")

            with open(cg_path, newline="") as inf:
                reader = csv.DictReader(inf, delimiter="\t")
                for row in reader:
                    cg = row["cell_group_name"]
                    cell = row["cell"]
                    if "." in cg:
                        lineage, sub = cg.split(".", 1)
                    else:
                        lineage, sub = cg, ""
                    is_tumor = lineage == "Tumor"
                    gene_cols = gene_values_for_cell_group(
                        genes_by_cg, cg, is_tumor, groups_in_pred
                    )

                    rec = {
                        "cell": cell,
                        "infercnv_sample_id": sample_id,
                        "infercnv_cell_group": cg,
                        "infercnv_lineage": lineage,
                        "infercnv_subcluster": sub,
                        "infercnv_is_tumor": "TRUE" if is_tumor else "FALSE",
                        "infercnv_tumor_clone": sub if is_tumor else "",
                        "infercnv_hmm_leiden_subclusters": "HMMi6_leiden_subclusters",
                    }
                    if include_cnv_call:
                        rec["infercnv_cnv_call"] = cnv_by_cell.get(cell, "")
                    rec.update(gene_cols)
                    w.writerow(rec)
                    n_written += 1

    print(f"Wrote {n_written} rows to {OUT_CSV}")
    print(f"Libraries: {len(cg_files)}")
    print(f"Cell groupings: {CG_NAME}")
    print(
        "Pred file: explicit INFERCNV_PRED_CNV_GENES or first of: "
        + ", ".join(PRED_CANDIDATES)
    )
    if n_missing_pred:
        print(f"WARN: {n_missing_pred} sample(s) had no matching pred_cnv_genes.dat in out/")
    if cnv_tsv:
        print(f"Per-cell CNV call TSV: {cnv_tsv} ({len(cnv_by_cell)} cell keys)")
    print(f"Genes: {','.join(GENES_OF_INTEREST)}")
    print(f"Neutral HMM state (missing subcluster/gene in pred): {NEUTRAL_HMM_STATE}")
    print(f"Gain (all genes) if HMM state >= {GAIN_MIN_HMM_STATE}")


if __name__ == "__main__":
    main()
