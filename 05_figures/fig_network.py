"""
fig_network.py — Co-expression network integrating WGCNA hub genes, DEGs, and
virulence-associated genes (Figure 5).

Nodes:
  fill colour → expression status (DEG up / DEG down / hub-only / co-expressed)
  border       → functional role (virulence phase, hub+DEG, hub-only)
  size         → hub: ∝ kME;  DEG: ∝ |log₂FC|;  other: fixed

Edges:
  Pearson |r| ≥ --corr-threshold, within the same WGCNA module only.
  Red = positive co-expression; blue = negative.

Layout:
  Each module is positioned independently (spring_layout), then placed at
  pre-defined centers and separated by a vectorised post-repulsion pass.
  The dimgrey module uses unweighted layout with stronger repulsion (k=1.4) to
  prevent the densely-connected cluster from collapsing.

Usage
-----
python fig_network.py \\
    --wgcna-dir    results/WGCNA \\
    --deseq-dir    results/DESeq2_results \\
    --dictionary   data/gene_dictionary.tsv \\
    --virulence-table data/virulence_table.csv \\
    --output       figures/fig_5_network.tiff

Input files expected in --wgcna-dir:
    module_assignments.csv   (from run_wgcna.py)
    kME.csv                  (from run_wgcna.py)
    expr_log.csv             (from run_wgcna.py — genes × samples)

Note on gene identifiers
------------------------
The network was built on concatenated 9a5c + Temecula1 IMG IDs stored in the
'gene' column of module_assignments.csv (e.g. "XF9a_12345|XFTem_67890").
Labels are converted to PD-locus-tag format (e.g. "PD0422") using the cross-
strain dictionary.
"""

import argparse
import glob
import os
import re
import sys
import warnings
from collections import defaultdict
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

warnings.filterwarnings("ignore")

try:
    from adjustText import adjust_text
    HAS_ADJUST_TEXT = True
except ImportError:
    HAS_ADJUST_TEXT = False
    print("[WARN] adjustText not installed — labels will not be auto-adjusted. "
          "Run: pip install adjustText")


# ---------------------------------------------------------------------------
# Visual constants
# ---------------------------------------------------------------------------

MOD_FILL = {
    "dimgrey":  "#707070",
    "darkgrey": "#a0a0a0",
    "silver":   "#c8c8e8",
}
COL_DEG_UP     = "#e05050"
COL_DEG_DOWN   = "#5088e0"
COL_DEFAULT    = "#bbbbbb"
BORDER_VIR_MOBILE  = "#FFD700"
BORDER_VIR_SESSILE = "#44bb44"
BORDER_HUB_DEG     = "#cc44cc"
BORDER_HUB_ONLY    = "#333333"
BORDER_DEFAULT     = "#999999"
COL_EDGE_POS   = "#cc3333"
COL_EDGE_NEG   = "#3355cc"

MODULE_LAYOUT_K = {"dimgrey": 1.4, "darkgrey": 0.35, "silver": 0.35}
MODULE_SCALE    = {"dimgrey": 0.50, "darkgrey": 0.27, "silver": 0.27}
MODULE_CENTERS  = {
    "dimgrey":  np.array([-0.2,   0.80]),
    "darkgrey": np.array([-1.05, -0.55]),
    "silver":   np.array([ 0.75, -1.45]),
}


# ---------------------------------------------------------------------------
# 1. Module assignments and hub genes
# ---------------------------------------------------------------------------

def load_wgcna(wgcna_dir: Path) -> tuple:
    mg_df  = pd.read_csv(wgcna_dir / "module_assignments.csv")
    kme_df = pd.read_csv(wgcna_dir / "kME.csv", index_col=0)
    gene_to_mod = dict(zip(mg_df["gene"], mg_df["module"]))

    # Map 9a5c IMG IDs to the full gene string (used to match DESeq2 gene IDs)
    img9_in_module: dict[str, str] = {}
    for _, row in mg_df.iterrows():
        m9 = re.search(r"(XF9a_\d+)", str(row["gene"]))
        if m9:
            img9_in_module[m9.group(1)] = row["gene"]

    hub_set: set  = set()
    hub_kme: dict = {}
    for mod in mg_df["module"].unique():
        me_col = f"ME{mod}"
        if me_col not in kme_df.columns:
            continue
        genes_mod = mg_df.loc[mg_df["module"] == mod, "gene"].tolist()
        valid     = [g for g in genes_mod if g in kme_df.index]
        if not valid:
            continue
        for g, v in kme_df.loc[valid, me_col].nlargest(10).items():
            hub_set.add(g)
            hub_kme[g] = v

    print(f"[INFO] {len(mg_df)} genes in {mg_df['module'].nunique()} modules; "
          f"{len(hub_set)} hub genes")
    return mg_df, kme_df, gene_to_mod, img9_in_module, hub_set, hub_kme


# ---------------------------------------------------------------------------
# 2. Gene ID → PD locus tag conversion
# ---------------------------------------------------------------------------

def build_id_maps(dict_path: Path) -> tuple:
    sep = "\t" if dict_path.suffix in {".tsv", ".txt"} else ","
    dct = pd.read_csv(dict_path, sep=sep)

    def split_tags(raw: str) -> list:
        return [t.strip().strip('"')
                for t in re.split(r'[\t\s;,]+', str(raw)) if t.strip()]

    xftem_to_pd: dict[str, str] = {}
    for _, row in dct.iterrows():
        imgT = str(row.get("Temecula1_IMG_ID", "")).strip()
        if not imgT or imgT == "nan":
            continue
        for tag in split_tags(row.get("Temecula1_old_locus_tags", "")):
            t = tag.upper()
            if re.match(r"PD_?\d+", t):
                xftem_to_pd[imgT] = t.replace("_", "")
                break

    def gene_to_pdlabel(gene_concat: str) -> str:
        mT = re.search(r"(XFTem_\d+)", gene_concat)
        if mT:
            pd_tag = xftem_to_pd.get(mT.group(1))
            if pd_tag:
                return pd_tag
        m9 = re.search(r"(XF9a_\d+)", gene_concat)
        return m9.group(1) if m9 else gene_concat

    return xftem_to_pd, gene_to_pdlabel


# ---------------------------------------------------------------------------
# 3. Virulence gene information
# ---------------------------------------------------------------------------

def load_virulence(vir_path: Path, xftem_to_pd: dict) -> dict:
    """Return dict: PD_tag → phase_label (1='mobile', 2='sessile')."""
    sep = "\t" if vir_path.suffix in {".tsv", ".txt"} else ","
    try:
        vir = pd.read_csv(vir_path, sep=sep, engine="python", encoding="latin-1")
    except Exception:
        vir = pd.read_csv(vir_path, sep=None, engine="python", encoding="latin-1")

    phase_col = next((c for c in vir.columns if "phase" in c.lower()), None)
    id_col    = next((c for c in vir.columns
                      if "gene_id" in c.lower() or "gene id" in c.lower()), None)
    if not phase_col or not id_col:
        print("[WARN] Could not locate phase/id columns in virulence table.")
        return {}

    phase_label = {1: "mobile", 2: "sessile"}
    return {
        str(row[id_col]).strip().upper(): phase_label[int(row[phase_col])]
        for _, row in vir.iterrows()
        if pd.notna(row[phase_col]) and int(row[phase_col]) in (1, 2)
    }


def virulence_info(
    gene_concat: str,
    xftem_to_pd: dict,
    vir_phase: dict,
) -> tuple:
    mT = re.search(r"(XFTem_\d+)", gene_concat)
    if not mT:
        return False, None
    pd_tag = xftem_to_pd.get(mT.group(1))
    if pd_tag and pd_tag in vir_phase:
        return True, vir_phase[pd_tag]
    return False, None


# ---------------------------------------------------------------------------
# 4. DEG summary across comparisons
# ---------------------------------------------------------------------------

def load_degs(
    deseq_dir: Path,
    img9_in_module: dict,
    padj_threshold: float,
    lfc_threshold: float,
) -> pd.DataFrame:
    records = []
    for comp_dir in deseq_dir.iterdir():
        if not comp_dir.is_dir():
            continue
        for csv_file in comp_dir.glob("*.csv"):
            if "GO_enrichment" in csv_file.name:
                continue
            try:
                df = pd.read_csv(csv_file)
                if not {"gene_id", "log2FoldChange", "padj"}.issubset(df.columns):
                    continue
                df["comparison"] = comp_dir.name
                records.append(df)
            except Exception as e:
                print(f"[WARN] {csv_file.name}: {e}")

    if not records:
        print("[WARN] No DESeq2 result files found.")
        return pd.DataFrame()

    all_degs = pd.concat(records, ignore_index=True).dropna(
        subset=["padj", "log2FoldChange"]
    )
    sig = all_degs[
        (all_degs["padj"] < padj_threshold) &
        (all_degs["log2FoldChange"].abs() >= lfc_threshold)
    ].copy()

    sig["gene_concat"] = sig["gene_id"].map(img9_in_module)
    sig = sig.dropna(subset=["gene_concat"])

    def summarize(group):
        n_up, n_down = (group["log2FoldChange"] > 0).sum(), (group["log2FoldChange"] < 0).sum()
        n = len(group)
        direction = "up" if n_up == n else ("down" if n_down == n else "mixed")
        return pd.Series({
            "log2FoldChange":      group["log2FoldChange"].mean(),
            "padj_representative": group["padj"].max(),
            "DEG_direction":       direction,
            "n_comparisons":       n,
            "n_up":                n_up,
            "n_down":              n_down,
            "comparisons":         "; ".join(sorted(group["comparison"].unique())),
        })

    grouped = sig.groupby("gene_concat").apply(summarize).reset_index()
    consistent = grouped[grouped["DEG_direction"] != "mixed"].copy()
    print(f"[INFO] {len(grouped)} DEGs mapped to WGCNA; "
          f"{len(consistent)} consistent (up/down), "
          f"{(grouped['DEG_direction']=='mixed').sum()} mixed (excluded)")
    return consistent.set_index("gene_concat")


# ---------------------------------------------------------------------------
# 5. Correlation graph
# ---------------------------------------------------------------------------

def build_graph(
    em: pd.DataFrame,            # samples × genes
    hub_set: set,
    deg_summary: pd.DataFrame,
    gene_to_mod: dict,
    corr_threshold: float,
    min_component: int,
) -> nx.Graph:
    relevant = list((hub_set | set(deg_summary.index)) & set(em.columns))
    print(f"[INFO] {len(relevant)} relevant genes for network")

    corr_mat = em[relevant].corr()

    edges = []
    for i in range(len(relevant)):
        for j in range(i + 1, len(relevant)):
            g1, g2 = relevant[i], relevant[j]
            if gene_to_mod.get(g1) != gene_to_mod.get(g2):
                continue
            r = corr_mat.loc[g1, g2]
            if abs(r) >= corr_threshold:
                edges.append((g1, g2, float(r)))
    print(f"[INFO] {len(edges)} edges (|r| ≥ {corr_threshold})")

    G = nx.Graph()
    for g in relevant:
        G.add_node(g)
    for g1, g2, r in edges:
        G.add_edge(g1, g2, weight=abs(r), r=r)

    # Remove small components and isolates
    for comp in [c for c in nx.connected_components(G) if len(c) < min_component]:
        G.remove_nodes_from(comp)
    G.remove_nodes_from(list(nx.isolates(G)))
    print(f"[INFO] After filter: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    return G


# ---------------------------------------------------------------------------
# 6. Layout (per-module spring + post-repulsion)
# ---------------------------------------------------------------------------

def post_repulsion(
    pos: dict,
    nodes_by_mod: dict,
    min_dist: float = 0.055,
    iterations: int = 80,
) -> dict:
    pos = {n: np.array(p, dtype=float) for n, p in pos.items()}
    for mod_nodes in nodes_by_mod.values():
        if len(mod_nodes) < 2:
            continue
        coords = np.array([pos[n] for n in mod_nodes])
        for _ in range(iterations):
            diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
            dist = np.linalg.norm(diff, axis=-1)
            np.fill_diagonal(dist, np.inf)
            mask = dist < min_dist
            if not mask.any():
                break
            safe_dist   = np.where(dist < 1e-9, 1e-9, dist)
            force_mag   = np.where(mask, (min_dist - dist) / safe_dist, 0.0)
            push        = force_mag[:, :, np.newaxis] * diff
            coords     += push.sum(axis=1) * 0.4
        for i, n in enumerate(mod_nodes):
            pos[n] = coords[i]
    return pos


def compute_layout(G: nx.Graph, gene_to_mod: dict) -> dict:
    pos: dict = {}
    for mod, center in MODULE_CENTERS.items():
        nodes_mod = [n for n in G.nodes() if gene_to_mod.get(n) == mod]
        if not nodes_mod:
            continue
        subG   = G.subgraph(nodes_mod)
        k_val  = MODULE_LAYOUT_K.get(mod, 0.35)
        n_iter = 600 if mod == "dimgrey" else 200
        weight = None if mod == "dimgrey" else "weight"
        sub_pos = nx.spring_layout(subG, seed=42, k=k_val,
                                   weight=weight, iterations=n_iter)
        coords = np.array(list(sub_pos.values()))
        scale  = MODULE_SCALE.get(mod, 0.32)
        coords = (coords - coords.mean(axis=0)) / (coords.std(axis=0) + 1e-9) * scale
        for node, coord in zip(sub_pos.keys(), coords):
            pos[node] = coord + center

    rng = np.random.default_rng(0)
    for node in G.nodes():
        if node not in pos:
            pos[node] = rng.standard_normal(2) * 0.05

    nodes_by_mod: dict = defaultdict(list)
    for n in G.nodes():
        nodes_by_mod[gene_to_mod.get(n, "unknown")].append(n)

    return post_repulsion(pos, nodes_by_mod)


# ---------------------------------------------------------------------------
# 7. Visual attributes
# ---------------------------------------------------------------------------

def node_attributes(
    G: nx.Graph,
    gene_to_mod: dict,
    hub_set: set,
    hub_kme: dict,
    deg_summary: pd.DataFrame,
    xftem_to_pd: dict,
    vir_phase: dict,
    white_bg: bool,
) -> tuple:
    fills, borders, bwidths, sizes = [], [], [], []
    for gene in G.nodes():
        mod      = gene_to_mod.get(gene, "unknown")
        is_hub   = gene in hub_set
        is_deg   = gene in deg_summary.index
        vir_flag, vir_ph = virulence_info(gene, xftem_to_pd, vir_phase)

        # Fill
        if is_hub and not is_deg:
            fill = MOD_FILL.get(mod, COL_DEFAULT)
        elif is_deg:
            fill = COL_DEG_UP if deg_summary.loc[gene, "log2FoldChange"] > 0 else COL_DEG_DOWN
        else:
            fill = COL_DEFAULT

        # Border
        if vir_flag:
            border = BORDER_VIR_MOBILE if vir_ph == "mobile" else BORDER_VIR_SESSILE
            bw     = 3.5
        elif is_hub and is_deg:
            border, bw = BORDER_HUB_DEG, 2.5
        elif is_hub:
            border = BORDER_HUB_ONLY if white_bg else "white"
            bw     = 1.8
        else:
            border, bw = BORDER_DEFAULT, 0.5

        # Size
        if is_hub:
            size = 180 + hub_kme.get(gene, 0) * 100
        elif is_deg:
            lfc  = deg_summary.loc[gene, "log2FoldChange"]
            size = 80 + min(abs(lfc) * 25, 70)
        elif vir_flag:
            size = 90
        else:
            size = 50
        size *= 1.5   # 50% larger for visibility

        fills.append(fill);  borders.append(border)
        bwidths.append(bw);  sizes.append(size)

    return fills, borders, bwidths, sizes


# ---------------------------------------------------------------------------
# 8. Draw
# ---------------------------------------------------------------------------

def draw(
    G: nx.Graph,
    pos: dict,
    fills, borders, bwidths, sizes,
    gene_to_mod: dict,
    hub_set: set,
    hub_kme: dict,
    deg_summary: pd.DataFrame,
    kme_df: pd.DataFrame,
    mg_df: pd.DataFrame,
    xftem_to_pd: dict,
    vir_phase: dict,
    gene_to_pdlabel,
    corr_threshold: float,
    white_bg: bool,
    output: Path,
):
    bg    = "white"   if white_bg else "#1a1a2e"
    tc    = "black"   if white_bg else "white"

    fig, ax = plt.subplots(figsize=(19, 16), facecolor=bg)
    ax.set_facecolor(bg)
    ax.axis("off")

    pos_edges = [(u, v) for u, v, d in G.edges(data=True) if d["r"] > 0]
    neg_edges = [(u, v) for u, v, d in G.edges(data=True) if d["r"] < 0]
    nx.draw_networkx_edges(G, pos, edgelist=pos_edges, ax=ax,
                           edge_color=COL_EDGE_POS, alpha=0.12, width=0.35)
    nx.draw_networkx_edges(G, pos, edgelist=neg_edges, ax=ax,
                           edge_color=COL_EDGE_NEG, alpha=0.12, width=0.35)
    nx.draw_networkx_nodes(G, pos, ax=ax,
                           node_color=fills, node_size=sizes,
                           edgecolors=borders, linewidths=bwidths, alpha=0.93)

    # Labels: top-10 hubs per module + virulence genes
    hub_labels: dict = {}
    for mod in mg_df["module"].unique():
        me_col = f"ME{mod}"
        if me_col not in kme_df.columns:
            continue
        genes_mod = [g for g in G.nodes() if gene_to_mod.get(g) == mod and g in hub_set]
        valid     = [g for g in genes_mod if g in kme_df.index]
        for g in kme_df.loc[valid, me_col].nlargest(10).index:
            hub_labels[g] = gene_to_pdlabel(g)

    vir_labels = {
        g: gene_to_pdlabel(g)
        for g in G.nodes() if virulence_info(g, xftem_to_pd, vir_phase)[0]
    }
    all_labels = {**hub_labels, **vir_labels}

    texts = []
    for gene, label in all_labels.items():
        x, y = pos[gene]
        t = ax.text(x, y, label, fontsize=5.0, color=tc,
                    fontweight="bold", ha="center", va="center")
        texts.append(t)

    if HAS_ADJUST_TEXT and texts:
        adjust_text(
            texts,
            x=[pos[g][0] for g in all_labels],
            y=[pos[g][1] for g in all_labels],
            ax=ax,
            expand=(1.4, 1.6),
            arrowprops=dict(arrowstyle="simple", color="#888888",
                            lw=0.5, shrinkA=5, shrinkB=3, mutation_scale=6),
            force_text=(0.8, 1.0),
            force_points=(0.3, 0.5),
            avoid_self=True,
        )

    # Module title labels
    for mod, center in MODULE_CENTERS.items():
        nodes_mod = [n for n in G.nodes() if gene_to_mod.get(n) == mod]
        if not nodes_mod:
            continue
        top_y = max(pos[n][1] for n in nodes_mod) + 0.10
        ax.text(center[0], top_y, f"Module: {mod}",
                color=MOD_FILL.get(mod, tc), fontsize=12, fontweight="bold",
                ha="center", transform=ax.transData,
                bbox=dict(boxstyle="round,pad=0.25",
                          facecolor=bg, edgecolor="none", alpha=0.85))

    # Legend
    legend_elements = [
        mpatches.Patch(color="none",         label="— Node fill —"),
        mpatches.Patch(color=COL_DEG_UP,     label="DEG up-regulated"),
        mpatches.Patch(color=COL_DEG_DOWN,   label="DEG down-regulated"),
        mpatches.Patch(color=MOD_FILL["dimgrey"],  label="Hub only — dimgrey module"),
        mpatches.Patch(color=MOD_FILL["darkgrey"], label="Hub only — darkgrey module"),
        mpatches.Patch(color=MOD_FILL["silver"],   label="Hub only — silver module"),
        mpatches.Patch(color=COL_DEFAULT,    label="Co-expressed (other)"),
        mpatches.Patch(color="none",         label=" "),
        mpatches.Patch(color="none",         label="— Node border —"),
        mpatches.Patch(facecolor="none", edgecolor=BORDER_VIR_MOBILE,  linewidth=3.0,
                       label="Virulence — mobile (gold)"),
        mpatches.Patch(facecolor="none", edgecolor=BORDER_VIR_SESSILE, linewidth=3.0,
                       label="Virulence — sessile (green)"),
        mpatches.Patch(facecolor="none", edgecolor=BORDER_HUB_DEG,     linewidth=2.5,
                       label="Hub + DEG (purple)"),
        mpatches.Patch(facecolor="none", edgecolor=BORDER_HUB_ONLY,    linewidth=1.8,
                       label="Hub only (dark border)"),
        mpatches.Patch(color="none",         label=" "),
        mpatches.Patch(color="none",         label="— Edges —"),
        Line2D([0], [0], color=COL_EDGE_POS, lw=2,
               label=f"Positive co-expression (r ≥ {corr_threshold})"),
        Line2D([0], [0], color=COL_EDGE_NEG, lw=2,
               label=f"Negative co-expression (r ≤ −{corr_threshold})"),
        mpatches.Patch(color="none",         label=" "),
        mpatches.Patch(color="none",         label="— Node size —"),
        mpatches.Patch(color="none",         label="Hub: proportional to kME"),
        mpatches.Patch(color="none",         label="DEG: proportional to |log₂FC|"),
    ]
    leg = ax.legend(handles=legend_elements, loc="lower left",
                    bbox_to_anchor=(0.01, 0.5), fontsize=8.0,
                    framealpha=1.0, facecolor="white" if white_bg else "#22223b",
                    edgecolor="#aaaaaa", labelcolor=tc,
                    title="Legend", title_fontsize=8.5)
    leg.get_title().set_color(tc)
    for handle, text in zip(leg.legend_handles, leg.get_texts()):
        lbl = text.get_text()
        if lbl.startswith("—") or lbl.strip() in ("", "Hub: proportional to kME",
                                                   "DEG: proportional to |log₂FC|"):
            handle.set_visible(False)

    ax.set_title("Co-expression network: WGCNA hubs × DEGs × Virulence genes",
                 color=tc, fontsize=13, pad=15)
    plt.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=600, bbox_inches="tight", facecolor=bg)
    plt.close(fig)
    print(f"[INFO] Network saved to {output}")


# ---------------------------------------------------------------------------
# 9. Export node table
# ---------------------------------------------------------------------------

def export_node_table(
    G: nx.Graph,
    gene_to_mod: dict,
    hub_set: set,
    hub_kme: dict,
    deg_summary: pd.DataFrame,
    xftem_to_pd: dict,
    vir_phase: dict,
    gene_to_pdlabel,
    output: Path,
):
    rows = []
    for gene in G.nodes():
        mod      = gene_to_mod.get(gene, "unknown")
        is_hub   = gene in hub_set
        is_deg   = gene in deg_summary.index
        vir_flag, vir_ph = virulence_info(gene, xftem_to_pd, vir_phase)
        if is_deg:
            r = deg_summary.loc[gene]
            lfc, padj_r, dirn = r["log2FoldChange"], r["padj_representative"], r["DEG_direction"]
            n_c, n_u, n_d, comps = int(r["n_comparisons"]), int(r["n_up"]), int(r["n_down"]), r["comparisons"]
        else:
            lfc = padj_r = dirn = comps = ""
            n_c = n_u = n_d = ""
        rows.append({
            "gene_concat":       gene,
            "pd_label":          gene_to_pdlabel(gene),
            "module":            mod,
            "is_hub":            is_hub,
            "hub_kME":           round(hub_kme[gene], 4) if is_hub else "",
            "is_DEG":            is_deg,
            "DEG_direction":     dirn,
            "log2FC_mean":       round(lfc, 4) if lfc != "" else "",
            "padj_worst":        f"{padj_r:.2e}" if padj_r != "" else "",
            "n_comparisons_sig": n_c,
            "n_up":              n_u,
            "n_down":            n_d,
            "comparisons":       comps,
            "is_virulence":      vir_flag,
            "virulence_phase":   vir_ph or "",
            "network_degree":    G.degree(gene),
        })
    df = (pd.DataFrame(rows)
          .sort_values(["module", "is_hub", "is_virulence"],
                       ascending=[True, False, False]))
    node_table_path = output.with_name(output.stem + "_node_table.csv")
    df.to_csv(node_table_path, index=False)
    print(f"[INFO] Node table saved to {node_table_path} ({len(df)} nodes)")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Co-expression network: WGCNA hubs × DEGs × virulence genes.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--wgcna-dir",       required=True, type=Path,
                   help="Directory with run_wgcna.py outputs "
                        "(module_assignments.csv, kME.csv, expr_log.csv).")
    p.add_argument("--deseq-dir",       required=True, type=Path)
    p.add_argument("--dictionary",      required=True, type=Path)
    p.add_argument("--virulence-table", required=True, type=Path)
    p.add_argument("--corr-threshold",  type=float, default=0.85)
    p.add_argument("--min-component",   type=int,   default=4)
    p.add_argument("--padj",            type=float, default=0.05)
    p.add_argument("--lfc",             type=float, default=1.0)
    p.add_argument("--white-bg",        action="store_true", default=True)
    p.add_argument("--dark-bg",         dest="white_bg", action="store_false")
    p.add_argument("--output",          required=True, type=Path)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    # 1. WGCNA data
    mg_df, kme_df, gene_to_mod, img9_in_module, hub_set, hub_kme = \
        load_wgcna(args.wgcna_dir)

    # 2. ID maps
    xftem_to_pd, gene_to_pdlabel = build_id_maps(args.dictionary)

    # 3. Virulence
    vir_phase = load_virulence(args.virulence_table, xftem_to_pd)
    print(f"[INFO] {len(vir_phase)} virulence genes (mobile + sessile)")

    # 4. DEGs — img9_in_module is defined BEFORE this call
    deg_summary = load_degs(args.wgcna_dir.parent / "DESeq2_results"
                            if args.deseq_dir is None else args.deseq_dir,
                            img9_in_module, args.padj, args.lfc)

    # 5. Expression matrix (genes × samples from run_wgcna.py → transpose)
    expr_path = args.wgcna_dir / "expr_log.csv"
    if not expr_path.exists():
        sys.exit(f"[ERROR] expr_log.csv not found in {args.wgcna_dir}")
    expr = pd.read_csv(expr_path, index_col=0)
    # expr_log.csv is genes × samples; need samples × genes for .corr()
    em = expr.T if expr.shape[0] > expr.shape[1] else expr

    # 6. Graph
    G = build_graph(em, hub_set, deg_summary, gene_to_mod,
                    args.corr_threshold, args.min_component)

    # 7. Layout
    pos = compute_layout(G, gene_to_mod)

    # 8. Visual attributes
    fills, borders, bwidths, sizes = node_attributes(
        G, gene_to_mod, hub_set, hub_kme, deg_summary,
        xftem_to_pd, vir_phase, args.white_bg,
    )

    # 9. Draw
    draw(G, pos, fills, borders, bwidths, sizes,
         gene_to_mod, hub_set, hub_kme, deg_summary, kme_df, mg_df,
         xftem_to_pd, vir_phase, gene_to_pdlabel,
         args.corr_threshold, args.white_bg, args.output)

    # 10. Node table
    export_node_table(G, gene_to_mod, hub_set, hub_kme, deg_summary,
                      xftem_to_pd, vir_phase, gene_to_pdlabel, args.output)


if __name__ == "__main__":
    main()
