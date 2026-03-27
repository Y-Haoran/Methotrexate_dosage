#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ase import Atoms
from ase.data import covalent_radii
from ase.data.colors import jmol_colors
from ase.io import read
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch


BASE_DIR = Path(__file__).resolve().parents[1]
FIGURE_DIR = BASE_DIR / "docs" / "figures"


def ensure_output_dir() -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)


def friendly_family_label(row: pd.Series) -> str:
    polymer = str(row["polymer_key"]).upper() if str(row["polymer_key"]).lower() != "hpc" else "HPC"
    cosolvent = str(row["cosolvent_system"]).replace("_", "/")
    return f"{polymer} + {cosolvent}"


def compact_family_label(polymer_key: str, cosolvent_system: str) -> str:
    polymer = polymer_key.upper() if polymer_key.lower() != "hpc" else "HPC"
    mapping = {
        "water": "H2O",
        "ethanol_water": "EtOH/H2O",
        "glycerol_water": "Gly/H2O",
    }
    cosolvent = mapping.get(cosolvent_system, cosolvent_system.replace("_", "/"))
    return f"{polymer}/{cosolvent}"


def render_methotrexate_family_ranking() -> None:
    source = BASE_DIR / "results" / "methotrexate_family_recommendation_v1" / "family_recommendations.csv"
    df = pd.read_csv(source)
    subsets = [
        ("neutral", "Neutral State", "#0b6e4f"),
        ("monoanion", "Monoanionic State", "#b35c1e"),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.4), dpi=220, sharex=False)
    fig.patch.set_facecolor("#f7f4ea")

    for ax, (state_id, title, color) in zip(axes, subsets):
        subset = df[df["state_id"] == state_id].nsmallest(999999, "state_rank").head(4).copy()
        subset["label"] = subset.apply(friendly_family_label, axis=1)
        subset = subset.sort_values("recommendation_score", ascending=True)

        ax.set_facecolor("#fffaf0")
        ax.barh(subset["label"], subset["recommendation_score"], color=color, edgecolor="#1f1f1f", linewidth=0.6)
        for _, row in subset.iterrows():
            ax.text(
                row["recommendation_score"] + 0.03,
                row["label"],
                f"{row['recommendation_score']:.2f}",
                va="center",
                ha="left",
                fontsize=9,
                color="#1f1f1f",
            )
        ax.set_title(title, fontsize=13, weight="bold", color="#1f1f1f")
        ax.set_xlabel("Recommendation Score", fontsize=10)
        ax.grid(axis="x", alpha=0.25, linestyle="--")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color("#555555")
        ax.spines["bottom"].set_color("#555555")
        ax.tick_params(labelsize=9)

    fig.suptitle("Methotrexate Example: Recommended Formulation Families", fontsize=16, weight="bold", y=0.98)
    fig.text(
        0.5,
        0.01,
        "Archived recommendation output from the current platform. Higher score means a stronger family-level prior plus descriptor match.",
        ha="center",
        fontsize=10,
        color="#444444",
    )
    fig.tight_layout(rect=(0.02, 0.05, 1, 0.94))
    fig.savefig(FIGURE_DIR / "methotrexate_family_ranking.png", bbox_inches="tight")
    plt.close(fig)


def render_paracetamol_relaxed_demo() -> None:
    source = BASE_DIR / "results" / "mechanistic_screen_relaxed_gpu_v3" / "mechanistic_screen_results_neutral.csv"
    df = pd.read_csv(source)
    subset = df[df["api_name"] == "paracetamol"].copy()
    subset = subset[subset["candidate_id"].isin(["screen_001", "screen_002", "screen_003", "screen_004"])].copy()
    subset["label"] = subset.apply(
        lambda row: f"{str(row['polymer_key']).upper() if str(row['polymer_key']).lower() != 'hpc' else 'HPC'} + {str(row['cosolvent_system']).replace('_', '/')}",
        axis=1,
    )
    subset = subset.sort_values("shortlist_rank_neutral", ascending=False)

    fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.4), dpi=220)
    fig.patch.set_facecolor("#f3f4f6")

    colors_pre = "#9ca3af"
    colors_post = "#0f766e"
    accent = "#7c2d12"

    axes[0].set_facecolor("#ffffff")
    axes[0].barh(subset["label"], subset["pre_screening_score_proxy"], color=colors_pre, label="Pre-relaxation")
    axes[0].barh(subset["label"], subset["screening_score_proxy"], color=colors_post, alpha=0.9, label="Post-relaxation")
    axes[0].set_title("Paracetamol Demo: Relaxation Changes the Ranking", fontsize=13, weight="bold")
    axes[0].set_xlabel("Score Proxy (higher is better)")
    axes[0].grid(axis="x", alpha=0.25, linestyle="--")
    axes[0].legend(frameon=False, loc="lower right")
    axes[0].spines["top"].set_visible(False)
    axes[0].spines["right"].set_visible(False)

    axes[1].set_facecolor("#ffffff")
    axes[1].barh(subset["label"], subset["score_delta_relaxation"], color=accent)
    for _, row in subset.iterrows():
        axes[1].text(
            row["score_delta_relaxation"] + 2.0,
            row["label"],
            f"+{row['score_delta_relaxation']:.1f}",
            va="center",
            ha="left",
            fontsize=9,
            color="#1f1f1f",
        )
    axes[1].set_title("Relaxation Shift", fontsize=13, weight="bold")
    axes[1].set_xlabel("Delta From Pre to Post")
    axes[1].grid(axis="x", alpha=0.25, linestyle="--")
    axes[1].spines["top"].set_visible(False)
    axes[1].spines["right"].set_visible(False)

    fig.suptitle("Mechanistic Screening Example: Archived Relaxed GPU Result", fontsize=16, weight="bold", y=0.98)
    fig.text(
        0.5,
        0.01,
        "This is not a mockup. It is generated from the tracked relaxed GPU archive and shows why post-relaxation ranking matters.",
        ha="center",
        fontsize=10,
        color="#444444",
    )
    fig.tight_layout(rect=(0.02, 0.05, 1, 0.94))
    fig.savefig(FIGURE_DIR / "paracetamol_relaxed_demo.png", bbox_inches="tight")
    plt.close(fig)


def render_platform_model_stack() -> None:
    fig = plt.figure(figsize=(14.0, 8.2), dpi=220)
    fig.patch.set_facecolor("#f8f5ef")

    canvas = fig.add_axes([0, 0, 1, 1])
    canvas.set_facecolor("#f8f5ef")
    canvas.axis("off")

    def draw_container(x: float, y: float, w: float, h: float, face: str, edge: str, title: str, subtitle: str) -> None:
        patch = FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.012,rounding_size=0.025",
            facecolor=face,
            edgecolor=edge,
            linewidth=1.8,
        )
        canvas.add_patch(patch)
        canvas.text(x + 0.02, y + h - 0.035, title, fontsize=14, weight="bold", color="#1f1f1f", va="top")
        canvas.text(x + 0.02, y + h - 0.078, subtitle, fontsize=8.9, color="#555555", va="top")

    def normalize_positions(positions: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
        centered = positions - positions.mean(axis=0, keepdims=True)
        scale = max(np.ptp(centered[:, axis]) for axis in range(3))
        scale = scale if scale > 0 else 1.0
        return centered / scale, centered, scale

    def draw_atoms(ax, atoms: Atoms, title: str, subtitle: str, elev: float, azim: float, offsets: list[np.ndarray] | None = None) -> None:
        if offsets is None:
            offsets = [np.zeros(3)]
        all_positions = []
        all_numbers = []
        all_original = []
        for offset in offsets:
            positions = atoms.get_positions() + offset[None, :]
            all_positions.append(positions)
            all_original.append(positions)
            all_numbers.extend(atoms.get_atomic_numbers().tolist())
        merged_positions = np.vstack(all_positions)
        normalized, original_centered, _ = normalize_positions(merged_positions)
        numbers = np.asarray(all_numbers, dtype=int)
        colors = jmol_colors[numbers]
        sizes = np.asarray([max(covalent_radii[number], 0.35) * 240 for number in numbers], dtype=float)

        for i in range(len(numbers)):
            for j in range(i + 1, len(numbers)):
                threshold = 1.18 * (covalent_radii[numbers[i]] + covalent_radii[numbers[j]])
                distance = np.linalg.norm(original_centered[i] - original_centered[j])
                if 0.1 < distance < threshold:
                    ax.plot(
                        [normalized[i, 0], normalized[j, 0]],
                        [normalized[i, 1], normalized[j, 1]],
                        [normalized[i, 2], normalized[j, 2]],
                        color="#7c7c7c",
                        linewidth=1.1,
                        alpha=0.9,
                    )

        ax.scatter(normalized[:, 0], normalized[:, 1], normalized[:, 2], s=sizes, c=colors, edgecolors="#333333", linewidths=0.25)
        ax.view_init(elev=elev, azim=azim)
        ax.set_axis_off()
        ax.set_box_aspect((1, 1, 1))
        ax.text2D(0.02, 0.94, title, transform=ax.transAxes, fontsize=11, weight="bold", color="#1f1f1f")
        ax.text2D(0.02, 0.86, subtitle, transform=ax.transAxes, fontsize=8.5, color="#555555")

    def draw_stage_badge(x: float, y: float, label: str, face: str, edge: str) -> None:
        badge = Circle((x, y), 0.018, facecolor=face, edgecolor=edge, linewidth=1.4)
        canvas.add_patch(badge)
        canvas.text(x, y, label, ha="center", va="center", fontsize=10.5, weight="bold", color="#1f1f1f")

    canvas.text(0.5, 0.968, "AI-Driven Platform: From Molecular Inputs To Formulation Shortlist", ha="center", va="top", fontsize=20.5, weight="bold", color="#1f1f1f")
    canvas.text(
        0.5,
        0.935,
        "Real 3D building blocks from the repo feed the recommendation layer and the checkpoint-based mechanistic screening layer.",
        ha="center",
        va="top",
        fontsize=9.8,
        color="#4b5563",
    )

    draw_container(0.03, 0.12, 0.26, 0.78, "#fff7ea", "#b7791f", "Input Chemistry", "Representative 3D structures used by the platform")
    draw_container(0.34, 0.12, 0.29, 0.78, "#eef6ff", "#2b6cb0", "Model Architecture", "Descriptor layer + mechanistic layer")
    draw_container(0.68, 0.12, 0.29, 0.78, "#f3f8f1", "#2f855a", "Platform Outputs", "What the user gets back")
    draw_stage_badge(0.05, 0.88, "1", "#fff3d6", "#b7791f")
    draw_stage_badge(0.36, 0.88, "2", "#e6f0ff", "#2b6cb0")
    draw_stage_badge(0.70, 0.88, "3", "#e8f5e9", "#2f855a")

    api_ax = fig.add_axes([0.058, 0.57, 0.19, 0.20], projection="3d")
    polymer_ax = fig.add_axes([0.058, 0.35, 0.19, 0.20], projection="3d")
    solvent_ax = fig.add_axes([0.058, 0.14, 0.19, 0.17], projection="3d")

    draw_atoms(api_ax, read(BASE_DIR / "SMILES_3D" / "api" / "methotrexate.xyz"), "API example", "methotrexate 3D structure", elev=18, azim=34)
    draw_atoms(polymer_ax, read(BASE_DIR / "SMILES_3D" / "polymer" / "pvp.xyz"), "Polymer example", "PVP fragment proxy", elev=20, azim=-48)

    water = read(BASE_DIR / "SMILES_3D" / "solvent" / "water.xyz")
    ethanol = read(BASE_DIR / "SMILES_3D" / "solvent" / "ethanol.xyz")
    glycerol = read(BASE_DIR / "SMILES_3D" / "solvent" / "glycerol.xyz")
    mixed = water.copy()
    mixed += ethanol.copy()
    mixed += glycerol.copy()
    draw_atoms(
        solvent_ax,
        mixed,
        "Solvent / co-solvent",
        "water + ethanol + glycerol",
        elev=18,
        azim=38,
    )

    model_ax = fig.add_axes([0.365, 0.17, 0.24, 0.68])
    model_ax.axis("off")

    model_boxes = [
        (0.05, 0.73, 0.90, 0.18, "#fff8e1", "#b7791f", "1. Descriptors", "RDKit\ncharge states"),
        (0.05, 0.46, 0.90, 0.18, "#eefbf3", "#2f855a", "2. Candidates", "family priors\nbuild systems"),
        (0.05, 0.19, 0.90, 0.20, "#edf2ff", "#2b6cb0", "3. Mechanics", "FAIRChem + ASE\nrelax + score"),
    ]
    for x, y, w, h, fc, ec, title, text in model_boxes:
        patch = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.04", facecolor=fc, edgecolor=ec, linewidth=1.8)
        model_ax.add_patch(patch)
        model_ax.text(x + 0.04, y + h - 0.04, title, fontsize=12, weight="bold", color="#1f1f1f", va="top")
        model_ax.text(x + 0.04, y + h - 0.11, text, fontsize=10.0, color="#374151", va="top", linespacing=1.35)
    model_ax.add_patch(FancyArrowPatch((0.50, 0.71), (0.50, 0.64), arrowstyle="Simple,head_length=8,head_width=8,tail_width=0.7", mutation_scale=1.0, linewidth=0.0, color="#7a8796", alpha=0.8))
    model_ax.add_patch(FancyArrowPatch((0.50, 0.44), (0.50, 0.40), arrowstyle="Simple,head_length=8,head_width=8,tail_width=0.7", mutation_scale=1.0, linewidth=0.0, color="#7a8796", alpha=0.8))

    ranking_ax = fig.add_axes([0.725, 0.515, 0.205, 0.225])
    ranking_ax.set_facecolor("#ffffff")
    ranking_df = pd.read_csv(BASE_DIR / "results" / "methotrexate_family_recommendation_v1" / "family_recommendations.csv")
    ranking_df = ranking_df[
        ((ranking_df["state_id"] == "neutral") & (ranking_df["state_rank"] <= 2))
        | ((ranking_df["state_id"] == "monoanion") & (ranking_df["state_rank"] <= 2))
    ].copy()
    ranking_df["label"] = ranking_df.apply(
        lambda row: f"{'N' if row['state_id'] == 'neutral' else 'I'} | {compact_family_label(str(row['polymer_key']), str(row['cosolvent_system']))}",
        axis=1,
    )
    ranking_df = ranking_df.sort_values("recommendation_score", ascending=False).reset_index(drop=True)
    positions = np.arange(len(ranking_df))
    bar_colors = ["#b35c1e" if state_id == "monoanion" else "#0b6e4f" for state_id in ranking_df["state_id"]]
    ranking_ax.barh(positions, ranking_df["recommendation_score"], color=bar_colors, edgecolor="#1f1f1f", linewidth=0.5)
    xmax = float(ranking_df["recommendation_score"].max()) * 1.18
    ranking_ax.set_xlim(0, xmax)
    ranking_ax.set_yticks([])
    ranking_ax.invert_yaxis()
    for idx, row in ranking_df.iterrows():
        ranking_ax.text(
            xmax * 0.03,
            positions[idx],
            str(row["label"]),
            va="center",
            ha="left",
            fontsize=8.2,
            color="#1f1f1f",
            bbox={"facecolor": "#ffffff", "edgecolor": "none", "alpha": 0.84, "pad": 0.12},
        )
        ranking_ax.text(
            float(row["recommendation_score"]) + xmax * 0.015,
            positions[idx],
            f"{float(row['recommendation_score']):.1f}",
            va="center",
            ha="left",
            fontsize=8.2,
            color="#1f1f1f",
        )
    ranking_ax.set_title("Shortlist", fontsize=10.8, weight="bold", pad=6)
    ranking_ax.grid(axis="x", alpha=0.22, linestyle="--")
    ranking_ax.tick_params(axis="x", labelsize=8.0)
    ranking_ax.spines["top"].set_visible(False)
    ranking_ax.spines["right"].set_visible(False)
    ranking_ax.spines["left"].set_visible(False)

    cluster_ax = fig.add_axes([0.725, 0.31, 0.205, 0.16])
    cluster_ax.set_facecolor("#ffffff")
    cluster_ax.axis("off")
    cluster_ax.set_title("Interaction map", fontsize=11.0, weight="bold", pad=4)
    nodes = {
        "API": (0.20, 0.58),
        "Polymer": (0.52, 0.78),
        "Solvent": (0.80, 0.56),
        "Ion": (0.50, 0.24),
    }
    edges = [
        ("API", "Polymer", "#0f766e", "bind"),
        ("API", "Solvent", "#b7791f", "compete"),
        ("API", "Ion", "#7c2d12", "ion"),
        ("Polymer", "Solvent", "#2563eb", "solvate"),
    ]
    for start, end, color, label in edges:
        x1, y1 = nodes[start]
        x2, y2 = nodes[end]
        cluster_ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-", linewidth=2.2, color=color, alpha=0.9))
        cluster_ax.text((x1 + x2) / 2, (y1 + y2) / 2 + 0.03, label, fontsize=7.1, color=color, ha="center")
    for name, (x, y) in nodes.items():
        circle = Circle((x, y), 0.08, facecolor="#f8fafc", edgecolor="#334155", linewidth=1.4)
        cluster_ax.add_patch(circle)
        cluster_ax.text(x, y, name, ha="center", va="center", fontsize=9.5, weight="bold", color="#1f1f1f")

    target_ax = fig.add_axes([0.725, 0.13, 0.205, 0.12])
    target_ax.set_facecolor("#ffffff")
    target_ax.axis("off")
    target_ax.text(0.02, 0.92, "Output bundle", fontsize=11.0, weight="bold", color="#1f1f1f", va="top")
    target_ax.text(
        0.02,
        0.69,
        "• ranked shortlist\n• interaction map\n• DFT / lab plan",
        fontsize=9.4,
        color="#374151",
        va="top",
        linespacing=1.5,
    )

    canvas.text(
        0.50,
        0.06,
        "Principle: the platform starts from real molecular structure, narrows plausible formulation families, then uses local mechanistic evidence to rank what should be tested first.",
        ha="center",
        va="center",
        fontsize=11,
        color="#444444",
    )

    fig.savefig(FIGURE_DIR / "platform_model_stack.png", bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    ensure_output_dir()
    render_methotrexate_family_ranking()
    render_paracetamol_relaxed_demo()
    render_platform_model_stack()
    print(str(FIGURE_DIR))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
