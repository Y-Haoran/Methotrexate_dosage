#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


BASE_DIR = Path(__file__).resolve().parents[1]
FIGURE_DIR = BASE_DIR / "docs" / "figures"


def ensure_output_dir() -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)


def friendly_family_label(row: pd.Series) -> str:
    polymer = str(row["polymer_key"]).upper() if str(row["polymer_key"]).lower() != "hpc" else "HPC"
    cosolvent = str(row["cosolvent_system"]).replace("_", "/")
    return f"{polymer} + {cosolvent}"


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
    fig, ax = plt.subplots(figsize=(13.0, 5.8), dpi=220)
    fig.patch.set_facecolor("#f8f5ef")
    ax.set_facecolor("#f8f5ef")
    ax.axis("off")

    boxes = [
        {
            "xy": (0.03, 0.24),
            "w": 0.26,
            "h": 0.52,
            "fc": "#fff3d6",
            "ec": "#8a5a00",
            "title": "Layer 1\nUser Input",
            "text": "API name\nAPI SMILES\ncontext\ncharge-state options",
        },
        {
            "xy": (0.37, 0.18),
            "w": 0.26,
            "h": 0.64,
            "fc": "#e8f4ea",
            "ec": "#1f6f4a",
            "title": "Layer 2\nCandidate Generation",
            "text": "RDKit descriptors\nrule-based API classes\nfamily prior table\ncandidate matrix generation",
        },
        {
            "xy": (0.71, 0.14),
            "w": 0.26,
            "h": 0.72,
            "fc": "#e9f0fb",
            "ec": "#205493",
            "title": "Layer 3\nMechanistic Screening",
            "text": "FAIRChem checkpoint\nASE + FIRE relaxation\ntask_name = omol\nreplicate scoring\naggregated shortlist",
        },
    ]

    for box in boxes:
        rect = plt.Rectangle(box["xy"], box["w"], box["h"], facecolor=box["fc"], edgecolor=box["ec"], linewidth=2.0)
        ax.add_patch(rect)
        x, y = box["xy"]
        ax.text(x + box["w"] / 2, y + box["h"] - 0.08, box["title"], ha="center", va="top", fontsize=16, weight="bold", color="#1f1f1f")
        ax.text(x + box["w"] / 2, y + box["h"] / 2 - 0.02, box["text"], ha="center", va="center", fontsize=12, color="#333333", linespacing=1.5)

    arrow_style = dict(arrowstyle="-|>", linewidth=2.2, color="#4b5563", shrinkA=0, shrinkB=0)
    ax.annotate("", xy=(0.37, 0.49), xytext=(0.29, 0.49), arrowprops=arrow_style)
    ax.annotate("", xy=(0.71, 0.49), xytext=(0.63, 0.49), arrowprops=arrow_style)

    ax.text(0.50, 0.94, "Platform Model Stack", ha="center", va="center", fontsize=20, weight="bold", color="#1f1f1f")
    ax.text(
        0.50,
        0.07,
        "Why this works: the platform first narrows plausible formulation families from structure, then ranks them with local mechanistic evidence.",
        ha="center",
        va="center",
        fontsize=11,
        color="#444444",
    )

    fig.tight_layout()
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
