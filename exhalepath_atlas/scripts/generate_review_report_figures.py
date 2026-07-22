#!/usr/bin/env python3
"""Generate review-report figures for the ExhalePath real-data evaluation."""

from __future__ import annotations

from pathlib import Path

import json
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "runs" / "review_report" / "figures"
ART = Path("/opt/cursor/artifacts/exhalepath_report/figures")
FIG.mkdir(parents=True, exist_ok=True)
ART.mkdir(parents=True, exist_ok=True)

# Visual system: deep ocean slate + warm amber (not purple / cream-terracotta / default AI look)
C = {
    "bg": "#0f1419",
    "panel": "#1a222c",
    "ink": "#e8eef4",
    "muted": "#8b9aab",
    "accent": "#e8a54b",
    "teal": "#3db8a0",
    "blue": "#5b8fd9",
    "rose": "#d4676a",
    "line": "#2a3544",
    "good": "#3db8a0",
    "warn": "#e8a54b",
    "bad": "#d4676a",
}


def style():
    plt.rcParams.update(
        {
            "figure.facecolor": C["bg"],
            "axes.facecolor": C["panel"],
            "axes.edgecolor": C["line"],
            "axes.labelcolor": C["ink"],
            "text.color": C["ink"],
            "xtick.color": C["muted"],
            "ytick.color": C["muted"],
            "grid.color": C["line"],
            "font.family": "DejaVu Sans",
            "axes.titleweight": "bold",
            "axes.titlesize": 13,
            "axes.labelsize": 11,
        }
    )


def save(fig, name: str):
    for dest in (FIG / name, ART / name):
        fig.savefig(dest, dpi=160, bbox_inches="tight", facecolor=C["bg"])
    plt.close(fig)
    print("wrote", name)


def fig_suite_scores():
    suites = [
        "Priority-10\ndirectional",
        "Public breath\ndirectional",
        "Literature\npublic breath",
        "Sci Data\nelev. direction",
        "Comorbidity\ndirectional",
        "Multisite\ncomposite",
        "Lit audit\ncase pass",
        "Vision\nfidelity",
        "Vision\nevidence-backed",
    ]
    vals = [100, 98.2, 100, 88.9, 83.5, 98.8, 100, 98.02, 97.28]
    colors = []
    for v in vals:
        if v >= 95:
            colors.append(C["teal"])
        elif v >= 85:
            colors.append(C["accent"])
        else:
            colors.append(C["rose"])

    fig, ax = plt.subplots(figsize=(11, 5.2))
    y = np.arange(len(suites))
    bars = ax.barh(y, vals, color=colors, height=0.68, edgecolor="none")
    ax.set_yticks(y)
    ax.set_yticklabels(suites, fontsize=10)
    ax.set_xlim(0, 112)
    ax.axvline(90, color=C["muted"], ls="--", lw=1, alpha=0.7)
    ax.text(90.5, len(suites) - 0.3, "90% gate", color=C["muted"], fontsize=9)
    ax.set_xlabel("Score (%)")
    ax.set_title("Real-data evaluation suite scores")
    for bar, v in zip(bars, vals):
        ax.text(v + 1.2, bar.get_y() + bar.get_height() / 2, f"{v:.1f}%", va="center", fontsize=9, color=C["ink"])
    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.35)
    save(fig, "01_suite_scores.png")


def fig_vision_breakdown():
    metrics = ["Direction", "Fold gates", "Top-k recall", "Pathway", "Cell states", "Sites"]
    vals = [99.0, 100.0, 99.3, 91.8, 100.0, 99.3]
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    angles = np.linspace(0, 2 * np.pi, len(metrics), endpoint=False).tolist()
    vals_c = vals + vals[:1]
    angles_c = angles + angles[:1]
    ax = plt.subplot(111, polar=True)
    ax.set_facecolor(C["panel"])
    ax.plot(angles_c, vals_c, color=C["teal"], lw=2.2)
    ax.fill(angles_c, vals_c, color=C["teal"], alpha=0.28)
    ax.set_xticks(angles)
    ax.set_xticklabels(metrics, fontsize=10, color=C["ink"])
    ax.set_ylim(0, 100)
    ax.set_yticks([50, 75, 100])
    ax.set_yticklabels(["50", "75", "100"], color=C["muted"], fontsize=8)
    ax.set_title("Vision suite metric means (50 diseases)", pad=18, color=C["ink"])
    fig = ax.figure
    fig.patch.set_facecolor(C["bg"])
    save(fig, "02_vision_radar.png")


def fig_vision_by_grade():
    grades = ["A\nmeasured", "B\ndirectional lit", "C\nheld-out lit", "D\nprior-approx"]
    vals = [99.9, 96.5, 97.7, 99.2]
    ns = [6, 20, 6, 18]
    fig, ax = plt.subplots(figsize=(8, 4.6))
    x = np.arange(len(grades))
    bars = ax.bar(x, vals, color=[C["teal"], C["blue"], C["accent"], "#6b7c8f"], width=0.62, edgecolor="none")
    ax.set_xticks(x)
    ax.set_xticklabels(grades)
    ax.set_ylim(90, 101.5)
    ax.set_ylabel("Mean composite (%)")
    ax.set_title("Vision fidelity by evidence grade")
    for i, (b, v, n) in enumerate(zip(bars, vals, ns)):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.25, f"{v:.1f}%\nn={n}", ha="center", fontsize=9)
    ax.grid(axis="y", alpha=0.35)
    note = ax.text(
        0.98,
        0.08,
        "Grade D GT aligns with atlas priors → consistency, not held-out proof",
        transform=ax.transAxes,
        ha="right",
        fontsize=8,
        color=C["muted"],
        style="italic",
    )
    save(fig, "03_vision_by_grade.png")


def fig_weakest_diseases():
    names = [
        "schizophrenia",
        "covid19",
        "depression",
        "pneumonia",
        "CKD",
        "epilepsy",
        "breast ca",
        "colon ca",
    ]
    comps = [85.0, 89.6, 91.5, 91.9, 91.9, 93.1, 93.5, 94.6]
    # stack: pathway as main soft axis for several
    pathway = [67, 50, 67, 50, 50, 67, 100, 67]
    fig, ax = plt.subplots(figsize=(9, 4.8))
    x = np.arange(len(names))
    ax.bar(x, comps, color=C["blue"], width=0.55, label="Composite %", edgecolor="none")
    ax.plot(x, pathway, color=C["accent"], marker="o", lw=2, label="Pathway axis %")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=25, ha="right")
    ax.set_ylim(40, 105)
    ax.set_ylabel("%")
    ax.set_title("Weakest vision profiles — composite vs pathway explainability")
    ax.legend(frameon=False, labelcolor=C["ink"])
    ax.grid(axis="y", alpha=0.35)
    save(fig, "04_weakest_diseases.png")


def fig_public_breath_cases():
    cases = [
        "sci asthma",
        "sci COPD",
        "sci bronch.",
        "lit T2D",
        "lit LUAD",
        "lit AD",
        "lit IBD",
    ]
    recall = [57.1, 0.0, 66.7, 100, 100, 100, 100]
    direction = [100, 100, 87.5, 100, 100, 100, 100]
    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    x = np.arange(len(cases))
    w = 0.36
    ax.bar(x - w / 2, recall, w, color=C["accent"], label="Elevated recall@15", edgecolor="none")
    ax.bar(x + w / 2, direction, w, color=C["teal"], label="Directional accuracy", edgecolor="none")
    ax.set_xticks(x)
    ax.set_xticklabels(cases, rotation=20, ha="right")
    ax.set_ylim(0, 115)
    ax.set_ylabel("%")
    ax.set_title("Public breath cases: Sci Data is harder than literature panels")
    ax.legend(frameon=False, labelcolor=C["ink"])
    ax.grid(axis="y", alpha=0.35)
    ax.annotate(
        "COPD elev set has 1 VOC;\ncorrect direction but not in top-15",
        xy=(1, 5),
        xytext=(1.6, 40),
        fontsize=8,
        color=C["muted"],
        arrowprops=dict(arrowstyle="->", color=C["muted"]),
    )
    save(fig, "05_public_breath_cases.png")


def fig_comorbidity():
    cases = ["MDD\nprimary", "MDD +\nobesity", "SZ\nprimary", "SZ +\nheart"]
    elev = [100, 100, 100, 100]
    supp = [100, 100, 50, 25]
    overall = [100, 100, 71.4, 62.5]
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    x = np.arange(len(cases))
    w = 0.25
    ax.bar(x - w, elev, w, color=C["teal"], label="Elevated direction", edgecolor="none")
    ax.bar(x, supp, w, color=C["rose"], label="Suppressed direction", edgecolor="none")
    ax.bar(x + w, overall, w, color=C["blue"], label="Overall direction", edgecolor="none")
    ax.set_xticks(x)
    ax.set_xticklabels(cases)
    ax.set_ylim(0, 120)
    ax.set_title("Comorbidity clinical: elevations strong; SZ suppressions weak")
    ax.legend(frameon=False, labelcolor=C["ink"], ncols=3, loc="upper center")
    ax.grid(axis="y", alpha=0.35)
    ax.text(
        0.98,
        0.05,
        "SZ wants ↓ acetone/methanol/TMA; model often predicts ↑",
        transform=ax.transAxes,
        ha="right",
        fontsize=8,
        color=C["muted"],
        style="italic",
    )
    save(fig, "06_comorbidity.png")


def fig_corpus_origins():
    labels = [
        "Literature\nmixed",
        "Literature\ndirectional",
        "Measured\ncohort relative",
        "Measured\ncohort",
        "Literature\nquantified",
    ]
    sizes = [64, 53, 35, 19, 18]
    colors = [C["blue"], C["teal"], C["accent"], "#c47a3a", "#6b7c8f"]
    fig, ax = plt.subplots(figsize=(7.2, 5))
    wedges, texts, autotexts = ax.pie(
        sizes,
        labels=labels,
        colors=colors,
        autopct=lambda p: f"{p:.0f}%",
        startangle=90,
        textprops={"color": C["ink"], "fontsize": 9},
        wedgeprops={"linewidth": 1, "edgecolor": C["bg"]},
        pctdistance=0.72,
    )
    for t in autotexts:
        t.set_color(C["bg"])
        t.set_fontweight("bold")
    ax.set_title("Training labels (n=189): real only — no synthetic y")
    save(fig, "07_corpus_origins.png")


def fig_calibrator_mae():
    vocs = [
        "ammonia",
        "heptanal",
        "ethanol",
        "pentane",
        "propanol",
        "nonanal",
        "2_butanone",
        "phenol",
        "benzene",
        "isopropanol",
        "octane",
        "ethylbenzene",
        "hexanal",
        "acetone",
        "2_pentanone",
        "isoprene",
    ]
    mae = [
        0.00003,
        0.0003,
        0.0035,
        0.0017,
        0.032,
        0.046,
        0.023,
        0.110,
        0.118,
        0.133,
        0.135,
        0.165,
        0.241,
        0.280,
        0.655,
        1.128,
    ]
    fig, ax = plt.subplots(figsize=(10, 4.8))
    colors = [C["teal"] if m < 0.15 else (C["accent"] if m < 0.4 else C["rose"]) for m in mae]
    ax.barh(np.arange(len(vocs)), mae, color=colors, edgecolor="none", height=0.7)
    ax.set_yticks(np.arange(len(vocs)))
    ax.set_yticklabels(vocs, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("MAE (log2 fold-change)")
    ax.set_title("Per-VOC calibrator fit error (lower is better)")
    ax.grid(axis="x", alpha=0.35)
    ax.text(
        0.98,
        0.05,
        "Most VOCs fit on small-n with no holdout → optimistic R² possible",
        transform=ax.transAxes,
        ha="right",
        fontsize=8,
        color=C["muted"],
        style="italic",
    )
    save(fig, "08_calibrator_mae.png")


def fig_pipeline_flow():
    fig, ax = plt.subplots(figsize=(11, 3.2))
    ax.set_xlim(0, 11)
    ax.set_ylim(0, 3)
    ax.axis("off")
    steps = [
        (0.3, "Real\ncorpus"),
        (2.1, "Train\ncalibrator"),
        (3.9, "Priority-10\n& public"),
        (5.7, "Comorbidity\n& multisite"),
        (7.5, "Audit &\nvision"),
        (9.3, "Completion\ngate"),
    ]
    for i, (x, label) in enumerate(steps):
        box = mpatches.FancyBboxPatch(
            (x, 1.0),
            1.5,
            1.2,
            boxstyle="round,pad=0.05,rounding_size=0.15",
            facecolor=C["panel"],
            edgecolor=C["teal"] if i < 5 else C["accent"],
            linewidth=1.8,
        )
        ax.add_patch(box)
        ax.text(x + 0.75, 1.6, label, ha="center", va="center", fontsize=9, color=C["ink"], fontweight="bold")
        if i < len(steps) - 1:
            ax.annotate(
                "",
                xy=(x + 1.65, 1.6),
                xytext=(x + 1.5, 1.6),
                arrowprops=dict(arrowstyle="->", color=C["muted"], lw=1.5),
            )
    ax.set_title("End-to-end real-data pipeline", color=C["ink"], pad=8)
    ax.text(5.5, 0.35, "Synthetic VOC labels blocked · y = MW / Sci Data / literature only", ha="center", color=C["muted"], fontsize=9)
    save(fig, "09_pipeline_flow.png")


def fig_diabetes_tempo():
    phases = ["Prediabetes", "Controlled\nT2D", "Poorly\ncontrolled", "T2D +\nobesity", "Ketotic"]
    acetone = [3.05, 12.06, 33.54, 33.54, 33.54]
    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    ax.plot(phases, acetone, color=C["accent"], marker="o", ms=9, lw=2.4)
    ax.fill_between(range(len(phases)), acetone, color=C["accent"], alpha=0.18)
    ax.axhline(1.8, color=C["teal"], ls="--", lw=1.2)
    ax.text(0.05, 2.1, "Literature gate ≥1.8×", color=C["teal"], fontsize=9)
    ax.set_ylabel("Acetone fold vs healthy")
    ax.set_title("Diabetes clinical tempo — breath acetone axis")
    ax.grid(axis="y", alpha=0.35)
    ax.text(
        0.98,
        0.08,
        "Physiology cap plateaus ~33× in later phases",
        transform=ax.transAxes,
        ha="right",
        fontsize=8,
        color=C["muted"],
        style="italic",
    )
    save(fig, "10_diabetes_tempo.png")


def main():
    style()
    fig_suite_scores()
    fig_vision_breakdown()
    fig_vision_by_grade()
    fig_weakest_diseases()
    fig_public_breath_cases()
    fig_comorbidity()
    fig_corpus_origins()
    fig_calibrator_mae()
    fig_pipeline_flow()
    fig_diabetes_tempo()
    print("done", FIG)


if __name__ == "__main__":
    main()
