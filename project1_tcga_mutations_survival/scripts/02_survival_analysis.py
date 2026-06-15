"""
Project 1: Mutation x Survival analysis for TCGA LUAD / BRCA / PAAD / COAD
============================================================================

For each cohort:
  - Loads clinical data and mutated-case lists produced by 01_download_data.py
  - Builds an overall-survival time/event table
      time  = days_to_death (if Dead) else days_to_last_follow_up (if Alive)
      event = 1 if Dead else 0
  - For a curated set of clinically relevant driver genes (chosen from the
    cohort's top-mutated-gene list), runs:
      * Kaplan-Meier curves (mutated vs wild-type) + log-rank test
      * A multivariable Cox proportional-hazards model adjusting for
        age at diagnosis, sex, and AJCC stage (when available)
  - Also runs the same comparison for TTN, the largest human gene, which is
    frequently in "top mutated gene" lists purely because of its size
    (a classic TCGA confounder) -- used as a negative-control comparison.

Outputs (per project), written to ../results/<project>/:
  - km_<GENE>.png            Kaplan-Meier plot
  - survival_summary.csv     log-rank p-values + Cox hazard ratios per gene
  - analysis_table.csv        merged per-case table used for the analysis
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
RESULTS_DIR = os.path.join(os.path.dirname(__file__), "..", "results")

# Curated driver genes per cohort (selected from the top-mutated-gene lists,
# excluding very large genes whose frequent mutation is mostly a passenger
# artifact of gene size). TTN is analyzed separately in every cohort as a
# size-driven negative control.
DRIVER_GENES = {
    "TCGA-LUAD": ["TP53"],
    "TCGA-BRCA": ["PIK3CA", "TP53", "CDH1", "GATA3"],
    "TCGA-PAAD": ["KRAS", "TP53", "SMAD4", "CDKN2A"],
    "TCGA-COAD": ["APC", "TP53", "KRAS"],
}
CONTROL_GENE = "TTN"


def build_survival_table(clin):
    df = clin.copy()
    df["event"] = (df["vital_status"] == "Dead").astype(int)
    df["time"] = df["days_to_death"].fillna(df["days_to_last_follow_up"])
    df = df.dropna(subset=["time"])
    df = df[df["time"] >= 0]
    df["age_years"] = df["age_at_diagnosis_days"] / 365.25

    # Simplify AJCC stage to an ordinal numeric covariate (I-IV)
    def stage_to_num(s):
        if not isinstance(s, str):
            return None
        s = s.upper().replace("STAGE", "").strip()
        for roman, num in [("IV", 4), ("III", 3), ("II", 2), ("I", 1)]:
            if s.startswith(roman):
                return num
        return None

    df["stage_num"] = df["ajcc_pathologic_stage"].apply(stage_to_num)
    df["sex_male"] = (df["gender"] == "male").astype(int)
    return df


def add_mutation_flag(df, proj_dir, gene):
    f = os.path.join(proj_dir, f"mutated_cases_{gene}.csv")
    mut_cases = set(pd.read_csv(f)["case_id"])
    df[f"mut_{gene}"] = df["case_id"].isin(mut_cases).astype(int)
    return df


def km_plot(df, gene, out_path, project_id):
    kmf = KaplanMeierFitter()
    fig, ax = plt.subplots(figsize=(6, 5))
    for label, mask, color in [
        (f"{gene} mutated", df[f"mut_{gene}"] == 1, "crimson"),
        (f"{gene} wild-type", df[f"mut_{gene}"] == 0, "steelblue"),
    ]:
        sub = df[mask]
        kmf.fit(sub["time"], sub["event"], label=f"{label} (n={len(sub)})")
        kmf.plot_survival_function(ax=ax, color=color, ci_show=True)

    lr = logrank_test(
        df.loc[df[f"mut_{gene}"] == 1, "time"],
        df.loc[df[f"mut_{gene}"] == 0, "time"],
        event_observed_A=df.loc[df[f"mut_{gene}"] == 1, "event"],
        event_observed_B=df.loc[df[f"mut_{gene}"] == 0, "event"],
    )
    ax.set_title(f"{project_id}: Overall survival by {gene} mutation status\nlog-rank p = {lr.p_value:.4g}")
    ax.set_xlabel("Days since diagnosis")
    ax.set_ylabel("Survival probability")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return lr.p_value


def cox_for_gene(df, gene):
    cols = ["time", "event", f"mut_{gene}", "age_years", "sex_male", "stage_num"]
    sub = df[cols].dropna()
    if sub[f"mut_{gene}"].nunique() < 2 or len(sub) < 20:
        return None
    cph = CoxPHFitter()
    try:
        cph.fit(sub, duration_col="time", event_col="event")
    except Exception as e:
        print(f"    Cox model failed for {gene}: {e}")
        return None
    summary = cph.summary
    row = summary.loc[f"mut_{gene}"]
    return {
        "gene": gene,
        "n": len(sub),
        "n_mutated": int(sub[f"mut_{gene}"].sum()),
        "HR": row["exp(coef)"],
        "HR_lower95": row["exp(coef) lower 95%"],
        "HR_upper95": row["exp(coef) upper 95%"],
        "cox_p": row["p"],
    }


def main():
    for project_id, genes in DRIVER_GENES.items():
        print(f"\n=== {project_id} ===")
        proj_data_dir = os.path.join(DATA_DIR, project_id)
        proj_results_dir = os.path.join(RESULTS_DIR, project_id)
        os.makedirs(proj_results_dir, exist_ok=True)

        clin = pd.read_csv(os.path.join(proj_data_dir, "clinical.csv"))
        df = build_survival_table(clin)
        print(f"  {len(df)} cases with valid survival time ({df['event'].sum()} deaths)")

        analyze_genes = genes + [CONTROL_GENE]
        summary_rows = []
        for gene in analyze_genes:
            mfile = os.path.join(proj_data_dir, f"mutated_cases_{gene}.csv")
            if not os.path.exists(mfile):
                print(f"  skipping {gene}: no mutation file")
                continue
            d = add_mutation_flag(df.copy(), proj_data_dir, gene)
            n_mut = d[f"mut_{gene}"].sum()
            n_wt = (d[f"mut_{gene}"] == 0).sum()
            if n_mut < 5 or n_wt < 5:
                print(f"  skipping {gene}: too few cases in one group (mut={n_mut}, wt={n_wt})")
                continue

            km_path = os.path.join(proj_results_dir, f"km_{gene}.png")
            p_logrank = km_plot(d, gene, km_path, project_id)
            cox_res = cox_for_gene(d, gene)

            row = {
                "gene": gene,
                "role": "driver_candidate" if gene != CONTROL_GENE else "size_control",
                "n_total": len(d),
                "n_mutated": int(n_mut),
                "n_wildtype": int(n_wt),
                "logrank_p": p_logrank,
            }
            if cox_res:
                row.update(
                    {
                        "cox_HR": cox_res["HR"],
                        "cox_HR_lower95": cox_res["HR_lower95"],
                        "cox_HR_upper95": cox_res["HR_upper95"],
                        "cox_p": cox_res["cox_p"],
                    }
                )
            summary_rows.append(row)
            print(
                f"  {gene}: n_mut={n_mut}, n_wt={n_wt}, logrank p={p_logrank:.4g}"
                + (f", Cox HR={cox_res['HR']:.2f} (p={cox_res['cox_p']:.4g})" if cox_res else ", Cox: n/a")
            )

            # keep the merged table from the last (largest) gene set for reference
            d.to_csv(os.path.join(proj_results_dir, "analysis_table.csv"), index=False)

        pd.DataFrame(summary_rows).to_csv(
            os.path.join(proj_results_dir, "survival_summary.csv"), index=False
        )

    print("\nDone.")


if __name__ == "__main__":
    main()
