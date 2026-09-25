#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ten-fold cross-validation model panel for the GEM perturbation-module genes.

This aligns the node11 bulk-support heatmap with common biomarker-screening
figures that report multiple algorithms evaluated by 10-fold CV.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.ensemble import ExtraTreesClassifier, GradientBoostingClassifier, RandomForestClassifier
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier


PROJECT_ROOT = Path(r"/path/to/project")
SCREEN = PROJECT_ROOT / "results" / "11_bulk_model"
TABLE = SCREEN / "tables"
TMP = SCREEN / "tmp"
SEED = 20260608


def model_zoo() -> dict[str, object]:
    return {
        "Logistic_L2": LogisticRegression(penalty="l2", solver="liblinear", random_state=SEED, max_iter=2000),
        "Logistic_L1": LogisticRegression(penalty="l1", solver="liblinear", random_state=SEED, max_iter=2000),
        "ElasticNet_logistic": LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=0.5, random_state=SEED, max_iter=5000),
        "Linear_SVM": SVC(kernel="linear", probability=True, random_state=SEED),
        "RBF_SVM": SVC(kernel="rbf", probability=True, random_state=SEED),
        "RandomForest": RandomForestClassifier(n_estimators=180, random_state=SEED, class_weight="balanced"),
        "ExtraTrees": ExtraTreesClassifier(n_estimators=180, random_state=SEED, class_weight="balanced"),
        "GradientBoosting": GradientBoostingClassifier(random_state=SEED),
        "DecisionTree": DecisionTreeClassifier(random_state=SEED, class_weight="balanced"),
        "KNN_3": KNeighborsClassifier(n_neighbors=3),
        "KNN_5": KNeighborsClassifier(n_neighbors=5),
        "GaussianNB": GaussianNB(),
    }


def main() -> None:
    df = pd.read_csv(TMP / "bulk_GEM_gene_matrix.csv")
    df["group"] = df["group"].replace({"ND": "CN", "NC": "CN", "Control": "CN", "AD": "AD", "UC": "AD"})
    y = (df["group"] == "AD").astype(int).to_numpy()
    X = df.drop(columns=["sample_id", "group"]).to_numpy(dtype=float)
    n_genes = X.shape[1]

    selectors = {
        "all": None,
        "top3": SelectKBest(f_classif, k=min(3, n_genes)),
        "top5": SelectKBest(f_classif, k=min(5, n_genes)),
        "top8": SelectKBest(f_classif, k=min(8, n_genes)),
        "top12": SelectKBest(f_classif, k=min(12, n_genes)),
        "top16": SelectKBest(f_classif, k=min(16, n_genes)),
    }
    cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=SEED)

    rows, pred_rows = [], []
    for model_name, estimator in model_zoo().items():
        for selector_name, selector in selectors.items():
            steps = [("scale", StandardScaler())]
            if selector is not None:
                steps.append(("select", selector))
            steps.append(("model", estimator))
            pipe = Pipeline(steps)
            pipeline_name = f"{selector_name}+{model_name}"
            try:
                pred = cross_val_predict(pipe, X, y, cv=cv, method="predict_proba", n_jobs=1)[:, 1]
                auc = roc_auc_score(y, pred)
                p = mannwhitneyu(pred[y == 1], pred[y == 0], alternative="two-sided").pvalue
                rows.append(
                    {
                        "model": model_name,
                        "selector": selector_name,
                        "pipeline": pipeline_name,
                        "cv": "10fold",
                        "auc": auc,
                        "p_value": p,
                    }
                )
                for sid, grp, score in zip(df["sample_id"], df["group"], pred):
                    pred_rows.append(
                        {
                            "sample_id": sid,
                            "group": grp,
                            "pipeline": pipeline_name,
                            "cv": "10fold",
                            "score": score,
                        }
                    )
            except Exception as exc:
                rows.append(
                    {
                        "model": model_name,
                        "selector": selector_name,
                        "pipeline": pipeline_name,
                        "cv": "10fold",
                        "auc": np.nan,
                        "p_value": np.nan,
                        "error": str(exc),
                    }
                )

    res = pd.DataFrame(rows).sort_values("auc", ascending=False)
    preds = pd.DataFrame(pred_rows)
    res.to_csv(TABLE / "GEM_bulk_multimodel_10fold_ROC_summary.csv", index=False)
    preds.to_csv(TABLE / "GEM_bulk_multimodel_10fold_predictions.csv", index=False)
    print(res.head(12).to_string(index=False))
    print(f"n_pipelines={res['pipeline'].nunique()} n_rows={len(res)}")


if __name__ == "__main__":
    main()
