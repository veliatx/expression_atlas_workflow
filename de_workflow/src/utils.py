"""
Utility functions for differential expression analysis workflow.

This module provides various statistical and data manipulation functions for RNA-seq analysis,
including:
- Design matrix releveling for DESeq2
- Likelihood ratio testing
- PCA variance analysis
- Mixed linear modeling for variance components
- TMM normalization
"""

import copy
from typing import List, Tuple

import anndata as ad
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.utils.parallel import parallel_backend

from scipy.stats import chi2
from scipy.stats.distributions import nbinom

import statsmodels.api as sm 
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

from pydeseq2.ds import DeseqStats
from pydeseq2.dds import DeseqDataSet
from pydeseq2.utils import build_design_matrix, nb_nll


def relevel_design(dds: DeseqDataSet, ref_level: List[str]) -> None:
    """Relevels pydeseq2 DeseqDataSet to level in ref_level. Rearranges coefficients to accomodate
    new reference level and rebuilds design matrix.

    Args:
        dds (DeseqDataset) pydeseq2 object to modify
        ref_level (List[str]) list of two elements ['condition','new_reference_level']
    """
    if ref_level[0] not in dds.obs.columns:
        raise ValueError(f"{ref_level[0]} condition not in design.")
    if ref_level[1] not in dds.obs[ref_level[0]].values:
        raise ValueError(f"{ref_level[1]} condition level not in {ref_level[0]}.")
    if not dds.ref_level:
        raise AttributeError(
            '%s define reference level "ref_level" for original design.'
        )

    if any(
        True if c.startswith(ref_level[0]) and c.endswith(ref_level[1]) else False
        for c in dds.obsm["design_matrix"].columns
    ):
        print(f"{ref_level[0]} already reference level for {ref_level[1]}")
        return

    design_matrix = build_design_matrix(
        metadata=dds.obs.copy(),
        design_factors=dds.design_factors,
        ref_level=ref_level,
    )

    if "LFC" not in dds.varm.keys():
        dds.deseq2()

    coef_df = dds.varm["LFC"].copy()

    refo = [c.split("_")[-1] for c in dds.varm["LFC"] if c.startswith(ref_level[0])][0]

    columns_to_relevel = [
        c for c in design_matrix.columns if c.startswith(ref_level[0])
    ]

    coef_df["intercept"] = (
        dds.varm["LFC"]["intercept"]
        + dds.varm["LFC"][f"{ref_level[0]}_{ref_level[1]}_vs_{refo}"]
    )

    for c in columns_to_relevel:
        ref, con = c.split(ref_level[0] + "_")[-1].split("_vs_")

        if f"{ref_level[0]}_{con}_vs_{ref}" in dds.varm["LFC"].columns:
            coef_df[c] = -1.0 * dds.varm["LFC"][f"{ref_level[0]}_{con}_vs_{ref}"]
        else:
            coef_df[c] = (
                dds.varm["LFC"][f"{ref_level[0]}_{ref}_vs_{refo}"]
                - dds.varm["LFC"][f"{ref_level[0]}_{con}_vs_{refo}"]
            )

    columns_drop = [c for c in dds.varm["LFC"] if c.startswith(ref_level[0])]

    coef_df.drop(columns_drop, axis=1, inplace=True)
    coef_df = coef_df[design_matrix.columns]

    dds.varm["LFC"] = coef_df.copy()
    dds.obsm["design_matrix"] = design_matrix.copy()
    dds.ref_level = ref_level

    print(f"dds releveled to {ref_level[0]}-{ref_level[1]}")


def likelihood_ratio_test(
    dds: DeseqDataSet, factors: List[str], alpha: float = 0.05
) -> pd.DataFrame:
    """Perform likelihood ratio test of full model against null model lacking factors described in arguments.

    Args:
        dds (DeseqDataSet) pydeseq2 object with full design
        factors (List[str]) factors to drop from full design
        alpha (float) alpha value for pvalue adjustment via bh-correction

    Returns:
        (pd.DataFrame) dataframe containg lr-statistic, pvalue, and padj of full model vs. null model
    """

    if any(True if c not in dds.design_factors else False for c in factors):
        raise ValueError(
            f"check to make sure all factors {','.join(factors)} in design factors."
        )

    # Calculate likelihood of fit model.
    mu_fit = np.stack(
        dds.varm["LFC"].apply(
            lambda x: dds.obsm["size_factors"]
            * np.exp(dds.obsm["design_matrix"].values @ x),
            axis=1,
        )
    )

    with parallel_backend("loky", inner_max_num_threads=1):
        res = Parallel(
            n_jobs=dds.n_processes,
            verbose=dds.joblib_verbosity,
            batch_size=dds.batch_size,
        )(
            delayed(nb_nll)(
                dds.X[:, i],
                mu_fit[i, :],
                dds.varm["dispersions"][i],
            )
            for i in range(dds.X.shape[1])
        )

    l_fit = 2.0 * np.array(res)

    # Calculate likelihood of null model.
    dds_null = copy.deepcopy(dds)

    dds_null.obsm["design_matrix"].drop(
        [
            c
            for c in dds_null.obsm["design_matrix"]
            if any(c.startswith(f) for f in factors)
        ],
        axis=1,
        inplace=True,
    )

    # Fit null model given estimated dispersion parameters from full model.
    dds_null.fit_LFC()

    mu_null = np.stack(
        dds_null.varm["LFC"].apply(
            lambda x: dds_null.obsm["size_factors"]
            * snp.exp(dds_null.obsm["design_matrix"].values @ x),
            axis=1,
        )
    )

    with parallel_backend("loky", inner_max_num_threads=1):
        res = Parallel(
            n_jobs=dds.n_processes,
            verbose=dds.joblib_verbosity,
            batch_size=dds.batch_size,
        )(
            delayed(nb_nll)(
                dds_null.X[:, i],
                mu_null[i, :],
                dds_null.varm["dispersions"][i],
            )
            for i in range(dds_null.X.shape[1])
        )

    l_null = 2.0 * np.array(res)

    # Perform LRT.
    lr_statistic = l_null - l_fit

    df = dds.obsm["design_matrix"].shape[1] - dds_null.obsm["design_matrix"].shape[1]

    p = chi2.sf(lr_statistic, df)

    padj = multipletests(p, alpha=alpha, method="fdr_bh")[1]

    return pd.DataFrame(
        np.array([lr_statistic, p, padj]).T,
        index=dds.varm["LFC"].index,
        columns=["lrstat", "pvalue", "padj"],
    )


def pca_anova_model(
    adata: ad.AnnData, var_cutoff: float = 0.6, skip_columns: List[str] = []
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Perform ANOVA analysis on PCA components to determine variance explained by different factors.

    Args:
        adata (ad.AnnData): AnnData object containing PCA results in obsm['X_pca']
        var_cutoff (float): Cumulative variance ratio cutoff to determine number of PCs to analyze
        skip_columns (List[str]): List of column names to exclude from the analysis

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: Two dataframes containing:
            1. Variance components for each factor and PC
            2. F-test p-values for each factor and PC
    """
    try:
        pc_n = np.where(np.cumsum(adata.uns["pca"]["variance_ratio"]) > var_cutoff)[0][
            0
        ]
    except:
        pc_n = 6

    columns = adata.obs.columns[
        (
            adata.obs.columns.str.startswith("sample_type")
            | adata.obs.columns.str.startswith("sample_condition")
            | adata.obs.columns.str.startswith("sample-condition")
            | adata.obs.columns.str.startswith("group")
        )
        & ~adata.obs.columns.isin(skip_columns)
    ]

    df = adata.obs.loc[:, columns]
    df = df.loc[:, (df.nunique() > 1) & ~adata.obs.isna().any(axis=0)]
    df.columns = df.columns.str.replace("-", "").str.replace("_", "")

    pc_df = pd.DataFrame(columns=df.columns.tolist() + ["residual", "PC"])
    f_df = pd.DataFrame(columns=df.columns.tolist() + ["PC"])

    for i in range(0, min(pc_n, 6)):

        df["PC"] = adata.obsm["X_pca"][:, i]

        try:
            lm = smf.ols(
                f'PC ~ {"+ ".join([f"C({c})" for c in df.columns if c != "PC"])}',
                data=df,
            ).fit()
            table = sm.stats.anova_lm(lm, type=2)

            re_var = []

            e_var = table.loc["Residual", "sum_sq"]
            re_var = table["sum_sq"][:-1].tolist()
            f = table["PR(>F)"][:-1].tolist()

        except:
            print(f"skipping pc{i+1}")

        pc_df.loc[len(pc_df)] = (np.array(re_var) / (sum(re_var) + e_var)).tolist() + [
            e_var / (sum(re_var) + e_var),
            i,
        ]
        f_df.loc[len(f_df)] = f + [i]

    return pc_df, f_df


def pca_variance_components_model(
    adata: ad.AnnData,
    var_cutoff: float = 0.6,
    skip_columns: List[str] = [],
) -> pd.DataFrame:
    """Fit mixed linear model to estimate variance components of PCA results.

    Args:
        adata (ad.AnnData): AnnData object containing PCA results in obsm['X_pca']
        var_cutoff (float): Cumulative variance ratio cutoff to determine number of PCs to analyze
        skip_columns (List[str]): List of column names to exclude from the analysis

    Returns:
        pd.DataFrame: Dataframe containing variance components for each factor and PC,
                     including residual variance and PC number
    """
    try:
        pc_n = np.where(np.cumsum(adata.uns["pca"]["variance_ratio"]) > var_cutoff)[0][
            0
        ]
    except:
        pc_n = 6

    columns = adata.obs.columns[
        (
            adata.obs.columns.str.startswith("sample_type")
            | adata.obs.columns.str.startswith("sample_condition")
            | adata.obs.columns.str.startswith("sample-condition")
            | adata.obs.columns.str.startswith("group")
        )
        & ~adata.obs.columns.isin(skip_columns)
    ]

    df = adata.obs.loc[:, columns]
    df = df.loc[:, (df.nunique() > 1) & ~adata.obs.isna().any(axis=0)]
    df.columns = df.columns.str.replace("-", "").str.replace("_", "")

    pc_df = pd.DataFrame(columns=df.columns.tolist() + ["residual", "PC"])

    for i in range(0, min(pc_n, 6)):

        df["PC"] = adata.obsm["X_pca"][:, i]

        model = sm.MixedLM.from_formula(
            "PC ~ 1",
            vc_formula={c: f"0 + C({c})" for c in df.columns if not c.startswith("PC")},
            groups=np.ones(df.shape[0]),
            data=df,
        )

        try:
            result = model.fit(method=["lbfgs"])

            re_var = []

            for c in df.columns:
                if c == "PC":
                    continue
                re = result.random_effects[1.0][
                    result.random_effects[1.0].index.str.startswith(c)
                ].values
                re_var.append(
                    (pd.get_dummies(df[c]).astype(int).values * re).sum(axis=1).var()
                )
            e_var = result.scale

        except:
            print(f"skipping pc{i+1}")

        pc_df.loc[len(pc_df)] = (np.array(re_var) / (sum(re_var) + e_var)).tolist() + [
            e_var / (sum(re_var) + e_var),
            i,
        ]

    return pc_df

def tmm_normalize(
    data: np.ndarray,
    trim_lfc=0.3,
    trim_mag=0.05,
) -> Tuple[np.ndarray, np.ndarray]:
    """TMM-normalization, reworked from conorm package.

    Args:
        data (np.ndarray): Expression matrix, rows as samples and columns as genes.
        trim_lfc (float): Cutoff for fold change.
        trim_mag (float): Cutoff for magnitude.
    Returns:
        (Tuple[np.ndarray,np.ndarray]): Normalized expression values, and normalization factors.
    """
    x = data.copy()

    lib_size = np.nansum(x, axis=1)
    mask = x == 0
    x[:, mask.all(axis=0)] = np.nan
    p75 = np.nanpercentile(x, 75, axis=1)
    ref = np.argmin(np.abs(p75 - p75.mean()))
    mask[:, mask[ref]] = True
    x[mask] = np.nan

    norm_x = x / lib_size[:, np.newaxis]
    log = np.log2(norm_x)
    m_g = log - log[ref]
    a_g = (log + log[ref]) / 2

    perc_m_g = np.nanquantile(
        m_g,
        [trim_lfc, 1 - trim_lfc],
        axis=1,
        method="nearest",
    )[..., np.newaxis]

    perc_a_g = np.nanquantile(
        a_g,
        [trim_mag, 1 - trim_mag],
        axis=1,
        method="nearest",
    )[..., np.newaxis]

    mask = mask | (m_g < perc_m_g[0]) | (m_g > perc_m_g[1])
    mask = mask | (a_g < perc_a_g[0]) | (a_g > perc_a_g[1])
    w_gk = (1 - norm_x) / x
    w_gk = 1 / (w_gk + w_gk[ref])

    w_gk[mask] = 0
    m_g[mask] = 0
    w_gk = w_gk / np.nansum(w_gk, axis=1)[:, np.newaxis]
    tmms = np.nansum(w_gk * m_g, axis=1)
    tmms = tmms - tmms.mean()
    tmms = np.exp2(tmms)

    return data / tmms[..., np.newaxis], tmms

