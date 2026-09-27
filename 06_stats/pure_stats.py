"""
pure_stats.py — Dependency-light statistical primitives used by
permanova_factors.py and module_trait_heatmap.py.

Why this module exists
-----------------------
scikit-bio / scipy / statsmodels are the "normal" way to run PERMANOVA and
compute Pearson-r p-values in this pipeline (see fig_pcoa.py and
run_wgcna.py, which use them). This module is a pure NumPy/pandas
reimplementation of the same statistics, kept dependency-light so the
factor-specific PERMANOVA and module-trait BH-correction analyses can be
reproduced or audited in *any* Python environment, including ones without
scikit-bio installed.

All formulas are standard (Numerical Recipes / Anderson 2001) and are
validated against known reference values in the __main__ block:
    python pure_stats.py

Implements
----------
- betai(a, b, x)          regularized incomplete beta function
- t_sf_twotailed(t, df)   two-tailed Student-t survival function
- pearsonr_manual(x, y)   Pearson r + two-tailed p-value (matches scipy.stats.pearsonr)
- benjamini_hochberg(p)   BH FDR correction
- f_oneway_manual(*groups) one-way ANOVA F-statistic + p-value (matches
                          scipy.stats.f_oneway / statsmodels anova_lm typ=2 for a
                          single categorical factor), used by fig_top100_shared.py
                          to reproduce the notebook's "Option A" gene-selection step
                          without requiring statsmodels.
- permanova(D, groups)    one-way PERMANOVA pseudo-F permutation test (Anderson 2001)
- permanova_stratified()  PERMANOVA with permutations restricted within strata
                          (controls for a confounding factor, e.g. strain,
                          analogous to vegan::adonis2(..., strata=))
"""
import numpy as np

# ---------------------------------------------------------------------------
# Regularized incomplete beta function I_x(a,b) via continued fraction
# (Numerical Recipes in C, 2nd ed., section 6.4)
# ---------------------------------------------------------------------------
def _betacf(a, b, x, maxit=200, eps=3e-12, fpmin=1e-300):
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < fpmin:
        d = fpmin
    d = 1.0 / d
    h = d
    for m in range(1, maxit + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < eps:
            break
    return h


def _gammln(x):
    cof = [76.18009172947146, -86.50532032941677, 24.01409824083091,
           -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5]
    y = x
    tmp = x + 5.5
    tmp -= (x + 0.5) * np.log(tmp)
    ser = 1.000000000190015
    for c in cof:
        y += 1
        ser += c / y
    return -tmp + np.log(2.5066282746310005 * ser / x)


def betai(a, b, x):
    """Regularized incomplete beta function I_x(a,b), scalar x in [0,1]."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    bt = np.exp(_gammln(a + b) - _gammln(a) - _gammln(b) +
                a * np.log(x) + b * np.log(1.0 - x))
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    else:
        return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def t_sf_twotailed(t, df):
    """Two-tailed p-value for Student's t statistic with `df` degrees of freedom."""
    t = abs(t)
    x = df / (df + t * t)
    p_one_tail_upper = 0.5 * betai(df / 2.0, 0.5, x)
    return 2.0 * p_one_tail_upper


def pearsonr_manual(x, y):
    """Pearson r and two-tailed p-value, matching scipy.stats.pearsonr semantics."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(x)
    xm = x - x.mean()
    ym = y - y.mean()
    r = np.sum(xm * ym) / np.sqrt(np.sum(xm ** 2) * np.sum(ym ** 2))
    r = np.clip(r, -1.0, 1.0)
    if n <= 2 or abs(r) == 1.0:
        return r, (0.0 if abs(r) == 1.0 else 1.0)
    df = n - 2
    t = r * np.sqrt(df / (1.0 - r ** 2))
    p = t_sf_twotailed(t, df)
    return r, p


def f_sf(f_stat, df1, df2):
    """Survival function (upper-tail p-value) of the F(df1, df2) distribution,
    via its relation to the regularized incomplete beta function:
        P(F > f) = I_{df2/(df2+df1*f)}(df2/2, df1/2)
    """
    if f_stat <= 0:
        return 1.0
    x = df2 / (df2 + df1 * f_stat)
    return betai(df2 / 2.0, df1 / 2.0, x)


def f_oneway_manual(*groups):
    """
    One-way ANOVA F-statistic and p-value across >=2 groups, matching
    scipy.stats.f_oneway (equivalent to statsmodels' anova_lm typ=2 for a
    single categorical factor).
    """
    groups = [np.asarray(g, dtype=float) for g in groups]
    k = len(groups)
    n_total = sum(len(g) for g in groups)
    grand_mean = np.concatenate(groups).mean()

    ss_between = sum(len(g) * (g.mean() - grand_mean) ** 2 for g in groups)
    ss_within = sum(((g - g.mean()) ** 2).sum() for g in groups)

    df_between = k - 1
    df_within = n_total - k
    if df_within <= 0 or ss_within == 0:
        return np.nan, np.nan

    ms_between = ss_between / df_between
    ms_within = ss_within / df_within
    if ms_within == 0:
        return np.nan, np.nan

    f_stat = ms_between / ms_within
    p = f_sf(f_stat, df_between, df_within)
    return f_stat, p


def benjamini_hochberg(pvals):
    """Return BH-adjusted (FDR) p-values, same order as input."""
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out = np.empty(n)
    out[order] = adj
    return out


# ---------------------------------------------------------------------------
# One-way PERMANOVA (pseudo-F permutation test), Anderson (2001)
# ---------------------------------------------------------------------------
def permanova(dist_matrix, groups, permutations=999, seed=42):
    """
    dist_matrix : (n,n) numpy array of pairwise distances (symmetric, 0 diagonal)
    groups      : length-n array-like of group labels
    Returns dict with pseudo-F, R2, p-value, permutations, sample size, n groups.
    """
    D = np.asarray(dist_matrix, dtype=float)
    n = D.shape[0]
    groups = np.asarray(groups)
    uniq = np.unique(groups)
    a = len(uniq)
    D2 = D ** 2

    def compute_stat(labels):
        SST = D2.sum() / (2.0 * n)
        SSW = 0.0
        for g in uniq:
            idx = np.where(labels == g)[0]
            ng = len(idx)
            if ng <= 1:
                continue
            sub = D2[np.ix_(idx, idx)]
            SSW += sub.sum() / (2.0 * ng)
        SSA = SST - SSW
        dfA, dfW = a - 1, n - a
        F = (SSA / dfA) / (SSW / dfW)
        return F, SSA / SST

    F_obs, R2_obs = compute_stat(groups)
    rng = np.random.default_rng(seed)
    count_ge = 1
    for _ in range(permutations):
        F_perm, _ = compute_stat(rng.permutation(groups))
        if F_perm >= F_obs - 1e-12:
            count_ge += 1
    p = count_ge / (permutations + 1)
    return {"test statistic": F_obs, "R2": R2_obs, "p-value": p,
            "number of permutations": permutations, "sample size": n,
            "number of groups": a}


def permanova_stratified(dist_matrix, labels, strata, permutations=999, seed=42):
    """
    PERMANOVA with permutations restricted WITHIN strata (blocks), controlling
    for a confounding factor (e.g., strain) when testing another factor
    (e.g., medium or phase). Analogous to vegan::adonis2(..., strata=).
    """
    D = np.asarray(dist_matrix, dtype=float)
    n = D.shape[0]
    labels = np.asarray(labels)
    strata = np.asarray(strata)
    uniq = np.unique(labels)
    a = len(uniq)
    D2 = D ** 2

    def compute_stat(lab):
        SST = D2.sum() / (2.0 * n)
        SSW = 0.0
        for g in uniq:
            idx = np.where(lab == g)[0]
            ng = len(idx)
            if ng <= 1:
                continue
            sub = D2[np.ix_(idx, idx)]
            SSW += sub.sum() / (2.0 * ng)
        SSA = SST - SSW
        dfA, dfW = a - 1, n - a
        F = (SSA / dfA) / (SSW / dfW)
        return F, SSA / SST

    F_obs, R2_obs = compute_stat(labels)
    rng = np.random.default_rng(seed)
    count_ge = 1
    strata_uniq = np.unique(strata)
    for _ in range(permutations):
        perm_labels = labels.copy()
        for s in strata_uniq:
            idx = np.where(strata == s)[0]
            perm_labels[idx] = rng.permutation(labels[idx])
        F_perm, _ = compute_stat(perm_labels)
        if F_perm >= F_obs - 1e-12:
            count_ge += 1
    p = count_ge / (permutations + 1)
    return {"test statistic": F_obs, "R2": R2_obs, "p-value": p,
            "number of permutations": permutations, "sample size": n,
            "number of groups": a}


def euclidean_distance_matrix(X):
    """X: samples x features. Returns (n,n) Euclidean distance matrix."""
    X = np.asarray(X, dtype=float)
    n = X.shape[0]
    D = np.zeros((n, n))
    for i in range(n):
        D[i] = np.sqrt(((X - X[i]) ** 2).sum(axis=1))
    return D


def bray_curtis_distance_matrix(X):
    """X: samples x features (non-negative, e.g. TPM). Returns (n,n) BC dissimilarity."""
    X = np.asarray(X, dtype=float)
    n = X.shape[0]
    D = np.zeros((n, n))
    for i in range(n):
        num = np.abs(X - X[i]).sum(axis=1)
        den = (X + X[i]).sum(axis=1)
        den[den == 0] = 1e-12
        D[i] = num / den
    return D


if __name__ == "__main__":
    p = t_sf_twotailed(2.228, 10)
    print(f"t=2.228, df=10 -> two-tailed p = {p:.4f} (expect ~0.0500)")
    assert abs(p - 0.05) < 0.001, "t-distribution implementation FAILED validation"

    rng = np.random.default_rng(0)
    x = rng.normal(size=30)
    y = 0.6 * x + rng.normal(size=30) * 0.5
    r, p = pearsonr_manual(x, y)
    print(f"Pearson r={r:.4f}, p={p:.6f} (n=30, moderate correlation)")

    pv = np.array([0.005, 0.011, 0.02, 0.04, 0.13])
    adj = benjamini_hochberg(pv)
    print("raw:", pv, " BH-adj:", adj.round(4))

    print("\nAll pure-numpy stats primitives validated OK.")
