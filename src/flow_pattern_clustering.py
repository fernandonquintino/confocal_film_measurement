import numpy as np
import pandas as pd

from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score
from sklearn.preprocessing import StandardScaler

from src.flow_pattern_criteria import DETECTED_FEATURES, MISSING_FEATURES, load_run

# =========================================================
# Feature space
# =========================================================

CLUSTER_FEATURES = MISSING_FEATURES + DETECTED_FEATURES

# Heavy-tailed features are log-transformed before scaling
LOG_FEATURES = ['bridge_freq_hz', 'bridge_mean_duration_s', 'gap_freq_hz', 'gap_mean_ms']


def feature_matrix(df, features=CLUSTER_FEATURES):
    """
    Standardised feature matrix (one row per run). Missing values (e.g. film
    statistics of runs without detected film) are filled with the column
    median; heavy-tailed features are log1p-transformed.
    """

    X = df[features].astype(float).copy()

    for column in features:
        if column in LOG_FEATURES:
            X[column] = np.log1p(X[column].clip(lower=0))

    X = X.fillna(X.median())

    return StandardScaler().fit_transform(X)


def kmeans_labels(X, k, seed=0):

    return KMeans(n_clusters=k, n_init=50, random_state=seed).fit_predict(X)


def ward_labels(X, k):

    return AgglomerativeClustering(n_clusters=k, linkage='ward').fit_predict(X)

# =========================================================
# Distribution space (t2 histograms)
# =========================================================

HIST_EDGES = np.append(np.arange(0, 8.0, 0.25), 8.1)   # last bin holds the 8.02 cap


def t2_histogram(df, include_missing=True, edges=HIST_EDGES):
    """
    Probability vector of t2: one bin per thickness class and, if
    ``include_missing``, an extra first bin with the undetected fraction.
    Samples flagged as mixer waves are excluded.
    """

    if 'mixer_wave' in df.columns:
        df = df.loc[~df['mixer_wave'].astype(bool)]

    t2 = df['t2'].to_numpy(dtype=float)
    counts, _ = np.histogram(t2[~np.isnan(t2)], bins=edges)

    if include_missing:
        counts = np.concatenate([[np.isnan(t2).sum()], counts])

    total = counts.sum()

    return counts / total if total > 0 else counts.astype(float)


def histogram_table(run_ids, base_path, trim=False, include_missing=True):
    """
    Matrix of t2 histograms (runs x bins).
    """

    return np.array([
        t2_histogram(load_run(run_id, base_path, trim=trim), include_missing)
        for run_id in run_ids
    ])


def js_distance(H):
    """
    Condensed Jensen-Shannon distance matrix (base 2, values in [0, 1]).
    """

    # pdist uses the natural log; dividing by sqrt(ln 2) converts to base 2
    return pdist(H, metric='jensenshannon') / np.sqrt(np.log(2))


def hierarchical_labels(distance, k, method='average'):

    return fcluster(linkage(distance, method=method), t=k, criterion='maxclust') - 1

# =========================================================
# Scores
# =========================================================

def purity(labels, reference):
    """
    Fraction of runs whose cluster majority label equals their own label.
    """

    df = pd.DataFrame({'cluster': labels, 'reference': reference})
    majority = df.groupby('cluster')['reference'].agg(lambda s: s.value_counts().index[0])

    return float((df['cluster'].map(majority) == df['reference']).mean())


def cluster_scores(labels, reference):
    """
    Agreement between a clustering and reference labels.
    """

    return {
        'ARI': adjusted_rand_score(reference, labels),
        'NMI': normalized_mutual_info_score(reference, labels),
        'purity': purity(labels, reference),
        'n_clusters': len(np.unique(labels)),
    }


def silhouette_curve(X=None, distance=None, method='kmeans', ks=range(2, 13)):
    """
    Silhouette score per number of clusters. Use ``X`` for k-means/Ward and
    a condensed ``distance`` for hierarchical clustering.
    """

    rows = []

    for k in ks:

        if method == 'kmeans':
            labels = kmeans_labels(X, k)
            score = silhouette_score(X, labels)

        elif method == 'ward':
            labels = ward_labels(X, k)
            score = silhouette_score(X, labels)

        else:
            labels = hierarchical_labels(distance, k, method=method)
            score = (
                silhouette_score(squareform(distance), labels, metric='precomputed')
                if len(np.unique(labels)) > 1 else np.nan
            )

        rows.append({'k': k, 'silhouette': score})

    return pd.DataFrame(rows)
