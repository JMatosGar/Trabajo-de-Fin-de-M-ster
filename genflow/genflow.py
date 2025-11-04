#Se importan los paquetes, simplificando el acceso a los mismos.
import sgkit as sg
import allel as al
import xarray as xr
import pandas as pd
import numpy as np

#Se importan los paquetes necesarios para las representaciones gráficas.
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.manifold import MDS
from sklearn.decomposition import PCA
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform

#Se construye una primera función que hace uso de los datos cargados directamente desde el dataset original o filtrado.
def ibs(dataset: xr.Dataset, sample_ids: list = None, threshold: float = None, plot: bool = True):
    """
    Calculates pairwise genetic similarity based on Identity by State between individuals.
    optionally allows the user to create a heatmap with the results.

    Parameters
    ----------
    dataset : xr.Dataset
        Dataset containing genotype information. Can be raw or filtered.
    sample_ids : list, default = None
        List of sample IDs to be used. If None is selected, all samples are used.

    threshold : float or None, optional
        Similarity threshold in [0, 1]. 
        Pairs with similarity <= threshold will be masked in the heatmap and excluded from the list. 
        If None, no filtering is applied.
    plot : bool, optional
        If True it will create a heatmap of the IBS matrix and returns it as Matplotlib Figure/Axes. 
        The function will always returns the numeric results.

    Returns
    -------
    ibs : np.ndarray with shape (n_samples, n_samples)
        Matrix of samples x samples containing pairwise genetic similarity between individuals.
    sample_pair : list of tuples
        List containing the tuples label_i, label_j and similarity for each i<j pair.
        If threshold is none, it will include every pair.
    fig : matplotlib.figure or None
        Heatmap created if the user established plot value as True.
    ax : matplotlib.axes or None
        Axes of the heatmat created if the user established plot value as True.

    Raises
    ------
    TypeError
        If `sample_ids` is not a list/tuple/ndarray or contains non-integer entries.
    IndexError
        If any sample index in `sample_ids` is out of range for the dataset.
    ValueError
        If the input genotype array cannot be found at dataset["calldata/GT"] or has incompatible dimensions.

    Examples
    --------
    >>> vcf = load_data("input.vcf", "output.zarr")
    >>> ds = qc_data(vcf, maf_thr = 0.01, var_thr = 0.95, sample_thr = 0.95, ploidy = 4, biallelic = False, min_var = 30)  
    >>> ibs_mat, pairs, fig, ax = ibs(ds, sample_ids=None, threshold=0.8, plot=True)
    >>> print(ibs_mat.shape)
    >>> print(pairs[:5])

    """

    #Se toman los valores del dataset.
    genome = dataset["calldata/GT"].values
    n_var = genome.shape[0]
    n_samples = genome.shape[1]

    #Se seleccionan las muestras indicadas por el usuario en caso de que este haya indicado alguna.
    if sample_ids is not None:

        #Se establecen unas condiciones para comprobar esta keyword.
        if not isinstance(sample_ids, (list, tuple, np.ndarray)):
            raise TypeError("sample_ids must be a list, tuple or NumPy array.")

        if not all(isinstance(i, (int, np.integer)) for i in sample_ids):
            raise TypeError("All elements in sample_ids must be integers.")

        if not all(0 <= i < n_samples for i in sample_ids): #Se verifica que todas las muestras están en el dataset.
            raise IndexError("There is at least one sample ID out of range")


        genome = genome[:, sample_ids, :] #Si todas las muestras son válidas se filtra el subset.
        n_sel = len(sample_ids)
        #Se lanza el error indicando el motivo por el que falla algún indice.


    else:
        n_sel = n_samples

    #Se toma el nombre de las muestras desde la columna sample seleccionando aquellas indicadas por el usuario.
    if "sample" in dataset:

        try:
            sample_names = [str(x) for x in dataset["sample"].values]
        except Exception:
            sample_names = None
    else:
        sample_names = None

    #En caso de que no puedan tomarse los nombres de las muestras se les asigna un nombre en base a su posición.
    if sample_names is None:
        labels = [str(i) for i in (sample_ids if sample_ids is not None else range(n_sel))]
    else:
        if sample_ids is None:
            labels = [str(i) for i in range(n_sel)]
        else:
            labels = [sample_names[i] for i in sample_ids]

    #Se genera el GenotypeArray de las muestras.
    gt = al.GenotypeArray(genome)
    nalt = gt.to_n_alt()

    #Se calcula la matriz IBS utilizando estos datos.
    mat = nalt.T #Se transpone la matriz
    ibs_mat = np.empty((n_sel, n_sel), dtype=float) #Se crea una matriz vacia de dimensiones igual al número de muestras seleccionadas.

    #Se itera sobre las muestras.
    for i in range(n_sel):
        ai = mat[i] #Se almacena el vector de la muestra i para evitar repeticiones internas.

        for j in range(i, n_sel):
            aj = mat[j] #Se almacena el vector de la muestra j para evitar repeticiones internas.

            valid = (ai != -1) & (aj != -1) #Se crea un vector booleano que comprueba que no haya valores de relleno.
            shared = np.count_nonzero((ai == aj) & valid) #Se cuentan las variantes donde los conteos alternativos son idénticos.
            total = np.count_nonzero(valid) #Se cuentan las variantes válidas.
            similarity = shared / float(total) if total > 0 else np.nan #Se obtiene la proporción de variantes compartidas.

            #Se actualiza la matriz para que mantenga simetría.
            ibs_mat[i, j] = similarity
            ibs_mat[j, i] = similarity

    #Se crea una lista que acumula tuples de pares de muestras con su respectiva similitud.
    sample_pair = []
    for i in range(n_sel):
        for j in range(i+1, n_sel):
            pair_sim = (ibs_mat[i, j])
            if threshold is None or pair_sim >= threshold: #se almacena la similitud tanto si establece un umbral como si no.
                sample_pair.append((labels[i], labels[j], pair_sim))

    #Se realiza un filtrado de la matriz en caso de que se haya establecido un threshold.
    if threshold is None:
        ibs_mat_filt = ibs_mat.copy()
    else:
        ibs_mat_filt = ibs_mat.copy()
        ibs_mat_filt[ibs_mat_filt <= threshold] = np.nan #Se establecen como NaN los puntos qeu no superan el umbral.

    #Se establecen las condiciones para generar el mapa de calor.
    fig = None
    ax = None

    if plot: #Solo se genera el mapa de calor si se especifica.
        heat_kwargs = dict(cmap="viridis", vmin=0.0, vmax=1.0, square=True)
        fig, ax = plt.subplots(figsize=(max(5, n_sel / 2.0), max(5, n_sel / 2.0)))
        sns.heatmap(ibs_mat_filt, xticklabels=labels, yticklabels=labels, ax = ax, **heat_kwargs)
        ax.set_title("IBS Matrix")
        ax.set_xlabel("Samples")
        ax.set_ylabel("Samples")
        fig.tight_layout()

    return ibs_mat, sample_pair, fig, ax


#Para evaluar mejor las características del flujo genético se pueden realizar unas visualizaciones sobre la matriz obtenida.
def ibs_analysis(ibs_mat: np.ndarray,
                 dist: str = "linear", cluster: str = "average",
                 mds_centers: int = 2, pcs_centers: int = 3,
                 n_clusters: int = None,
                 seed: int = 10137):
    """
    Perform exploratory analyses on an IBS similarity matrix: convert to distance,
    compute MDS and PCA (double-centred similarity), and produce hierarchical clustering.

    Parameters
    ----------
    ibs_mat : numpy.ndarray
        Symmetric pairwise similarity matrix (n_samples x n_samples). Values should
        be floats in [0, 1] or NaN for missing comparisons. The function will validate
        that the array is 2-dimensional and square.
    dist : {'linear', 'sqrt', 'angular'}, optional
        Method to convert similarity to distance:
        - 'linear'    : d = 1 - s
        - 'sqrt'      : d = sqrt(max(0, 1 - s))
        - 'angular'   : d = arccos(clipped_s) / pi
        Default is 'linear'.
    cluster : str, optional
        Linkage method used by scipy.cluster.hierarchy.linkage for hierarchical
        clustering. Default is 'average'.
    mds_centers : int, optional
        Number of dimensions for MDS embedding. Default is 2.
    pcs_centers : int, optional
        Number of principal components to compute from double-centred IBS matrix.
        Default is 3.
    n_clusters : int or None, optional
        If provided, cut the hierarchical tree to produce this number of flat clusters
        using fcluster (criterion='maxclust'). If None, no flat clustering is returned.
    seed : int, optional
        Random seed forwarded to sklearn.MDS for reproducible embeddings. Default 10137.

    Returns
    -------
    dict
    Dictionary with the following keys:
        - "ibs_mat" : numpy.ndarray
            Input IBS similarity matrix as generated by ibs function.
        - "dist_mat" : numpy.ndarray
            Symmetric distance matrix.
        - "mds_coords" : numpy.ndarray
            Coordinates from MDS with shape (n_samples, mds_centers).
        - "pca_coords" : numpy.ndarray
            Coordinates from PCA with shape (n_samples, min(pcs_centers, n_samples)).
        - "variance_ratio" : list of floats
            PCA explained variance ratio for the returned components.
        - "Z" : ndarray
            Linkage matrix.
        - "clusters" : ndarray or None
            Flat cluster labels if `n_clusters` was specified, otherwise None.
        - "fig_mds" : (matplotlib.figure, matplotlib.axes)
            Figure and axes for the MDS scatter plot.
        - "fig_pca" : (matplotlib.figure, matplotlib.axes)
            Figure and axes for the PCA scatter plot.
        - "fig_den" : (matplotlib.figure, matplotlib.axes)
            Figure and axes for the hierarchical clustering dendrogram.

    Raises
    ------
    ValueError
        If `ibs_mat` is not a squared 2D array
        if `dist` is not one of the accepted methods
        if there are no valid pairs for clustering after NaN handling
        if numerical assumptions are incorrect.
    TypeError
        If argument types are incorrect.

    Examples
    --------
    >>> vcf = load_data("input.vcf", "output.zarr")
    >>> ds = qc_data(vcf, maf_thr = 0.01, var_thr = 0.95, sample_thr = 0.95, ploidy = 4, biallelic = False, min_var = 30)  
    >>> ibs_mat, pairs, fig, ax = ibs(ds, sample_ids=None, threshold=0.8, plot=True)
    >>> ibs = def ibs_analysis(ibs_mat: np.ndarray, dist: "linear", cluster: "average",
                 mds_centers = 2, pcs_centers = 3, n_clusters = 4, seed = 435)
"""

                   
    #Se comprueba que la matriz tiene la forma y estructura correcta.
    ibs_mat = np.asarray(ibs_mat, dtype=float)
    if ibs_mat.ndim != 2 or ibs_mat.shape[0] != ibs_mat.shape[1]:
        raise ValueError("IBS matrix must be a 2D array of (n samples x n samples).")

    #Se toman las etiquetas desde la matriz.
    shape = ibs_mat.shape[0]
    labels = [str(i) for i in range(shape)]

    #Se seleccionan los valores NaN.
    nan_mask = np.isnan(ibs_mat)

    #Se convierte la matriz de similitud en una matriz de distancias.
    dist_method = ("linear", "sqrt", "angular") #Se crea una lista con las formas de calcular las distancias.

    if dist not in dist_method: #Si el método seleccionado es valido se devuelve un mensaje.
        raise ValueError(f"Method to calculate distances is not permitted. Available methods are: {dist_method}.")

    #Se crea una copia de la matriz original para evitar modificarla.
    ibs = ibs_mat.copy()

    #En caso de seleccionar un método válido se comprueba cual se ha seleccionado y se calcula la distancia de la manera correspondiente.
    else:
        if dist == "linear":
            dist_mat = 1.0 - ibs

        elif dist == "sqrt":
            dist_mat = np.sqrt(np.maximum(0.0, 1.0 - ibs))

        elif dist == "angular":
            clip = np.clip(ibs, -1.0, 1.0)
            dist_mat = np.arccos(clip) / np.pi

        #Se filtra la matriz de distancias y se establecen las características necesarias para poder realizar el análisis de la mism.
        dist_mat[nan_mask] = np.nan
        np.fill_diagonal(dist_mat, 0.0)
        dist_mat = (dist_mat + dist_mat.T) / 2.0
        dist_mat = np.maximum(dist_mat, 0.0)

    #Primero se realiza un escalado multidimensional (MDS).
    mds_dist = dist_mat.copy()

    #Se sustituyen los NaN por el valor máximo.
    if np.isnan(mds_dist).any():
        max = np.nanmax(mds_dist)
        mds_dist[np.isnan(mds_dist)] = max

    mds = MDS(n_components=mds_centers, dissimilarity="precomputed", random_state=seed)
    mds_coord = mds.fit_transform(mds_dist)

    #Se representa graficamente el MDS.
    fig_mds, ax_,mds = plt.subplots(figsize=(max(5, 10), max(4, 8)))
    ax_mds.scatter(mds_coord[:, 0], mds_coord[: 1], s = 50, c = "C0")

    for i, label in enumerate(labels):
        ax_mds.text(mds_coord[i, 0], mds_coord[i, 1], label, fontsize=10, ha="center", va="center")

    ax_mds.set_title("MDS")
    ax_mds.set_xlabel("MDS1")
    ax_mds.set_ylabel("MDS2")
    fig_mds.tight_layout()

    #A continuación se calcula el PCA con centrado doble.
    pca_dist = dist_mat.copy()

    if np.isnan(pca_dist).any():
        row_mean = np.nanmean(pca_dist, axis=1, keepdims=True)
        ind = np.where(np.isnan(pca_dist))
        for i, j in zip(*ind):
            pca_dist[i, j] = row_mean[i] if not np.isnan(row_mean[i]) else 0.0

    J = np.eye(pca_dist.shape[0]) - np.ones((pca_dist.shape[0], pca_dist.shape[0])) / pca_dist.shape[0]
    pca_dist = J @ pca_dist @ J

    pca = PCA(n_components=min(pcs_centers, n), random_state=seed)
    pca_coord = pca.fit_transform(pca_dist)
    variance_ratio = pca.explained_variance_ratio_.tolist()

    #Se diseña la representación gráfica.
    fig_pca, ax_pca = plt.subplots(figsize=(max(5, 10), max(4, 8)))
    ax_pca.scatter(pca_coord[:, 0], pca_coord[:, 1], s=50, c="C1")
    for i, label in enumerate(labels):
        ax_pca.text(pca_coord[i, 0], pca_coord[i, 1], label, fontsize=10, ha="center", va="center")

    ax_pca.set_title("PCA")
    ax_pca.set_xlabel("PC1")
    ax_pca.set_ylabel("PC2")
    fig_pca.tight_layout()

    #Finalmente se realiza un clustering jerárquico.
    clust_dist = dist_mat.copy()

    if np.isnan(clust_dist).any():
        fill = np.nanmax(clust_dist)
        clust_dist[np.isnan(clust_dist)] = fill

    dvec = squareform(clust_dist, checks=False)
    Z = sch.linkage(dvec, method=cluster)
    clusters = None
    if n_clusters is not None:
        clusters = sch.fcluster(Z, n_clusters, criterion="maxclust")

    #Se utiliza el clustering jerarquico para construir un dendograma.
    fig_den, ax_den = plt.subplots(figsize=(max(5, 10), max(4, 8)))
    sch.dendrogram(Z, labels=labels, leaf_rotation=90, ax=ax_den)
    ax_den.set_title("Dendrogram of hierarchical clustering")
    fig_den.tight_layout()

    results = {
        "ibs_mat": ibs_mat,
        "dist_,mat": dist_mat,
        "mds_coord": mds_coord,
        "pca_coord": pca_coord,
        "variance_ratio": variance_ratio,
        "clusters": clusters,
        "Z": Z,
        "fig_mds": (fig_mds, ax_mds),
        "fig_pca": (fig_pca, ax_pca),
        "fig_den": (fig_den, ax_den)}

    return
