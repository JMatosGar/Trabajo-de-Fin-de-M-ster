#Se importan los paquetes, simplificando el acceso a los mismos.
import sgkit as sg
import allel as al
import xarray as xr
import pandas as pd
import numpy as np

#Este análisis requiere de múltiples pasos por lo que es preferible separar distintas funciones para su ejecución:

#Primero se construye una función que permite estimar el valor de Ne (población efectiva)
def temporal_Ne(dataset_t0: xr.Dataset, dataset_t1: xr.Dataset,
                t0_sample: list = None, t1_sample: list = None,
                generations: int = 5, ploidy: int = 2,
                maf: float = 0.01, loci_min: int = 10, ci: tuple[float, float] = (2.5, 97.5),
                boot_rep: int = 5000, seed: int = 10137) -> tuple(pd.DataFrame, dict, pd.DataFrame):
    """
    Estimates the efective population size (Ne) based on temporal allel frequency changes.

    Parameters
    ----------
    dataset_t0 : xarray.Dataset
        sgkit Dataset containing genomic of the poblation at time 0. Can be filtered or not.
        Must contain at least "calldata/GT" and "sample" coordinates.
    dataset_t1 : xarray.Dataset
        sgkit Dataset containing genomic of the poblation at time 1. Can be filtered or not.
        Must contain at least "calldata/GT" and "sample" coordinates.
    t0_sample : optional, list of str
        Sample identifiers present in the dataset corresponding to time 0.
    t1_sample : optional, list of str
        Sample identifiers present in the dataset corresponding to time 1.
    generations : int, default = 5
        Number of generations separating times 0 and 1.
        Default value is 5 but can be specified.
    ploidy : int, default = 2
        Number of chromosomes inherited by the descendants.
        Default value is 2 (diploid) but can be specified.
    maf : float, default = 0.01
        Minimum allel frequency threshold.
    loci_min : int, default = 10
        Minimum number of loci requiered to calculate Ne.
    boot_rep : int, default = 5000
        Number of bootstrap replicates over loci. Enhance confidence intervals.
        Default value is 5000 but can be specified.
    seed : int, default = 10137
        Seed for reproducibility.

    Returns
    -------
    Ne_df : pandas.DataFrame
        Single-row DataFrame that summarizes the estimates of Ne.
        Columns:
        - Ne_estimate : median per-locus Ne.
        - median_boot : bootstrap median of the estimator.
        - low_CI : 5 percentile from bootstrap.
        - up_CI : 95 percentile from bootstrap.
        - n_loci : number of loci contributing to estimation.

    Ne_LD : dictionary
        Dictionary containing arrays of quantities per locus. Designed to be used by ne_ld function.
        Keys:
        - "Ne_locus":
        - "pbar":
        - "delta":

    variants : pandas.DataFrame
        Single-row DataFrame that summarizes the samples and variants counts after the filtering steps:
        - t0:
        - t1:
        - Common:
        - n_samples_t0:
        - n_samples_t1:

    Raises
    ------
      ValueError:
        if ploidy is not an integer with a value of 1 or above.
        if generations is not an integer with a value of 1 or above
        if maf is not a positive value.
        if ci is not a tuple with two values between 0 and 100 and former value is above the latter.
        if "calldata/GT" is not in dataset_t0 or dataset_t1.
        if dataset_t0 or dataset_t1 does not contain samples or variants.
        if dataset_t0 or dataset_t1 does not contain genomic data.
        if common_id is empty.
        if there are not enough loci to estimate Ne.
        if there are not enough loci with valid allel frequencies.
        if there are not enough valid per-locus Ne estimates.

    Example:
    --------
      >>> vcf_t0 = load_data("input_t0.vcf", "output_t0.zarr")
      >>> vcf_t1 = load_data("input_t1.vcf", "output_t1.zarr")
      >>> ne_df, ne_ld, variants = temporal_Ne(vcf_t0, vcf_t1,
                t0_sample = None, t1_sample = None,
                generations = 10, ploidy = 4,
                maf = 0.05, loci_min = 25, ci: = (0.5, 99.5),
                boot_rep = 1000, seed = 435)

    """
#De manera similar a otros módulos se comprueban los datos introducidos.
    if not (isinstance(ploidy, int) and ploidy >= 2):
        raise ValueError(f"Ploidy ({ploidy}) must be an integer with a value of 1 or above")

    if not (isinstance(generations, int) and generations >= 1):
        raise ValueError(f"Generations ({generations}) must be an integer with a value of 1 or above")

    if not (isinstance(maf, (int, float)) and maf > 0):
        raise ValueError(f"MAF ({maf}) must be a positive value.")

    if not (isinstance(ci, (list, tuple)) and len(ci) == 2 and 0 <= ci[0] < ci[1] <= 100):
        raise ValueError("ci must be a tuple with two values (lower | upper) between 0 and 100")

    if not (isinstance(boot_rep, int) and boot_rep >0):
        raise ValueError(f"Boot_rep ({boot_rep}) must be an integer above 0")

    if "calldata/GT" not in dataset_t0:
        raise ValueError("Time 0 dataset does not contain genomic data")
    if "calldata/GT" not in dataset_t1:
        raise ValueError("Time 1 dataset does not contain genomic data")

    n_samples_t0 = len(dataset_t0.dims["samples"])
    n_var_t0 = len(dataset_t0.dims["variants"])
    if n_samples_t0 == 0 or n_var_t0 == 0:
        raise ValueError("Time 0 dataset does not contain samples or variants")

    n_samples_t1 = len(dataset_t1.dims["samples"])
    n_var_t1 = len(dataset_t1.dims["variants"])
    if n_samples_t1 == 0 or n_var_t1 == 0:
        raise ValueError("Time 1 dataset does not contain samples or variants")

#Se utiliza una semilla para crear un generador aleatorio para el bootstrap.
    rng = np.random.default_rng(seed)

#Se identifican las variantes dentro del dataset a tiempo 0.
    if "variants/ID" in dataset_t0:
        id_t0 = dataset_t0["variants/ID"].values.astype(str)
    elif all(i in dataset_t0 for i in ("variants/POS","variants/REF","variants/ALT")):
        pos = dataset_t0["variants/POS"].values.astype(str)
        ref = dataset_t0["variants/REF"].values.astype(str)
        alt = dataset_t0["variants/ALT"].values.astype(str)
        id_t0 = np.char.add(np.char.add(pos, "_"), np.char.add(ref, ("_" + alt)))
    else:
        raise ValueError("Time 1 dataset must contain variants/ID or POS+REF+ALT")

#Se identifican las variantes dentro del dataset a tiempo 1.
    if "variants/ID" in dataset_t1:
        id_t1 = dataset_t1["variants/ID"].values.astype(str)
    elif all(i in dataset_t1 for i in ("variants/POS","variants/REF","variants/ALT")):
        pos = dataset_t1["variants/POS"].values.astype(str)
        ref = dataset_t1["variants/REF"].values.astype(str)
        alt = dataset_t1["variants/ALT"].values.astype(str)
        id_t1 = np.char.add(np.char.add(pos, "_"), np.char.add(ref, ("_" + alt)))
    else:
        raise ValueError("Time 1 dataset must contain variants/ID or POS+REF+ALT")


#Se identifican las variantes comunes entre los dos conjuntos.
    set_t1 = set(id_t1) #Se crea un conjunto para llevar a cabo la busqueda.
    common_mask = np.array([i in set_t1 for i in id_t0], dtype=bool) #Se crea un array booleano que devuelve true cuando hay coincidencia.
    common_id = id_t0[common_mask] #Se extraen los ID del set a tiempo 0 que se encuentran en el set a tiempo 1.

    #Se establece una condición para comprobar que se ha seleccionado al menos una variante.
    if common_id.size == 0:
        raise ValueError("There are no common variants between time 0 dataset and time 1 dataset.")

    vars_t0 = np.where(common_mask)[0] #Se obtienen los indices de las variantes comunes que se encuentran en el set a tiempo 0.

    idx_t1 = {j:i for i, j in enumerate(id_t1)} #Se crea un diccionario que registra el indice en el set a tiempo 1.
    vars_t1 = np.array([idx_t1[j] for j in common_id]) #Se crea un array que toma el indice en el set a tiempo 1 en el mismo orden que el set a tiempo 0.

    #Se crean dos subsets que contienen esta selección de variantes ordenada.
    ds_t0 = dataset_t0.isel(variants=vars_t0)
    ds_t1 = dataset_t1.isel(variants=vars_t1)

    #Adicionalmente se extraen las muestras que contienen las variantes seleccionadas si estas no se han instroducido con los datos.
    if t0_sample is None:
        t0_id = np.arange(ds_t0.dims["samples"], dtype = int)
    else:
        t0_id = np.array(t0_sample, dtype = int)

    if t1_sample is None:
        t1_id = np.arange(ds_t1.dims["samples"], dtype = int)
    else:
        t1_id = np.array(t1_sample, dtype = int)

    #A modo de control se puede crear un diccionario que registra las variantes existente en cada set.
    variants = pd.DataFrame({"t0": int(dataset_t0.dims.get("variants", 0)),
                "t1": int(dataset_t1.dims.get("variants", 0)),
                "Common": int(ds_t0.dims.get("variants", 0)),
                "n_samples_t0": int(len(t0_id)),
                "n_samples_t1": int(len(t1_id))})

#Una vez seleccionadas y registradas las variantes a comparar, se puede proceder a analizar la Deriva genética.
    #Primero se extraen los datos genómicos en ambos tiempos.
    g_t0 = al.GenotypeArray(ds_t0["calldata/GT"].values)
    allel_count_t0 = al.count_alleles(g_t0.subset(samples=t0_id))

    g_t1 = al.GenotypeArray(ds_t1["calldata/GT"].values)
    allel_count_t1 = al.count_alleles(g_t1.subset(samples=t1_id))

    #Se calculan las frecuencias alélicas del alelo alternativo en cada variante en los dos tiempos.
    with np.errstate(divide="ignore", invalid="ignore"): #Se evita que se muestren warnings.
        #Se crean matrices de conteo de alelos en ambos tiempos
        allel_mat_t0 = (allel_count_t0[:, 1] / allel_count_t0.sum(axis=1)).astype(float)
        allel_mat_t1 = (allel_count_t1[:, 1] / allel_count_t1.sum(axis=1)).astype(float)

    #Se crea un filtro para eliminar valores infinitos y NaN.
    mask_inf = np.isfinite(allel_mat_t0) & np.isfinite(allel_mat_t1)

    #Se establece un error en caso de que no hayan loci validos.
    if mask_inf.sum() == 0:
        raise ValueError("There are no loci with valid allel frequencies.")

    else:
        #Se aplica este filtro sobre las matrices.
        p_t0 = allel_mat_t0[mask_inf]
        p_t1 = allel_mat_t1[mask_inf]

    #Se realiza un nuevo filtrado del maf (minimum allel frequency) sobre la media entre los datos entre los dos tiempos.
    allel_mat_mean = (p_t0 + p_t1) / 2
    mask_maf = (allel_mat_mean >= maf) & (allel_mat_mean <= 1 - maf)

    #Se establece un error en caso de que no hayan loci validos.
    if mask_maf.sum() == 0:
        raise ValueError("There are not enough loci with valid allel frequencies.")

    else:
        #Se aplica este filtro sobre las matrices.
        p_t0 = p_t0[mask_maf]
        p_t1 = p_t1[mask_maf]
        allel_mat_mean = allel_mat_mean[mask_maf]

    #Se establece una condición que permite comprobar que hay suficientes loci para realizar el cálculo.
    if p_t0.size < loci_min or p_t1.size < loci_min:
        raise ValueError("There are not enough loci to estimate Ne")

    #Se ajustan las muestras según el número de cromosomas y se calcula la varianza esperada por muestreo.
    n_t0 = ploidy * len(t0_id)
    n_t1 = ploidy * len(t1_id)
    var_sample = allel_mat_mean * (1 - allel_mat_mean) * ((1.0 /  n_t0) + (1 / n_t1))

        #A continuación se calcula el cambio de frecuencias alélicas para cada muestra.
    delta = p_t1 - p_t0

    #Se calcula la varianza observada de este cambio de frecuencias y se ajusta con la varianza del muestreo.
    var_delta = delta ** 2
    var_delta = np.maximum(var_delta - var_sample, 0.0) #Si el resultado fuese negativo se pasa a cero para evitar errores.

    #Se calcula el valor de Ne por locus utilizando la ecuación 𝑁𝑒 = 𝑝(1−𝑝) / 2*(Var(Δ𝑝)/𝑡)
    with np.errstate(divide="ignore", invalid="ignore"): #Se ignoran los warnings producidos al dividir por cero.
      Ne_locus = (allel_mat_mean * (1 - allel_mat_mean)) / (2 * (var_delta / generations))
      Ne_locus = Ne_locus[np.isfinite(Ne_locus) & (Ne_locus > 0)] #Se eliminan los valores negativos y los infinitos.
      if Ne_locus.size == 0:
          raise ValueError("There are not enough valid per-locus Ne estimates")

    #Se crea un diccionario que contiene los valores necesarios para determinar la deriva genética.
    #Este diccionario requiere de la extracción de una serie de parámetros que serán utilizados posteriormente.

    #Se obtiene una matriz 2D de muestras frente a variantes.
    try:
        nalt = g_t0.subset(samples=t0_id).to_n_alt()
    except Exception:
        nalt = al.GenotypeArray(ds_t0["calldata/GT"].values).subset(samples=t0_id).to_n_alt()


    #Se añade una comprobación de los datos.
    if nalt.ndim != 2:
        raise ValueError("Nalt must be a 2D array (variants x samples).")

    #Se obtienen las posiciones alineadas a los loci filtrados.
    final_pos = None
    if "variants/POS" in ds_t0:
        all_pos = ds_t0["variants/POS"].value.astype(int)
        try:
            final_pos = all_pos[mask_inf]
            final_pos = final_pos[mask_maf]
        except Exception:
            final_pos = None

    #Se añade una comprobación de los datos.
    if final_pos.ndim != 1 or final_pos.shape[0] != nalt.shape[0]:
        raise ValueError("Positions must be a 1D array aligned with nalt.")

    #Se obtienen los ids asociados a los loci filtrados.
    final_var = None
    try:
        all_id = id_t0
        final_var = all_id[mask_inf]
        final_var = final_var[mask_maf].astype(str)
    except Exception:
        final_var = None

    #Se obtiene el número de muestras empleadas.
    if (isinstance(nalt, np.ndarray) or hasattr(nalt, "shape")):
        n_sample = int(nalt.shape[1])
    else:
        n_sample = None

    #Se cpnstruye el diccionario con los datos necesarios.
    Ne_LD = {"Ne_locus": np.asarray(Ne_locus),
             "pbar": np.asarray(allel_mat_mean),
             "delta": np.asarray(delta),
             "nalt": np.asarray(nalt),
             "positions": np.asarray(final_pos) if final_pos is not None else None,
             "variant_ids": np.asarray(final_var) if final_var is not None else None,
             "n_sample": n_sample,
             "sample_id_t0": np.asarray(t0_id, dtype=int),
             "sample_id_t1": np.asarray(t1_id, dtype=int)}

    #Se calcula la mediana.
    Ne_median = float(np.median(Ne_locus))

    #Se aplica un bootstrap para simular diversas distribuciones de los datos con lo que se obtiene un intervalo de confianza.
    boot = rng.integers(0, Ne_locus.size, size=(boot_rep, Ne_locus.size))
    boot_est = np.median(Ne_locus[boot], axis=1)
    low, med, up = np.percentile(boot_est, (ci[0], 50, ci[1]))

    #Finalmente se crea un Data Frame que contiene un registro de los valores obtenidos.
    Ne_df = pd.DataFrame({
        "Ne_estimate": [Ne_median],
        "median_boot": [float(med)],
        "low_CI": [float(low)],
        "up_CI": [float(up)],
        "n_loci": [int(Ne_locus.size)]})

    return Ne_df, Ne_LD, variants

#Una vez completada esta función se devuelven tres conjuntos de datos, dos de los cuales permiten una evaluación.
#El diccionario puede ser utilizado para calcular el desequilibrio de ligamiento que permite determinar la deriva genética.

def ne_ld(
    ne_dict: dict,
    method: str = "moderate_density",
    window: int = None,
    pairs: int = None,
    ploidy: int = 2,
    boot_rep: int = 500,
    ci: tuple[float, float] = (2.5, 97.5),
    seed: int = 10137) -> pd.DataFrame:
    """
    Estimates Ne from Linkage Disequilibrium (LD) using windowed pairing.

    Parameters
    ----------
    ne_dict : dict
        Dictionary containing quantities per locus. Requiered by this funcion to work properly.
    method : str, default = "moderate_density"
        Method used to calculate LD.
        Default value is "moderate_density" but can be specified.
    window_size : int, default = None
        Window size used to calculate LD.
        Default value is None but can be specified.
    pairs : int, default = None
        Number of pairs used to calculate LD.
        Default value is None but can be specified.
    ploidy : int, default = 2
        Number of chromosomes inherited by the descendants.
    boot_rep : int, default = 500
        Number of bootstrap replicates over loci. Enhance confidence intervals.
        Default value is 500 but can be specified.
    ci : tuple[float, float], default = (2.5, 97.5)
        Confidence interval.
        Default value is (2.5, 97.5) but can be specified.
    seed : int, default = 10137
        Seed for reproducibility.

    Returns
    -------
    pd.DataFrame
        Single-row DataFrame containing the estimate and summary statistics:
        - Ne : Estimate of effective population size (1 / (3 * r2_mean)).
        - Ne_median_boot : Median of bootstrapped Ne estimates.
        - Ne_low_CI :  Lower percentile of bootstrapped Ne according to `ci`.
        - Ne_up_CI : Upper percentile of bootstrapped Ne according to `ci`.
        - method :  Name of the method used.
        are provided while using a custom method name).
        - window : Window size in base pairs used.
        - pairs : Number of pairs to be analysed.
        - n_variants : Number of variants.
        - n_samples : Number of samples used.
        - n_pairs : Number of candidate pairs.
        - n_pairs_used : Number of pairs  used to compute r2 after filtering values.
        - r2_mean : Mean of r^2 values after sample-size bias correction.
        - bias : Bias term applied to r^2: 1 / (ploidy * n_sample).

    Example:
    --------
      >>> vcf_t0 = load_data("input_t0.vcf", "output_t0.zarr")
      >>> vcf_t1 = load_data("input_t1.vcf", "output_t1.zarr")
      >>> ne_df, ne_ld, variants = temporal_Ne(vcf_t0, vcf_t1,
                t0_sample = None, t1_sample = None,
                generations = 10, ploidy = 4,
                maf = 0.05, loci_min = 25, ci: = (0.5, 99.5),
                boot_rep = 1000, seed = 435)
      >>> ne_ld = ne_ld( ne_ld, method = "high_density",
                window = 50000, pairs = 4000,
                ploidy = 4, boot_rep = 100, ci: = (5, 95), seed: int = 435)
                
    """
    #Se crea un diccionario que acumula los métodos predefinidos y diseñados por el usuario de forma similar a módulos previos.
    methods = {
        "high_density": {"window": 50000, "pairs": 10000},
        "moderate_density": {"window": 250000, "pairs": 50000},
        "long_region": {"window": 1000000, "pairs": 200000}}

   #Se crean unas normas para hacer uso de este diccionario.
   if window is None and pairs is None:
        if method not in methods:
            raise ValueError(f"Method '{method}' is not valid. Accepted methods are: {list(methods.keys())}")
        window = methods[method]["window"]
        pairs = methods[method]["pairs"]

    elif window is not None and pairs is not None:
        if method in methods:
            raise ValueError("A predefined method cannot be selected when custom values are provided")
        else:
            methods[method] = {"window": window, "pairs": pairs}
            print(f"New method '{method}' has been registered:\n  • window = {window}\n  • pairs   = {pairs}")
            window = methods[method]["window"]
            pairs = methods[method]["pairs"]

    elif window is None or pairs is None:
        raise ValueError(f"Both window and pairs must be provided. Otherwise, use a predefined method: {list(methods.keys())}")

    #Se establecen los controles necesarios para las keywords utilizadas de manera similar a las funciones previas.
    if not (isinstance(ploidy, int) and ploidy >= 2):
        raise ValueError(f"Ploidy ({ploidy}) must be an integer with a value of 2 or above")

    if not (isinstance(boot_rep, int) and boot_rep >0):
        raise ValueError(f"Boot_rep ({boot_rep}) must be an integer above 0")

    if not (isinstance(ci, (list, tuple)) and len(ci) == 2 and 0 <= ci[0] < ci[1] <= 100):
        raise ValueError("ci must be a tuple with two values (lower | upper) between 0 and 100")

    #Se establece el rng.
    rng = np.random.default_rng(seed)

    #El primer paso consiste en tomar los valores almacenados en el diccionario de entrada.
    nalt = np.asarray(ne_dict["nalt"])
    pos = np.asarray(ne_dict["postitions"])
    n_sample = int(ne_dict["n_sample"])

    #Se calcula la mitad del tamaño de la ventana.
    half_window = int(window // 2)

    #Se crean dos listas vacias para los identificadores:
    i_id, j_id = [], []

    #Se itera por cada variante y se definen los limites de la ventana genómica centrado la variante i.
    for i in range(nalt.shape[0]):
        left = pos[i] - half_window
        right = pos[i] + half_window

    #Se buscan las variantes dentro de la ventana genómica.
        j_var = np.where((pos >= left) & (pos <= right))[0]

        #Se itera por cada variante dentro de ña ventana centrada en i.
        for j in j_var:
            if j <= i:
                continue #Se evita emparejar las variantes consigo mismas o con las anteriores.
            i_id.append(i)
            j_id.append(j)
            if len(i_id) >= pairs:
                break #Se detiene la busqueda cuando se encuentra el número de pares máximo establecido por el usuario.
        if len(i_id) >= pairs:
            break #se detiene el bucle externo cuando se alcanza este número.

    #Se comprueba que haya al menos un par de variantes en la ventana.
    if len(i_id) == 0:
        raise ValueError("There are not enough variants in the window. Use a bigger window size.")

    #Se conviertes estas listas en arrays de numpy para poder calcular la r2.
    i_id = np.array(i_id)
    j_id = np.array(j_id)

    #Se calcula el valor de r2 del desequilibrio de ligamientos para cada par de variantes.
    r2 = []
    for i, j in zip(i_id, j_id):
        try:
            r = al.rogers_huff_r2(nalt[i], nalt[j])
            r2.append(float(r ** 2))
        except Exception:
            continue

    #Se convierte la lista en un array de numpy y se conservan aquellos positivos.
    r2 = np.array(r2)
    r2 = r2[np.isfinite(r2) & (r2 > 0)]

    #Se comprueba que hay al menos un valor de r2.
    if r2.size == 0:
        raise ValueError("There are no valid r2 values.")

    #Se corrige el sesgo de muestreo.
    bias = 1 / (ploidy * n_sample)
    r2 = r2 - bias
    r2 = r2[np.isfinite(r2) & (r2 > 0)]

    if r2.size == 0:
        raise ValueError("There are no valid r2 values after bias correction.")

    #Se obtiene la media de los valores corregidos.
    r2_mean = float(np.mean(r2))

    #Se obtiene Ne a partir de la fórmula usando LD: 1/(3*r2).
    Ne = 1 / (3 * r2_mean)

    #Se aplica el bootstrap para comprobar estos valores de r2 de manera similar a la función previa.
    boot = rng.integers(0, r2.size, size=(boot_rep, r2.size))
    boot_est = np.mean(r2[boot], axis=1)

    #Se recalcula Ne para cada media del bootstrap.
    boot_ne = np.where(boot_est > 0, 1 / (3 * boot_est), np.nan) #Se sustituyen la smedias negativas por NaN.
    boot_ne = boot_ne[np.isfinite(boot_ne)] #Se eliminan los valores NaN.

    #Se comprueba que haya al menos un resultado válido.
    if boot_ne.size == 0:
        raise ValueError("There are no valid Ne estimates after bootstrap.")
    else:
      low, med, up = np.percentile(boot_ne, (ci[0], 50, ci[1]))

    #Se crea un Data Frame de pandas para almacenar los resultados.
    return pd.DataFrame({
        "Ne": Ne,
        "Ne_median_boot": med,
        "Ne_low_CI": low,
        "Ne_up_CI": up,
        "method": method,
        "window": window,
        "pairs": pairs,
        "n_variants": nalt.shape[0],
        "n_samples": n_sample,
        "n_pairs": len(i_id),
        "n_pairs_used": int(r2.size),
        "r2_mean": r2_mean
        "bias": bias,}, index=[0])
