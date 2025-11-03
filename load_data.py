#Se importan los paquetes, simplificando el acceso a los mismos.
import sgkit as sg
import bio2zarr.vcf as v2z
import bio2zarr.tskit as t2z

#Se importan los paquetes adicionales necesarios para el filtrado de los datos.
import allel as al
import xarray as xr
import pandas as pd
import numpy as np

#Se construye una función que permita cargar los datos.

def load_data(input: str, output: str):
    """
    Loads genomic data and converts it to ZARR format.

    Parameters:
    -----------
    input: str
        Path to the input file. Support .vcf, .tskit and .vcf/.tskit.gz

    output: str
        Path to the output file. The file will be converted to ZARR format and will be created if it doesn't exist.

    Returns:
    --------
    sgkit Dataset
        Genomic dataset stored in ZARR and with Xarray structure.

    Raises:
    -------
    ValueError
        If the file type is not supported.

    Example:
    --------
    >>> vcf = load_data("input.vcf", "output.zarr")

    """
#Se establece una serie de elif conditions que permiten comprobar el tipo de fichero para hacer uso de la función correspondiente
    if input.endswith(".vcf") or input.endswith(".vcf.gz"): #Cuando se introduce un vcf o vcf.gz se utiliza v2z.convert().
        try:
            v2z.convert(input, output)
        except Exception as error:
            print(f"VCF file could not be converted: {error}")

    elif input.endswith(".tskit") or input.endswith(".tskit.gz"): #Cuando se introduce un tskit o tskit.gz se utiliza t2z.convert().
        try:
            t2z.convert(input, output)
        except Exception as error:
            print(f"TSKIT file could not be converted: {error}")

    elif input.endswith(".PLINK") or input.endswith(".PLINK.gz"): #Cuando se introduce un PLINK se lanza una excepción única.
        raise ValueError("PLINK files are not currently supported by this program")

    else:
        raise ValueError("The file type is not supported")

#A continuación se cargan los datos para que puedan ser utilizados por el usuario
    return sg.load_dataset(output)



#Se añade una función que permite relizar un preprocesado de los datos cargados para realizar un filtrado de los mismos.
def qc_data(dataset: xr.Dataset, maf_thr: float = 0.05, var_thr: float = 0.90, sample_thr: float = None, ploidy: int = 2, biallelic: bool = True, min_var: int = 10):
    """
    Preprocess genomic data in VCF ZARR format and creates a subset based on a series of personalised filters.
    Filters are to be introduced by the user and can be changed if needed.

    Parameters:
    -----------
    dataset: xarray.Dataset
        Dataset in xarray format as returned from sg.load_dataset().
        Must contain variables similar to those used by the RoH module:
            - "sample"
            - "calldata/GT"
            - "variants/POS"

    maf_thr: float, default = 0.05
        Minimum Alelic Frequency that has serves to filter the data.
        Default value is 0.05 but can be specified.

    var_thr: float, default = 0.05
        Minimum per Variant call rate requiered to retain a variant in order to filter the data.
        Default value is 0.90 but can be specified.

    sample_thr: float, default = None
        Minimum per Sample call rate requiered to filter the data.
        Default value is None but can be specified.
        If value is 0 or None, this filter is skipped.

    ploidy: int, default = 2
        Number of chromosomes inherited by the descendants.
        Default value is 2 (diploid) but can be specified.

    biallelic: bool, default = True
        Forces the filter of any multiallelic variant.
        False value allows the use of multiallelic variants.

    min_var: int, default = 10
        Minimum number of variants that allows a proper analysis.
        Default value is 10 but can be specified based on the characteristics of the genome and organism.

    Returns:
    --------
    Filtered sgkit Dataset
        Genomic dataset with Xarray structure filtered according to the user's parameters.

    Raises:
    -------
    RuntimeError
        If the genomic dataset could not be processed.

    ValueError
        If values are out of range.
        If ploidy value is negative or zero.
        If filtered dataset is empty.
        If filtered variants number is lower than designed threshold.

    Example:
    --------
    Using standarized parameters:
      >>> vcf = load_data("input.vcf", "output.zarr")
      >>> ds = qc_data(vcf)

    Using custom parameters:
      >>> vcf = load_data("input.vcf", "output.zarr")
      >>> ds = qc_data(vcf, maf_thr = 0.01, var_thr = 0.95, sample_thr = 0.95, ploidy = 4, biallelic = False, min_var = 30)
    """
#Antes de comenzar se comprueban los valores introducidos por el usuario.
    if not (0.0 <= maf_thr <= 1.0):
        raise ValueError(f"Minimum Allel Frecuency threshold ({maf_thr}) must be between 0 and 1")
    if not (0.0 <= var_thr <= 1.0):
        raise ValueError(f"Variant Call rate threshold ({var_thr}) must be between 0 and 1")
    if sample_thr is not None and not (0.0 <= sample_thr <= 1.0):
        raise ValueError(f"Sample Call rate ({sample_thr}) must be between 0 and 1 or None")
    if not (isinstance(ploidy, int) and ploidy >= 1):
        raise ValueError(f"Ploidy ({ploidy}) must be an integer with a value of 1 or above")
    if not isinstance(min_var, int) or min_var < 1:
        raise ValueError(f"Minimum number of variants requiered ({min_var}) must be an integer with a value of 1 or above")

#Se crea una copia del dataset para evitar una posible modificación de los datos
    ds = dataset.copy()

#Se replica la selección de datos de manera similar a los módulos previos.
    try:
#Se contean los alelos para cada variante y para cada muestra, lo que permite calcular las frecuencias alélicas así como las estadísticas de calidad.
        ds = sg.count_call_alleles(ds)
        ds = sg.call_allele_frequencies(ds)
        ds = sg.sample_heterozygosity(ds)

    except Exception as e:
        raise RuntimeError(f"The genomic dataset could not be processed: {e}")

#Se toma el número de variantes y muestras encontradas
    n_var = int(ds.dims["variants"])
    n_sam = int(ds.dims["samples"])

#Se realiza un filtrado de los datos en base a las características estipuladas por el usuario:
    try:
        #Se establece una condición que permite desactivar el filtrado de las muestras dejando como 0 o None el parámetro.
        if sample_thr is not None and sample_thr > 0:
            keep_samples = np.ones(n_sam, dtype = bool) #Se crea un array de booleanos.
            sample_call = ds["call_allele_counts"].sum(dim = ("variants", "alleles")).values.astype(float) #Se suman en un vector los conteos alélicos para cada variante y alelo.
            sample_call_rate = sample_call / (ploidy * n_var) #Se transforman los datos de alelos según el tipo de ploidía del organismo.

            #Se filtran las muestras cuyo conteo sea superior al umbral establecido.
            filt_samples = sample_call_rate >= sample_thr

            #Se comprueba que las muestras hayan pasado el filtrado.
            if filt_samples.sum() == 0:
                raise ValueError("No sample has passed the filter. Lowering the threshold might be requiered.")

            #Se seleccionan solo las muestras que han pasado el filtro.
            else:
                sample_ID = ds["sample"].values[filt_samples]
                ds = ds.sel(sample = sample_ID)

                #Se actualiza el listado de muestras.
                n_sam = int(ds.dims["samples"])
        else:
          filt_samples = np.ones(int(ds.dims["samples"]), dtype=bool)

        #A continuación se comprueban los alelos replicando la estructura previa.

        allel_freq = ds["call_allele_frequencies"] #Se extraen los valore de frecuencia, alelos y muestras.
        allel_count = ds["call_allele_counts"].sum(dim="alleles") #Se seleccionan los conteos alélicos.
        allel_call = allel_count > 0 #Se construye una matriz booleana que indica aquellos conteos con al menos una llamada no missing.
        variant_rate = (allel_call.sum(dim="samples") / int(ds.dims["samples"])).values.astype(float) #Cuantifica el ratio de muestras llamadas por cada variante

        #Se filtran las variantes cuyo conteo sea superior al umbral establecido.
        filt_variant = variant_rate >= var_thr

        #Se obtiene la frecuencia media de cada alelo para cada variante.
        allel_freq_mean = allel_freq.mean(dim="samples").values.astype(float)

        #Si las variantes son bialélicas se toma el valor de p del alelo alternativo.
        if allel_freq_mean.ndim == 2 and allel_freq_mean.shape[1] >= 2:
            max_af = allel_freq_mean.max(axis=1)
            maf = 1.0 - max_af

        #En caso de que las variantes no sean bialélicas se obtiene la frecuencia del alelo más frecuente y se utiliza para el cálculo.
        else:
            maf = np.zeros(allel_freq_mean.shape[0], dtype=float)

        #Se seleccionan unicamente las muestras que superen el umbral.
        filt_maf = maf >= maf_thr


        #Se comprueba si se ha establecido bialelismo obligatorio o no:
        if biallelic:
            if allel_freq_mean.ndim == 2:
                nonzero_allel = (allel_freq_mean > 0).sum(axis=1)
            else:
                nonzero_allel = np.zeros(allel_freq_mean.shape[0], dtype=int)

            bial = nonzero_allel == 2

        else:
            bial = np.ones_like(filt_maf, dtype=bool)

        #Finalmente se establecen las "mascaras", seleccionando unicamente las variables que cumplen con los requisitos que se han establecido.
        filt_data = filt_variant & filt_maf & bial
        post_filt = int(np.asarray(filt_data).sum())

        if post_filt == 0:
            raise ValueError("No variant has passed the filter. Lowering the threshold might be requiered.")

        if post_filt < min_var:
            raise ValueError(f"Number of variants selected {post_filt} is lower than the expected threshold {min_var}.")

        filt_sample_id = np.where(filt_data)[0]
        ds = ds.isel(variants=filt_sample_id)

        #Se crean dos diccionarios que contienen las máscaras y la variación en el número de muestras y variantes.
        masks = {"sample_mask": np.asarray(filt_samples, dtype=bool),
                 "var_mask": np.asarray(filt_variant, dtype=bool),
                 "maf_mask": np.asarray(filt_maf, dtype=bool),
                 "biallelic": np.asarray(bial, dtype=bool),}

        variation = {"n_samples": n_sam,
                     "n_samples_filt": int(np.asarray(filt_samples).sum()),
                     "n_var": n_var,
                     "n_var_filt": int(ds.dims["variants"])}

        return ds, masks, variation

    except Exception as e:
        raise RuntimeError(f"The genomic dataset could not be filtered: {e}")
