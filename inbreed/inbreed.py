#Se importan los paquetes, simplificando el acceso a los mismos.
import sgkit as sg
import xarray as xr
import pandas as pd
import numpy as np

#Para calcular la consanguinidad se crea una función que contiene las funciones y cálculos necesarios.
def inbreed(dataset: xr.Dataset, method: str = "diploid", half_parent: bool = False) -> pd.DataFrame:
    """
    Loads genomic data in VCF ZARR format and calculates the values of Heterozygosity, Homozygosity and Inbreeding per sample.
    Allows the calculation of the pedigree values provided that parent values are included in the Xarray dataset.

    Parameters:
    -----------
    dataset: xarray.Dataset
        Dataset in xarray format as returned from sg.load_dataset().
        Must contain:
            - Sample ID
            - Heterozygosity per sample ("sample_qc.het_rate")
            - Allelic frequency per variant ("allele_frequency")
        Can contain:
            - Sample pedigree ("parent")

    method: str, default = "diploid"
        Method used by the pedigree inbreeding calculation to evaluate different chromosomic characteristics.
        Default value is "diploid" (diploid organisms).
        "Hamilton-Kerr" value is used to work with utopolyploids and mixed-ploids.



    half_parent: boolean, default = False
        Raises or not a ValueError when individuals are provided with just one parent.
        Default value is False (Raises ValueError)
        True value:
            - Using "diploid" value -- Unrecorded parents are considered to be unique founders unrelated to others.
            - Using "Hamilton-Kerr" value -- Tau and Lambda parameters are assumed to be equal to those of the known parent.

    Returns:
    --------
    pd.DataFrame
        Dataframe indexed by sample containing:
        - "Heterozygosity" (H0,i)
        - "Homozygosity" (1 - H0,i)
        - "Inbreeding" (Fi = 1 - H0,i/HE)
        - "Pedigree" (if parent values are provided) (F = Σ (1/2) n1 + n2 + 1 (1 + F*A))

    Raises:
    -------
    RuntimeError
        If the genomic dataset could not be processed.

    ValueError
        If the expected heterozygosity is 0 or above 1.
        If the input method is not "diploid" or "Hamilton-Kerr".

    Example:
    --------
    >>> vcf = load_data("input.vcf", "output.zarr")
    >>> df = inbreed(vcf)
    >>> print(df)

    >>> vcf = load_data("input.vcf", "output.zarr")
    >>> df = inbreed(vcf, "diploid")  || inbreed(vcf)
    >>> print(df)

    >>> vcf = load_data("input.vcf", "output.zarr")
    >>> df = inbreed(vcf, "Hamilton-Kerr")
    >>> print(df)
    """

#El primer paso es extraer los datos necesarios para poder determinar la heterocigosidad.
    try:
#Se contean los alelos para cada variante y para cada muestra, lo que permite calcular las frecuencias alélicas así como las estadísticas de calidad.
        Data = sg.count_call_alleles(dataset)

#Se obtiene la tasa de heterocigosidad observada de cada individuo(H0,i).
        Data = sg.sample_heterozygosity(Data)

#Se calculan las frecuencias alélicas para obtener la frecuencia del alelo alternativo (p).
        Data = sg.call_allele_frequencies(Data)

    except Exception as e:
        raise RuntimeError("The genomic dataset could not be processed")
        #Se establece un RuntimeError utilizando la excepción para que se devuelva el error original acompañando a este mensaje.

#A continuación se extrae el valor de heterocigosidad observada y se obtiene la homocigosidad para incluirla posteriormente en la tabla de salida.
    H0 = Data["sample_heterozygosity"] #Se selecciona la heterocigosidad de cada muestra.
    Hom0 = 1 - H0 #Se obtiene la homocigosidad como la diferencia entre la proporción de heterocigosidad y el total (1).

#Se obtiene la frecuencia del alelo alternativo (p) en cada posición utilizando las metricas extraidas previamente.
    p = Data["call_allele_frequencies"] #Para evitar asumir bialelismo se toma la matriz de frecuencias alélicas extraidas.

#Se calcula la heterocigosidad esperada por variante para cada muestra según la fórmula para loci multialélicos (HE = 1 - [Σ_i p_i,j^2])
    he = 1 - (p**2).sum(dim = "alleles")
    HEm = he.mean(dim = "variants", skipna = True) #Se obtiene la media de todas las variantes.
    HE = HEm.() #Se extrae el valor del Xarray y se transforma en float para evitar errores de calculo posteriores.

#Finalmente se calcula la consanguinidad de cada muestra utilizando la fórmula Fi = 1 - H0,i/HE.
#Se establece la condición de que la heterocigosidad esperada se encuentre entre 0 y 1.
    if (0 < HE <= 1):
        Fi = 1 - H0/HE
    else:
        raise ValueError(f"Expected heterozygosity ({HE}) cannot be 0 or above 1")

#Se crea un dataframe que contiene los datos obtenidos en este paso y que se ordenen según el número de muestra.
    df = pd.DataFrame({
        "Heterozygosity": H0.values,
        "Homozygosity": Hom0.values,
        "Inbreeding": Fi.values},
                      index = Data["sample"].values)

#Se establece el nombre del indice del dataframe, en este caso correspondiente a las muestras.
    df.index.name = "sample"

#En caso de disponer del pedigree en formato válido (Pandas Data frame), se incorpora el cálculo de la consanguinidad según el pedigree.
    if "parent" in dataset.coords or "parent" in dataset.data_vars:

#Se crea una serie de condciciones que permiten hacer uso de los distintos métodos incluidos en la función de inbreeding con pedigree.
         #Se crea una condición que compruebe que el método a utilizar es correcto.
         if method not in ["diploid", "Hamilton-Kerr"]:
            raise ValueError(f"Selected method {method} is not valid. Accepted methods are:\n  • diploid\n  • Hamilton-Kerr")

        #Si se utiliza el método diploide.
        if method == "diploid":
            ped_ds = sg.pedigree_inbreeding(
                Data,
                method = method, #Se aplica el método que haya sido seleccionado.
                parent = "parent",
                allow_half_founders = half_parent, #Se comprueba el permiso del usuario frente a parentales desconocidos.
                merge = False) #Se evita que se sobreescriba el dataset original.

        #Si se utiliza el método de Hamilton-Kerr.
        if method == "Hamilton-Kerr":

            #Se comprueba que están las variables requeridas en el dataset cargado.
            if "stat_Hamilton_Kerr_tau" not in Data or "stat_Hamilton_Kerr_lambda" not in Data:
                raise ValueError("Tau and Lambda parameters not found in dataset")

            ped_ds = sg.pedigree_inbreeding(
                Data,
                method = method, #Se aplica el método que haya sido seleccionado.
                parent = "parent",
                allow_half_founders = half_parent, #Se comprueba el permiso del usuario frente a parentales desconocidos.
                merge = False) #Se evita que se sobreescriba el dataset original.


#Se crea un data frame que contiene el valor del pedigree.
        ped_df = pd.DataFrame({"Pedigree": ped_ds["pedigree_inbreeding"].values},
                              index=Data["sample"].values)
        ped_df.index.name = "sample"

        df = df.join(ped_df, how = "left") #En caso de que haya una variable parents, se fusionan ambos data frames para incluir el valor del pedigree calculado.

    return df
