#Se importan los paquetes, simplificando el acceso a los mismos.
import sgkit as sg
import allel as al
import xarray as xr
import pandas as pd
import numpy as np


#Se construye la función que permite calcular las Runs of Homozygosity (RoH).
def roh(dataset: xr.Dataset, mode: str = "standard", window: int = None, size: int = None):
    """
    Loads genomic data in VCF ZARR format and calculates the runs of homozygosity using a Hidden Markov Model of Poisson.
    Sizes are to be introduced by the user and can be changed if needed.

    Parameters:
    -----------
    dataset: xarray.Dataset
        Dataset in xarray format as returned from sg.load_dataset().
        Must contain:
            - Sample ID
            - Heterozygosity per sample ("sample_qc.het_rate")
            - Allelic frequency per variant ("allele_frequency")

    mode: str, default = "standard"
        Predefined values to be used by the Poisson HMM.
        Default value is "standard" (size = 25000 | min_size = 2500) but can be specified.

    window: int, optional
        Minimun size of the stretch to be considered a RoH.
        Default value is 2500 but can be specified. Numbers must be separated by a lower dash (underscore). Requires a size value to be used.

    size: int, optional
        Window size used by the Poisson HMM method to detect any homozygous stretches.
        Default value is 25000 but can be specified. Numbers must be separated by a lower dash (underscore). Requieres a window value to be used.

    Returns:
    --------
    pd.DataFrame
        Dataframe indexed by sample containing:
        - "nRoH" (number of RoHs)
        - "Genome size" (size in pb of the genome analysed)
        - "Length" (pb of the RoH)
        - "Proportion of RoH"(proportion of the genome occupied by RoH)

    Raises:
    -------
    RuntimeError
        If the genomic dataset could not be processed.

    ValueError
        If selected mode is not valid.
        If a mode is selected and sizes are manually included.
        If one of the sizes is missing.
        If one or both of the sizes exceed the genome size.

    Example:
    --------
    Using standarized parameters:
      >>> vcf = load_data("input.vcf", "output.zarr")
      >>> df = roh(vcf)
      >>> print(df)

    Using custom parameters:
      >>> vcf = load_data("input.vcf", "output.zarr")
      >>> df = roh(vcf, window = 100000, size = 3000)
      >>> print(df)

    """
#Se crea un diccionario que permita almacenar diversos perfiles personalizados según el interés del usuario.
    Method = {
        "standard": {"window": 25000, "size": 2500}}

#Se hace uso de un elif condition para comprobar si el usuario ha indicado uno de estos valores predefinidos o si ha incluido unos propios.
    if window is None and size is None: #En caso de que no se haya incluido ningún tamaño.
        if mode not in Method: #Si el método incluido no se encuentra en el diccionario creado anteriormente
            raise ValueError(f"Mode {mode} is not valid. Accepted modes are: {list(Method.keys())}") #Se lanza un error que indica que este modo no es valido.

#En caso de que no se indiquen los tamaños y el modo indicado sea correcto, se aplican los valores correspondientes.
        window = Method[mode]["window"]
        size = Method[mode]["size"]


    elif window is not None and size is not None: #En caso de que se indiquen los tamaños deseados.
        if mode in Method: #Si se indica uno de los metodos predefinidos se lanza un mensaje de incompatibilidad
            raise ValueError("A predefined mode cannot be selected when indicating desired values")

        else: #En caso de que el método indicado no se encuentre en el diccionario en el momento crea una clave en el diccionario personalizada para el usuario.
            Method[mode] = {"size": size, "window": window}
            print(f"New mode {mode}, has been registered:\n  • window      = {window}\n  • size  = {size}\n")

            window = Method[mode]["window"]
            size = Method[mode]["size"]

    elif window is None or size is None: #Esta condición comprueba que se han indicado las dos medidas requeridas.
            raise ValueError(f"Both sizes are requiered. Else a predefined mode can be used: {list(Method.keys())}")

#Una vez comprobados los tamaños y métodos incluidos se replica la selección de datos de manera similar al caso del módulo de consanguinidad.
    try:
#Se contean los alelos para cada variante y para cada muestra, lo que permite calcular las frecuencias alélicas así como las estadísticas de calidad.
        Data = sg.count_call_alleles(dataset)

#Se obtiene la tasa de heterocigosidad observada de cada individuo(H0,i).
        Data = sg.sample_heterozygosity(Data)

#Se calculan las frecuencias alélicas para obtener la frecuencia del alelo alternativo (p).
        Data = sg.call_allele_frequencies(Data)

#Se obtienen los valores requeridos para determinar las RoH.
        Genotype = al.GenotypeArray(Data["calldata/GT"].values) #Se crea una matriz con los datos.
        Positions = Data["variants/POS"].values #Se obtiene la posición genómica de cada variante.
        Samples = Data ["sample"].values #Se obtiene el número de cada muestra.

#Se obtiene el tamaño total del genoma a analizar haciendo uso de las posiciones inicial y final.
        Gen_size = int(Positions.max() - Positions.min())

#Se establece un error en caso de que uno de los tamaños indicados exceda el tamaño del genoma de trabajo.
        if window > Gen_size:
            raise ValueError(f"Indicated size ({window}) exceed genomic size ({Gen_size})")

        if size > Gen_size:
            raise ValueError(f"Indicated min_size ({size}) exceed genomic size ({Gen_size})")

    except Exception as e:
        raise RuntimeError(f"The genomic dataset could not be processed: {e}")
        #Se establece un RuntimeError utilizando la excepción para que se devuelva el error original acompañando a este mensaje.

#Se crea un diccionario vacio para almacenar los datos de los RoH.
    ROH = {}

#Se hace uso de un For Loop para iterar por cada posición de cada muestra
    for i, sample in enumerate(Samples):
        n = Genotype[:, i].to_n_alt() #Se obtienen los recuentos del alelo alternativo.

#Se utiliza el HMM de Markov para determinar los RoH.
        df_roh, froh = allel.roh_poissonhmm(
            n, Positions,
            phet_roh=0.001,             #Probabilidad de heterocigosis en el RoH. Depende de la tasa de mutación "de novo" y de la tasa de error del genotipo.
            phet_nonroh=(0.0025, 0.01), #Probabilidad de heterocigosis fuera del RoH. Depende de la diversidad nucleotídica de la población y de las tasas de mutación y de error del genotipo.
            transition=0.001,           #Probabilidad de transición entre estados. Depende del tamaño de la ventana seleccionada.
            window_size=window,         #Tamaño de la ventana para considerar el RoH.
            min_roh=size,               #Longitud mínima que permite considerar un RoH válido.
            contig_size=Gen_size        #Tamaño del fragmento con el que se está trabajando
        )

#Se actualiza el diccionario creado anteriormente extrayendo las metricas obtenidas.
        ROH[sample] = {
            "nRoH": len(df_roh),
            "Genome size (pb)": Gen_size,
            "Length (pb)": df_roh["length"].sum(),
            "Proportion of RoH": froh}

#Se crea un dataframe con estos datos.
    roh_df = pd.DataFrame.from_dict(ROH, orient = "index")
    roh_df.index.name = "sample"

    return roh_df
