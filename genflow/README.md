Este módulo hace uso de la información que se ha preparado previamente haciendo uso de la función load_data().

Estas funciones se ha diseñado con el objetivo de simplificar la determinación del flujo genético entre los individuos de una población. De forma habitual este estudio requeriría de numerosos datos de diversa índole pero en este caso se ha optado por construir un proxy del cálculo del Identity by Descendant (IBD) que consiste en obtener el valor de Identical by State(IBS) que es un indicador de lo similes que son dos secuencias. 

Este método no asegura que la similitud sea heredada pero puede ser útil en algunos datasets o como primera aproximación para estudios posteriores. Esta función simplifica el procesado de los datos y permite al usuario obtener una serie de datos relevantes para el proceso de conservación de la especie, aunque estos deben analizarse en conjunto com el resto de pruebas, tal y como se indica en el manual de usuario.

La primera función permite hacer uso del dataset cargado para llevar a cabo una determunación de la similitud entre los genomas de estos individuos. Los resultados de esta se devuelven como matriz de pares y como heatmap, lo que permite analizar visualmente los datos. Además se han incluido la posibilidad de filtrar unicamente los pares con una similitud superior a un umbral establecido por el usuario, de manera que se centre unicamente en aquellos puntos realmente relevantes.

La segunda función de este módulo consiste en una ampliación de esta primera función, permitiendo realizar unos análisis complementarios que aumentan la información generada a cambio de una mayor carga computacional. En este caso se han incluido un análisis de componentes principales (PCA), un análisis multidimensional (MDS) y un agrupamiento jerárquico que permite crear un árbol de distancias que muestra la cercanía entre los individuos.

En este módulo no se han utilizado funciones nuevas de los paquetes sgkit o scikit-allel, pero si se han incluido otras funciones de matplotlib, sklearn y scipy que permiten llevar a cabo estos cálculos.
