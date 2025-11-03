# Trabajo de fin de Máster
**1. _Introducción:_**

Repositorio de las funciones diseñadas para generar un software de análisis genómico orientado a la conservación de especies amenazadas.

Este repositorio contiene las funciones diseñadas para llevar a cabo un análisis del genome de los individuos que componen una o más poblaciones de especies que se encuentran en peligro de extinción, dentro de cualquiera de los niveles establecidos por la Unión Internacional por la Conservación de la Naturaleza (IUCN).

Las funciones diseñadas para este fin pretenden estar diseñadas de forma que puedan ser utilizadas incluso por equipos de investigación que no cuenten con una formación especializada en bioinformática. Esto significa que, de forma general, las funciones serán más simples y no tendran el alcance al que se podría llegar creando un paquete más profesional, ofreciendo a cambio una reducción del coste computacional y una simplificación de la salida de datos.

Para crear estas funciones se han utilizado los paquetes sgkit(1) y scikit-allel(2) que incluyen una gran variedad de funciones diseñadas para trabajar sobre la información genómica generada tras la secuenciación de los genomas. Estos paquetes habilitan el trabajo con ZARR (3) que es un proyecto diseñado para el almacenamiento en la nube de arrays multidimensionales de gran tamaño, como los que se generan al analizar genomas, así como el posterior trabajo con los mismos.

Siguiendo las indicaciones incluidas en los GitHub de los creadores de estos paquetes, la instalación debe realizarse haciendo uso del comando [pip install] seguido del nombre del paquete correspondiente. En cuanto a ZARR, será necesario incluir bio2zarr(4) que se encuentra dentro de sgkit y que permite realizar la conversión de una serie de tipos de fichero a este. Además es necesario instalar otros paquetes que permiten llevar a cabo el trabajo con los datos y llevar a cabo las funciones como pueden ser Xarray(5) y pomegranate(6).

**2. Contenido:_**

Este repositorio se encuentra dividido en carpetas que contienen los distintos módulos diseñados para llevar a cabo las distintas funciones además de ficheros README que contienen una breve explicación del funcionamiento, requisitos y salida de datos. A continuación se listan las distintas carpetas y su contenido:

  -  load_data.py: Contiene un módulo de carga de datos a partir de ficheros vcf o tslink. Se adjunta un README con una explicación detallada y el enlace a las funciones empleadas.
  -  inbreed.py: Contiene un módulo que hace uso de los ficheros cargados para calcular la consanguinidad de los individuos analizados. Incluye la posibilidad de calcular consanguinidad según pedigree. Se adjunta un README con una explicación detallada y el enlace a las funciones empleadas.
  -  roh.py: Contiene un módulo que hace uso de los ficheros cargados para determinar el tamaño y proporción del genoma ocupado por Runs of Homozygocity. Se adjunta un README con una explicación detallada y el enlace a las funciones empleadas.
  -  gendrift.py: Contiene un módulo que hace uso de los ficheros cargados para inferir el efecto de la deriva genética sobre los individuos. Este efecto se infiere haciendo uso del tamaño de población efectiva y el tamaño de población efectiva de acuerdo al desequilibrio de ligamientos (LD). Se adjunta un README con una explicación detallada y el enlace a las funciones empleadas.
  -  genflow.py: Contiene un módulo que hace uso de los ficheros cargados para inferir el flujo genético existente entre los individuos analizados. Se infiere a partir del análisis de la identidad por estado (IBS) de los genomas de los individuos analizados. Se adjunta un README con una explicación detallada y el enlace a las funciones empleadas.


**2. _Referencias:_**

  -  (1) Paquete sgkit -> https://sgkit-dev.github.io/sgkit/latest/index.html
  -  (2) Paquete scikit-allel -> https://github.com/cggh/scikit-allel/blob/master/docs/index.rst
  -  (3) Proyecto ZARR -> https://zarr.dev/
  -  (4) Funciones bio2zarr -> https://github.com/sgkit-dev/bio2zarr/tree/main/bio2zarr
  -  (5) Paquete Xarray -> https://pypi.org/project/xarray/
  -  (6) Paquete pomegranate -> https://pypi.org/project/pomegranate/

