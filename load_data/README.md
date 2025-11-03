**1.  _Carga y procesado de los datos._**

Este módulo contiene dos funciones: Una función está diseñada para cargar facilmente los datos siguiendo las indicaciones incluidas en los paquetes sgkit y scikit-allel mientras que la otra permite realizar un filtrado de los datos a decisión del usuario.

Siguiendo las indicaciones incluidas en los githubs de ambos creadores de estos paquetes, la instalación debe realizarse haciendo uso del comando [pip install] seguido del nombre del paquete correspondiente. Además será imprescindible cargar los paquetes en la interfaz de trabajo para lo cual se simplificará el acceso a los mismos haciendo uso de las siglas más adecuadas.

Además de los paquetes previamente mencionados, se ha incluido el paquete Xarray que aparece en el desglose de los paquetes de utilidad utilizados por sgkit y scikit-allel ya que estos trabajan preferentemente utilizando este tipo de estructura. Los paquetes de numpy y pandas deben incluirse para la acceder a algunas de sus funciones de interés como pueden ser en este caso .values() y la creación de data frames.

Las funciones incluidas en los paquetes scikit-allel y sgkit diseñadas para llevar a cabo un análisis análisis genómico permiten trabajar con numerosos formatos de datos, pero debido a la cantidad de información, así como la complejidad de la misma, se ha delegado en la conversión a ZARR, que genera bloques de información que permiten trabajar con una gran cantidad de datos sin requerir de una igualmente grande capacidad computacional.

En este módulo la carga de los datos se lleva a cabo haciendo uso de una función que comprueba el formato del fichero de entrada y hace uso de la función más apropiada para cada uno. Se ha incorporado una excepción única al formato PLINK ya que tanto sgkit como bio2zarr permiten trabajar con este tipo de formato y convertirlo a VCF ZARR pero no hay actualmente (5 septiembre) funciones para esto en python, aunque si están disponibles en Linux.

La función convert() de bio2zarr.vcf/.tskit permite realizar la conversión de estos formatos a VCF ZARR en un único paso, para lo cual solo es imprescindible incorporar el "path" a los ficheros de entrada y salida, el cual será creado en caso de no existir en el momento de utilizar el módulo.

Una vez realizada la conversión, el fichero .zarr puede ser cargado directamente utilizando la función load_dataset() de sgkit, lo que permite almacenar en una variable la información contenida en el fichero.

El módulo diseñado para el filtrado de los datos tiene una estructura considerablemente más compleja ya que requiere una serie de palabras clave que pemiten llevar a cabo el filtrado según las características propuestas por el usuario. El funcionamiento general consiste en la creación de filtros (masks) que reducen el tamaño del dataset según las características especificadas, pudiendo devolver un mensaje en cualquiera de los pasos si el filtrado ha sido demasiado estricto.

Dado que la función requiere de la entrada de numerosos datos por parte del usuario, antes de comenzar a trabajar con los datos se comprueban aquellas palabras clave que puedan inducir errores en el procesado de los datos, estipulando unos límites lógicos pero que devuelvan al usuario un mensaje claro sobre el fallo que puede haber ocurrido. Finalmente se devuelve el dataset filtrado además de dos diccionarios que contienen los filtros empleados y la cantidad de muestras y variantes que han quedado fuera de los mismos.

Para llevar a cabo el filtrado de los datos no se han empleado funciones de los paquetes instalados ya que no están incluidas, por lo que se ha llevado a cabo mediante la implementación de condiciones sobre una copia de los datos originales. 


___________________________________
2. Referencias:
*  [bio2zarr.vcf.convert()](https://sgkit-dev.github.io/bio2zarr/vcf2zarr/python_api.html)
*  [bio2zarr.tskit.convert()](https://sgkit-dev.github.io/bio2zarr/tskit2zarr/python_api.html)
*  [load_dataset()](https://sgkit-dev.github.io/sgkit/latest/generated/sgkit.load_dataset.html#sgkit.load_dataset)

No se referencian aquí las funciones de Numpy empleadas ya que estas pertenecen a la colección de funciones de python que no se encuentran relacionadas directamente con el trabajo.
