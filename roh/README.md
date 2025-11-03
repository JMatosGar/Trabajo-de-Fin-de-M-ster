Este módulo parte de la información que se ha preparado anteriormente haciendo uso de la función load_data().

Esta función se ha diseñado con el objetivo de ofrecer al usuario la posibilidad de profundizar en la información ofrecida por la consanguinidad. Siendo este un módulo "optativo", se ha mantenido separado de la función de consanguinidad para conservar el objetivo de no generar una gran carga computacional en el lado del usuario.

Con esta función se determinan el número de Runs of Homozygosity, así como su tamaño total y la proporción que ocupan en el genoma de forma que el usuario puede inferir una relación entre eventos, naturales o artificiales, que pueden haber provocado cuellos de botella en la(s) poblacion(es).

El paquete sgkit no ofrece actualmente funciones que permiten determinar directamente estos valores, pero el paquete scikit-allele si contiene una función que permite realizar un Hidden Markov Model de Poisson que permite identificar estas RoHs. En todo caso, el paquete sgkit será necesario para obtener los valores requeridos por esta función a partir del dataset.

___________________________________
#Funciones:
*  [allel.GenotypeArray()](https://scikit-allel.readthedocs.io/en/stable/model/ndarray.html)
*  [allel.roh_poissonhmm()](https://scikit-allel.readthedocs.io/en/stable/stats/roh.html#)
