Este módulo parte de la información que se ha preparado anteriormente haciendo uso de la función load_data().

Esta función se ha diseñado con el objetivo de ofrecer al usuario una visión clara y directa del efecto de la consanguinidad en los individuos de la especie de trabajo. Para esto, se ha hecho uso de las funciones incluidas en el paquete sgkit que permiten obtener los parámetros de interés que se pueden encontrar en el fichero vcf.zarr como pueden ser las frecuencias alélicas y los valores de heterocigosidad observada entre otros.

En base a estos parámetros la función puede calcular un valor de consanguinidad que oscila entre 0 y 1 y lo representa en una tabla acompañado de los valores de heterocigosidad y homocigosidad observadas, los cuales permiten realizar un análisis más profundo del efecto que tiene la consanguinidad sobre la(s) poblacion(es). El cálculo de la consanguinidad haciendo uso de la heterocigosidad observada y esperada (Fi = 1 - H0/HE) se ha llevado a cabo a mano para poder extraer estos valores y poder así construir la tabla objetivo.

El paquete sgkit incluye además una función que permite realizar el cálculo de la consanguinidad con los datos de pedigree (F = Σ[(1/2) n1 + n2 + 1 (1 + F*A)]) que es un método que permite obtener el valor de consanguinidad haciendo uso de valores parentales. Este método requiere que el dataset de partida incluya entre sus variables información relacionada con los parentales de cada individuo, incluyendo el número de cromosomas y la proporción de Identity by Descendant (IBD) que aporta a la descendencia, aunque estos valores solo son utilizados en uno de los métodos disponibles.

En conjunto, las funciones permiten generar un Data Frame que contiene valores de heterocigosidad, homocigosidad, consanguinidad y consanguinidad según pedigree, siendo esta última columna dependiente de la información incluida en el dataset de origen de forma que la propia función se adapta a la información ofrecida por el usuario y devuelve unos resultados adaptados a esto.

En principio los resultados ofrecidos se pueden interpretar directamente en si mismos, pero es importante realizar una relación entre los diversos valores para identificar así los individuos/poblaciones en mayor riesgo, lo que permite gestionar adecuadamente el proceso de conservación de la especie.

El workflow que sigue esta función es el siguiente:
1. Se cargan los datos y se introducen las keywords de interés (tipo de análisis y características del mismo) siempre y cuando se hayan incluido los datos parentales.
2. El programa extrae los valores de interés a partir de la información almacenada en el dataset.
3. En base a estos valores se obtiene la homocigosidad observada y la heterocigosidad esperada. Se calcula la consanguinidad y se genera un Data Frame con los valores obtenidos.
4. En caso de que existan datos de los parentales se lleva a cabo un segundo cálculo de la consanguinidad según pedigree haciendo uso de las funciones internas del paquete. El resultado se genera como Data Frame y se combina con el previo.
___________________________________
**Funciones:**
*  [sgkit.count_call_alleles()](https://sgkit-dev.github.io/sgkit/latest/generated/sgkit.count_call_alleles.html)
*  [sgkit.observed_heterozygosity()](https://sgkit-dev.github.io/sgkit/latest/generated/sgkit.observed_heterozygosity.html)
*  [sgkit.call_allele_frequencies()](https://sgkit-dev.github.io/sgkit/latest/generated/sgkit.call_allele_frequencies.html)
*  [sgkit.pedigree_inbreeding()](https://sgkit-dev.github.io/sgkit/latest/generated/sgkit.pedigree_inbreeding.html)
