Este módulo hace uso de la información que se ha preparado previamente haciendo uso de la función load_data().

Estas funciones se ha diseñado con el objetivo de simplificar la determinación del estado de deriva genética en los individuos de una misma población. De forma habitual este estudio requeriría de numerosos datos de diversa índole pero es posible realizar una inferencia del efecto de la deriva sobre la salud de la población mediante el cálculo del tamaño de población efectiva, lo que facilita el acceso a esta información.

En este módulo se han incluido dos funciones que permiten obtener una visión general haciendo uso de dos datasets a tiempos diferentes. La primera función realiza un cálculo más "superficial" mientras que la segunda hace uso de los valores extraidos en la primera para generar un dataframe que contiene los estadísticos que permiten evaluar esta característica.

Estas funciones requieren de una serie de palabras clave que permiten al usuario ajustar el cálculo a las particularidades de su especie de estudio. En el manual de usuario se incluyen algunos parámetros utilizados de manera general en este tipo de estudios pero no se limita su uso a estos. En el caso de la segunda función se ha incluido una estructura que permite almacenar modos personalizados.

Este módulo emplea, directa o indirectamente funciones de sgkit y scikit-allel que se han empleado en los módulos previos por lo que no serán anotadas en este apartado.
___________________________________
**Funciones:**  
  -  [allel.rogers_huff_r](https://scikit-allel.readthedocs.io/en/stable/stats/ld.html)
