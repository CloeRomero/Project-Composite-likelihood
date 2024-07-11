# Project-Composite-likelihood

En este repositorio se encuentra los códigos que realicé en el marco de mi trabajo de memoria.

En ellos implementé el método de la verosimilitud compuesta por pares de bloques (CL) propuesta por Eidsvik et al. en su paper ![Estimation and Prediction in Spatial Models With Block Composite Likelihoods]{https://www.researchgate.net/publication/262576230_Estimation_and_Prediction_in_Spatial_Models_With_Block_Composite_Likelihoods}, junto con el Algoritmo de Fisher que propone.

El método fue probado con dos enfoques distintos para la selección de la forma de los bloques, el enfoque por filas y el enfoque por columnas, propuestos por Rahul Mukerjee en su paper ![Improving upon the effective sample size based on Godambe information for block likelihood inference]{https://link.springer.com/article/10.1007/s00180-023-01328-6}.

A continuación, una breve descripción de los archivos.

- `codigoCL.R`: En éste archivo se encuentra la implementación del método la verosimilitud compuesta por pares de bloques para realizar estimaciones de los parámetros de un modelo de covarianza de Mátern para datos espaciales en una grilla regular. Esto para la configuración de bloques según el enfoque por filas

- `funciones.R`: En éste archivo se encuentran las funciones necesarias para la implementación del método, tales como: la función que clasifica los datos por bloques según el enfoque por filas, la que crea función de convarianza, la que invierte estas funciones de covarianzas y la función que contiene el algoritmo de Fisher.

- `codigoCL_col.R`: Al igual que el archivo `codigoCL.R`, se encuentra la implementación del método CL, pero esta vez adaptado a la configuración de bloques según el enfoque columnas.

- `funciones_col.R`: Éste archivo es el análogo a `funciones.R`, pero adaptado a la configuración de bloques según el enfoque columna. Varía la función que clasifica los datos por bloque y levemente el Algoritmo de Fisher.
