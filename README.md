# Letalidad vial por localidad en Bogotá: estimación bayesiana de la probabilidad de siniestro con muertos (2015–2020)

Proyecto Integrador 2 - Estadística Bayesiana (ECI Julio Garavito, 2026-2).
Modelos Beta-Binomial y Gamma-Poisson, con pruebas de bondad de ajuste, verificación en JAGS y modelo jerárquico.

## Informe principal
- `Analisis PRI2IT Siniestros Viales.Rmd`: análisis completo y reproducible (supuestos, prior histórica, Beta-Binomial por localidad, cuarentena de 2020, sensibilidad, validación predictiva, Gamma-Poisson, intervalos HPD, funciones de pérdida, ROPE y factor de Bayes, JAGS y extensión jerárquica).
- `docs/informe_siniestros_viales.html`: render anterior, **desactualizado** respecto al Rmd actual; volver a tejer antes de usarlo como informe final.
- `borrador_usme_JS.Rmd`: borrador previo centrado en Usme; queda como referencia.

## Datos
"Siniestros Viales Consolidados Bogotá D.C." (Secretaría Distrital de Movilidad, CC BY 4.0).
Cobertura real: **2015-01-01 a 2020-12-31** (196.152 siniestros; 422.416 actores viales).

1. Descargar el xlsx desde <https://datosabiertos.bogota.gov.co/dataset/siniestros-viales-consolidados-bogota-d-c>
   y guardarlo en `data/raw/siniestros_viales_consolidados_bogota_dc.xlsx` (no se versiona por su tamaño, ~49 MB).
2. Abrir el Rmd en RStudio y ejecutar o hacer Knit. Leer el xlsx tarda un par de minutos.

## Requisitos
- R con los paquetes `readxl`, `dplyr`, `tidyr`, `ggplot2`, `rjags` y `coda`.
- **JAGS** instalado en el sistema (<https://mcmc-jags.sourceforge.io>); `rjags` no funciona sin él.

## Scripts auxiliares
`scripts/00_importar.R`, `scripts/01_eda.R` y `scripts/02_eda_completo.R` corresponden al EDA exploratorio inicial. `scripts/03_beta_binomial_usme.R` y las salidas `10_usme_*` conservan un análisis piloto anterior (prior Bogotá 2015–2019 y Usme 2020); **no son los resultados del diseño final**. El análisis vigente está en el Rmd principal. `scripts/04_figuras_presentacion.R` genera figuras para Canva con el periodo y prior del diseño final.

## Hallazgos principales
- Proporción de siniestros con muertos: 1,53 % en Bogotá. No es constante entre años ni entre localidades.
- Usme (3,2 %), Ciudad Bolívar (2,6 %) y Rafael Uribe Uribe (2,3 %) están claramente por encima de Bogotá (1,40 % en 2018-2020); Kennedy es compatible con la ciudad.
- En la cuarentena de 2020 hubo menos siniestros pero más mortales (1,89 % frente a 1,40 %).
- Usme se eligió por tener la mayor proporción observada; el modelo jerárquico corrige ese sesgo de selección y mantiene la conclusión.
