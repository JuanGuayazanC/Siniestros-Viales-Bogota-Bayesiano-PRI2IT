# ANÁLISIS PILOTO HISTÓRICO — no corresponde al modelo final del proyecto.
# Modelo Beta-Binomial para la proporción de siniestros fatales en Usme
#   Prior:         Bogotá 2015-2019 (varias calibraciones, ver más abajo)
#   Verosimilitud: Usme 2020  ->  X | p ~ Binomial(n, p)
#   Posterior:     p | X ~ Beta(a + x, b + n - x)
# Requiere: scripts/01_eda.R ejecutado (data/processed/siniestros_limpio.rds)
suppressMessages({library(dplyr); library(tidyr); library(ggplot2)})
set.seed(2026)
dir.create("outputs/tables",  showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)
theme_set(theme_minimal(base_size = 11))

LOCALIDAD <- "Usme"
ROPE_HALF <- 0.0025   # ±0,25 pp alrededor de p0: supuesto de "diferencia práctica nula"

s <- readRDS("data/processed/siniestros_limpio.rds")

# ---- Datos ------------------------------------------------------------------
hist_bog <- s %>% filter(anio <= 2019)
x_b <- sum(hist_bog$fatal); n_b <- nrow(hist_bog)            # Bogotá 2015-2019
p_b <- x_b / n_b
x <- s %>% filter(localidad == LOCALIDAD, anio == 2020) %>% summarise(sum(fatal)) %>% pull()
n <- s %>% filter(localidad == LOCALIDAD, anio == 2020) %>% nrow()
p0_2020 <- mean(s$fatal[s$anio == 2020])                     # Bogotá 2020 (referencia del contraste)
cat(sprintf("Bogotá 2015-19: %d/%d = %.5f | Bogotá 2020: %.5f | %s 2020: %d/%d = %.5f\n",
            x_b, n_b, p_b, p0_2020, LOCALIDAD, x, n, x / n))

# ---- Calibración del prior informativo --------------------------------------
# Tomar Beta(1 + x_b, 1 + n_b - x_b) usa 173.000 siniestros como "tamaño muestral previo": asume que
# TODAS las localidades tienen la misma p y aplasta a los datos de Usme (n = 597). Alternativa:
# media = p_b y concentración m = a + b estimada con la dispersión ENTRE localidades en 2015-2019
# (máxima verosimilitud marginal Beta-Binomial; Sumapaz se excluye por tener n = 5).
lbb <- function(x, n, a, b) lchoose(n, x) + lbeta(x + a, n - x + b) - lbeta(a, b)   # log m(x)
loc19 <- hist_bog %>% filter(localidad != "Sumapaz") %>%
  group_by(localidad) %>% summarise(n = n(), x = sum(fatal), .groups = "drop")
nll <- function(logm) -sum(lbb(loc19$x, loc19$n, p_b * exp(logm), (1 - p_b) * exp(logm)))
m_hat <- exp(optimize(nll, c(log(5), log(1e6)))$minimum)
cat(sprintf("Concentración entre localidades m = a + b = %.0f  (prior calibrado: Beta(%.1f, %.1f))\n",
            m_hat, p_b * m_hat, (1 - p_b) * m_hat))

priors <- tibble(
  prior = c("Calibrado (entre localidades)", "Bogotá 2015-19 completo", "Jeffreys Beta(0.5, 0.5)", "Uniforme Beta(1, 1)"),
  a     = c(p_b * m_hat, 1 + x_b, 0.5, 1),
  b     = c((1 - p_b) * m_hat, 1 + n_b - x_b, 0.5, 1))

# ---- Funciones de inferencia ------------------------------------------------
hpd_beta <- function(a, b, level = 0.95) {              # unimodal si a, b > 1; si no, cae a colas iguales
  if (a <= 1 || b <= 1) return(qbeta(c((1 - level) / 2, 1 - (1 - level) / 2), a, b))
  f <- function(p) qbeta(p + level, a, b) - qbeta(p, a, b)
  p <- optimize(f, c(0, 1 - level))$minimum
  qbeta(c(p, p + level), a, b)
}
resumen <- function(a, b, p0) {
  A <- a + x; B <- b + n - x
  h <- hpd_beta(A, B)
  rope <- pbeta(p0 + ROPE_HALF, A, B) - pbeta(p0 - ROPE_HALF, A, B)
  # Factor de Bayes: H0: p = p0  vs  H1: p ~ Beta(a, b)  (prior de H1 = el prior del modelo)
  logBF01 <- dbinom(x, n, p0, log = TRUE) - lbb(x, n, a, b)
  tibble(media = A / (A + B), mediana = qbeta(.5, A, B), moda = ifelse(A > 1 & B > 1, (A - 1) / (A + B - 2), NA),
         ic_lo = qbeta(.025, A, B), ic_hi = qbeta(.975, A, B), hpd_lo = h[1], hpd_hi = h[2],
         P_p_mayor_p0 = 1 - pbeta(p0, A, B), P_rope = rope,
         BF01 = exp(logBF01))
}

# P(H0|datos) con probabilidades previas 1/2-1/2 es BF01 / (1 + BF01)
tabla <- bind_rows(lapply(c(Bogota_2020 = p0_2020, Bogota_2015_19 = p_b), function(p0) {
  bind_rows(lapply(seq_len(nrow(priors)), function(i) cbind(prior = priors$prior[i], a = priors$a[i], b = priors$b[i],
                                                           resumen(priors$a[i], priors$b[i], p0)))) %>% mutate(p0 = p0)
}), .id = "referencia_p0") %>%
  mutate(P_H0_dado_datos = BF01 / (1 + BF01), across(where(is.numeric), ~ signif(.x, 5)))
write.csv(tabla, "outputs/tables/10_usme_posterior.csv", row.names = FALSE)
print(as.data.frame(tabla %>% filter(referencia_p0 == "Bogota_2020") %>%
  select(prior, a, b, media, hpd_lo, hpd_hi, P_p_mayor_p0, P_rope, BF01, P_H0_dado_datos)), row.names = FALSE)

# ---- Verificación por Monte Carlo (muestreo de la posterior) ----------------
i <- 1                                                     # prior calibrado
A <- priors$a[i] + x; B <- priors$b[i] + n - x
draws <- rbeta(1e5, A, B)
cat(sprintf("\nMonte Carlo (prior calibrado, 1e5 sim): media=%.5f, P(p>p0)=%.4f  [analítico: %.5f, %.4f]\n",
            mean(draws), mean(draws > p0_2020), A / (A + B), 1 - pbeta(p0_2020, A, B)))
# Predictiva posterior del número de fatales en n siniestros
pred <- rbinom(1e5, n, draws)
cat(sprintf("Predictiva posterior fatales en n=%d: media=%.1f, IC95%%=[%d, %d]; observado: %d\n",
            n, mean(pred), quantile(pred, .025), quantile(pred, .975), x))

# ---- Gráficos ----------------------------------------------------------------
grid <- seq(0, 0.06, length.out = 800)
dens <- bind_rows(lapply(seq_len(nrow(priors)), function(k) {
  A <- priors$a[k] + x; B <- priors$b[k] + n - x
  bind_rows(tibble(p = grid, d = dbeta(grid, priors$a[k], priors$b[k]), tipo = "Prior", prior = priors$prior[k]),
            tibble(p = grid, d = dbeta(grid, A, B), tipo = "Posterior", prior = priors$prior[k]))
})) %>% mutate(prior = factor(prior, priors$prior), tipo = factor(tipo, c("Prior", "Posterior")))
g1 <- ggplot(dens, aes(p, d, colour = tipo, linetype = tipo)) + geom_line(linewidth = .8) +
  geom_vline(xintercept = p0_2020, colour = "grey30", linetype = 3) +
  geom_vline(xintercept = x / n, colour = "firebrick", linetype = 2) +
  facet_wrap(~prior, scales = "free_y") + scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Proporción de siniestros fatales p", y = "Densidad", colour = NULL, linetype = NULL,
       title = sprintf("%s: prior y posterior según la elección del prior", LOCALIDAD),
       subtitle = "Punteada gris: Bogotá 2020 (p0). Roja discontinua: proporción observada en la localidad (2020)")
ggsave("outputs/figures/10_usme_prior_posterior.png", g1, width = 10, height = 6, dpi = 150)

A <- priors$a[1] + x; B <- priors$b[1] + n - x
g2 <- ggplot(tibble(p = draws), aes(p)) + geom_histogram(aes(y = after_stat(density)), bins = 60, fill = "#9ecae1", colour = NA) +
  stat_function(fun = dbeta, args = list(shape1 = A, shape2 = B), colour = "#08519c") +
  geom_vline(xintercept = p0_2020, linetype = 3) +
  labs(x = "p", y = "Densidad", title = "Posterior (prior calibrado): muestreo Monte Carlo vs densidad analítica")
ggsave("outputs/figures/10_usme_montecarlo.png", g2, width = 7, height = 4.5, dpi = 150)
cat("\nListo: outputs/tables/10_usme_posterior.csv y figuras 10_usme_*.png\n")
