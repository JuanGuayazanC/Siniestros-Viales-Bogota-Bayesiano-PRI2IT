# Figuras para la presentación (paleta de la plantilla de Canva: gris #434343, naranja #df632f, amarillo #f1c930)
# Requiere data/processed/siniestros_limpio.rds (scripts/01_eda.R) y JAGS instalado.
suppressMessages({library(dplyr); library(ggplot2); library(rjags); library(coda)})
set.seed(2026)
dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)

GRIS <- "#434343"; NARANJA <- "#df632f"; AMARILLO <- "#f1c930"; AZUL <- "#4f6266"; FONDO <- "#f3f3f3"
tema <- theme_minimal(base_size = 15) +
  theme(plot.background = element_rect(fill = NA, colour = NA), panel.background = element_rect(fill = NA, colour = NA),
        text = element_text(colour = GRIS), axis.text = element_text(colour = GRIS), panel.grid.minor = element_blank(),
        legend.position = "none")

s <- readRDS("data/processed/siniestros_limpio.rds") %>%
  mutate(cuarentena = FECHA >= as.Date("2020-03-25") & FECHA <= as.Date("2020-08-31"))
hist_ <- s %>% filter(anio <= 2017)
rec   <- s %>% filter(anio >= 2018, !cuarentena)
p0    <- mean(rec$fatal)

# Prior histórica por localidad (método de los momentos, sin Sumapaz)
th <- hist_ %>% filter(CODIGO_LOCALIDAD != 20) %>% group_by(CODIGO_LOCALIDAD) %>% summarise(n = n(), y = sum(fatal), p = y / n, .groups = "drop")
m_h <- sum(th$y) / sum(th$n); v_obs <- sum((th$p - m_h)^2) / (nrow(th) - 1); v_bin <- mean(m_h * (1 - m_h) / th$n)
peso <- m_h * (1 - m_h) / (v_obs - v_bin) - 1; aH <- m_h * peso; bH <- (1 - m_h) * peso
cat(sprintf("Prior Beta(%.2f; %.2f), peso %.0f, p0 = %.4f\n", aH, bH, peso, p0))

tab <- rec %>% filter(CODIGO_LOCALIDAD != 20) %>% group_by(localidad) %>% summarise(n = n(), y = sum(fatal), .groups = "drop") %>%
  mutate(a = aH + y, b = bH + n - y, media = a / (a + b), li = qbeta(.025, a, b), ls = qbeta(.975, a, b))

# 1. Letalidad por año (y cuarentena)
grupos <- bind_rows(
  s %>% filter(anio <= 2019) %>% mutate(g = as.character(anio)),
  s %>% filter(anio == 2020, !cuarentena) %>% mutate(g = "2020\nsin cuar."),
  s %>% filter(cuarentena) %>% mutate(g = "Cuarentena\n2020")) %>%
  group_by(g) %>% summarise(n = n(), y = sum(fatal), .groups = "drop") %>%
  mutate(p = y / n, li = qbeta(.025, y + 1, n - y + 1), ls = qbeta(.975, y + 1, n - y + 1),
         g = factor(g, c(as.character(2015:2019), "2020\nsin cuar.", "Cuarentena\n2020")), cua = grepl("Cuarentena", g))
f1 <- ggplot(grupos, aes(g, p, colour = cua)) + geom_errorbar(aes(ymin = li, ymax = ls), width = .2, linewidth = 1) +
  geom_point(size = 4) + scale_colour_manual(values = c(`FALSE` = AZUL, `TRUE` = NARANJA)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = .1, decimal.mark = ","), limits = c(0.010, 0.025)) +
  labs(x = NULL, y = "Siniestros con muertos") + tema
ggsave("outputs/figures/pres_01_letalidad_anio.png", f1, width = 8, height = 4.95, dpi = 200, bg = "transparent")

# 2. Localidades: posterior (prior histórica) frente a Bogotá
tab <- tab %>% mutate(v = ifelse(li > p0, "por encima", ifelse(ls < p0, "por debajo", "compatible")))
cols_v <- c("por encima" = NARANJA, compatible = "#7c898e", "por debajo" = AZUL)
f2 <- ggplot(tab, aes(media, reorder(localidad, media), colour = v)) + geom_vline(xintercept = p0, linetype = "dashed", colour = GRIS) +
  geom_errorbarh(aes(xmin = li, xmax = ls), height = 0, linewidth = 1) + geom_point(size = 3) +
  scale_colour_manual(values = cols_v) + scale_x_continuous(labels = scales::percent_format(accuracy = .1, decimal.mark = ",")) +
  labs(x = "Probabilidad de siniestro con muertos (IC 95 %)", y = NULL, caption = "Línea discontinua: Bogotá (1,40 %)") + tema
ggsave("outputs/figures/pres_02_localidades.png", f2, width = 7.4, height = 7.0, dpi = 200, bg = "transparent")

# 3. Contrastes frente al resto de la ciudad (diferencia en puntos porcentuales)
NS <- 1e5; ty <- sum(rec$fatal); tn <- nrow(rec)
con <- do.call(rbind, lapply(seq_len(nrow(tab)), function(i) {
  yj <- tab$y[i]; nj <- tab$n[i]
  d <- rbeta(NS, aH + yj, bH + nj - yj) - rbeta(NS, 1 + ty - yj, 1 + (tn - nj) - (ty - yj))
  data.frame(localidad = tab$localidad[i], dif = 100 * mean(d), lo = 100 * quantile(d, .025), hi = 100 * quantile(d, .975), P = mean(d > 0))
})) %>% mutate(v = ifelse(P > .975, "por encima", ifelse(P < .025, "por debajo", "compatible")))
f3 <- ggplot(con, aes(dif, reorder(localidad, dif), colour = v)) + geom_vline(xintercept = 0, linetype = "dashed", colour = GRIS) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 1) + geom_point(size = 3) +
  scale_colour_manual(values = cols_v) +
  labs(x = "Diferencia con el resto de Bogotá (pp, IC 95 %)", y = NULL) + tema
ggsave("outputs/figures/pres_03_contrastes.png", f3, width = 6.7, height = 8.06, dpi = 200, bg = "transparent")

# 4. Trazas de JAGS para Usme (3 cadenas, valores iniciales dispersos)
usm <- rec %>% filter(CODIGO_LOCALIDAD == 5)
modelo <- "model { y ~ dbin(p, n); p ~ dbeta(a0, b0) }"
ini <- lapply(1:3, function(i) list(p = c(0.005, 0.03, 0.30)[i], .RNG.name = "base::Mersenne-Twister", .RNG.seed = 2026 + i))
m <- jags.model(textConnection(modelo), data = list(y = sum(usm$fatal), n = nrow(usm), a0 = aH, b0 = bH), inits = ini,
                n.chains = 3, n.adapt = 2000, quiet = TRUE)
update(m, 4000, progress.bar = "none")
mc <- coda.samples(m, "p", n.iter = 40000, progress.bar = "none")
cat("R-hat:", gelman.diag(mc)$psrf[1, 1], "| n_eff:", round(effectiveSize(mc)), "\n")
tr <- do.call(rbind, lapply(1:3, function(i) data.frame(iter = 1:1500, p = as.numeric(mc[[i]][1:1500, "p"]), cadena = factor(i))))
f4 <- ggplot(tr, aes(iter, p, colour = cadena)) + geom_line(alpha = .8, linewidth = .35) +
  scale_colour_manual(values = c(NARANJA, AMARILLO, AZUL)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1, decimal.mark = ",")) +
  labs(x = "Iteración", y = "p (Usme)") + tema
ggsave("outputs/figures/pres_04_trazas_jags.png", f4, width = 6.54, height = 3.95, dpi = 220, bg = "transparent")
cat("Figuras listas en outputs/figures/pres_*.png\n")
