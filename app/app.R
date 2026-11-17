# App Shiny local: letalidad vial por localidad en Bogotá (Beta-Binomial, 2018-2020 sin cuarentena)
# Replica el cálculo de scripts/04_figuras_presentacion.R (prior por momentos, sin Sumapaz).
# Cómo correrla (desde la raíz del proyecto):  shiny::runApp("app")

suppressMessages({library(shiny); library(dplyr); library(ggplot2)})

ruta_datos <- c("data/processed/siniestros_limpio.rds", "../data/processed/siniestros_limpio.rds")
ruta_datos <- ruta_datos[file.exists(ruta_datos)][1]
if (is.na(ruta_datos)) stop("No encuentro data/processed/siniestros_limpio.rds. Ejecuta scripts/00_importar.R y 01_eda.R.")

# ---- Datos y calibración (idéntico a scripts/04) ----------------------------
s <- readRDS(ruta_datos) %>%
  mutate(cuarentena = FECHA >= as.Date("2020-03-25") & FECHA <= as.Date("2020-08-31"))
hist_ <- s %>% filter(anio <= 2017)
rec   <- s %>% filter(anio >= 2018, !cuarentena)
p0    <- mean(rec$fatal)

th <- hist_ %>% filter(CODIGO_LOCALIDAD != 20) %>% group_by(CODIGO_LOCALIDAD) %>%
  summarise(n = n(), y = sum(fatal), p = y / n, .groups = "drop")
m_h   <- sum(th$y) / sum(th$n)
v_obs <- sum((th$p - m_h)^2) / (nrow(th) - 1)
v_bin <- mean(m_h * (1 - m_h) / th$n)
peso  <- m_h * (1 - m_h) / (v_obs - v_bin) - 1
aH    <- m_h * peso
bH    <- (1 - m_h) * peso

tab <- rec %>% filter(CODIGO_LOCALIDAD != 20) %>% group_by(localidad) %>%
  summarise(n = n(), y = sum(fatal), .groups = "drop") %>%
  mutate(a = aH + y, b = bH + n - y, media = a / (a + b),
         li = qbeta(.025, a, b), ls = qbeta(.975, a, b))

# ---- Reglas de decisión -----------------------------------------------------
K_OPC <- c(1.25, 1.5, 2)     # k: intervenir si p > k * p0
R_OPC <- c(2, 3, 5)          # c2/c1: cuánto cuesta omitir frente a intervenir

pct <- function(x, d = 2) paste0(formatC(100 * x, format = "f", digits = d, decimal.mark = ","), " %")

prob_exceso <- function(loc, k) {
  f <- tab[tab$localidad == loc, ]
  1 - pbeta(k * p0, f$a, f$b)
}
umbral <- function(r) 1 / (1 + r)          # c1 / (c1 + c2)
intervenir <- function(loc, k, r) prob_exceso(loc, k) > umbral(r)

grilla <- expand.grid(k = K_OPC, r = R_OPC)
conteo_escenarios <- function(loc) sum(mapply(intervenir, loc, grilla$k, grilla$r))

# ---- Interfaz ---------------------------------------------------------------
ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-family: 'Segoe UI', sans-serif; color: #434343; }
    .caja { background:#f3f3f3; border-radius:10px; padding:14px 18px; margin-bottom:12px; }
    .grande { font-size: 1.6em; font-weight: 600; color: #df632f; }
    table.tab { border-collapse: collapse; width: 100%; }
    table.tab th, table.tab td { border-bottom: 1px solid #ddd; padding: 6px 8px; text-align: left; }
    .si { color: #df632f; font-weight: 600; } .no { color: #7c898e; }
  "))),
  titlePanel("Letalidad vial por localidad en Bogotá (2018–2020)"),
  sidebarLayout(
    sidebarPanel(width = 4,
      selectInput("loc", "Localidad", choices = sort(tab$localidad), selected = "Usme"),
      selectInput("k", "k: qué tan exigente es el riesgo (intervenir si p > k · p0)",
                  choices = setNames(K_OPC, paste0("k = ", formatC(K_OPC, format = "f", digits = 2, decimal.mark = ",")
                                                   , " (umbral ", pct(K_OPC * p0), ")")), selected = 1.5),
      selectInput("r", "c₂/c₁: cuánto más cuesta omitir que intervenir",
                  choices = setNames(R_OPC, paste0("c₂/c₁ = ", R_OPC)), selected = 3),
      helpText(paste0("Bogotá (p0, 2018–2020 sin cuarentena): ", pct(p0)), br(),
               "Prior histórico: Beta(", formatC(aH, format = "f", digits = 2, decimal.mark = ","), "; ",
               formatC(bH, format = "f", digits = 2, decimal.mark = ","), "). Datos de 2015–2017.")
    ),
    mainPanel(width = 8,
      div(class = "caja", uiOutput("resumen")),
      plotOutput("densidad", height = 320),
      h4("Los 9 escenarios para esta localidad"),
      uiOutput("matriz"),
      h4("Todas las localidades con el escenario elegido"),
      uiOutput("todas")
    )
  )
)

# ---- Servidor ---------------------------------------------------------------
server <- function(input, output, session) {

  fila <- reactive({ tab[tab$localidad == input$loc, ] })

  output$resumen <- renderUI({
    f <- fila(); k <- as.numeric(input$k); r <- as.numeric(input$r)
    pe <- prob_exceso(input$loc, k)
    decision <- if (pe > umbral(r)) "Intervenir" else "No intervenir"
    tagList(
      div(paste0(input$loc, ": probabilidad estimada ", pct(f$media)),
          span(style = "font-size:.8em; font-weight:normal;",
               paste0("   (ICr 95 %: ", pct(f$li), " – ", pct(f$ls), ")"))),
      p(paste0("P(p > ", pct(k * p0), " | datos) = ", formatC(pe, format = "f", digits = 3, decimal.mark = ","),
               ".  Umbral de decisión = ", formatC(umbral(r), format = "f", digits = 2, decimal.mark = ","), ".")),
      div(class = "grande", decision)
    )
  })

  output$densidad <- renderPlot({
    f <- fila(); k <- as.numeric(input$k)
    xmax <- max(0.06, f$ls * 1.3)
    g <- seq(0, xmax, length.out = 600)
    dd <- bind_rows(
      data.frame(p = g, d = dbeta(g, aH, bH), tipo = "Prior"),
      data.frame(p = g, d = dbeta(g, f$a, f$b), tipo = "Posterior"))
    ggplot(dd, aes(p, d, colour = tipo, linetype = tipo)) +
      geom_line(linewidth = 1) +
      geom_vline(xintercept = p0, linetype = "dotted", colour = "#434343") +
      geom_vline(xintercept = k * p0, linetype = "dashed", colour = "#df632f") +
      annotate("text", x = p0, y = Inf, label = " Bogotá", hjust = 0, vjust = 1.5, colour = "#434343") +
      annotate("text", x = k * p0, y = Inf, label = paste0(" k·p0"), hjust = 0, vjust = 3, colour = "#df632f") +
      scale_colour_manual(values = c(Prior = "#7c898e", Posterior = "#df632f")) +
      scale_linetype_manual(values = c(Prior = "dotted", Posterior = "solid")) +
      scale_x_continuous(labels = function(x) paste0(formatC(100 * x, format = "f", digits = 1, decimal.mark = ","), " %")) +
      labs(x = "Probabilidad de siniestro con muertos (p)", y = "Densidad", colour = NULL, linetype = NULL) +
      theme_minimal(base_size = 14) + theme(legend.position = "top")
  })

  output$matriz <- renderUI({
    filas <- lapply(R_OPC, function(r) {
      tags$tr(tags$th(paste0("c₂/c₁ = ", r)),
              lapply(K_OPC, function(k) {
                si <- intervenir(input$loc, k, r)
                tags$td(class = if (si) "si" else "no", if (si) "Intervenir" else "No")
              }))
    })
    tagList(
      tags$table(class = "tab",
                 tags$thead(tags$tr(tags$th(""), lapply(K_OPC, function(k) tags$th(paste0("k = ", formatC(k, format = "f", digits = 2, decimal.mark = ",")))))),
                 tags$tbody(filas)),
      p(paste0("Intervenir en ", conteo_escenarios(input$loc), " de 9 escenarios."))
    )
  })

  output$todas <- renderUI({
    k <- as.numeric(input$k); r <- as.numeric(input$r)
    d <- tab %>% arrange(desc(media)) %>%
      mutate(P = mapply(prob_exceso, localidad, k),
             decision = ifelse(P > umbral(r), "Intervenir", "No"))
    tags$table(class = "tab",
      tags$thead(tags$tr(tags$th("Localidad"), tags$th("Media"), tags$th("ICr 95 %"),
                         tags$th(paste0("P(p > k·p0)")), tags$th("Decisión"))),
      tags$tbody(lapply(seq_len(nrow(d)), function(i) {
        tags$tr(tags$td(d$localidad[i]), tags$td(pct(d$media[i])),
                tags$td(paste0(pct(d$li[i]), " – ", pct(d$ls[i]))),
                tags$td(formatC(d$P[i], format = "f", digits = 3, decimal.mark = ",")),
                tags$td(class = if (d$decision[i] == "Intervenir") "si" else "no", d$decision[i]))
      })))
  })
}

shinyApp(ui, server)
