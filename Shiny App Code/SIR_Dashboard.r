# app.R
library(shiny)
library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(bslib)

# ---- Helpers ----
sir_ode <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dS <- -beta * S * I
    dI <-  beta * S * I - gamma * I
    dR <-  gamma * I
    list(c(dS, dI, dR))
  })
}

run_sir <- function(S0, I0, R0, beta, gamma, days, N) {
  sum0 <- S0 + I0 + R0
  if (!is.finite(sum0) || sum0 <= 0) sum0 <- 1
  S0n <- S0 / sum0; I0n <- I0 / sum0; R0n <- R0 / sum0
  
  y0 <- c(S = S0n, I = I0n, R = R0n)
  times <- seq(0, days, by = 0.1)
  parms <- c(beta = beta, gamma = gamma)
  
  out <- ode(y = y0, times = times, func = sir_ode, parms = parms, method = "rk4")
  out <- as.data.frame(out)
  
  out |>
    mutate(S_count = S * N, I_count = I * N, R_count = R * N)
}

nice_number <- function(x) {
  if (!is.finite(x)) return("")
  number(x, accuracy = ifelse(abs(x) >= 1000, 1, 0.01), scale_cut = cut_short_scale())
}

# ---- UI ----
ui <- page_fluid(  # bslib v5 layout
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  title = "SIR Modelling Dashboard",
  layout_sidebar(
    sidebar = sidebar(
      h4("Initial Conditions"),
      sliderInput("S0", "Initial susceptible (proportion)", min = 0, max = 1, value = 0.99, step = 0.001),
      sliderInput("I0", "Initial infectious (proportion)",  min = 0, max = 1, value = 0.01, step = 0.001),
      sliderInput("R0_init", "Initial recovered (proportion)", min = 0, max = 1, value = 0.00, step = 0.001),
      div(style = "margin-top:-8px; color:#666;", textOutput("sum_check")),
      tags$hr(),
      h4("Parameters"),
      sliderInput("beta",  HTML("&beta;: transmission rate (per day)"), min = 0.05, max = 1.5, value = 0.30, step = 0.01),
      sliderInput("gamma", HTML("&gamma;: recovery rate (per day)"),    min = 0.05, max = 1.0, value = 0.10, step = 0.01),
      sliderInput("days", "Simulation horizon (days)", min = 30, max = 360, value = 160, step = 5),
      numericInput("N", "Population size (for counts)", value = 100000, min = 1, step = 1000),
      checkboxInput("logy", "Log scale (y-axis)", value = FALSE),
      actionButton("reset", "Reset defaults"),
      tags$hr(),
      downloadButton("dl_csv", "Download CSV")
    ),
    card(
      card_header(h3("S, I, R over time", class = "m-0")),
      card_body(plotOutput("sir_plot", height = "480px"))
    ),
    # ---- HERE: Key Metrics + Notes cards ----
    layout_columns(
      col_widths = c(6, 6),
      card(
        card_header(h4("Key Metrics", class = "m-0")),
        card_body(
          fluidRow(
            column(6, div(class = "border rounded p-3 mb-2",
                          strong("Basic reproduction number (R0)"),
                          div(style = "font-size:1.4rem;", textOutput("R0txt")))),
            column(6, div(class = "border rounded p-3 mb-2",
                          strong("Mean infectious period"),
                          div(style = "font-size:1.4rem;", textOutput("infper"))))
          ),
          fluidRow(
            column(6, div(class = "border rounded p-3 mb-2",
                          strong("Peak infectious (count)"),
                          div(style = "font-size:1.4rem;", textOutput("ipeak")))),
            column(6, div(class = "border rounded p-3 mb-2",
                          strong("Time of peak (days)"),
                          div(style = "font-size:1.4rem;", textOutput("tpeak"))))
          )
        )
      ),
      card(
        card_header(h4("Notes", class = "m-0")),
        card_body(
          tags$ul(
            tags$li(HTML("ODE system: &dot;S = −&beta;SI, &nbsp; &dot;I = &beta;SI − &gamma;I, &nbsp; &dot;R = &gamma;I.")),
            tags$li("Proportions are normalized internally if they do not sum to 1."),
            tags$li(HTML("R0 = &beta; / &gamma;. Time units are days."))
          )
        )
      )
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {
  
  observeEvent(input$reset, {
    updateSliderInput(session, "S0", value = 0.99)
    updateSliderInput(session, "I0", value = 0.01)
    updateSliderInput(session, "R0_init", value = 0.00)
    updateSliderInput(session, "beta", value = 0.30)
    updateSliderInput(session, "gamma", value = 0.10)
    updateSliderInput(session, "days", value = 160)
    updateCheckboxInput(session, "logy", value = FALSE)
    updateNumericInput(session, "N", value = 100000)
  })
  
  output$sum_check <- renderText({
    s <- input$S0 + input$I0 + input$R0_init
    if (abs(s - 1) < 1e-6) sprintf("Sum: %.3f (OK)", s)
    else sprintf("Sum: %.3f — will be normalized.", s)
  })
  
  sol <- reactive({
    run_sir(
      S0 = input$S0, I0 = input$I0, R0 = input$R0_init,
      beta = input$beta, gamma = input$gamma,
      days = input$days, N = input$N
    )
  })
  
  output$R0txt  <- renderText(sprintf("%.3f", input$beta / input$gamma))
  output$infper <- renderText(sprintf("%.2f days", 1 / input$gamma))
  
  observe({
    df <- sol()
    i_idx <- which.max(df$I_count)
    output$ipeak <- renderText(nice_number(df$I_count[i_idx]))
    output$tpeak <- renderText(sprintf("%.1f", df$time[i_idx]))
  })
  
  output$sir_plot <- renderPlot({
    df <- sol() |>
      select(time, S = S_count, I = I_count, R = R_count) |>
      pivot_longer(-time, names_to = "Compartment", values_to = "Count")
    
    lab_fun <- label_number(scale_cut = cut_short_scale())
    
    ggplot(df, aes(x = time, y = Count, color = Compartment)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = c(S = "#1f77b4", I = "#d62728", R = "#2ca02c")) +
      labs(x = "Time (days)", y = "Number of individuals",
           title = sprintf("SIR dynamics (N = %s, R0 = %.2f)",
                           nice_number(input$N), input$beta / input$gamma)) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "top",
            plot.title.position = "plot",
            plot.title = element_text(face = "bold")) +
      { if (input$logy) scale_y_log10(labels = lab_fun)
        else scale_y_continuous(labels = lab_fun) }
  })
  
  output$dl_csv <- downloadHandler(
    filename = function() {
      paste0("sir_solution_N", input$N, "_R0_", sprintf("%.2f", input$beta / input$gamma), ".csv")
    },
    content = function(file) {
      df <- sol() |>
        select(time, S_prop = S, I_prop = I, R_prop = R, S_count, I_count, R_count)
      write.csv(df, file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
