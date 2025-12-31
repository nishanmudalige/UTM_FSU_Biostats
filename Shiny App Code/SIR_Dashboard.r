# app.R
# Deterministic SIR Dashboard (counts)

library(shiny)
library(bslib)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(deSolve)

# -----------------------------
# Deterministic SIR (counts)
# -----------------------------
sir_ode_counts <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dS <- -beta * S * I / N
    dI <-  beta * S * I / N - gamma * I
    dR <-  gamma * I
    list(c(dS, dI, dR))
  })
}

run_deterministic <- function(S0, I0, R0, beta, gamma, days, dt, N) {
  times <- seq(0, days, by = dt)
  state <- c(S = S0, I = I0, R = R0)
  parms <- c(beta = beta, gamma = gamma, N = N)
  
  out <- ode(
    y     = state,
    times = times,
    func  = sir_ode_counts,
    parms = parms,
    method = "lsoda"
  )
  
  as.data.frame(out) |>
    as_tibble() |>
    mutate(across(c(S, I, R), ~ pmax(.x, 0)))
}

# -----------------------------
# UI
# -----------------------------
ui <- page_fillable(
  theme = bs_theme(version = 5),
  
  layout_sidebar(
    sidebar = sidebar(
      width = 320,
      
      card(
        card_header("Inputs"),
        sliderInput("N", "Population (N)", min = 100, max = 500000, value = 10000, step = 100),
        sliderInput("I0", "Initial infected (I0)", min = 1, max = 5000, value = 25, step = 1),
        sliderInput("R0_init", "Initial recovered (R0)", min = 0, max = 5000, value = 0, step = 1),
        
        sliderInput("beta", "Transmission rate (β)", min = 0.01, max = 10, value = 0.40, step = 0.01),
        sliderInput("gamma", "Recovery rate (γ)", min = 0.01, max = 10, value = 0.10, step = 0.01),
        
        sliderInput("days", "Days", min = 10, max = 365, value = 200, step = 1),
        sliderInput("dt", "Time step (dt)", min = 0.1, max = 2, value = 0.25, step = 0.05)
      )
    ),
    
    # Main
    card(
      full_screen = TRUE,
      card_header("Deterministic SIR (Counts)"),
      
      plotOutput("sir_plot", height = "520px"),
      
      # Key metrics panel (spans the full width of the plot)
      card(
        class = "mt-3",
        card_header("Key metrics"),
        layout_columns(
          col_widths = c(4, 4, 4),
          
          card(
            card_header("Reproduction"),
            uiOutput("metric_r0"),
            uiOutput("metric_rt0")
          ),
          
          card(
            card_header("Peak infection"),
            uiOutput("metric_peak_i"),
            uiOutput("metric_t_peak")
          ),
          
          card(
            card_header("Final size"),
            uiOutput("metric_final_size"),
            uiOutput("metric_final_recovered")
          )
        )
      )
    )
  )
)

# -----------------------------
# Server
# -----------------------------
server <- function(input, output, session) {
  
  sim <- reactive({
    req(input$N, input$I0, input$R0_init, input$beta, input$gamma, input$days, input$dt)
    
    N  <- input$N
    I0 <- input$I0
    R0 <- input$R0_init
    S0 <- max(N - I0 - R0, 0)
    
    run_deterministic(
      S0 = S0, I0 = I0, R0 = R0,
      beta = input$beta, gamma = input$gamma,
      days = input$days, dt = input$dt, N = N
    )
  })
  
  metrics <- reactive({
    df <- sim()
    N  <- input$N
    
    peak_row <- df |> slice(which.max(I))
    
    list(
      R0_basic = input$beta / input$gamma,
      Rt0 = (input$beta / input$gamma) * (df$S[1] / N),
      
      peak_I = peak_row$I,
      t_peak = peak_row$time,
      
      final_size = N - df$S[nrow(df)],
      final_never_infected = df$S[nrow(df)]
    )
  })
  
  output$sir_plot <- renderPlot({
    df <- sim()
    
    long <- df |>
      pivot_longer(cols = c(S, I, R), names_to = "Compartment", values_to = "Count") |>
      mutate(Compartment = factor(Compartment, levels = c("S", "I", "R")))
    
    ggplot(long, aes(x = time, y = Count, color = Compartment)) +
      geom_line(linewidth = 1.2) +
      scale_y_continuous(labels = comma) +
      labs(x = "Time", y = "Count", color = NULL) +
      theme_minimal(base_size = 13) +
      theme(
        legend.position = "top",
        panel.grid.minor = element_blank()
      )
  })
  
  output$metric_r0 <- renderUI({
    m <- metrics()
    div(
      tags$div(style = "font-size: 28px; font-weight: 700;", sprintf("%.3f", m$R0_basic)),
      tags$div(style = "opacity: 0.8;", "R\u2080 = \u03B2 / \u03B3")
    )
  })
  
  output$metric_rt0 <- renderUI({
    m <- metrics()
    div(
      tags$div(style = "font-size: 22px; font-weight: 650;", sprintf("%.3f", m$Rt0)),
      tags$div(style = "opacity: 0.8;", "R\u209C at t = 0 (accounts for S\u2080/N)")
    )
  })
  
  output$metric_peak_i <- renderUI({
    m <- metrics()
    div(
      tags$div(style = "font-size: 28px; font-weight: 700;", comma(round(m$peak_I))),
      tags$div(style = "opacity: 0.8;", "Maximum I(t)")
    )
  })
  
  output$metric_t_peak <- renderUI({
    m <- metrics()
    div(
      tags$div(style = "font-size: 22px; font-weight: 650;", sprintf("%.2f days", m$t_peak)),
      tags$div(style = "opacity: 0.8;", "Time of peak infection")
    )
  })
  
  output$metric_final_size <- renderUI({
    m <- metrics()
    div(
      tags$div(style = "font-size: 28px; font-weight: 700;", comma(round(m$final_size))),
      tags$div(style = "opacity: 0.8;", "Total infected over the outbreak (N - S\u2091\u2099\u2090\u2097)")
    )
  })
  
  output$metric_final_recovered <- renderUI({
    m <- metrics()
    div(
      tags$div(style = "font-size: 22px; font-weight: 650;", comma(round(m$final_never_infected))),
      tags$div(style = "opacity: 0.8;", "Total never infected")
    )
  })
}

shinyApp(ui, server)
