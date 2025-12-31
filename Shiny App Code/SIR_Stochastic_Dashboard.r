# app.R
# Changes:
# - Default for "Emphasize the most recently generated path" is now FALSE
# - Legend order forced to S, I, R via factor levels + scale_color_discrete(drop=FALSE)
# - Legend key lines made thicker via guides(... override.aes ...)
# - Layer controls moved below plot (above key statistics)
# - "Emphasize..." and "Log y-axis" moved into the same layer-controls panel, to the right of Deterministic (ODE)
# - Removed the empty panel to the right of the plot

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
  y0 <- c(S = S0, I = I0, R = R0)
  parms <- c(beta = beta, gamma = gamma, N = N)
  out <- ode(y = y0, times = times, func = sir_ode_counts, parms = parms, method = "rk4")
  as.data.frame(out)
}

# -----------------------------
# Stochastic SIR (Gillespie SSA)
# -----------------------------
gillespie_sir <- function(S0, I0, R0, beta, gamma, days, N, max_events = 5e6) {
  t <- 0
  S <- as.integer(S0); I <- as.integer(I0); R <- as.integer(R0)
  
  times <- numeric(1); Ss <- integer(1); Is <- integer(1); Rs <- integer(1)
  times[1] <- 0; Ss[1] <- S; Is[1] <- I; Rs[1] <- R
  
  ev <- 0L
  while (t < days && I > 0L && ev < max_events) {
    rate_inf <- beta * S * I / N
    rate_rec <- gamma * I
    rate_sum <- rate_inf + rate_rec
    if (!is.finite(rate_sum) || rate_sum <= 0) break
    
    dt <- rexp(1, rate_sum)
    t_new <- t + dt
    if (t_new > days) break
    
    u <- runif(1)
    if (u < rate_inf / rate_sum) {
      if (S > 0L) { S <- S - 1L; I <- I + 1L }
    } else {
      if (I > 0L) { I <- I - 1L; R <- R + 1L }
    }
    
    t <- t_new
    ev <- ev + 1L
    
    times <- c(times, t)
    Ss    <- c(Ss, S)
    Is    <- c(Is, I)
    Rs    <- c(Rs, R)
  }
  
  tibble(time = times, S = Ss, I = Is, R = Rs)
}

# Evaluate a step path on a common grid (LOCF)
step_to_grid <- function(path_df, grid_times) {
  idx <- findInterval(grid_times, path_df$time, left.open = FALSE, rightmost.closed = TRUE)
  idx[idx == 0] <- 1
  tibble(
    time = grid_times,
    S = path_df$S[idx],
    I = path_df$I[idx],
    R = path_df$R[idx]
  )
}

value_card_small <- function(title_ui, id, value_font = "1.10rem") {
  div(
    class = "border rounded p-3",
    div(style="font-weight:700; margin-bottom:2px; font-size:0.85rem; color:#555;", title_ui),
    div(style=paste0("font-size:", value_font, "; line-height:1.1;"), textOutput(id))
  )
}

# -----------------------------
# UI
# -----------------------------
ui <- page_fluid(
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  title = "Stochastic SIR Dashboard (Accumulating Paths)",
  
  tags$style(HTML("
    .plot-and-controls {
      display: flex;
      gap: 12px;
      align-items: stretch;
      width: 100%;
    }
    .plot-pane { flex: 1 1 auto; min-width: 0; }

    .controls-block-title { font-weight:700; margin-top:10px; margin-bottom:4px; }
    .tight-top { margin-top: 10px; }

    .actionbar {
      display: flex;
      align-items: flex-end;
      flex-wrap: wrap;
      gap: 0;
    }
    .actionbar .form-group { margin-bottom: 0; }
    .actionbar .vsep {
      width: 1px;
      height: 40px;
      background: #d0d0d0;
      margin: 0 14px;
    }
    .actionbar .nblock {
      display: flex;
      align-items: flex-end;
      gap: 10px;
    }
    .actionbar .nblock label {
      font-size: 0.8rem;
      color: #666;
      margin-bottom: 4px;
    }
  ")),
  
  layout_sidebar(
    sidebar = sidebar(
      numericInput("N", "Population size (N)", value = 1000, min = 10, step = 10),
      numericInput("I0", "Initial infectious (count)", value = 20, min = 1, step = 1),
      numericInput("R0", "Initial recovered (count)", value = 0, min = 0, step = 1),
      div(style="color:#666; font-size:0.9rem;", textOutput("S0_txt")),
      tags$hr(),
      
      h4("Rates"),
      sliderInput("beta",  HTML("&beta; (transmission per day)"), min = 0.05, max = 2.5, value = 0.30, step = 0.01),
      sliderInput("gamma", HTML("&gamma; (recovery per day)"),    min = 0.02, max = 1.5, value = 0.10, step = 0.01),
      tags$hr(),
      
      h4("Time grid"),
      sliderInput("days", "Horizon (days)", min = 30, max = 365, value = 160, step = 5),
      sliderInput("dt", "Plot time step (days)", min = 0.1, max = 2, value = 0.5, step = 0.1),
      numericInput("seed", "Base random seed", value = sample(1:100000, 1), min = 1, step = 1)
    ),
    
    card(
      card_body(
        
        # --- Plot (full width; no right-side empty panel) ---
        div(
          class = "plot-and-controls",
          div(
            class = "plot-pane",
            card(card_body(
              plotOutput("path_plot", height = "680px")
            ))
          )
        ),
        
        # --- Layer toggle controls panel (below plot, above key statistics) ---
        div(class = "tight-top",
            card(
              card_body(
                layout_columns(
                  col_widths = c(3, 3, 3, 3),
                  
                  div(
                    div(class="controls-block-title", "Sample paths"),
                    checkboxGroupInput(
                      "show_paths", label = NULL,
                      choices = c("Susceptible (S)"="S", "Infectious (I)"="I", "Recovered (R)"="R"),
                      selected = c("S","I","R")
                    )
                  ),
                  
                  div(
                    div(class="controls-block-title", "Running mean"),
                    checkboxGroupInput(
                      "show_mean", label = NULL,
                      choices = c("Susceptible (S)"="S", "Infectious (I)"="I", "Recovered (R)"="R")
                      # selected = c("S","I","R"),
                    )
                  ),
                  
                  div(
                    div(class="controls-block-title", "Deterministic (ODE)"),
                    checkboxGroupInput(
                      "show_det", label = NULL,
                      choices = c("Susceptible (S)"="S", "Infectious (I)"="I", "Recovered (R)"="R")
                      # selected = c("S","I","R")
                    )
                  ),
                  
                  div(
                    div(class="controls-block-title", "Visual Options"),
                    checkboxInput(
                      "highlight_latest",
                      "Most recent path",
                      value = FALSE
                    ),
                    checkboxInput("logy", "Log y-axis", value = FALSE)
                  )
                )
              )
            )
        ),
        
        # --- Key stats ---
        div(class="tight-top",
            layout_columns(
              col_widths = c(2, 2, 3, 3, 2),
              
              value_card_small("Paths generated", "npaths_txt"),
              
              div(
                class = "border rounded p-3",
                div(style="font-weight:700; margin-bottom:2px; font-size:0.85rem; color:#555;",
                    HTML("R<sub>0</sub> = &beta;/&gamma;")),
                div(style="font-size:1.10rem; line-height:1.1; font-weight:800;",
                    textOutput("R0_txt"))
              ),
              
              value_card_small("Estimated extinction (I hits 0)", "ext_txt"),
              value_card_small("Mean peak I (stochastic)", "peak_txt"),
              value_card_small("Time of mean peak I", "tpeak_txt")
            )
        ),
        
        # --- Actions ---
        div(class="tight-top",
            card(
              class = "mt-2",
              card_body(
                div(
                  class = "actionbar",
                  actionButton("draw", "Generate a new stochastic path", class = "btn-primary"),
                  div(class = "vsep"),
                  div(
                    class = "nblock",
                    numericInput("batch_n", "N paths", value = 10, min = 1, step = 1, width = "110px"),
                    actionButton("drawN", "Generate N stochastic paths", class = "btn-outline-primary")
                  ),
                  div(class = "vsep"),
                  actionButton("reset", "Reset", class = "btn-outline-secondary")
                )
              )
            )
        )
      )
    )
  )
)

# -----------------------------
# SERVER
# -----------------------------
server <- function(input, output, session) {
  
  S0_count <- reactive({
    S0 <- input$N - input$I0 - input$R0
    as.integer(max(0, S0))
  })
  
  rv <- reactiveValues(
    n = 0L,
    all_paths = NULL,
    latest_path_id = NA_integer_
  )
  
  output$npaths_txt <- renderText(as.character(rv$n))
  
  observeEvent(input$reset, {
    rv$n <- 0L
    rv$all_paths <- NULL
    rv$latest_path_id <- NA_integer_
  }, ignoreInit = TRUE)
  
  add_paths <- function(k) {
    validate(
      need(input$N > 0, "N must be positive."),
      need(input$I0 >= 1, "I0 must be at least 1."),
      need(input$R0 >= 0, "R0 must be nonnegative."),
      need(S0_count() + input$I0 + input$R0 == input$N, "Initial counts must sum to N (S0 implied)."),
      need(is.finite(k) && k >= 1, "Batch size must be at least 1."),
      need(k <= 500, "Batch size too large (cap is 500 to keep the app responsive).")
    )
    
    grid_times <- seq(0, input$days, by = input$dt)
    new_ids <- (rv$n + 1L):(rv$n + as.integer(k))
    
    new_df <- bind_rows(lapply(new_ids, function(id) {
      set.seed(input$seed + id)
      raw_path <- gillespie_sir(
        S0 = S0_count(), I0 = input$I0, R0 = input$R0,
        beta = input$beta, gamma = input$gamma,
        days = input$days, N = input$N
      )
      step_to_grid(raw_path, grid_times) |>
        mutate(path_id = id)
    }))
    
    rv$all_paths <- bind_rows(rv$all_paths, new_df)
    rv$n <- rv$n + as.integer(k)
    rv$latest_path_id <- rv$n
  }
  
  observeEvent(input$draw,  { add_paths(1) }, ignoreInit = TRUE)
  observeEvent(input$drawN, { add_paths(as.integer(input$batch_n)) }, ignoreInit = TRUE)
  
  observeEvent(TRUE, {
    session$sendInputMessage("draw", list(value = 1))
  }, once = TRUE)
  
  output$R0_txt <- renderText(sprintf("%.3f", input$beta / input$gamma))
  
  output$ext_txt <- renderText({
    req(rv$all_paths)
    max_t <- max(rv$all_paths$time)
    fin <- rv$all_paths |>
      group_by(path_id) |>
      filter(time == max_t) |>
      summarise(I_end = first(I), .groups = "drop")
    percent(mean(fin$I_end == 0), accuracy = 0.1)
  })
  
  output$peak_txt <- renderText({
    req(rv$all_paths)
    peaks <- rv$all_paths |>
      group_by(path_id) |>
      summarise(peakI = max(I), .groups = "drop")
    format(round(mean(peaks$peakI), 1), big.mark = ",")
  })
  
  output$tpeak_txt <- renderText({
    req(rv$all_paths)
    mean_I <- rv$all_paths |>
      group_by(time) |>
      summarise(meanI = mean(I), .groups = "drop")
    sprintf("%.1f", mean_I$time[which.max(mean_I$meanI)])
  })
  
  output$path_plot <- renderPlot({
    req(rv$all_paths)
    
    paths_comp <- input$show_paths
    mean_comp  <- input$show_mean
    det_comp   <- input$show_det
    
    validate(
      need(length(paths_comp) + length(mean_comp) + length(det_comp) >= 1,
           "Turn on at least one layer/compartment.")
    )
    
    comp_levels <- c("S", "I", "R")
    
    all_long_paths <- NULL
    if (length(paths_comp) > 0) {
      all_long_paths <- rv$all_paths |>
        select(time, path_id, all_of(paths_comp)) |>
        pivot_longer(cols = all_of(paths_comp), names_to = "Compartment", values_to = "Count") |>
        mutate(Compartment = factor(Compartment, levels = comp_levels))
    }
    
    mean_long <- NULL
    if (length(mean_comp) > 0) {
      mean_long <- rv$all_paths |>
        select(time, path_id, all_of(mean_comp)) |>
        pivot_longer(cols = all_of(mean_comp), names_to = "Compartment", values_to = "Count") |>
        group_by(time, Compartment) |>
        summarise(Count = mean(Count), .groups = "drop") |>
        mutate(Compartment = factor(Compartment, levels = comp_levels))
    }
    
    det_long <- NULL
    if (length(det_comp) > 0) {
      det <- run_deterministic(
        S0 = S0_count(), I0 = input$I0, R0 = input$R0,
        beta = input$beta, gamma = input$gamma,
        days = input$days, dt = input$dt, N = input$N
      )
      det_long <- det |>
        select(time, all_of(det_comp)) |>
        pivot_longer(cols = all_of(det_comp), names_to = "Compartment", values_to = "Count") |>
        mutate(Compartment = factor(Compartment, levels = comp_levels))
    }
    
    hist_long <- NULL
    latest_long <- NULL
    if (!is.null(all_long_paths) && is.finite(rv$latest_path_id)) {
      hist_long   <- all_long_paths |> filter(path_id != rv$latest_path_id)
      latest_long <- all_long_paths |> filter(path_id == rv$latest_path_id)
    } else if (!is.null(all_long_paths)) {
      hist_long <- all_long_paths
    }
    
    alpha_hist <- 0.16
    lw_hist    <- 0.45
    
    p <- ggplot() +
      labs(
        x = "Time (days)",
        y = "Count",
        title = "Stochastic SIR sample paths with running mean and deterministic overlay",
        subtitle = sprintf("Paths generated: %d  |  Latest path: %s",
                           rv$n, ifelse(is.finite(rv$latest_path_id), rv$latest_path_id, "NA")),
        color = NULL
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "bottom"
      ) +
      scale_color_discrete(drop = FALSE) +
      guides(color = guide_legend(override.aes = list(linewidth = 2.2)))
    
    if (!is.null(hist_long) && nrow(hist_long) > 0) {
      p <- p +
        geom_line(
          data = hist_long,
          aes(x = time, y = Count, group = interaction(path_id, Compartment), color = Compartment),
          linewidth = lw_hist,
          alpha = alpha_hist
        )
    }
    
    if (!is.null(latest_long) && nrow(latest_long) > 0) {
      if (isTRUE(input$highlight_latest)) {
        p <- p +
          geom_line(
            data = latest_long,
            aes(x = time, y = Count, group = Compartment, color = Compartment),
            linewidth = 1.25,
            alpha = 0.95
          )
      } else {
        p <- p +
          geom_line(
            data = latest_long,
            aes(x = time, y = Count, group = interaction(path_id, Compartment), color = Compartment),
            linewidth = lw_hist,
            alpha = alpha_hist
          )
      }
    }
    
    if (!is.null(mean_long)) {
      p <- p +
        geom_line(
          data = mean_long,
          aes(x = time, y = Count, group = Compartment, color = Compartment),
          linewidth = 1.8,
          alpha = 0.95
        )
    }
    
    if (!is.null(det_long)) {
      p <- p +
        geom_line(
          data = det_long,
          aes(x = time, y = Count, group = Compartment, color = Compartment),
          linewidth = 1.2,
          linetype = "dashed",
          alpha = 0.95
        )
    }
    
    if (isTRUE(input$logy)) {
      p + scale_y_log10(labels = label_number(scale_cut = cut_short_scale()))
    } else {
      p + scale_y_continuous(labels = label_number(scale_cut = cut_short_scale()))
    }
  })
}

shinyApp(ui, server)
