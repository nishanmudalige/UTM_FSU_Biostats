# app.R
# install.packages(c("shiny", "deSolve", "ggplot2", "dplyr"))

library(shiny)
library(deSolve)
library(ggplot2)
library(dplyr)
library(grid)   # for unit()

##---------------------------
## SIR with demography (fractions)
##---------------------------
sir_demographic <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dS <- mu - beta * S * I - mu * S
    dI <- beta * S * I - (gamma + mu) * I
    dR <- gamma * I - mu * R
    list(c(dS, dI, dR))
  })
}

simulate_SIR <- function(S0, I0, R0, parms, tmax = 500, n_steps = 1000) {
  state <- c(S = S0, I = I0, R = R0)
  times <- seq(0, tmax, length.out = n_steps)
  as.data.frame(ode(
    y     = state,
    times = times,
    func  = sir_demographic,
    parms = parms
  ))
}

## Helper to format complex eigenvalues nicely
format_complex <- function(z, digits = 3) {
  re <- Re(z); im <- Im(z)
  if (abs(im) < 1e-6) {
    return(sprintf(paste0("%.", digits, "f"), re))
  }
  sprintf(
    paste0("%.", digits, "f %s %.", digits, "fi"),
    re,
    ifelse(im >= 0, "+", "-"),
    abs(im)
  )
}

##---------------------------
## UI
##---------------------------
ui <- fluidPage(
  withMathJax(),
  titlePanel("SIR Phase Plane: Vector Field + Trajectories"),
  sidebarLayout(
    sidebarPanel(
      h4("Parameters"),
      sliderInput("beta",  "β (transmission)",
                  min = 0.1, max = 4,  value = 1.5, step = 0.1),
      sliderInput("gamma", "γ (recovery)",
                  min = 1/20, max = 1, value = 1/5, step = 0.01),
      sliderInput("mu",    "μ (birth/death)",
                  min = 0,   max = 0.1, value = 1/70, step = 0.001),
      
      hr(),
      h4("Simulation"),
      sliderInput("tmax", "Time horizon",
                  min = 50, max = 1500, value = 500, step = 50),
      sliderInput("grid_n", "Vector field grid density",
                  min = 10, max = 35, value = 20),
      
      hr(),
      h4("Initial conditions mode"),
      radioButtons(
        "init_mode", NULL,
        choices = c(
          "Random region"                    = "random",
          "Manual S values (per trajectory)" = "manual"
        ),
        selected = "random"
      ),
      
      # ---- Random mode controls ----
      conditionalPanel(
        condition = "input.init_mode == 'random'",
        sliderInput("n_traj", "Number of trajectories",
                    min = 1, max = 12, value = 6),
        br(),
        strong("Random S–I region"),
        sliderInput("S_min", "Min S", min = 0, max = 0.9,
                    value = 0.1, step = 0.05),
        sliderInput("S_max", "Max S", min = 0.1, max = 1,
                    value = 0.9, step = 0.05),
        sliderInput("I_min", "Min I", min = 0.001, max = 0.5,
                    value = 0.01, step = 0.005),
        sliderInput("I_max", "Max I", min = 0.005, max = 0.5,
                    value = 0.1, step = 0.005)
      ),
      
      # ---- Manual mode controls (individualized S0) ----
      conditionalPanel(
        condition = "input.init_mode == 'manual'",
        helpText("Set individual S0 values for each trajectory."),
        numericInput("n_manual", "Number of manual trajectories",
                     value = 4, min = 1, max = 10, step = 1),
        sliderInput("I_manual", "Initial I (same for all S0)",
                    min = 0.001, max = 0.5,
                    value = 0.01, step = 0.005),
        uiOutput("manual_S_ui")  # dynamic S0 sliders per trajectory
      ),
      
      hr(),
      strong("Basic reproduction number R0:"),
      verbatimTextOutput("R0_txt", placeholder = TRUE)
    ),
    
    mainPanel(
      plotOutput("phasePlot", height = "600px"),
      br(),
      h4("Jacobian at Endemic Equilibrium and Local Stability"),
      uiOutput("jacobian_info")
    )
  )
)

##---------------------------
## SERVER
##---------------------------
server <- function(input, output, session) {
  
  parms_reactive <- reactive({
    list(beta = input$beta, gamma = input$gamma, mu = input$mu)
  })
  
  output$R0_txt <- renderText({
    p <- parms_reactive()
    R0 <- p$beta / (p$gamma + p$mu)
    sprintf("R0 = %.3f", R0)
  })
  
  ## Dynamic sliders for manual S0 values
  output$manual_S_ui <- renderUI({
    req(input$init_mode == "manual")
    n <- input$n_manual
    sliders <- lapply(seq_len(n), function(j) {
      sliderInput(
        inputId = paste0("S_manual_", j),
        label   = paste("S0 for trajectory", j),
        min = 0.0, max = 0.99,
        value = j / (n + 1),
        step = 0.01
      )
    })
    do.call(tagList, sliders)
  })
  
  ## Reactive initial conditions (random vs manual)
  init_grid_reactive <- reactive({
    if (input$init_mode == "random") {
      # ---- Random region mode ----
      n  <- input$n_traj
      S0 <- runif(n, min = input$S_min, max = input$S_max)
      I0 <- runif(n, min = input$I_min, max = input$I_max)
      
      # enforce S + I <= 1
      I0 <- pmin(I0, 1 - S0 - 1e-3)
      I0[I0 < 0] <- 1e-3
      
      data.frame(S0 = S0, I0 = I0)
      
    } else {
      # ---- Manual S mode with individualized S0 ----
      n <- input$n_manual
      S0 <- numeric(0)
      for (j in seq_len(n)) {
        val <- input[[paste0("S_manual_", j)]]
        if (!is.null(val)) S0 <- c(S0, val)
      }
      S0 <- S0[S0 > 0 & S0 < 1]
      if (length(S0) == 0) {
        S0 <- c(0.2, 0.4, 0.6, 0.8)
      }
      
      I0 <- rep(input$I_manual, length(S0))
      
      # enforce S + I <= 1
      I0 <- pmin(I0, 1 - S0 - 1e-3)
      I0[I0 < 0] <- 1e-3
      
      data.frame(S0 = S0, I0 = I0)
    }
  })
  
  ##---------------------------
  ## Jacobian / eigenvalues at endemic equilibrium
  ##---------------------------
  jacobian_info_reactive <- reactive({
    p <- parms_reactive()
    beta  <- p$beta
    gamma <- p$gamma
    mu    <- p$mu
    
    R0 <- beta / (gamma + mu)
    
    if (R0 <= 1) {
      return(list(
        has_endemic = FALSE,
        R0          = R0
      ))
    }
    
    # Endemic equilibrium in S–I plane
    S_star <- (gamma + mu) / beta     # = 1/R0
    I_star <- mu * (R0 - 1) / beta
    
    # Jacobian of (S, I) subsystem:
    # dS/dt = mu - beta S I - mu S
    # dI/dt = beta S I - (gamma + mu) I
    J <- matrix(
      c(
        -beta * I_star - mu,   -beta * S_star,
        beta * I_star,        beta * S_star - (gamma + mu)
      ),
      nrow = 2, byrow = TRUE
    )
    
    eig <- eigen(J)
    vals <- eig$values
    
    # Classification
    re1 <- Re(vals[1]); re2 <- Re(vals[2])
    im1 <- Im(vals[1]); im2 <- Im(vals[2])
    
    eps <- 1e-6
    classification <- ""
    
    if (any(c(re1, re2) > eps)) {
      # At least one eigenvalue with positive real part
      if (abs(im1) > eps || abs(im2) > eps) {
        classification <- "Unstable spiral (focus)"
      } else if (re1 * re2 < -eps) {
        classification <- "Saddle (unstable)"
      } else {
        classification <- "Unstable node"
      }
    } else {
      # All real parts <= 0 -> stable or center
      if (abs(im1) > eps || abs(im2) > eps) {
        classification <- "Stable spiral (focus)"
      } else if (re1 * re2 < -eps) {
        classification <- "Saddle (unstable)"
      } else {
        classification <- "Stable node"
      }
    }
    
    list(
      has_endemic   = TRUE,
      R0            = R0,
      S_star        = S_star,
      I_star        = I_star,
      J             = J,
      eigenvalues   = vals,
      classification = classification
    )
  })
  
  output$jacobian_info <- renderUI({
    info <- jacobian_info_reactive()
    
    if (!info$has_endemic) {
      txt <- paste0(
        "<p>",
        "For the chosen parameters, \\(R_0 \\le 1\\), ",
        "so there is <strong>no endemic equilibrium</strong>.<br>",
        "The only equilibrium is the disease-free state ",
        "\\((S^*, I^*, R^*) = (1, 0, 0)\\), which is locally asymptotically ",
        "stable in this regime.",
        "</p>"
      )
      return(withMathJax(HTML(txt)))
    }
    
    J  <- info$J
    ev <- info$eigenvalues
    
    html <- paste0(
      "<p>",
      "For the chosen parameters, \\(R_0 = ",
      sprintf("%.3f", info$R0),
      " &gt; 1\\), so an endemic equilibrium exists.",
      "</p>",
      "<p>",
      "Endemic equilibrium (in the \\((S, I)\\) plane):<br>",
      "\\(S^* = ", sprintf("%.3f", info$S_star),
      ",\\ I^* = ", sprintf("%.3f", info$I_star), "\\)",
      "</p>",
      "<p>",
      "Jacobian of the \\((S, I)\\) subsystem at the endemic equilibrium:<br>",
      "\\[ J(S^*, I^*) = ",
      "\\begin{pmatrix}",
      sprintf("%.3f", J[1, 1]), " & ", sprintf("%.3f", J[1, 2]), " \\\\ ",
      sprintf("%.3f", J[2, 1]), " & ", sprintf("%.3f", J[2, 2]),
      "\\end{pmatrix} \\]",
      "</p>",
      "<p>",
      "Eigenvalues of \\(J(S^*, I^*)\\):<br>",
      "\\(\\lambda_1 = ", format_complex(ev[1]),
      ",\\ \\lambda_2 = ", format_complex(ev[2]), "\\)",
      "</p>",
      "<p>",
      "Classification of the endemic equilibrium: ",
      "<strong>", info$classification, "</strong>.",
      "</p>"
    )
    
    withMathJax(HTML(html))
  })
  
  ##---------------------------
  ## Phase-plane plot
  ##---------------------------
  output$phasePlot <- renderPlot({
    parms <- parms_reactive()
    beta  <- parms$beta
    gamma <- parms$gamma
    mu    <- parms$mu
    
    #---------------- Vector field ----------------
    S_vals <- seq(0.01, 0.99, length.out = input$grid_n)
    I_vals <- seq(0.01, 0.99, length.out = input$grid_n)
    
    field <- expand.grid(S = S_vals, I = I_vals) |>
      mutate(
        dS = mu - beta * S * I - mu * S,
        dI = beta * S * I - (gamma + mu) * I
      )
    
    max_len <- max(sqrt(field$dS^2 + field$dI^2))
    scale_factor <- 0.1 / max_len
    field <- field |>
      mutate(
        u = dS * scale_factor,
        v = dI * scale_factor
      )
    
    #---------------- Trajectories ----------------
    init_grid <- init_grid_reactive()
    
    traj_list <- lapply(seq_len(nrow(init_grid)), function(i) {
      S0 <- init_grid$S0[i]
      I0 <- init_grid$I0[i]
      R0 <- 1 - S0 - I0
      df <- simulate_SIR(
        S0, I0, R0,
        parms   = parms,
        tmax    = input$tmax,
        n_steps = 800
      )
      df$traj_id <- i
      df
    })
    traj <- bind_rows(traj_list)
    
    start_pts <- traj |>
      group_by(traj_id) |>
      slice(1)
    
    ggplot() +
      geom_segment(
        data = field,
        aes(x = S, y = I, xend = S + u, yend = I + v),
        arrow = arrow(length = unit(0.08, "cm")),
        linewidth = 0.3,
        alpha = 0.7
      ) +
      geom_path(
        data = traj,
        aes(x = S, y = I, group = factor(traj_id), colour = factor(traj_id)),
        linewidth = 0.9,
        alpha = 0.9
      ) +
      geom_point(
        data = start_pts,
        aes(x = S, y = I, colour = factor(traj_id)),
        size = 2.3
      ) +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      labs(
        x = "S (susceptible fraction)",
        y = "I (infected fraction)",
        colour = "Trajectory",
        title = "SIR Phase Plane: Vector Field + Trajectories"
      ) +
      theme_minimal()
  })
}

shinyApp(ui = ui, server = server)
