library(shiny)
library(tidyverse)
library(DT)
library(splines)
library(glmnet)

# ----------------------------
# Helpers
# ----------------------------
mae <- function(y, yhat) mean(abs(y - yhat))

# Natural spline feature builder used in spline-stepwise formula
build_spline_step_features <- function(df) {
  df %>%
    mutate(
      logShucked = log(Shucked.Weight + 1),
      logShell   = log(Shell.Weight + 1),
      Length2    = Length^2,
      Height2    = Height^2,
      Shucked2   = Shucked.Weight^2,
      Shell2     = Shell.Weight^2,
      Len_Diam   = Length * Diameter,
      Ht_Wt      = Height * Weight,
      Shell_Shucked = Shell.Weight * Shucked.Weight
    )
}

# ----------------------------
# UI
# ----------------------------
ui <- navbarPage(
  "Abalone Age Explorer",
  tabPanel(
    "Data",
    sidebarLayout(
      sidebarPanel(
        checkboxGroupInput("sex_filter", "Sex", choices = c("F", "I", "M"),
                           selected = c("F","I","M")),
        sliderInput("age_filter", "Age (train only)", min = 1, max = 29, value = c(1, 29))
      ),
      mainPanel(
        DTOutput("data_table"),
        br(),
        verbatimTextOutput("quick_summary")
      )
    )
  ),
  
  tabPanel(
    "Relationships",
    sidebarLayout(
      sidebarPanel(
        selectInput("xvar", "X variable",
                    choices = c("Length","Diameter","Height","Weight",
                                "Shucked.Weight","Viscera.Weight","Shell.Weight"),
                    selected = "Shell.Weight"),
        checkboxInput("color_by_sex", "Color by Sex", TRUE),
        checkboxInput("show_lm", "Show linear fit (lm)", TRUE),
        checkboxInput("show_loess", "Show LOESS smoother", TRUE),
        sliderInput("alpha_pts", "Point transparency", min = 0.05, max = 1, value = 0.25)
      ),
      mainPanel(
        plotOutput("rel_plot", height = "520px")
      )
    )
  ),
  
  tabPanel(
    "Correlations",
    sidebarLayout(
      sidebarPanel(
        helpText("Correlation heatmap uses numeric variables only.")
      ),
      mainPanel(
        plotOutput("corr_heat", height = "520px"),
        br(),
        plotOutput("corr_age_bar", height = "420px")
      )
    )
  ),
  
  tabPanel(
    "Model Comparison",
    sidebarLayout(
      sidebarPanel(
        numericInput("seed", "Random seed", value = 123, min = 1),
        sliderInput("train_prop", "Training proportion", min = 0.6, max = 0.9, value = 0.8, step = 0.05),
        actionButton("fit_models", "Fit / Refit models"),
        hr(),
        selectInput("resid_model", "Residual plot model",
                    choices = c("Baseline", "Elastic Net", "Spline-stepwise"),
                    selected = "Spline-stepwise")
      ),
      mainPanel(
        plotOutput("mae_plot", height = "420px"),
        br(),
        plotOutput("resid_plot", height = "420px"),
        br(),
        DTOutput("mae_table")
      )
    )
  ),
  
  tabPanel(
    "Predict Age",
    sidebarLayout(
      sidebarPanel(
        selectInput("p_sex", "Sex", choices = c("F","I","M"), selected = "F"),
        numericInput("p_length", "Length", value = 1.3, min = 0),
        numericInput("p_diam", "Diameter", value = 1.0, min = 0),
        numericInput("p_height", "Height", value = 0.35, min = 0),
        numericInput("p_weight", "Weight", value = 23, min = 0),
        numericInput("p_shucked", "Shucked.Weight", value = 10, min = 0),
        numericInput("p_viscera", "Viscera.Weight", value = 5, min = 0),
        numericInput("p_shell", "Shell.Weight", value = 6.7, min = 0),
        hr(),
        selectInput("p_model", "Model", choices = c("Baseline","Elastic Net","Spline-stepwise"),
                    selected = "Spline-stepwise"),
        actionButton("predict_btn", "Predict")
      ),
      mainPanel(
        h3(textOutput("pred_out")),
        verbatimTextOutput("pred_note")
      )
    )
  )
)

# ----------------------------
# Server
# ----------------------------
server <- function(input, output, session) {
  
  train <- reactive({
    df <- read.csv("train2.csv")
    df$Sex <- factor(df$Sex)
    df
  })
  
  # Filtered data for EDA tabs
  train_f <- reactive({
    df <- train() %>%
      filter(Sex %in% input$sex_filter) %>%
      filter(Age >= input$age_filter[1], Age <= input$age_filter[2])
    df
  })
  
  output$data_table <- renderDT({
    datatable(train_f(), options = list(pageLength = 10, scrollX = TRUE))
  })
  
  output$quick_summary <- renderPrint({
    df <- train_f()
    list(
      n = nrow(df),
      age_summary = summary(df$Age),
      sex_counts = table(df$Sex)
    )
  })
  
  output$rel_plot <- renderPlot({
    df <- train_f()
    x <- input$xvar
    
    p <- ggplot(df, aes_string(x = x, y = "Age")) +
      geom_point(alpha = input$alpha_pts)
    
    if (isTRUE(input$color_by_sex)) {
      p <- ggplot(df, aes_string(x = x, y = "Age", color = "Sex")) +
        geom_point(alpha = input$alpha_pts)
    }
    
    if (isTRUE(input$show_lm)) {
      p <- p + geom_smooth(method = "lm", se = FALSE)
    }
    if (isTRUE(input$show_loess)) {
      p <- p + geom_smooth(method = "loess", se = FALSE)
    }
    
    p + labs(title = paste("Age vs", x))
  })
  
  output$corr_heat <- renderPlot({
    df <- train_f() %>%
      select(Age, Length, Diameter, Height, Weight, Shucked.Weight, Viscera.Weight, Shell.Weight)
    
    cm <- cor(df, use = "complete.obs")
    cor_long <- as.data.frame(as.table(cm))
    names(cor_long) <- c("Var1","Var2","Corr")
    
    ggplot(cor_long, aes(Var1, Var2, fill = Corr)) +
      geom_tile() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Correlation Heatmap (Numeric Variables)", x = "", y = "")
  })
  
  output$corr_age_bar <- renderPlot({
    df <- train_f() %>%
      select(Age, Length, Diameter, Height, Weight, Shucked.Weight, Viscera.Weight, Shell.Weight)
    cm <- cor(df, use = "complete.obs")
    age_cors <- sort(cm["Age", -1], decreasing = TRUE)
    tib <- tibble(var = names(age_cors), corr = as.numeric(age_cors))
    
    ggplot(tib, aes(reorder(var, corr), corr)) +
      geom_col() +
      coord_flip() +
      labs(title = "Correlation with Age (sorted)", x = "Variable", y = "Correlation")
  })
  
  # Store fitted models + MAEs
  fitted <- reactiveValues(
    split = NULL,
    baseline = NULL,
    enet = NULL,
    spline_step = NULL,
    mae_tbl = NULL
  )
  
  observeEvent(input$fit_models, {
    df <- train()
    
    set.seed(input$seed)
    n <- nrow(df)
    idx <- sample(seq_len(n), size = floor(input$train_prop * n))
    tr <- df[idx, ]
    te <- df[-idx, ]
    
    # ---------------- Baseline ----------------
    baseline_mod <- lm(
      Age ~ Sex + Length + Diameter + Height + Weight +
        Shucked.Weight + Viscera.Weight + Shell.Weight,
      data = tr
    )
    pred_base <- predict(baseline_mod, newdata = te)
    mae_base <- mae(te$Age, pred_base)
    
    # ---------------- Elastic Net ----------------
    # Model matrix (handles factor Sex)
    x_tr <- model.matrix(
      Age ~ Sex + Length + Diameter + Height + Weight +
        Shucked.Weight + Viscera.Weight + Shell.Weight,
      data = tr
    )[,-1]
    y_tr <- tr$Age
    
    x_te <- model.matrix(
      Age ~ Sex + Length + Diameter + Height + Weight +
        Shucked.Weight + Viscera.Weight + Shell.Weight,
      data = te
    )[,-1]
    
    # alpha=0.5 = Elastic Net balance; you can tune this later if you want
    cvfit <- cv.glmnet(x_tr, y_tr, alpha = 0.5)
    pred_en <- as.numeric(predict(cvfit, newx = x_te, s = "lambda.min"))
    mae_en <- mae(te$Age, pred_en)
    
    # ---------------- Spline-stepwise (your chosen spec) ----------------
    tr_tf <- build_spline_step_features(tr)
    te_tf <- build_spline_step_features(te)
    
    big <- lm(
      Age ~ Sex +
        ns(Diameter, df = 4) +
        ns(Weight, df = 4) +
        ns(Shell.Weight, df = 4) +
        logShucked + logShell +
        Length2 + Height2 + Shucked2 + Shell2 +
        Len_Diam + Ht_Wt + Shell_Shucked +
        Shucked.Weight + Viscera.Weight,
      data = tr_tf
    )
    # If you want *exact* stepwise selection every run, uncomment:
    # big <- step(big, direction = "both", trace = FALSE)
    
    pred_ss <- predict(big, newdata = te_tf)
    mae_ss <- mae(te$Age, pred_ss)
    
    mae_tbl <- tibble(
      Model = c("Baseline", "Elastic Net", "Spline-stepwise"),
      MAE = c(mae_base, mae_en, mae_ss)
    ) %>% arrange(MAE)
    
    fitted$split <- list(tr = tr, te = te, tr_tf = tr_tf, te_tf = te_tf)
    fitted$baseline <- baseline_mod
    fitted$enet <- cvfit
    fitted$spline_step <- big
    fitted$mae_tbl <- mae_tbl
  }, ignoreInit = FALSE)
  
  output$mae_plot <- renderPlot({
    req(fitted$mae_tbl)
    ggplot(fitted$mae_tbl, aes(reorder(Model, MAE), MAE)) +
      geom_col() +
      coord_flip() +
      labs(title = "Test MAE by Model", x = "Model", y = "MAE (years)")
  })
  
  output$mae_table <- renderDT({
    req(fitted$mae_tbl)
    datatable(fitted$mae_tbl, rownames = FALSE)
  })
  
  output$resid_plot <- renderPlot({
    req(fitted$split)
    te <- fitted$split$te
    te_tf <- fitted$split$te_tf
    
    if (input$resid_model == "Baseline") {
      pred <- predict(fitted$baseline, newdata = te)
      res <- te$Age - pred
      fit <- pred
      ttl <- "Residuals vs Fitted — Baseline"
    } else if (input$resid_model == "Elastic Net") {
      x_te <- model.matrix(
        Age ~ Sex + Length + Diameter + Height + Weight +
          Shucked.Weight + Viscera.Weight + Shell.Weight,
        data = te
      )[,-1]
      pred <- as.numeric(predict(fitted$enet, newx = x_te, s = "lambda.min"))
      res <- te$Age - pred
      fit <- pred
      ttl <- "Residuals vs Fitted — Elastic Net"
    } else {
      pred <- predict(fitted$spline_step, newdata = te_tf)
      res <- te$Age - pred
      fit <- pred
      ttl <- "Residuals vs Fitted — Spline-stepwise"
    }
    
    ggplot(tibble(fit = fit, res = res), aes(fit, res)) +
      geom_point(alpha = 0.25) +
      geom_hline(yintercept = 0) +
      labs(title = ttl, x = "Fitted Age", y = "Residual")
  })
  
  observeEvent(input$predict_btn, {
    req(fitted$baseline, fitted$enet, fitted$spline_step)
    
    new_df <- tibble(
      Sex = factor(input$p_sex, levels = levels(train()$Sex)),
      Length = input$p_length,
      Diameter = input$p_diam,
      Height = input$p_height,
      Weight = input$p_weight,
      Shucked.Weight = input$p_shucked,
      Viscera.Weight = input$p_viscera,
      Shell.Weight = input$p_shell
    )
    
    if (input$p_model == "Baseline") {
      yhat <- as.numeric(predict(fitted$baseline, newdata = new_df))
    } else if (input$p_model == "Elastic Net") {
      x_new <- model.matrix(
        ~ Sex + Length + Diameter + Height + Weight +
          Shucked.Weight + Viscera.Weight + Shell.Weight,
        data = new_df
      )[,-1]
      yhat <- as.numeric(predict(fitted$enet, newx = x_new, s = "lambda.min"))
    } else {
      new_tf <- build_spline_step_features(new_df)
      yhat <- as.numeric(predict(fitted$spline_step, newdata = new_tf))
    }
    
    output$pred_out <- renderText(sprintf("Predicted age: %.2f years", yhat))
    output$pred_note <- renderPrint({
      list(
        model_used = input$p_model,
        note = "Prediction is based on the model fitted in the Model Comparison tab (same session)."
      )
    })
  })
}

shinyApp(ui, server)
