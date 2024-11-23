# Load required libraries
library(shiny)
library(ggplot2)
library(shinythemes)

# Function to check if a number is prime
isPrime <- function(n) {
  if (n <= 1) return(FALSE)
  if (n <= 3) return(TRUE)
  if (n %% 2 == 0 || n %% 3 == 0) return(FALSE)
  i <- 5
  while (i * i <= n) {
    if (n %% i == 0 || n %% (i + 2) == 0) return(FALSE)
    i <- i + 6
  }
  return(TRUE)
}

# Function to generate points on the elliptic curve
generatePoints <- function(a, b, p) {
  pts <- data.frame(x = integer(), y = integer())
  for (x in 0:(p - 1)) {
    rhs <- (x^3 + a * x + b) %% p
    for (y in 0:(p - 1)) {
      lhs <- (y^2) %% p
      if (lhs == rhs) {
        pts <- rbind(pts, data.frame(x = x, y = y))
      }
    }
  }
  return(pts)
}

# Extended GCD function for modular inverse
gcd_extended <- function(a, b) {
  if (b == 0) {
    return(list(gcd = a, x = 1, y = 0))
  } else {
    g <- gcd_extended(b, a %% b)
    x1 <- g$y
    y1 <- g$x - floor(a / b) * g$y
    return(list(gcd = g$gcd, x = x1, y = y1))
  }
}

# Function to compute modular inverse
modinv <- function(a, p) {
  g <- gcd_extended(a %% p, p)
  if (g$gcd != 1) {
    return(NULL)  # No inverse exists
  } else {
    return((g$x %% p + p) %% p)
  }
}

# Function for point addition on the elliptic curve
pointAdd <- function(P, Q, a, p) {
  if (is.null(P)) return(Q)
  if (is.null(Q)) return(P)
  if (P$x == Q$x && P$y == (-Q$y %% p)) {
    # P + (-P) = O (point at infinity)
    return(NULL)
  }
  if (P$x == Q$x && P$y == Q$y) {
    # Point doubling
    s_num <- (3 * P$x^2 + a) %% p
    s_den <- (2 * P$y) %% p
  } else {
    # Point addition
    s_num <- (Q$y - P$y) %% p
    s_den <- (Q$x - P$x) %% p
  }
  s_inv <- modinv(s_den, p)
  if (is.null(s_inv)) {
    # Denominator inverse does not exist
    return(NULL)
  }
  s <- (s_num * s_inv) %% p
  x_r <- (s^2 - P$x - Q$x) %% p
  y_r <- (s * (P$x - x_r) - P$y) %% p
  return(list(x = x_r %% p, y = y_r %% p))
}

# UI code
ui <- fluidPage(
  theme = shinytheme("cerulean"),
  titlePanel("Elliptic Curve Point Generator"),
  tags$head(
    tags$style(HTML("
      body {
        background-color: #f5f5f5;
      }
      .shiny-input-panel {
        background-color: #fff;
        padding: 20px;
        border-radius: 5px;
      }
      #ellipticCurvePlot {
        border: 1px solid #ddd;
        background-color: #fff;
      }
      .btn {
        margin-top: 10px;
      }
      #currentPointInfo {
        font-size: 16px;
        font-weight: bold;
      }
    "))
  ),
  sidebarLayout(
    sidebarPanel(
      numericInput("a", "Parameter a:", value = 2, min = 0),
      numericInput("b", "Parameter b:", value = 3, min = 0),
      numericInput("p", "Prime p:", value = 17, min = 2),
      hr(),
      withMathJax(
        uiOutput("equationText")
      ),
      actionButton("showPoints", "Select Curve", class = "btn-primary"),
      uiOutput("curveInfo")
    ),
    mainPanel(
      uiOutput("generatorPointUI"),
      tags$br(),
      tags$br(),
      verbatimTextOutput("currentPointInfo"),
      plotOutput("ellipticCurvePlot", height = "600px")
    )
  )
)

# Server code
server <- function(input, output, session) {
  output$equationText <- renderUI({
    req(input$a, input$b, input$p)
    withMathJax(sprintf("$$ y^2 \\equiv x^3 + %dx + %d \\mod %d $$", input$a, input$b, input$p))
  })
  
  points <- reactiveVal()
  initialG <- reactiveVal()
  currentPoint <- reactiveVal()
  kValue <- reactiveVal(0)
  pathPoints <- reactiveVal(data.frame())
  pathLines <- reactiveVal(data.frame())
  
  observeEvent(input$reset, {
    initialG(NULL)
    currentPoint(NULL)
    kValue(0)
    pathPoints(data.frame())
    pathLines(data.frame())
    
    # Re-render the plot to show only the base points
    output$ellipticCurvePlot <- renderPlot({
      if (!is.null(points())) {
        plotPoints(points(), input$a, input$b, input$p)
      }
    })
    
    # Clear generator point selection and current point info
    output$generatorPointUI <- renderUI({
      if (!is.null(points())) {
        tagList(
          h3("Point Generation"),
          selectInput("generatorPoint", "Select Generator Point G:",
                      choices = sapply(1:nrow(points()), function(i) paste0("(", points()$x[i], ",", points()$y[i], ")"))),
          actionButton("startGenerating", "Start/Reset Generating", class = "btn-success"),
          actionButton("nextPoint", "Next Point", class = "btn-info")
        )
      }
    })
    
    output$currentPointInfo <- renderText({ NULL })
  })
  
  
  observeEvent(input$showPoints, {
    a <- input$a
    b <- input$b
    p <- input$p
    
    # Validate that p is prime
    if (!isPrime(p)) {
      showNotification("p must be a prime number.", type = "error")
      return()
    }
    
    # Check discriminant
    discriminant <- (4 * a^3 + 27 * b^2) %% p
    if (discriminant == 0) {
      showNotification("The discriminant is zero. The curve is singular.", type = "error")
      return()
    }
    
    # Generate points
    pts <- generatePoints(a, b, p)
    if (nrow(pts) == 0) {
      showNotification("No valid points found for given parameters.", type = "error")
      return()
    }
    points(pts)
    
    # Display curve information
    output$curveInfo <- renderUI({
      tagList(
        h4("Curve Information:"),
        p(sprintf("Discriminant (Δ): %d", discriminant)),
        p(sprintf("Characteristic (p): %d", p)),
        p(sprintf("Number of Points on Curve: %d", nrow(pts)))
      )
    })
    
    # Update the plot
    output$ellipticCurvePlot <- renderPlot({
      plotPoints(pts, a, b, p)
    })
    
    # Update generator point selection
    output$generatorPointUI <- renderUI({
      tagList(
        h3("Point Generation"),
        selectInput("generatorPoint", "Select Generator Point G:",
                    choices = sapply(1:nrow(pts), function(i) paste0("(", pts$x[i], ",", pts$y[i], ")"))),
        actionButton("startGenerating", "Start/Reset Generating", class = "btn-success"),
        actionButton("nextPoint", "Next Point", class = "btn-info")
      )
    })
  })
  
  observeEvent(input$startGenerating, {
    req(input$generatorPoint)
    gStr <- input$generatorPoint
    gStr <- gsub("[()]", "", gStr)
    coords <- as.numeric(strsplit(gStr, ",")[[1]])
    G <- list(x = coords[1], y = coords[2])
    
    initialG(G)
    currentPoint(G)
    kValue(1)
    pathPoints(data.frame(x = G$x, y = G$y))
    pathLines(data.frame())
    
    # Update the plot
    output$ellipticCurvePlot <- renderPlot({
      plotPointsWithPath(points(), pathPoints(), pathLines(), initialG())
    })
    
    # Update current point info
    output$currentPointInfo <- renderText({
      sprintf("k = %d\nCurrent Point: (%d, %d)", kValue(), currentPoint()$x, currentPoint()$y)
    })
  })
  
  observeEvent(input$nextPoint, {
    req(currentPoint())
    a <- input$a
    p <- input$p
    G <- initialG()
    P <- currentPoint()
    
    H <- pointAdd(P, G, a, p)
    
    if (is.null(H)) {
      showNotification("Point at infinity reached.", type = "warning")
      return()
    }
    
    kValue(kValue() + 1)
    currentPoint(H)
    
    # Update path points and lines
    pathPoints(rbind(pathPoints(), data.frame(x = H$x, y = H$y)))
    pathLines(rbind(pathLines(), data.frame(
      x = P$x, y = P$y,
      xend = H$x, yend = H$y
    )))
    
    # Update the plot
    output$ellipticCurvePlot <- renderPlot({
      plotPointsWithPath(points(), pathPoints(), pathLines(), initialG())
    })
    
    # Update current point info
    output$currentPointInfo <- renderText({
      sprintf("k = %d\nCurrent Point: (%d, %d)", kValue(), currentPoint()$x, currentPoint()$y)
    })
  })
  
  # Function to plot points using ggplot2
  plotPoints <- function(pts, a, b, p) {
    ggplot(pts, aes(x = x, y = y)) +
      geom_point(size = 3) +
      ggtitle(sprintf("Elliptic Curve over F_%d", p)) +
      xlab("x") + ylab("y") +
      theme_minimal()
  }
  
  # Function to plot points with highlighted path
  plotPointsWithPath <- function(pts, pathPts, pathLns, G) {
    p <- ggplot(pts, aes(x = x, y = y)) +
      geom_point(size = 3, color = "grey") +
      xlab("x") + ylab("y") +
      theme_minimal()
    
    if (!is.null(G)) {
      p <- p + annotate("point", x = G$x, y = G$y, color = "red", size = 4)
    }
    
    if (nrow(pathPts) > 0) {
      p <- p + geom_point(data = pathPts, aes(x = x, y = y), color = "blue", size = 4)
    }
    
    if (nrow(pathLns) > 0) {
      p <- p + geom_segment(data = pathLns, aes(x = x, y = y, xend = xend, yend = yend),
                            linetype = "dashed", color = "green")
    }
    
    p
  }
  
}

# Run the Shiny app
shinyApp(ui = ui, server = server)
