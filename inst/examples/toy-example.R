library(ABCmaps)

# Store results beside this example file by default.
source_files <- vapply(sys.frames(), function(frame) {
  if (is.null(frame$ofile)) "" else frame$ofile
}, character(1))
source_files <- source_files[basename(source_files) == "toy-example.R"]
file_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(source_files) > 0) {
  example_dir <- dirname(normalizePath(source_files[length(source_files)]))
} else if (length(file_argument) > 0) {
  example_dir <- dirname(normalizePath(sub("^--file=", "", file_argument[1])))
} else {
  example_dir <- getwd()
}

# Toy model and paper settings.
yobs <- 1
sigma <- 0.2
if (!exists("N", inherits = FALSE)) N <- 10000
if (!exists("pilot_N", inherits = FALSE)) pilot_N <- 5000
if (!exists("num_rep", inherits = FALSE)) num_rep <- 10
if (!exists("output_dir", inherits = FALSE)) output_dir <- example_dir

n_L <- 1
n_H <- 1
alpha <- 0.7
alpha_L <- 0.7
epsilon_T <- 0.05
gamma <- 0.5

f_L <- function(theta) 4 * theta^2
f_H <- function(theta) 4 * theta^2 + 0.3 * cos(5 * pi * theta)

Dis_fun <- function(theta, n_sim, type = "high") {
  theta <- as.numeric(theta)
  mu <- if (type == "low") f_L(theta) else f_H(theta)
  Y <- matrix(rnorm(length(theta) * n_sim, 0, sigma),
              length(theta), n_sim) + mu
  (Y - yobs)^2
}

Prior <- list(
  sample = function(n) runif(n, -2, 2),
  density = function(theta) dunif(theta, -2, 2)
)

# Estimate the LF threshold and pilot FPR.
set.seed(20260812)
pilot_time <- system.time({
  pilot <- estimate_pilot_fpr(
    Theta = Prior$sample(pilot_N),
    Dis_fun = Dis_fun,
    n_L = n_L,
    n_H = n_H,
    kappa = 0.05,
    epsilon_target = epsilon_T,
    a_L = 0.01,
    cost_mode = "fixed",
    c_L = 20,
    c_H = 500
  )
})[["elapsed"]]
epsilon_L_T <- pilot$summary$epsilon_L

# Run 10 paired ABC-ASMC and MAPS repetitions.
runs <- vector("list", 2 * num_rep)
for (i in seq_len(num_rep)) {
  seed <- 20260806 + 10000 + i

  set.seed(seed)
  asmc_time <- system.time({
    asmc_fit <- ABC_ASMC(
      Prior, Dis_fun,
      N = N, n_H = n_H, alpha = alpha,
      epsilon_T = epsilon_T, gamma = gamma
    )
  })[["elapsed"]]

  set.seed(seed)
  maps_time <- system.time({
    maps_fit <- MAPS(
      Prior, Dis_fun,
      N = N, n_L = n_L, n_H = n_H,
      alpha = alpha, alpha_L = alpha_L,
      epsilon_T = epsilon_T, epsilon_L_T = epsilon_L_T,
      gamma = gamma
    )
  })[["elapsed"]]

  runs[[2 * i - 1]] <- list(
    method = "ABC-ASMC",
    theta = asmc_fit$Theta[, 1],
    weights = asmc_fit$weights,
    time = asmc_time,
    n_high = asmc_fit$n_high[length(asmc_fit$n_high)]
  )
  runs[[2 * i]] <- list(
    method = "MAPS",
    theta = maps_fit$Theta[, 1],
    weights = maps_fit$weights,
    time = maps_time + pilot_time,
    n_high = maps_fit$n_high[length(maps_fit$n_high)]
  )
  message("Finished repetition ", i, "/", num_rep)
}

# Calculate ISE against the exact HF-ABC density.
grid <- seq(-2, 2, length.out = 1201)
radius <- sqrt(epsilon_T)
truth <- pnorm(yobs + radius, f_H(grid), sigma) -
  pnorm(yobs - radius, f_H(grid), sigma)
truth <- truth / sum(diff(grid) * (truth[-1] + truth[-length(truth)]) / 2)

comparison <- do.call(rbind, lapply(seq_along(runs), function(i) {
  fit <- runs[[i]]
  estimate <- density(fit$theta, weights = fit$weights, bw = 0.03,
                      from = -2, to = 2, n = length(grid))$y
  error <- (estimate - truth)^2
  data.frame(
    repetition = ceiling(i / 2),
    method = fit$method,
    time = fit$time,
    ISE = sum(diff(grid) * (error[-1] + error[-length(error)]) / 2),
    n_high = fit$n_high
  )
}))

summary <- aggregate(cbind(ISE, time, n_high) ~ method, comparison, mean)

# Save the numerical results and the requested comparison plot.
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(comparison, file.path(output_dir, "toy-comparison.csv"),
          row.names = FALSE)
write.csv(summary, file.path(output_dir, "toy-summary.csv"), row.names = FALSE)

draw_plot <- function() {
  colors <- c("ABC-ASMC" = "#0072B2", "MAPS" = "#D55E00")
  shapes <- c("ABC-ASMC" = 16, "MAPS" = 17)
  plot(comparison$time, comparison$ISE,
       xlab = "Wall-clock time (seconds)", ylab = "ISE",
       col = colors[comparison$method], pch = shapes[comparison$method],
       cex = 1.2)
  legend("topright", names(colors), col = colors, pch = shapes, bty = "n")
}

pdf(file.path(output_dir, "toy-comparison.pdf"), width = 7, height = 6)
draw_plot()
dev.off()

png(file.path(output_dir, "toy-comparison.png"),
    width = 1400, height = 1200, res = 180)
draw_plot()
dev.off()

print(pilot$summary[c("epsilon_L", "fpr", "fpr_max")])
print(summary)
