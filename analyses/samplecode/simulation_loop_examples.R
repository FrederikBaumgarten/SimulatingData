## ===========================================================================
##
##  SIMULATE BEFORE YOU SAMPLE
##  A three-step loop for designing ecological studies
##
##  Supplementary R script -- worked examples
##  Baumgarten, Amrhein & Wolkovich
##
## ---------------------------------------------------------------------------
##
##  HOW TO USE THIS SCRIPT
##
##  Run it top to bottom, or step through section by section. Every section is
##  self-contained after Section 0 and follows the same shape:
##
##      (a) a comment block explaining WHAT WE ARE ASKING
##      (b) the simulation code
##      (c) a plot
##      (d) a comment block explaining WHAT IT SHOWS
##
##  Base R only. No packages are required except for the optional hierarchical
##  example in Section 5, which uses lme4 and is skipped if it is not installed.
##
##  The three steps of the simulation loop:
##      Step 1  Translate a verbal hypothesis into an explicit generative model
##      Step 2  Simulate and visualise the data that model implies
##      Step 3  Iterate design choices against a precision target
##
## ===========================================================================


## ---------------------------------------------------------------------------
## SECTION 0. Setup
## ---------------------------------------------------------------------------

set.seed(1)

## Set to TRUE to write each figure to a PDF in the working directory.
SAVE_PLOTS <- FALSE

## Number of replicate studies per design. Lower this to 200 for a quick pass;
## raise it to 2000 for smoother curves in the final figures.
R_REPS <- 500

## Small helper so every figure is either drawn to screen or written to file.
open_plot <- function(name, width = 7, height = 5) {
  if (SAVE_PLOTS) pdf(paste0(name, ".pdf"), width = width, height = height)
}
close_plot <- function() {
  if (SAVE_PLOTS) dev.off()
}

## Consistent look for all figures.
par(bty = "l", las = 1)



## ===========================================================================
## SECTION 1. LINEAR EXAMPLE -- plant growth and nitrogen
## ===========================================================================
##
## Research question: what is the effect of nitrogen fertilisation on sunflower
## biomass? We plan to grow plants across a gradient of nitrogen concentrations
## and harvest biomass at the end of the season.
##
## ---------------------------------------------------------------------------


## ---------------------------------------------------------------------------
## 1.1  STEP 1 -- From a verbal hypothesis to a generative model
## ---------------------------------------------------------------------------
##
## The verbal hypothesis is "plants grow more when nitrogen increases". That
## sentence is not yet something we can simulate from. To simulate, we have to
## commit to three things the sentence leaves open:
##
##      - the SHAPE of the relationship (here: linear)
##      - the SIZE of the effect        (here: the slope b)
##      - the CONSISTENCY of the effect (here: the residual SD sigma)
##
## The model:
##
##      y_i = a + b * x_i + e_i,     e_i ~ Normal(0, sigma)
##
## Note what has already happened: a vague hypothesis has become three numbers
## we have to defend. This is the main work of Step 1, and it is usually where
## you discover the hypothesis was less well specified than it felt.

a     <- 30   # intercept: expected biomass (g) at x = 0
b     <-  7   # slope: expected biomass gain (g) per unit nitrogen  <- the estimand
sigma <-  5   # residual SD (g): process variation + measurement error

## These are not "truth". They are assumptions, made explicit so they can be
## argued with. Everything below is conditional on them, which is exactly the
## point: we are asking "IF the world looks roughly like this, what will my
## planned study be able to tell me?"


## ---------------------------------------------------------------------------
## 1.2  STEP 2 -- Simulate the data the model implies, and look at it
## ---------------------------------------------------------------------------
##
## We simulate y THROUGH the model structure rather than drawing y directly
## from some marginal distribution: generate the predictor under the planned
## design, push it through the mean function, then add noise.

n <- 50
x <- rnorm(n, mean = 10, sd = 4)              # planned nitrogen gradient
y <- a + b * x + rnorm(n, mean = 0, sd = sigma)
fake <- data.frame(x = x, y = y)

open_plot("fig_linear_fakedata")
plot(fake$x, fake$y,
     xlab = "Nitrogen concentration",
     ylab = "Biomass (g)",
     main = "Step 2: the hypothesis, made visible",
     pch = 16, col = "grey30")
abline(a, b, lwd = 2)
legend("topleft", legend = "assumed true relationship", lwd = 2, bty = "n")
close_plot()

## WHAT IT SHOWS
## This is the hypothesis in the same visual form as a real results figure --
## effect size and scatter together. The useful question at this point is not
## statistical but biological: does this look like data I could plausibly
## collect from this system? If the scatter looks far too tight, or the
## biomass values are impossible, the assumptions are wrong and it is far
## cheaper to find out now than after a field season.
##
## In an experiment with fixed treatment levels rather than a gradient, only
## the line generating x changes, e.g.
##      x <- rep(c(0, 5, 10, 15, 20), each = 10)
## The rest of the loop is identical.


## ---------------------------------------------------------------------------
## 1.3  Sanity check -- can the planned analysis recover what we put in?
## ---------------------------------------------------------------------------

fit_lm <- lm(y ~ x, data = fake)
summary(fit_lm)
confint(fit_lm)

est_b <- unname(coef(fit_lm)["x"])       # our estimate of the slope
ci_b  <- unname(confint(fit_lm)["x", ])  # its 95% interval: c(lower, upper)

open_plot("fig_linear_fit", width = 10, height = 5)
par(mfrow = c(1, 2), bty = "l", las = 1)

## LEFT: the data, with the fitted and the true line.
plot(fake$x, fake$y,
     xlab = "Nitrogen concentration", ylab = "Biomass (g)",
     main = "Fitted vs. generating model",
     pch = 16, col = "grey30")
abline(fit_lm, lwd = 2, col = "firebrick")
abline(a, b, lwd = 2, lty = 2)
legend("topleft", legend = c("fitted", "true"),
       lwd = 2, lty = c(1, 2), col = c("firebrick", "black"), bty = "n")

## RIGHT: the estimate and its interval -- the quantity the whole rest of this
## script is about. Everything from Section 1.4 onwards is just asking how the
## LENGTH of this red bar responds to design choices.
plot(1, type = "n",
     xlim = range(c(ci_b, b)) + c(-0.6, 0.6), ylim = c(0.5, 1.5),
     yaxt = "n", ylab = "",
     xlab = "Slope b (g biomass per unit nitrogen)",
     main = "The estimate and its 95% interval")
abline(v = b, lwd = 2, lty = 2)
arrows(ci_b[1], 1, ci_b[2], 1, code = 3, angle = 90, length = 0.08,
       lwd = 3, col = "firebrick")
points(est_b, 1, pch = 16, cex = 1.6, col = "firebrick")
text(b, 1.40, "true b = 7", pos = 4, cex = 0.85)
text(mean(ci_b), 0.75, sprintf("width = %.2f", ci_b[2] - ci_b[1]), cex = 0.85)
par(mfrow = c(1, 1))
close_plot()

## WHAT IT SHOWS
## The estimate should land near b = 7. This is a check on the code, not on
## the science: if the analysis cannot recover the parameters from data we
## generated ourselves under ideal conditions, something is wrong with the
## model, the code, or both -- and we would never have known from the empirical
## fit alone.
##
## It is tempting to look at the p-value here. We deliberately do not. The
## quantity we carry forward is the WIDTH of the interval for b -- the red bar
## in the right-hand panel. A narrow bar means the study pinned the nitrogen
## effect down; a wide bar means it is compatible with many different effect
## sizes, whatever the p-value says.
##
## Note that this is ONE study. Run it again with a different seed and the bar
## shifts and changes length. That variability is exactly what Section 1.4
## quantifies.


## ---------------------------------------------------------------------------
## 1.4  STEP 3a -- How does precision depend on sample size?
## ---------------------------------------------------------------------------
##
## Now we play the planned study through many times and ask what we typically
## get. One "run" = generate a dataset, fit the planned model, record the
## width of the 95% interval for the parameter we care about.

one_run <- function(n, a, b, sigma, mean_x = 10, sd_x = 4) {
  x  <- rnorm(n, mean = mean_x, sd = sd_x)
  y  <- a + b * x + rnorm(n, 0, sigma)
  ci <- confint(lm(y ~ x))["x", ]
  unname(ci[2] - ci[1])
}

set.seed(2)

## Sample sizes to try. Deliberately dense at the low end, where the curve is
## steepest and where most real ecological studies actually sit.
ns <- c(8, 10, 12, 15, 20, 25, 30, 40, 50, 65, 80, 100, 125, 160, 200)

width_by_n <- sapply(ns, function(nn) {
  w <- replicate(R_REPS, one_run(nn, a, b, sigma))
  c(mean = mean(w), q80 = unname(quantile(w, 0.8)))
})

open_plot("fig_precision_vs_n")
plot(ns, width_by_n["mean", ], type = "n", ylim = c(0, max(width_by_n)),
     xaxt = "n",
     xlab = "Sample size (n)",
     ylab = "95% interval width for slope b",
     main = "Design-precision curve")
abline(v = ns, col = "grey93")                       # guide line at every point
lines(ns, width_by_n["mean", ], type = "b", pch = 16)
lines(ns, width_by_n["q80",  ], type = "b", pch = 1, lty = 2)
axis(1, at = ns, labels = FALSE, tcl = -0.3)         # tick under every point
axis(1, at = c(10, 20, 30, 50, 80, 125, 200))        # labels only where they fit
legend("topright", legend = c("mean width", "80th percentile of width"),
       pch = c(16, 1), lty = c(1, 2), bty = "n")
close_plot()

## WHAT IT SHOWS
## Precision improves with n, but with sharply diminishing returns -- interval
## width falls roughly as 1/sqrt(n), so halving the interval costs four times
## the sample size. Planning on that curve is a different exercise from asking
## whether a test will clear a threshold.


## ---------------------------------------------------------------------------
## 1.4b  Why there are TWO lines on that plot
## ---------------------------------------------------------------------------
##
## "The precision of a design with n = 50" is not a single number. Every study
## you could run under that design gives a slightly different interval, so the
## design produces a whole DISTRIBUTION of possible widths. Let us look at it.

set.seed(21)
w50 <- replicate(2000, one_run(50, a, b, sigma))

open_plot("fig_width_distribution")
hist(w50, breaks = 40, col = "grey85", border = "white",
     xlab = "95% interval width for b, in a single study with n = 50",
     main = "One design, many possible outcomes")
abline(v = mean(w50),           lwd = 2, col = "firebrick")
abline(v = quantile(w50, 0.8),  lwd = 2, lty = 2, col = "steelblue")
legend("topright",
       legend = c(sprintf("mean = %.2f", mean(w50)),
                  sprintf("80th percentile = %.2f", quantile(w50, 0.8))),
       lwd = 2, lty = c(1, 2), col = c("firebrick", "steelblue"), bty = "n")
close_plot()

## WHAT IT SHOWS
## This is the point the two lines in the previous figure were making.
##
## The MEAN width is a poor thing to design against. By construction, about
## half of the studies you might actually run come out WIDER than the mean. You
## only get to run the study once -- so planning to the mean gives you roughly
## a coin flip's chance of missing the precision you promised in your proposal.
##
## The 80th PERCENTILE answers a more useful question: "what width am I
## reasonably confident of achieving?" Design to it and 80% of the studies you
## could run would meet or beat your target. In the clinical trials literature
## this idea is called ASSURANCE. The 80% is a convention, not a law -- pick the
## percentile that matches how badly you can afford to be disappointed.
##
## The gap between the two matters most at small n, where the distribution of
## widths is most skewed. That is precisely where ecological studies live.


## ---------------------------------------------------------------------------
## 1.5  Validating the simulation machinery against the analytic result
## ---------------------------------------------------------------------------
##
## For an ordinary linear model there is a closed form, so we can check that
## the simulation reproduces it. This matters: it is the only case in this
## script where we can verify the machinery against a known answer, and having
## done so we can trust the same machinery in the cases below where no formula
## exists.
##
## The standard error of a regression slope, and the interval width it implies:
##
##      se(b)  =  sigma / (sd_x * sqrt(n - 1))
##      width  =  2 * t(0.975, n - 2) * se(b)
##
## where -- as a reminder, since these symbols were set a while ago:
##
##      sigma   residual SD of y around the line (we assumed 5 g)
##      sd_x    spread of the predictor values the design covers (we used 4)
##      n       number of observations
##      t(...)  the 97.5th percentile of a t distribution on n - 2 degrees of
##              freedom: the multiplier that turns a standard error into the
##              half-width of a 95% interval (roughly 1.96 once n is large)
##
## Read the formula as a statement about design: precision improves if the
## noise shrinks, if the predictor is spread more widely, or if the sample
## grows -- and only as sqrt(n) in that last case. Those three levers are
## exactly what Sections 1.4, 1.6 and 1.7 explore one at a time.

analytic_width <- function(n, sigma, sd_x = 4) {
  2 * qt(0.975, df = n - 2) * sigma / (sd_x * sqrt(n - 1))
}

cbind(n         = ns,
      simulated = round(width_by_n["mean", ], 3),
      analytic  = round(analytic_width(ns, sigma), 3))

open_plot("fig_analytic_check")
plot(ns, width_by_n["mean", ], type = "p", pch = 16,
     xlab = "Sample size (n)", ylab = "95% interval width for slope b",
     main = "Simulation reproduces the closed-form result")
curve(analytic_width(x, sigma), from = min(ns), to = max(ns), add = TRUE, lwd = 2)
legend("topright", legend = c("simulated", "analytic"),
       pch = c(16, NA), lwd = c(NA, 2), bty = "n")
close_plot()

## WHAT IT SHOWS
## The points sit on the curve. In this simple case simulation was not
## necessary -- and it is worth being honest that a formula would have done.
## The reason to simulate is that the formula runs out almost immediately:
## it does not survive non-linear mean functions (Section 2), structurally
## wrong models (Section 3), confounding (Section 4), or nested designs
## (Section 5). The simulation approach handles all of them unchanged.


## ---------------------------------------------------------------------------
## 1.6  STEP 3b -- Precision also depends on WHERE you sample
## ---------------------------------------------------------------------------
##
## Reminder of what one_run() does, because we are now turning a different dial.
## It simulates one complete study from
##
##      x ~ Normal(mean_x, sd_x)          the nitrogen gradient the design covers
##      y  = a + b*x + e,  e ~ Normal(0, sigma)
##
## fits lm(y ~ x), and returns the WIDTH of the 95% interval for b. Its
## arguments in order are (n, a, b, sigma, mean_x, sd_x).
##
## Below we hold everything fixed at n = 50, a = 30, b = 7, sigma = 5 and vary
## ONLY sd_x -- how widely the nitrogen values are spread. Note that sd_x is a
## DESIGN choice (which plots you fertilise at which rates), not a fact about
## nature, so it is genuinely yours to change.

set.seed(3)
sd_x_values <- c(1, 2, 4, 6, 8)

width_by_range <- sapply(sd_x_values, function(s) {
  mean(replicate(R_REPS, one_run(50, a, b, sigma, mean_x = 10, sd_x = s)))
})

open_plot("fig_precision_vs_range")
plot(sd_x_values, width_by_range, type = "b", pch = 16,
     xlab = "SD of nitrogen values (spread of the predictor)",
     ylab = "95% interval width for slope b",
     main = "Same n = 50, very different precision")
close_plot()

## WHAT IT SHOWS
## At a FIXED sample size, spreading the predictor out shrinks the interval
## dramatically. A study with a narrow nitrogen gradient and n = 50 can be far
## less informative than one with a wide gradient and n = 20. This lever is
## free -- it costs nothing extra to sample across a wider range -- and it is
## invisible to a conventional power calculation that treats n as the only
## design variable.


## ---------------------------------------------------------------------------
## 1.7  STEP 3c -- And on how noisy the measurements are
## ---------------------------------------------------------------------------
##
## Third dial. Same model as before,
##
##      y = a + b*x + e,     e ~ Normal(0, sigma)
##
## but now we vary sigma -- the SD of the scatter around the line. Two things
## live in sigma and it is worth keeping them apart in your head:
##   - PROCESS variation: plants genuinely differ, soils are patchy
##   - MEASUREMENT error: your balance, your proxy, your observer
## You cannot do much about the first, but the second is a design choice, which
## is what makes this dial actionable.

set.seed(4)
sigma_hi <- 10   # e.g. a cruder biomass proxy, or more heterogeneous plots

width_hi <- sapply(ns, function(nn) mean(replicate(R_REPS, one_run(nn, a, b, sigma_hi))))

open_plot("fig_precision_vs_noise")
plot(ns, width_by_n["mean", ], type = "b", pch = 16, ylim = c(0, max(width_hi)),
     xlab = "Sample size (n)", ylab = "95% interval width for slope b",
     main = "Measurement quality trades off against sample size")
lines(ns, width_hi, type = "b", pch = 1, lty = 2)
legend("topright", legend = c(paste("sigma =", sigma), paste("sigma =", sigma_hi)),
       pch = c(16, 1), lty = c(1, 2), bty = "n")
close_plot()

## WHAT IT SHOWS
## Doubling the noise costs about a fourfold increase in sample size to recover
## the same precision. Framed as a design decision: if a better instrument, a
## repeated measurement, or a more careful protocol can halve sigma, that may
## be far cheaper than quadrupling the number of plants.


## ---------------------------------------------------------------------------
## 1.8  Turning a scientific requirement into a sample size
## ---------------------------------------------------------------------------
##
## This is the step that makes the whole exercise decision-relevant. Instead of
## "how many plants do I need for significance?", ask:
##
##      How precisely do I need to know b for the answer to be USEFUL?
##
## Suppose a fertiliser recommendation requires the nitrogen effect to be known
## to within +/- 0.25 g biomass per unit nitrogen. That is a target interval
## width of 0.5. We can read the required n straight off the design-precision
## curve -- using the 80th percentile, not the mean.
##
## The general rule, since it is easy to lose in the arithmetic:
##
##      if you need the effect to within +/- k,  target an interval width of 2k
##
## because a 95% interval runs k below the estimate and k above it. Here k =
## 0.25 and so the target width is 0.5. Everything else stays as before: model
## y = a + b*x + e with a = 30, b = 7, sigma = 5, and x spread with sd_x = 4.

target_width <- 0.5

set.seed(5)
n_grid <- seq(20, 300, by = 20)
q80_by_n <- sapply(n_grid, function(nn) {
  quantile(replicate(R_REPS, one_run(nn, a, b, sigma)), 0.8)
})

n_required <- n_grid[which(q80_by_n <= target_width)[1]]
cat("Sample size needed to achieve width <=", target_width,
    "in 80% of studies: n =", n_required, "\n")

open_plot("fig_design_target")
plot(n_grid, q80_by_n, type = "b", pch = 16,
     xlab = "Sample size (n)", ylab = "80th percentile of 95% interval width",
     main = "Reading a design target off the curve")
abline(h = target_width, lty = 2, col = "firebrick")
abline(v = n_required, lty = 3)
text(max(n_grid), target_width, "required precision", pos = 3, adj = 1,
     col = "firebrick", cex = 0.8)
close_plot()

## WHAT IT SHOWS
## A defensible sample size, derived from what the science needs rather than
## from convention or from a power target. The number is conditional on the
## assumed sigma and predictor range -- which is a feature, not a limitation:
## it makes the basis of the design auditable, and reviewers can argue with the
## assumptions rather than with the number.


## ---------------------------------------------------------------------------
## 1.9  Why under-powered designs mislead even when they "work"
## ---------------------------------------------------------------------------
##
## Everything above concerns precision. This section concerns what happens to
## INTERPRETATION when precision is poor -- the Type S (sign) and Type M
## (magnitude) errors of Gelman & Carlin (2014).
##
## Consider a more realistic ecological situation than our flagship example: a
## small true effect, noisy measurements, and a modest sample.

b_small     <- 0.5
sigma_noisy <- 15
n_small     <- 20

set.seed(6)
sim <- replicate(4000, {
  x   <- rnorm(n_small, 10, 4)
  y   <- a + b_small * x + rnorm(n_small, 0, sigma_noisy)
  s   <- summary(lm(y ~ x))$coefficients["x", ]
  c(est = unname(s[1]), p = unname(s[4]))
})
est <- sim["est", ]
p   <- sim["p", ]

signif_runs <- p < 0.05
power       <- mean(signif_runs)
type_s      <- mean(sign(est[signif_runs]) != sign(b_small))
type_m      <- mean(abs(est[signif_runs])) / abs(b_small)

cat(sprintf("Power                                  : %.1f%%\n", 100 * power))
cat(sprintf("Type S (wrong sign | significant)      : %.1f%%\n", 100 * type_s))
cat(sprintf("Type M (exaggeration | significant)    : %.1fx\n", type_m))

open_plot("fig_type_sm")
hist(est, breaks = 50, col = "grey85", border = "white",
     xlab = "Estimated slope b", main = "What this design would report")
hist(est[signif_runs], breaks = 50, col = "firebrick", border = "white", add = TRUE)
abline(v = b_small, lwd = 2, lty = 2)
legend("topright",
       legend = c("all studies", "'significant' studies", "true b"),
       fill = c("grey85", "firebrick", NA), border = NA,
       lty = c(NA, NA, 2), lwd = c(NA, NA, 2), bty = "n")
close_plot()

## WHAT IT SHOWS
## The red distribution -- the studies that would get published -- barely
## overlaps the true value. Conditional on clearing p < 0.05, this design
## reports an effect several times too large, and a non-trivial fraction of the
## time with the wrong sign. The design does not merely fail to detect the
## effect; it actively manufactures misleading literature.
##
## This cannot be diagnosed after the fact from a single dataset. It is only
## visible before data collection, from the design itself.



## ===========================================================================
## SECTION 2. NON-LINEAR EXAMPLE -- a saturating response
## ===========================================================================
##
## Many ecological responses saturate. Here biomass rises with nitrogen and
## then levels off as other resources become limiting:
##
##      y_i = Vmax * x_i / (K + x_i) + e_i,     e_i ~ Normal(0, sigma)
##
## Vmax = asymptotic biomass; K = half-saturation constant (the nitrogen value
## at which expected biomass is Vmax / 2).
##
## In the linear case, design was mostly about n. Here it is mostly about WHERE
## the predictor is sampled -- and no closed-form power formula is available.
## ---------------------------------------------------------------------------

Vmax    <- 120
K       <- 8
sigma_n <- 10

mu_mm <- function(x, Vmax, K) Vmax * x / (K + x)


## ---------------------------------------------------------------------------
## 2.1  STEP 2 -- Two candidate designs, same sample size
## ---------------------------------------------------------------------------

set.seed(7)
n <- 50

x_narrow <- runif(n, 0, 12)    # Design A: low-to-moderate nitrogen only
x_wide   <- runif(n, 0, 40)    # Design B: extends past the bend

fake_narrow <- data.frame(x = x_narrow, y = mu_mm(x_narrow, Vmax, K) + rnorm(n, 0, sigma_n))
fake_wide   <- data.frame(x = x_wide,   y = mu_mm(x_wide,   Vmax, K) + rnorm(n, 0, sigma_n))

open_plot("fig_mm_designs", width = 10, height = 5)
par(mfrow = c(1, 2), bty = "l", las = 1)
plot(fake_narrow, xlim = c(0, 40), ylim = c(0, 140), pch = 16, col = "grey30",
     xlab = "Nitrogen", ylab = "Biomass (g)", main = "Design A: x in [0, 12]")
curve(mu_mm(x, Vmax, K), 0, 40, add = TRUE, lwd = 2)
plot(fake_wide,   xlim = c(0, 40), ylim = c(0, 140), pch = 16, col = "grey30",
     xlab = "Nitrogen", ylab = "Biomass (g)", main = "Design B: x in [0, 40]")
curve(mu_mm(x, Vmax, K), 0, 40, add = TRUE, lwd = 2)
par(mfrow = c(1, 1))
close_plot()

## WHAT IT SHOWS
## Over the narrow range the saturating curve is nearly a straight line. The
## data contain almost no information about the asymptote, because the design
## never went anywhere near it. Same n, same model, same noise -- but Design A
## cannot see the thing the model is about.


## ---------------------------------------------------------------------------
## 2.2  STEP 3 -- Fit both designs, and track convergence honestly
## ---------------------------------------------------------------------------
##
## Reminder of the model we are fitting, since it is the first non-linear one
## in this script:
##
##      y  =  Vmax * x / (K + x)  +  e,      e ~ Normal(0, sigma)
##
##      Vmax  the asymptote: biomass approached as nitrogen becomes unlimiting
##            (we set 120 g)
##      K     the half-saturation constant: the nitrogen value at which expected
##            biomass is Vmax/2, i.e. how FAST the curve bends (we set 8)
##      sigma residual SD as before (we set 10 g)
##
## Vmax and K are not interchangeable descriptions of "how much nitrogen
## helps": one is about the ceiling, the other about how quickly you approach
## it. The whole point of Section 2 is that a design can pin down one without
## the other -- or neither.
##
## Non-linear least squares can fail to converge. It is tempting to set
## warnOnly = TRUE and move on, but that returns non-converged fits as if they
## were estimates and quietly contaminates everything downstream. Instead we
## record convergence status and report the failure rate -- which is itself one
## of the most useful things a design simulation can tell you.

fit_mm <- function(dat) {
  fit <- try(nls(y ~ Vmax * x / (K + x),
                 data    = dat,
                 start   = list(Vmax = max(dat$y), K = median(dat$x)),
                 control = nls.control(maxiter = 200)),
             silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  fit
}

one_run_mm <- function(n, x_max, Vmax, K, sigma) {
  x   <- runif(n, 0, x_max)
  dat <- data.frame(x = x, y = mu_mm(x, Vmax, K) + rnorm(n, 0, sigma))
  fit <- fit_mm(dat)
  if (is.null(fit)) return(c(Vmax_hat = NA, K_hat = NA))
  co <- coef(fit)
  c(Vmax_hat = unname(co["Vmax"]), K_hat = unname(co["K"]))
}

set.seed(8)
estA <- as.data.frame(t(replicate(R_REPS, one_run_mm(50, 12, Vmax, K, sigma_n))))
estB <- as.data.frame(t(replicate(R_REPS, one_run_mm(50, 40, Vmax, K, sigma_n))))

cat(sprintf("Design A (narrow): %.1f%% of runs failed to converge\n",
            100 * mean(is.na(estA$K_hat))))
cat(sprintf("Design B (wide)  : %.1f%% of runs failed to converge\n",
            100 * mean(is.na(estB$K_hat))))

estA <- estA[complete.cases(estA), ]
estB <- estB[complete.cases(estB), ]

## NOTE: dropping failed fits is itself a selection step, of exactly the kind
## this paper warns about elsewhere. The surviving estimates from Design A are
## a biased subset -- the runs that happened to be informative. Always report
## how many were dropped alongside the results.


## ---------------------------------------------------------------------------
## 2.3  The cloud of estimates -- seeing a parameter trade-off
## ---------------------------------------------------------------------------

open_plot("fig_mm_cloud", width = 10, height = 5)
par(mfrow = c(1, 2), bty = "l", las = 1)
lim_K <- c(0, 40); lim_V <- c(0, 400)
plot(estA$K_hat, estA$Vmax_hat, xlim = lim_K, ylim = lim_V,
     pch = 16, col = rgb(0, 0, 0, 0.25),
     xlab = "K estimate", ylab = "Vmax estimate",
     main = "Design A: strong trade-off")
abline(v = K, h = Vmax, lty = 2, col = "firebrick")
plot(estB$K_hat, estB$Vmax_hat, xlim = lim_K, ylim = lim_V,
     pch = 16, col = rgb(0, 0, 0, 0.25),
     xlab = "K estimate", ylab = "Vmax estimate",
     main = "Design B: identified")
abline(v = K, h = Vmax, lty = 2, col = "firebrick")
par(mfrow = c(1, 1))
close_plot()

## WHAT IT SHOWS
## Under Design A the estimates smear along a ridge: many (K, Vmax) pairs
## produce almost identical curves over x in [0, 12], so the data cannot tell
## them apart. Under Design B the cloud collapses around the true values.
##
## This is a design problem, not an analysis problem. No cleverer statistical
## method recovers information the sampling design never collected.


## ---------------------------------------------------------------------------
## 2.4  More data, or better-placed data?
## ---------------------------------------------------------------------------

set.seed(9)
ns_mm <- c(20, 50, 100, 200)

spread_K <- function(n, x_max) {
  est <- t(replicate(R_REPS, one_run_mm(n, x_max, Vmax, K, sigma_n)))
  est <- est[complete.cases(est), , drop = FALSE]
  IQR(est[, "K_hat"])
}

iqr_narrow <- sapply(ns_mm, spread_K, x_max = 12)
iqr_wide   <- sapply(ns_mm, spread_K, x_max = 40)

open_plot("fig_mm_n_vs_range")
plot(ns_mm, iqr_narrow, type = "b", pch = 16, ylim = c(0, max(iqr_narrow)),
     xlab = "Sample size (n)", ylab = "IQR of K estimates",
     main = "Where you sample beats how much you sample")
lines(ns_mm, iqr_wide, type = "b", pch = 1, lty = 2)
legend("topright", legend = c("x in [0, 12]", "x in [0, 40]"),
       pch = c(16, 1), lty = c(1, 2), bty = "n")
close_plot()

## WHAT IT SHOWS
## Quadrupling n under the narrow design buys less than simply extending the
## nitrogen gradient at the original n. If you only have budget for one change,
## the simulation tells you which one to make -- before you spend it.



## ===========================================================================
## SECTION 3. MISSPECIFICATION -- a precise answer to the wrong question
## ===========================================================================
##
## So far the fitted model always matched the generating model. Real analyses
## are not so lucky. Here we fit the WRONG structural model -- a straight line
## -- to data generated by the saturating process, using the narrow design
## where the curve looks almost linear.
## ---------------------------------------------------------------------------

fit_wrong <- lm(y ~ x, data = fake_narrow)

cat(sprintf("R-squared of the wrong model : %.3f\n", summary(fit_wrong)$r.squared))
cat(sprintf("95%% CI for the slope         : [%.2f, %.2f]\n",
            confint(fit_wrong)["x", 1], confint(fit_wrong)["x", 2]))
cat(sprintf("Predicted biomass at x = 40  : %.0f g\n",
            predict(fit_wrong, newdata = data.frame(x = 40))))
cat(sprintf("True expected biomass at 40  : %.0f g\n", mu_mm(40, Vmax, K)))

open_plot("fig_misspecification")
plot(fake_narrow$x, fake_narrow$y, xlim = c(0, 45), ylim = c(0, 320),
     pch = 16, col = "grey30",
     xlab = "Nitrogen", ylab = "Biomass (g)",
     main = "A well-fitting model, extrapolated")
curve(mu_mm(x, Vmax, K), 0, 45, add = TRUE, lwd = 2)
abline(fit_wrong, lwd = 2, lty = 2, col = "firebrick")
rect(0, 0, 12, 320, col = rgb(0, 0, 0, 0.05), border = NA)
text(6, 300, "sampled", cex = 0.8)
text(30, 300, "extrapolated", cex = 0.8)
legend("bottomright", legend = c("true (saturating)", "fitted (linear)"),
       lwd = 2, lty = c(1, 2), col = c("black", "firebrick"), bty = "n")
close_plot()

## WHAT IT SHOWS
## Every conventional diagnostic looks fine. R-squared is respectable, the
## slope is estimated to within a couple of percent, nothing in the output
## flags a problem -- and the prediction at x = 40 is wrong by more than a
## factor of two.
##
## This is the single most important result in the script. PRECISION IS NOT
## CORRECTNESS. Narrowing an interval around a parameter of a wrong model just
## makes you more confident about the wrong thing, and collecting more data
## makes it worse rather than better. Designing to a precision target is only
## meaningful alongside a defensible model structure -- which is why Step 1 of
## the loop is where it is.



## ===========================================================================
## SECTION 4. CONFOUNDING -- precision is not validity either
## ===========================================================================
##
## An observational design. Soil moisture drives both nitrogen availability and
## biomass, so nitrogen and biomass are associated for reasons that have
## nothing to do with a causal effect of nitrogen. We have a large sample:
## n = 500.
##
## The generative model now has THREE variables rather than two, which is the
## whole point -- so, written out:
##
##      moisture ~ Normal(0, 1)                        the confounder
##      x        = 10 + 3*moisture + noise             N tracks moisture
##      y        = 30 + 7*x + 15*moisture + noise      both affect biomass
##
## The 7 is the true causal effect of nitrogen: what biomass would change by if
## we intervened on nitrogen and left moisture alone. The 15 is moisture's own
## effect. Because moisture pushes x and y in the same direction, a model that
## omits it will attribute part of moisture's effect to nitrogen.
## ---------------------------------------------------------------------------

set.seed(10)
n_obs    <- 500
moisture <- rnorm(n_obs, 0, 1)
x_obs    <- 10 + 3 * moisture + rnorm(n_obs, 0, 1)          # N tracks moisture
y_obs    <- 30 + 7 * x_obs + 15 * moisture + rnorm(n_obs, 0, 5)
obs      <- data.frame(x = x_obs, y = y_obs, moisture = moisture)

naive <- lm(y ~ x, data = obs)                # moisture unmeasured
adj   <- lm(y ~ x + moisture, data = obs)     # moisture measured and included

cat("True causal effect of nitrogen : 7.00\n")
cat(sprintf("Naive estimate    : %.2f  95%% CI [%.2f, %.2f]\n",
            coef(naive)["x"], confint(naive)["x", 1], confint(naive)["x", 2]))
cat(sprintf("Adjusted estimate : %.2f  95%% CI [%.2f, %.2f]\n",
            coef(adj)["x"],   confint(adj)["x", 1],   confint(adj)["x", 2]))

open_plot("fig_confounding")
plot(1, type = "n", xlim = c(6, 13), ylim = c(0.5, 2.5), yaxt = "n",
     xlab = "Estimated effect of nitrogen on biomass", ylab = "",
     main = "n = 500: narrow, confident, and wrong")
abline(v = 7, lwd = 2, lty = 2)
arrows(confint(naive)["x", 1], 2, confint(naive)["x", 2], 2,
       code = 3, angle = 90, length = 0.05, lwd = 3, col = "firebrick")
arrows(confint(adj)["x", 1], 1, confint(adj)["x", 2], 1,
       code = 3, angle = 90, length = 0.05, lwd = 3, col = "steelblue")
points(coef(naive)["x"], 2, pch = 16, cex = 1.5, col = "firebrick")
points(coef(adj)["x"],   1, pch = 16, cex = 1.5, col = "steelblue")
axis(2, at = c(2, 1), labels = c("naive", "adjusted"), las = 1)
text(7, 2.4, "true effect", pos = 4, cex = 0.8)
close_plot()

## WHAT IT SHOWS
## The naive estimate is about 65% too large, and its interval is so narrow
## that it excludes the truth comfortably. More data would shrink that interval
## further while leaving the bias untouched -- large n buys precision about a
## biased quantity, not accuracy.
##
## The simulation loop cannot fix this. What it can do is make it visible at
## the design stage, where you still have the option to measure the confounder,
## randomise, or state plainly that the study estimates an association rather
## than an effect. Simulating a confounder you suspect exists is a cheap way to
## find out how much it would matter if you ignored it.



## ===========================================================================
## SECTION 5. HIERARCHICAL DESIGN -- how many sites, how many plots?
## ===========================================================================
##
## The case where simulation stops being optional. Plots are nested within
## sites, sites differ in baseline biomass and in their response to nitrogen,
## and the design question is how to split a fixed budget of plots between
## "more sites" and "more plots per site". There is no useful closed-form
## answer.
##
## The generative model, spelled out -- this is the first one with more than one
## source of variation, so it is worth reading slowly:
##
##      for each site s:    a_s ~ Normal(30, sd_site)   the site's own intercept
##                          b_s ~ Normal(b,  sd_slope)  the site's own N effect
##      for each plot i in site s:
##                          y_i = a_s + b_s * x_i + e_i,  e_i ~ Normal(0, sigma)
##
##      b         the AVERAGE nitrogen effect across sites -- our estimand (7)
##      sd_slope  how much that effect genuinely differs between sites (2)
##      sd_site   how much baseline biomass differs between sites (10)
##      sigma     plot-to-plot noise within a site (5)
##
## The key quantity is sd_slope. If the nitrogen effect were identical at every
## site, plots would be interchangeable and only the total would matter. Because
## it varies, each site gives you what amounts to ONE noisy observation of the
## effect -- and the number of sites starts to dominate everything else.
##
## Requires lme4. Skipped automatically if not installed.
## ---------------------------------------------------------------------------

if (!requireNamespace("lme4", quietly = TRUE)) {
  message("Section 5 skipped: install.packages('lme4') to run the hierarchical example.")
} else {

  sim_hier <- function(n_site, n_plot, b = 7, sigma = 5,
                       sd_site = 10, sd_slope = 2) {
    site   <- rep(seq_len(n_site), each = n_plot)
    a_site <- rnorm(n_site, 30, sd_site)[site]     # site-level intercepts
    b_site <- rnorm(n_site, b,  sd_slope)[site]    # site-level slopes
    x      <- rnorm(n_site * n_plot, 10, 4)
    y      <- a_site + b_site * x + rnorm(n_site * n_plot, 0, sigma)
    data.frame(y = y, x = x, site = factor(site))
  }

  ci_width_hier <- function(n_site, n_plot) {
    d <- sim_hier(n_site, n_plot)
    m <- try(suppressMessages(suppressWarnings(
           lme4::lmer(y ~ x + (x | site), data = d,
                      control = lme4::lmerControl(calc.derivs = FALSE)))),
             silent = TRUE)
    if (inherits(m, "try-error")) return(NA_real_)
    ## NOTE: do NOT write sqrt(diag(vcov(m))) here. vcov() on an lmerMod
    ## returns a Matrix-package S4 object, and base diag() fails on it with
    ## "long vectors not supported yet". Take the standard error from the
    ## fixed-effects coefficient table instead, which is a plain matrix.
    se <- summary(m)$coefficients["x", "Std. Error"]
    unname(2 * 1.96 * se)
  }

  ## A fixed budget of 120 plots, split five different ways.
  designs <- data.frame(n_site = c( 4,  6, 12, 20, 30),
                        n_plot = c(30, 20, 10,  6,  4))
  designs$total <- designs$n_site * designs$n_plot   # 120 in every row

  set.seed(11)
  reps_h <- max(50, R_REPS %/% 5)   # mixed models are slower; use fewer reps
  res <- mapply(function(s, p) {
    w <- replicate(reps_h, ci_width_hier(s, p))
    c(width = mean(w, na.rm = TRUE), fail_rate = mean(is.na(w)))
  }, designs$n_site, designs$n_plot)

  designs$width     <- res["width", ]
  designs$fail_rate <- res["fail_rate", ]
  print(designs)

  open_plot("fig_hierarchical")
  bp <- barplot(designs$width,
                names.arg = paste0(designs$n_site, " sites\n",
                                   designs$n_plot, " plots"),
                ylab = "95% interval width for the average nitrogen effect",
                main = "Same 120 plots, five ways to spend them",
                col = "grey80", border = NA,
                ylim = c(0, max(designs$width) * 1.15))
  text(bp, designs$width, sprintf("%.2f", designs$width), pos = 3, cex = 0.85)
  close_plot()

  ## WHAT IT SHOWS
  ## The same total effort produces very different precision depending on how
  ## it is allocated -- roughly a two-fold difference in interval width between
  ## the best and the worst allocation, at identical cost.
  ##
  ## Why: when the nitrogen effect genuinely varies between sites, adding more
  ## plots INSIDE a site mostly buys a better estimate of that site's own slope,
  ## which you already knew reasonably well. It does nothing about the
  ## site-to-site variation -- and that is what limits your knowledge of the
  ## AVERAGE effect. The effective sample size for the average is closer to the
  ## number of sites than to the number of plots.
  ##
  ## Watch the fail_rate column too. With only 4 sites, a random-slope model is
  ## often singular or fails to converge: the design is asking the data for
  ## something they cannot supply. That failure rate is itself a design
  ## diagnostic, and worth reporting.
  ##
  ## No formula gives you this. Simulating the design does, in a few lines.
}



## ===========================================================================
## SECTION 6. What to report
## ===========================================================================
##
## If you use this loop to justify a design, report enough for someone else to
## disagree with you productively:
##
##   1. The generative model -- the equation for the mean function and for the
##      stochastic component, on the response scale.
##   2. The assumed parameter values, with their justification (literature,
##      pilot data, or explicitly labelled expert judgement).
##   3. What the error term represents: process variation, measurement error,
##      or both; and whether hierarchy or autocorrelation was included.
##   4. The designs explored and the performance criterion (here: the 80th
##      percentile of the 95% interval width for the focal parameter).
##   5. The resulting design target and how it maps to a scientific or
##      management requirement.
##   6. This script, a seed, and sessionInfo().
##
## ---------------------------------------------------------------------------

sessionInfo()

## ===========================================================================
## END
## ===========================================================================
