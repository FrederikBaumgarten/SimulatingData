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
##  Run it top to bottom, or step through section by section. Each section
##  follows the same shape:
##      (a) a comment block explaining WHAT WE ARE ASKING
##      (b) the simulation code
##      (c) a figure
##      (d) a comment block explaining WHAT IT SHOWS
##
##  Base R only, except the optional hierarchical example in Section 5, which
##  uses lme4 and is skipped if it is not installed.
##
##  The three steps of the simulation loop:
##      Step 1  Translate a verbal hypothesis into an explicit generative model
##      Step 2  Simulate and visualise the data that model implies
##      Step 3  Iterate design choices against a precision target
##
##
## ===========================================================================


## ---------------------------------------------------------------------------
## SECTION 0. Setup
## ---------------------------------------------------------------------------

set.seed(1)

SAVE_PLOTS <- FALSE                                   # TRUE to write PDFs
FIG_DIR    <- "~/github/SimulatingData/analyses/figures/"
R_REPS     <- 500                                     # 100 to draft, 2000 final

## Figures are built at the size they will be printed. \textwidth in the
## manuscript is 15.59 cm = 6.14 inches; include at width=\textwidth and do
## not scale, or the type drops below the journal minimum.
TEXTWIDTH_IN <- 6.14

## Shared panel style: closed box, inward ticks, tight margins.
panel_par <- function(mfrow = c(1, 1), mar = c(3.9, 3.9, 2.0, 0.9)) {
  par(mfrow = mfrow, bty = "o", las = 1, mar = mar,
      mgp = c(2.3, 0.6, 0), tcl = 0.3,          # tcl > 0 = ticks point inwards
      cex.main = 1.0, font.main = 1, cex.lab = 1.0, cex.axis = 0.9)
}

open_plot <- function(name, height = 4.2, width = TEXTWIDTH_IN,
                      mfrow = c(1, 1), mar = c(3.9, 3.9, 2.0, 0.9)) {
  if (SAVE_PLOTS) {
    pdf(file = file.path(FIG_DIR, paste0(name, ".pdf")),
        width = width, height = height, pointsize = 10, useDingbats = FALSE)
  }
  panel_par(mfrow, mar)
}
close_plot <- function() if (SAVE_PLOTS) dev.off()

## Draw the closing box after every panel (some plot types overwrite it).
fin <- function(label = NULL) {
  box(bty = "o")
  if (!is.null(label)) mtext(label, side = 3, adj = 0, line = 0.5, font = 2)
}

panel_par()



## ===========================================================================
## SECTION 1. LINEAR EXAMPLE -- plant growth and nitrogen
## ===========================================================================
##
## Research question: what is the effect of nitrogen fertilisation on sunflower
## biomass? We plan to grow plants across a nitrogen gradient and harvest
## biomass at the end of the season.
## ---------------------------------------------------------------------------


## ---------------------------------------------------------------------------
## 1.1  STEP 1 -- From a verbal hypothesis to a generative model
## ---------------------------------------------------------------------------
##
## "Plants grow more when nitrogen increases" cannot be simulated from. To
## simulate we must commit to three things it leaves open:
##      - the SHAPE of the relationship   (here: linear)
##      - the SIZE of the effect          (the slope b)
##      - the CONSISTENCY of the effect   (the residual SD sigma)
##
##      y_i = a + b * x_i + e_i,     e_i ~ Normal(0, sigma)

a     <- 30   # intercept: expected biomass (g) at zero nitrogen
b     <-  7   # slope: biomass gain (g) per unit nitrogen   <- the estimand
sigma <-  5   # residual SD (g): process variation + measurement error

## These are assumptions, not truth. Everything below is conditional on them.


## ---------------------------------------------------------------------------
## 1.2  STEP 2 -- Simulate the data the model implies, and look at it
## ---------------------------------------------------------------------------

n    <- 50
x    <- rnorm(n, mean = 10, sd = 4)                  # planned nitrogen gradient
y    <- a + b * x + rnorm(n, mean = 0, sd = sigma)
fake <- data.frame(x = x, y = y)

fit_lm <- lm(y ~ x, data = fake)
est_b  <- unname(coef(fit_lm)["x"])
ci_b   <- unname(confint(fit_lm)["x", ])

open_plot("fig_step2", height = 2.6, mfrow = c(1, 3))

## (a) the hypothesis, made visible
plot(fake$x, fake$y, pch = 16, col = "grey30",
     xlab = "Nitrogen", ylab = "Biomass (g)", main = "(a) Simulated data")
abline(a, b, lwd = 2)
fin()

## (b) does the planned analysis recover what we put in?
plot(fake$x, fake$y, pch = 16, col = "grey30",
     xlab = "Nitrogen", ylab = "Biomass (g)", main = "(b) Fitted and true line")
abline(fit_lm, lwd = 2, col = "firebrick")
abline(a, b, lwd = 2, lty = 2)
legend("topleft", legend = c("fitted", "true"), lwd = 2, lty = c(1, 2),
       col = c("firebrick", "black"), bty = "n", cex = 0.85, inset = 0.02)
fin()

## (c) the quantity the rest of the paper is about
plot(1, type = "n", xlim = range(c(ci_b, b)) + c(-0.7, 0.7), ylim = c(0.4, 1.6),
     yaxt = "n", ylab = "", xlab = "Slope b (g per unit N)",
     main = "(c) Estimate and interval")
abline(v = b, lwd = 2, lty = 2)
arrows(ci_b[1], 1, ci_b[2], 1, code = 3, angle = 90, length = 0.06,
       lwd = 3, col = "firebrick")
points(est_b, 1, pch = 16, cex = 1.5, col = "firebrick")
text(b, 1.45, "true b", pos = 4, cex = 0.85)
text(mean(ci_b), 0.62, sprintf("width = %.2f", diff(ci_b)), cex = 0.85)
fin()

close_plot()

## WHAT IT SHOWS
## (a) is the hypothesis in the same visual form as a real results figure:
## effect size and scatter together. The question to ask is biological, not
## statistical: could this system produce these data?
##
## (b) is a check on the code, not on the science. If a correct analysis
## cannot recover the parameters from data we generated ourselves under ideal
## conditions, something is wrong, and the empirical fit would never have told
## us. (A blunter title for this panel would be "recovery check"; we keep the
## descriptive one because the panel simply shows two lines.)
##
## (c) is the quantity we design against: the WIDTH of the interval for b.
## Run this again with a different seed and the bar moves and changes length.
## Section 1.4 quantifies that.


## ---------------------------------------------------------------------------
## 1.3  STEP 3 -- Three design levers, and one validity check
## ---------------------------------------------------------------------------
##
## One "run" = generate a dataset under the planned design, fit the planned
## model, record the width of the 95% interval for b.

one_run <- function(n, a, b, sigma, mean_x = 10, sd_x = 4) {
  x  <- rnorm(n, mean = mean_x, sd = sd_x)
  y  <- a + b * x + rnorm(n, 0, sigma)
  ci <- confint(lm(y ~ x))["x", ]
  unname(ci[2] - ci[1])
}

set.seed(2)
ns <- c(8, 10, 12, 15, 20, 25, 30, 40, 50, 65, 80, 100, 125, 160, 200)

width_by_n <- sapply(ns, function(nn) {
  w <- replicate(R_REPS, one_run(nn, a, b, sigma))
  c(mean = mean(w), q80 = unname(quantile(w, 0.8)))
})

w50 <- replicate(2000, one_run(50, a, b, sigma))     # one design, many outcomes

analytic_width <- function(n, sigma, sd_x = 4)       # closed form, for checking
  2 * qt(0.975, df = n - 2) * sigma / (sd_x * sqrt(n - 1))

sd_x_values    <- c(1, 2, 4, 6, 8)                   # lever 2: predictor spread
width_by_range <- sapply(sd_x_values, function(s)
  mean(replicate(R_REPS, one_run(50, a, b, sigma, sd_x = s))))

sigma_hi <- 10                                       # lever 3: measurement noise
width_hi <- sapply(ns, function(nn)
  mean(replicate(R_REPS, one_run(nn, a, b, sigma_hi))))

axis_ns <- c(10, 20, 30, 50, 80, 125, 200)           # labelled ticks

## The closed-form check is reported at the console rather than plotted: it
## validates the machinery but adds little to the manuscript figure.
cat("\nsimulated vs closed-form interval width\n")
print(cbind(n = ns,
            simulated = round(width_by_n["mean", ], 3),
            formula   = round(analytic_width(ns, sigma), 3)))

open_plot("fig_precision", height = 4.6, mfrow = c(2, 2))

## (a) design-precision curve
plot(ns, width_by_n["mean", ], type = "n", ylim = c(0, max(width_by_n)),
     xaxt = "n", xlab = "Sample size (n)", ylab = "95% interval width",
     main = "(a) Precision vs. sample size")
lines(ns, width_by_n["mean", ], type = "b", pch = 16, cex = 0.8)
lines(ns, width_by_n["q80",  ], type = "b", pch = 1,  cex = 0.8, lty = 2)
axis(1, at = ns, labels = FALSE, tcl = 0.2)
axis(1, at = axis_ns)
legend("topright", legend = c("mean", "80th pct"), pch = c(16, 1),
       lty = c(1, 2), bty = "n", cex = 0.85, inset = 0.02)
fin()

## (b) the width is itself a random variable
hist(w50, breaks = 40, col = "grey85", border = "white",
     xlab = "95% interval width", main = "(b) Width varies between studies")
abline(v = mean(w50),          lwd = 2, col = "firebrick")
abline(v = quantile(w50, 0.8), lwd = 2, lty = 2, col = "steelblue")
mtext("n = 50", side = 3, line = -1.3, adj = 0.95, cex = 0.75)
legend("topright", legend = c("mean", "80th pct"), lwd = 2, lty = c(1, 2),
       col = c("firebrick", "steelblue"), bg = "white", box.col = "grey80",
       cex = 0.8, inset = c(0.02, 0.12))
fin()

## (c) lever 2: where you sample
plot(sd_x_values, width_by_range, type = "b", pch = 16, cex = 0.8,
     xlab = "SD of nitrogen values", ylab = "95% interval width",
     main = "(c) Precision vs. spread of x")
mtext("n = 50", side = 3, line = -1.3, adj = 0.95, cex = 0.75)
fin()

## (d) lever 3: how well you measure
plot(ns, width_by_n["mean", ], type = "n", ylim = c(0, max(width_hi)),
     xaxt = "n", xlab = "Sample size (n)", ylab = "95% interval width",
     main = "(d) Precision vs. noise")
lines(ns, width_by_n["mean", ], type = "b", pch = 16, cex = 0.8)
lines(ns, width_hi,             type = "b", pch = 1,  cex = 0.8, lty = 2)
axis(1, at = ns, labels = FALSE, tcl = 0.2)
axis(1, at = c(10, 20, 50, 100, 200))
legend("topright", legend = c(bquote(sigma == .(sigma)),
                              bquote(sigma == .(sigma_hi))),
       pch = c(16, 1), lty = c(1, 2), bty = "n", cex = 0.85, inset = 0.02)
fin()

close_plot()

## WHAT IT SHOWS
## (a) Precision improves with n, but as 1/sqrt(n): halving the interval costs
## four times the sample. Two lines, because of (b).
##
## (b) A design does not have a precision; it has a DISTRIBUTION of precisions.
## About half of the studies you might run come out wider than the mean, and
## you only run the study once. Design against a high quantile instead. In the
## clinical trials literature this is called ASSURANCE.
##
## (c) At FIXED n, spreading the predictor out shrinks the interval sharply.
## This lever costs nothing and is invisible to a power calculation.
##
## For this model there is also a formula, se(b) = sigma / (sd_x * sqrt(n-1)),
## and the simulation reproduces it (printed above, not plotted). Simulation
## was not needed here; we check the machinery where it CAN be checked so we
## can trust it in Sections 2 to 4, where no formula exists.
##
## (d) Doubling the noise costs about four times the sample size. Halving sigma
## with a better instrument may be far cheaper than quadrupling the plants.


## ---------------------------------------------------------------------------
## 1.4  Turning a scientific requirement into a sample size
## ---------------------------------------------------------------------------
##
## If a claim requires the effect to within +/- k, target an interval of width
## 2k. Suppose a fertiliser recommendation needs b to within +/- 0.25 g per
## unit nitrogen: target width 0.5, read off at the 80th percentile.

target_width <- 0.5
set.seed(5)
n_grid   <- seq(20, 300, by = 20)
q80_by_n <- sapply(n_grid, function(nn)
  quantile(replicate(R_REPS, one_run(nn, a, b, sigma)), 0.8))
n_required <- n_grid[which(q80_by_n <= target_width)[1]]
cat(sprintf("n needed for width <= %.2f in 80%% of studies: %d\n",
            target_width, n_required))

open_plot("fig_design_target", height = 3.4, width = 4.2)
plot(n_grid, q80_by_n, type = "b", pch = 16, cex = 0.8,
     xlab = "Sample size (n)", ylab = "Interval width (80th pct)",
     main = "Choosing n from a precision target")
abline(h = target_width, lty = 2, col = "firebrick")
abline(v = n_required,   lty = 3)
text(mean(range(n_grid)), target_width, "required precision",
     pos = 3, col = "firebrick", cex = 0.85)
text(n_required, max(q80_by_n) * 0.9, paste("n =", n_required),
     pos = 4, cex = 0.85)
fin()
close_plot()

## WHAT IT SHOWS
## A defensible sample size, derived from what the science needs. It is
## conditional on the assumed sigma and predictor range, which is a feature: a
## reviewer can argue with sigma, a statement about the system, rather than
## with n, which is merely its consequence.


## ---------------------------------------------------------------------------
## 1.5  What an imprecise design does to the published literature
## ---------------------------------------------------------------------------
##
## A more typical ecological situation: small true effect, noisy measurement,
## modest sample.

b_small <- 0.5; sigma_noisy <- 15; n_small <- 20

set.seed(6)
sim <- replicate(4000, {
  x <- rnorm(n_small, 10, 4)
  y <- a + b_small * x + rnorm(n_small, 0, sigma_noisy)
  s <- summary(lm(y ~ x))$coefficients["x", ]
  c(est = unname(s[1]), p = unname(s[4]))
})
sig    <- sim["p", ] < 0.05
power  <- mean(sig)
type_s <- mean(sign(sim["est", sig]) != sign(b_small))
type_m <- mean(abs(sim["est", sig])) / abs(b_small)

cat(sprintf("Power %.0f%%   |   wrong sign %.0f%%   |   exaggeration %.1fx\n",
            100 * power, 100 * type_s, type_m))

open_plot("fig_type_sm", height = 3.4, width = 4.6)
hist(sim["est", ], breaks = 50, col = "grey85", border = "white",
     xlab = "Estimated slope b", main = "What this design would publish")
hist(sim["est", sig], breaks = 50, col = "firebrick", border = "white",
     add = TRUE)
abline(v = b_small, lwd = 2, lty = 2)
legend("topleft", legend = c("all studies", "reach p < 0.05", "true b"),
       fill = c("grey85", "firebrick", NA), border = NA,
       lty = c(NA, NA, 2), lwd = c(NA, NA, 2), bg = "white",
       box.col = "grey80", cex = 0.8, inset = 0.02)
fin()
close_plot()

## WHAT IT SHOWS
## The dark distribution is what would be published. It barely overlaps the
## truth. Conditional on p < 0.05 this design reports an effect several times
## too large and sometimes of the wrong sign. The design does not merely fail
## to detect the effect; it manufactures misleading results. None of this is
## visible from inside a single study, but all of it is visible beforehand.



## ===========================================================================
## SECTION 2. NON-LINEAR EXAMPLE -- a saturating response
## ===========================================================================
##
## Biomass rises with nitrogen, then levels off as other resources limit:
##
##      y_i = Vmax * x_i / (K + x_i) + e_i,     e_i ~ Normal(0, sigma)
##
##  Vmax  the CEILING: biomass approached as nitrogen becomes unlimiting.
##  K     the HALF-SATURATION CONSTANT: the nitrogen value at which biomass
##        reaches Vmax/2. K controls HOW FAST the curve rises, not how high.
##        Small K = rises steeply and flattens early. Large K = rises slowly
##        and needs much more nitrogen to approach the same ceiling.
##
## Vmax and K answer different biological questions, and a design can pin down
## one, the other, both, or neither. That is what this section is about.
## ---------------------------------------------------------------------------

Vmax    <- 120
K       <-   8
sigma_n <-  10

mu_mm <- function(x, Vmax, K) Vmax * x / (K + x)


## ---------------------------------------------------------------------------
## 2.1  STEP 2 -- Two designs, same sample size, different nitrogen range
## ---------------------------------------------------------------------------

set.seed(7)
n <- 50
x_narrow <- runif(n, 0, 12)   # Design A: stops short of saturation
x_wide   <- runif(n, 0, 40)   # Design B: extends past the bend

fake_narrow <- data.frame(x = x_narrow,
                          y = mu_mm(x_narrow, Vmax, K) + rnorm(n, 0, sigma_n))
fake_wide   <- data.frame(x = x_wide,
                          y = mu_mm(x_wide,   Vmax, K) + rnorm(n, 0, sigma_n))

open_plot("fig_mm_designs", height = 3.0, mfrow = c(1, 2))
plot(fake_narrow, xlim = c(0, 40), ylim = c(0, 140), pch = 16, col = "grey30",
     xlab = "Nitrogen", ylab = "Biomass (g)", main = "(a) Narrow range")
curve(mu_mm(x, Vmax, K), 0, 40, add = TRUE, lwd = 2)
fin()
plot(fake_wide,   xlim = c(0, 40), ylim = c(0, 140), pch = 16, col = "grey30",
     xlab = "Nitrogen", ylab = "Biomass (g)", main = "(b) Wide range")
curve(mu_mm(x, Vmax, K), 0, 40, add = TRUE, lwd = 2)
fin()
close_plot()

## WHAT IT SHOWS
## Over the narrow range the saturating curve is nearly a straight line. The
## data carry almost no information about the ceiling, because the design
## never went near it. Same n, same model, same noise.


## ---------------------------------------------------------------------------
## 2.2  STEP 3 -- Fit, track convergence, and see the parameter trade-off
## ---------------------------------------------------------------------------
##
## nls can fail to converge. Setting warnOnly = TRUE returns non-converged
## fits as if they were estimates, which quietly contaminates everything
## downstream. Record convergence instead: the failure rate is a design
## diagnostic in its own right.

fit_mm <- function(dat) {
  fit <- try(nls(y ~ Vmax * x / (K + x), data = dat,
                 start = list(Vmax = max(dat$y), K = median(dat$x)),
                 control = nls.control(maxiter = 200)), silent = TRUE)
  if (inherits(fit, "try-error")) NULL else fit
}

one_run_mm <- function(n, x_max) {
  x   <- runif(n, 0, x_max)
  dat <- data.frame(x = x, y = mu_mm(x, Vmax, K) + rnorm(n, 0, sigma_n))
  fit <- fit_mm(dat)
  if (is.null(fit)) return(c(Vmax_hat = NA, K_hat = NA))
  co <- coef(fit)
  c(Vmax_hat = unname(co["Vmax"]), K_hat = unname(co["K"]))
}

set.seed(8)
estA <- as.data.frame(t(replicate(R_REPS, one_run_mm(50, 12))))   # narrow
estB <- as.data.frame(t(replicate(R_REPS, one_run_mm(50, 40))))   # wide
failA <- mean(is.na(estA$K_hat)); failB <- mean(is.na(estB$K_hat))
cat(sprintf("convergence failures: narrow %.0f%%, wide %.0f%%\n",
            100 * failA, 100 * failB))
estA <- estA[complete.cases(estA), ]
estB <- estB[complete.cases(estB), ]

## NOTE: dropping failed fits is itself a selection step, of the kind this
## paper warns against elsewhere. Always report how many were dropped.

open_plot("fig_mm_cloud", height = 2.9, mfrow = c(1, 3))

lim_K <- c(0, 40); lim_V <- c(0, 400)

## (a)+(b) where the estimates land
plot(estA$K_hat, estA$Vmax_hat, xlim = lim_K, ylim = lim_V, pch = 16,
     col = rgb(0, 0, 0, 0.25), xlab = "K estimate", ylab = "Vmax estimate",
     main = "(a) Narrow: K, Vmax trade off")
abline(v = K, h = Vmax, lty = 2, col = "firebrick"); fin()

plot(estB$K_hat, estB$Vmax_hat, xlim = lim_K, ylim = lim_V, pch = 16,
     col = rgb(0, 0, 0, 0.25), xlab = "K estimate", ylab = "Vmax estimate",
     main = "(b) Wide: both identified")
abline(v = K, h = Vmax, lty = 2, col = "firebrick"); fin()

## (c) WHY they trade off: those different (K, Vmax) pairs draw the same curve
##     over the sampled range, and only separate beyond it.
plot(0, type = "n", xlim = c(0, 40), ylim = c(0, 220),
     xlab = "Nitrogen", ylab = "Biomass (g)",
     main = "(c) Curves coincide where sampled")
rect(0, 0, 12, 220, col = "grey92", border = NA)
show <- unique(round(seq(1, nrow(estA), length.out = 25)))   # a spread of fits
for (i in show)
  curve(mu_mm(x, estA$Vmax_hat[i], estA$K_hat[i]), 0, 40, add = TRUE,
        col = rgb(0, 0, 0, 0.25))
curve(mu_mm(x, Vmax, K), 0, 40, add = TRUE, lwd = 2.5, col = "firebrick")
text(6, 205, "sampled", cex = 0.8)
text(28, 205, "extrapolated", cex = 0.8)
fin()

close_plot()

## WHAT IT SHOWS
## (a) Under the narrow design the estimates smear along a ridge: a large K
## paired with a large Vmax fits the sampled data just as well as a small K
## with a small Vmax. Neither parameter is pinned down.
##
## (b) Under the wide design the cloud collapses onto the true values.
##
## (c) is the explanation. Each grey line is the curve implied by one narrow-
## design fit. Inside the shaded region, where the data are, they are
## indistinguishable. Outside it they fan out enormously. The data cannot
## choose between them because the design never visited the conditions that
## separate them.
##
## This is a design problem, not an analysis problem. The wide design does not
## analyse better; it simply visits the conditions where the ceiling shows.


## ---------------------------------------------------------------------------
## 2.3  Where to spend the next observation
## ---------------------------------------------------------------------------
##
## Section 1.3(d) already showed that spreading x out improves precision for a
## LINEAR slope. The new question here is a budget comparison in a case where
## the parameter is barely identifiable at all: given a fixed amount of money,
## is it better to buy more plants, or to buy a wider nitrogen gradient?

set.seed(9)
ns_mm <- c(10, 20, 35, 50, 75, 100, 150, 200)

spread_K <- function(n, x_max) {
  est <- t(replicate(R_REPS, one_run_mm(n, x_max)))
  IQR(est[complete.cases(est), "K_hat"])
}
iqr_narrow <- sapply(ns_mm, spread_K, x_max = 12)
iqr_wide   <- sapply(ns_mm, spread_K, x_max = 40)

open_plot("fig_mm_tradeoff", height = 3.4, width = 4.6)
plot(ns_mm, iqr_narrow, type = "b", pch = 16, cex = 0.8, log = "x",
     ylim = c(0, max(iqr_narrow)), xaxt = "n", xlab = "Sample size (n)",
     ylab = "Spread of K estimates (IQR)",
     main = "Range beats sample size")
axis(1, at = c(10, 20, 50, 100, 200))
lines(ns_mm, iqr_wide, type = "b", pch = 1, cex = 0.8, lty = 2)
axis(1, at = ns_mm, labels = FALSE, tcl = 0.2)
legend("topright", legend = c("x in [0, 12]", "x in [0, 40]"),
       pch = c(16, 1), lty = c(1, 2), bg = "white", box.col = "grey80",
       cex = 0.85, inset = 0.02)
fin()
close_plot()

## WHAT IT SHOWS
## The comparison the budget actually faces. Quadrupling the sample size under
## the narrow design leaves K less well determined than n = 10 does under the
## wide design. Twenty times the plants lose to a wider gradient.
##
## The point is stronger than in Section 1.3(d). There, a wider range simply
## improved precision. Here it decides whether the parameter can be estimated
## at all, because the narrow design never reaches the conditions that
## distinguish K from Vmax.



## ===========================================================================
## SECTION 3. WHEN PRECISION IS NOT ENOUGH
## ===========================================================================
##
## Two ways a design can be precise and still wrong: the model has the wrong
## shape, or the parameter is not the one we mean.
## ---------------------------------------------------------------------------

## --- (a) a model with the wrong shape --------------------------------------
fit_wrong <- lm(y ~ x, data = fake_narrow)
r2_wrong  <- summary(fit_wrong)$r.squared
ci_wrong  <- confint(fit_wrong)["x", ]
pred40    <- unname(predict(fit_wrong, newdata = data.frame(x = 40)))
true40    <- mu_mm(40, Vmax, K)

## For comparison: how well does the TRUE curve fit the same data? It is not
## fitted, so this is the ceiling any model could reach here. If the wrong
## model comes close to it, R^2 cannot tell the two apart.
ss_tot   <- sum((fake_narrow$y - mean(fake_narrow$y))^2)
r2_curve <- 1 - sum((fake_narrow$y - mu_mm(fake_narrow$x, Vmax, K))^2) / ss_tot

cat(sprintf("straight line: R2 = %.2f | true curve: R2 = %.2f\n", r2_wrong, r2_curve))
cat(sprintf("line predicts %.0f g at x = 40, truth %.0f g; interval width %.2f\n",
            pred40, true40, diff(ci_wrong)))

## --- (b,c) a confounded observational design -------------------------------
set.seed(10)
n_obs    <- 500
moisture <- rnorm(n_obs, 0, 1)                            # the confounder
x_obs    <- 10 + 3 * moisture + rnorm(n_obs, 0, 1)        # N tracks moisture
y_obs    <- 30 + 7 * x_obs + 15 * moisture + rnorm(n_obs, 0, 5)
obs      <- data.frame(x = x_obs, y = y_obs, moisture = moisture)

naive <- lm(y ~ x, data = obs)
adj   <- lm(y ~ x + moisture, data = obs)
cat(sprintf("true 7.00 | naive %.2f [%.2f, %.2f] | adjusted %.2f\n",
            coef(naive)["x"], confint(naive)["x", 1], confint(naive)["x", 2],
            coef(adj)["x"]))

n_seq <- c(100, 250, 500, 1000, 2500, 5000, 10000)
naive_by_n <- t(sapply(n_seq, function(nn) {
  m  <- rnorm(nn, 0, 1)
  xx <- 10 + 3 * m + rnorm(nn, 0, 1)
  yy <- 30 + 7 * xx + 15 * m + rnorm(nn, 0, 5)
  f  <- lm(yy ~ xx)
  c(est = unname(coef(f)[2]), lo = confint(f)[2, 1], hi = confint(f)[2, 2])
}))

open_plot("fig_limits", height = 2.9, mfrow = c(1, 3))

## (a) the wrong shape
plot(fake_narrow$x, fake_narrow$y, xlim = c(0, 45),
     ylim = c(0, max(pred40, true40) * 1.3), pch = 16, col = "grey30",
     xlab = "Nitrogen", ylab = "Biomass (g)", main = "(a) Wrong shape")
rect(0, 0, 12, max(pred40, true40) * 1.3, col = "grey92", border = NA)
points(fake_narrow$x, fake_narrow$y, pch = 16, col = "grey30")
curve(mu_mm(x, Vmax, K), 0, 45, add = TRUE, lwd = 2)
abline(fit_wrong, lwd = 2, lty = 2, col = "firebrick")
legend("topleft", legend = c("true", "fitted line"), lwd = 2, lty = c(1, 2),
       col = c("black", "firebrick"), bty = "n", cex = 0.8, inset = 0.02)
## Label each R^2 with the object it belongs to, colour-matched to the lines.
text(45, 0.16 * max(pred40, true40) * 1.3,
     bquote("curve:" ~ R^2 == .(sprintf("%.2f", r2_curve))),
     col = "black", pos = 2, cex = 0.78)
text(45, 0.06 * max(pred40, true40) * 1.3,
     bquote("line:" ~ R^2 == .(sprintf("%.2f", r2_wrong))),
     col = "firebrick", pos = 2, cex = 0.78)
fin()

## (b) confounding at n = 500
plot(1, type = "n", xlim = c(6, 13), ylim = c(0.5, 2.5), yaxt = "n",
     xlab = "Estimated effect of N", ylab = "", main = "(b) Confounded, n = 500")
abline(v = 7, lwd = 2, lty = 2)
arrows(confint(naive)["x", 1], 2, confint(naive)["x", 2], 2, code = 3,
       angle = 90, length = 0.04, lwd = 3, col = "firebrick")
arrows(confint(adj)["x", 1], 1, confint(adj)["x", 2], 1, code = 3,
       angle = 90, length = 0.04, lwd = 3, col = "steelblue")
points(coef(naive)["x"], 2, pch = 16, cex = 1.3, col = "firebrick")
points(coef(adj)["x"],   1, pch = 16, cex = 1.3, col = "steelblue")
axis(2, at = c(2, 1), labels = c("naive", "adjusted"))
text(7, 2.4, "truth", pos = 4, cex = 0.8)
fin()

## (c) it gets worse, not better, with n
plot(n_seq, naive_by_n[, "est"], log = "x", type = "n", ylim = c(6, 13),
     xlab = "Sample size (n)", ylab = "Estimated effect of N",
     main = "(c) More data, same bias")
abline(h = 7, lwd = 2, lty = 2)
arrows(n_seq, naive_by_n[, "lo"], n_seq, naive_by_n[, "hi"], code = 3,
       angle = 90, length = 0.025, lwd = 2, col = "firebrick")
points(n_seq, naive_by_n[, "est"], pch = 16, cex = 0.8, col = "firebrick")
text(min(n_seq), 7, "truth", pos = 3, cex = 0.8)
fin()

close_plot()

## WHAT IT SHOWS
## (a) Every conventional diagnostic looks fine: R2 around 0.75, a slope
## estimated to a couple of per cent. Extrapolated to x = 40 the fitted line
## is wrong by more than a factor of two.
##
## The two R2 values are the sharper point. The wrong straight line reaches
## 0.75; the TRUE curve, which is the best any model could do here, reaches
## only 0.79. Over this range R2 barely separates the correct model from a
## badly wrong one, so a good fit statistic is not evidence that the shape is
## right. Precision is not correctness, and neither is goodness of fit.
##
## (b) The naive analysis is precise and excludes the truth. The residuals are
## well behaved and the model IS correctly specified as a straight line. The
## problem is the estimand, not the fit.
##
## (c) The bias does not depend on n. Larger samples tighten the interval
## around the wrong value, so the most confident studies are the ones whose
## bias is hardest to notice.



## ===========================================================================
## SECTION 4. HIERARCHICAL DESIGN -- sites versus plots per site
## ===========================================================================
##
## The case where simulation stops being optional. Plots sit within sites,
## sites differ in baseline biomass AND in their nitrogen response, and a
## fixed budget of 120 plots can be split many ways.
##
##      for each site s:  a_s ~ Normal(30, sd_site)     site intercept
##                        b_s ~ Normal(b,  sd_slope)    site nitrogen effect
##      for each plot i:  y_i = a_s + b_s * x_i + e_i,  e_i ~ Normal(0, sigma)
##
## The key quantity is sd_slope. If every site responded identically, plots
## would be interchangeable and only the total would matter. Because sites
## differ, each site is close to ONE noisy observation of the effect, so the
## number of sites is what limits knowledge of the AVERAGE effect.
##
## Is more sites therefore always better? Over the range below, yes: precision
## improves monotonically as effort moves from plots into sites. The limit is
## not a turning point in the interval but a failure of the model: with fewer
## than about four plots per site there is too little within-site information
## to separate a site's own slope from noise, the random-slope model becomes
## singular, and the fit stops being interpretable at all. Those designs are
## therefore reported by their failure rate rather than by a width.
##
## Requires lme4; skipped if not installed.
## ---------------------------------------------------------------------------

if (!requireNamespace("lme4", quietly = TRUE)) {
  message("Section 4 skipped: install.packages('lme4') to run it.")
} else {

  sim_hier <- function(n_site, n_plot, b = 7, sigma = 5,
                       sd_site = 10, sd_slope = 2) {
    site   <- rep(seq_len(n_site), each = n_plot)
    a_site <- rnorm(n_site, 30, sd_site)[site]
    b_site <- rnorm(n_site, b,  sd_slope)[site]
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
    ## NOTE: not sqrt(diag(vcov(m))). vcov() on an lmerMod is a Matrix S4
    ## object and base diag() fails on it. Use the coefficient table.
    unname(2 * 1.96 * summary(m)$coefficients["x", "Std. Error"])
  }

  ## A fixed budget of 120 plots, allocated five ways. Allocations with fewer
  ## than four plots per site are omitted: there is then too little within-site
  ## information to separate a site's slope from noise, the random-slope model
  ## is singular, and a width averaged over the runs that happened to converge
  ## would flatter a design that does not work. That limit is discussed in the
  ## text rather than plotted.
  designs <- data.frame(n_site = c( 4,  6, 12, 20, 30),
                        n_plot = c(30, 20, 10,  6,  4))
  designs$total <- designs$n_site * designs$n_plot     # 120 throughout

  set.seed(11)
  reps_h <- max(50, R_REPS %/% 5)
  res <- mapply(function(s, p) {
    w <- replicate(reps_h, ci_width_hier(s, p))
    c(width = mean(w, na.rm = TRUE), fail = mean(is.na(w)))
  }, designs$n_site, designs$n_plot)
  designs$width <- res["width", ]
  designs$fail  <- res["fail",  ]
  print(designs)

  open_plot("fig_hierarchical", height = 3.5, width = 5.0,
            mar = c(4.6, 4.2, 2.0, 0.9))
  bp <- barplot(designs$width, col = "grey80", border = NA,
                ylim = c(0, max(designs$width, na.rm = TRUE) * 1.20),
                ylab = "95% interval width",
                main = "Sites vs. plots per site")
  text(bp, designs$width, sprintf("%.2f", designs$width), pos = 3, cex = 0.8)
  axis(1, at = bp, labels = paste0(designs$n_site, "×", designs$n_plot),
       tick = FALSE, line = -0.4, cex.axis = 0.85)
  mtext("sites × plots per site   (120 plots throughout)", side = 1,
        line = 2.4, cex = 0.85)
  fin()
  close_plot()

  ## WHAT IT SHOWS
  ## The same 120 plots, spent seven ways. Across the designs that can be
  ## fitted, precision improves steadily as effort moves from plots into
  ## sites: roughly a two-fold difference in interval width for identical
  ## cost. When sites differ in their response, the number of sites is what
  ## limits knowledge of the average effect, and adding plots within a site
  ## mostly refines a slope you already knew.
  ##
  ## The trend does not continue indefinitely, but the limit appears as a model
  ## that will not fit rather than as a wider interval. Below about four plots
  ## per site there is no within-site information left to separate a site's own
  ## slope from residual noise, the random-slope model becomes singular, and
  ## the estimate stops being interpretable. Those allocations are left out of
  ## the figure for that reason; try adding 60 x 2 to the designs above and
  ## watch the failure rate rather than the width.
  ##
  ## Read the failure rate as a design diagnostic, not a nuisance. A model
  ## that will not converge is saying the design asks the data for something
  ## they cannot supply, and it says so before any field work.
}



## ===========================================================================
## SECTION 5. What to report
## ===========================================================================
##
##   1. The generative model: f(x; theta) and the stochastic component.
##   2. The assumed parameter values, and their justification.
##   3. What the error term represents, and whether hierarchy is included.
##   4. The design variables explored, and which were held fixed.
##   5. The performance criterion and target: interval width at a stated
##      quantile, and the requirement it derives from.
##   6. Failures: the proportion of fits that did not converge.
##   7. Code, seed and software versions, archived with a DOI.
## ---------------------------------------------------------------------------

sessionInfo()

## ===========================================================================
## END
## ===========================================================================
