# --- Setup ---
set.seed(42)
# If needed:
# install.packages("rstanarm")
library(rstanarm)

# --- Design & parameters ---
n_per_group <- 20          # samples per group

a <- 15           # baseline mean (Control)
b <- 5                  # treatment effect (difference in means: Treatment - Control)
sigma <- 2             # common SD of measurement noise

# --- Simulate data for two treatment levels ---
treatment <- factor(rep(c("Control", "Treatment"), each = n_per_group),
                    levels = c("Control","Treatment"))
y <- a + (treatment == "Treatment") * b + rnorm(2 * n_per_group, 0, sigma)

fake <- data.frame(treatment, y)

# visual check 
par(mfrow = c(1,1))
boxplot(y ~ treatment, data = fake, main = "Outcome by Treatment", ylab = "y")
stripchart(y ~ treatment, data = fake, vertical = TRUE, method = "jitter",
           add = TRUE, pch = 20)

# Plot with fitted means (from Bayesian model below)
# (We'll fit first, then add means to a second panel.)
fit_bayes <- stan_glm(y ~ treatment, data = fake, refresh = 0)  # Gaussian model
print(fit_bayes, digits = 2)

# Extract posterior draws and summaries
post <- as.matrix(fit_bayes)
delta_post <- post[, "treatmentTreatment"]        # Treatment - Control
prob_gt0 <- mean(delta_post > 0)
ci95 <- quantile(delta_post, c(0.025, 0.975))

# Fitted means from posterior mean coefficients
a_hat <- coef(fit_bayes)[1]                        # mean(Control)
tau_hat <- coef(fit_bayes)[2]                      # mean(Treatment) - mean(Control)

plot(jitter(as.numeric(treatment), amount = 0.08), y, xaxt = "n",
     xlab = "", ylab = "y", main = "Data with Fitted Group Means", pch = 20)
axis(1, at = c(1,2), labels = levels(treatment))
abline(h = a_hat, lwd = 3)                         # fitted mean: Control
abline(h = a_hat + tau_hat, lwd = 3, lty = 2)      # fitted mean: Treatment
legend("topleft", legend = c("Control mean", "Treatment mean"),
       lwd = c(3,3), lty = c(1,2), bty = "n")

par(mfrow = c(1,1))

cat(sprintf("\nBayesian difference (Treatment - Control):\n  mean = %.2f,  95%% CI = [%.2f, %.2f],  Pr(diff > 0) = %.3f\n",
            mean(delta_post), ci95[1], ci95[2], prob_gt0))

# --- Frequentist check: Welch's t-test (no equal-variance assumption) ---
tt <- t.test(y ~ treatment, var.equal = TRUE)
print(tt)
fake

# Equivalent linear model; get CI and robust SEs (HC3)
m <- lm(y ~ treatment, data = fake)
summary(m)




# --- Optional: effect size (Cohen's d) ---
yC <- y[treatment == "Control"]; yT <- y[treatment == "Treatment"]
s_pooled <- sqrt(((length(yC)-1)*var(yC) + (length(yT)-1)*var(yT)) /
                   (length(yC) + length(yT) - 2))
d <- (mean(yT) - mean(yC)) / s_pooled
cat(sprintf("\nCohen's d ≈ %.2f\n", d))

# --- Basic diagnostics worth a peek ---
par(mfrow = c(1,2))
plot(fitted(fit_bayes), resid(fit_bayes), main = "Residuals vs Fitted",
     xlab = "Fitted", ylab = "Residuals")
abline(h = 0, lty = 2)
qqnorm(resid(fit_bayes)); qqline(resid(fit_bayes))
par(mfrow = c(1,1))
