# ---------------- Setup ----------------
set.seed(42)
# install.packages("emmeans")   # run once if needed
library(emmeans)

# ---------------- Simulate count data (overdispersed by design) -------------
n <- 60
treatment <- factor(rep(c("Control","Drought"), each = n), levels = c("Control","Drought"))
mu_control <- 30
Delta_true <- -10   # absolute difference wanted (Drought - Control)
mu_drought <- max(0.1, mu_control + Delta_true)
theta_true <- 2                        # overdispersion (NB parameter)

# generate NB counts, but we'll first (mis)fit Poisson on purpose to see dispersion
buds <- rnbinom(2*n, size = theta_true, mu = ifelse(treatment=="Control", mu_control, mu_drought))
dat <- data.frame(treatment, buds)

# ---------------- Base-R plot: boxplot + stripchart -------------------------
boxplot(buds ~ treatment, data = dat, main = "Buds by Treatment", ylab = "Count")
stripchart(buds ~ treatment, data = dat, vertical = TRUE, method = "jitter",
           pch = 20, col = "gray40", add = TRUE)

# ---------------- GLMs with explicit family -------------------------------
# Poisson GLM (equidispersion)
m_pois <- glm(buds ~ treatment, family = poisson, data = dat)
summary(m_pois)
# Quasi-Poisson (adjusts SE for overdispersion; same mean model as Poisson)
m_qpois <- glm(buds ~ treatment, family = quasipoisson, data = dat)
summary(m_qpois)
# ---------------- Absolute difference on response scale --------------------
# Use emmeans to get group means on response scale, then take their difference
emm_pois  <- emmeans(m_pois,  ~ treatment, type = "response")
emm_pois
emm_qpois <- emmeans(m_qpois, ~ treatment, type = "response")
emm_qpois







# ---------------- Quick overdispersion check (Poisson) ---------------------
overdisp <- sum(residuals(m_pois, type = "pearson")^2) / df.residual(m_pois)
cat(sprintf("\nPoisson overdispersion estimate (Pearson χ^2 / df): %.2f\n", overdisp))

# ---------------- Optional: Negative Binomial (if overdispersion >> 1) -----
# install.packages("MASS")
# library(MASS)
# m_nb <- MASS::glm.nb(buds ~ treatment, data = dat)
# summary(contrast(emmeans(m_nb, ~ treatment, type = "response"), "revpairwise"), infer = c(TRUE, TRUE))
