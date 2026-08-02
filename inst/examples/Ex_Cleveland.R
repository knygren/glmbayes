############################### Start of Cleveland dataset example ####################

data("Cleveland")
head(Cleveland)
summary(Cleveland)

# OpenCL-accelerated Bayesian logistic regression example
# This only runs off CRAN when OpenCL is available
if (identical(Sys.getenv("NOT_CRAN"), "true") && has_opencl()) {
  ps <- Prior_Setup(
    hd ~ age + sex + cp + trestbps + chol +
      fbs + restecg + thalach + exang + oldpeak + slope + ca + thal,
    family = binomial(logit),
    data = Cleveland
  )

  fit <- glmb(
    hd ~ age + sex + cp + trestbps + chol +
      fbs + restecg + thalach + exang + oldpeak + slope + ca + thal,
    family       = binomial(link = "logit"),
    pfamily      = dNormal(mu = ps$mu, Sigma = ps$Sigma),
    data         = Cleveland,
    n            = 1000,
    Gridtype     = 2,
    use_parallel = TRUE,
    use_opencl   = TRUE,
    verbose      = FALSE
  )
  summary(fit)
} else {
  message("Skipping OpenCL example (CRAN check or OpenCL not built).")
}
###############################################################################
## End of Cleveland dataset example
###############################################################################
