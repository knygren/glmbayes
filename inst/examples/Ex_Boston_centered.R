############################### Boston_centered dataset example ####################

data("Boston_centered")
head(Boston_centered)
summary(Boston_centered)

## Predictors are mean-centered (column means ~0)
predictors <- setdiff(names(Boston_centered), "medv")
colMeans(Boston_centered[predictors])

form <- medv ~
  crim + zn +
  indus + chas + nox + age + dis + rad + tax + ptratio + black + lstat + rm

lm.boston <- lm(form, data = Boston_centered, x = TRUE, y = TRUE)
summary(lm.boston)

ps <- Prior_Setup(form, gaussian(), data = Boston_centered)

lmb.boston <- lmb(
  form,
  data = Boston_centered,
  pfamily = dNormal(
    mu = ps$mu,
    Sigma = ps$Sigma,
    dispersion = ps$dispersion
  )
)
summary(lmb.boston)

lmb.boston_v2 <- lmb(
  form,
  data = Boston_centered,
  pfamily = dNormal_Gamma(
    mu = ps$mu,
    Sigma_0 = ps$Sigma_0,
    shape = ps$shape,
    rate = ps$rate
  )
)
summary(lmb.boston_v2)

## Independent Normal-Gamma (OpenCL path when available off CRAN)
if (identical(Sys.getenv("NOT_CRAN"), "true") && has_opencl()) {
  lmb.boston_v3 <- lmb(
    n = 1000L,
    form,
    data = Boston_centered,
    pfamily = dIndependent_Normal_Gamma(
      ps$mu,
      ps$Sigma,
      shape = ps$shape_ING,
      rate = ps$rate
    ),
    use_parallel = TRUE,
    use_opencl = TRUE,
    verbose = FALSE
  )
  summary(lmb.boston_v3)
} else {
  message("Skipping OpenCL example (CRAN check or OpenCL not built).")
}
###############################################################################
## End of Boston_centered dataset example
###############################################################################
