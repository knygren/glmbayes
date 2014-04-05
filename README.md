# glmbayes

The `glmbayes` package produces iid samples for Bayesian Genereralized Linear Models and is intended as a Bayesian version of the `glm` function for classical models.

##  Details

 Estimation can be performed using three main functions. For models with fixed dispersion parameters, the `rglmb` function is the workhorse function and comes with a minimialistic interface for the input and output. It is also suitable for use as part of block Gibbs sampling procedures. The `glmb` function is essentially a wrapper function for the `rglmb` function that provides an interface closer to that of `glm` function. The `rglmbdisp` function can be leveraged in order to produce samples for the dispersion parameters associated with the gaussian and Gamma link functions. Most methods defined for the output of the `glm` function are also defined for the `rglmb` and `glmb`, and `rglmbdisp` functions (see their respective documentation for details).

For the regression parameters, multivariate normal priors are assumed. Simulation for the gaussian family with the identify link function is performed using standard procedures for multivariate normal densities. For all other families and link functions, simulation is performed using the likelihood subgradient approach of Nygren and Nygren (2006). This approach involves the construction of an enveloping function for the full posterior density followed by accept-reject sampling. For models that are approximately multivariate normal, the expected number of draws required per acceptance are bounded from above as noted in Nygren and Nygren (2006).

Currently implemented models include the gaussian (identity link), poisson/quasipoisson (log link), binomial/quasibinomial (logit, probit, and cloglog links), and Gamma (log link) families. These models all have log-concave likelihood functions that allow us to leverage the likelihood-subgradient approach for the iid sampling. Models that fail to have log-concave likelihood functions are not implemented. Our demos (viewable by entering the `demo()` command) provides examples of each of these families and links.

The current implementation requires separate use of the `rglmbdisp` function in order to generate samples for dispersion parameters (gaussian, Gammma, quasipoisson, quasibinomaial families). Our demos include examples of the joint use of the `rglmb` and `rglmbdisp` to produce samples for both regression and dispersion parameters using two-block Gibbs samplers. As these two-block Gibbs samplers likely are geometrically ergodic, future implementations may incorporate these two-block Gibbs samplers into the `rglmb` and `glmb` functions by leveraging theoretical bounds om convergence rates derived using Rosenthal (1996) type drift and minorization conditions.

The `rglmb` function can also be used in Block-Gibbs sampling implementations for Hierarchical Bayesian models. The demos associated with this package contains examples of such models.
