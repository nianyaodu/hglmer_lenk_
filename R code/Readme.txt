Hglmer_Lenk

* ordered_dirichlet.R : Generate random sample from ordered dirichlet distribution
* hglmer_linear.R : Fit finite mixtures of glm with random effects in case of linear regression
* hglmer_linear_simulation.R : Run simulation for hglmer_linear under setting of Section 5.1
* hglmer_logistic.R : Fit finite mixtures of glm with random effects in case of logistic regression
* hglmer_logistic_simulation.R : Run simulation for hglmer_logistic under setting of Section 5.1

Notes

(1) The 'hglmer_linear.R' file implements a function to compute the log marginal likelihood for the corresponding model. A preliminary analysis shows that the results of marginal likelihoods are largely consistent with those reported in the paper, although the inferred number of clusters K is not stable.

(2) As an alternative, we consider the following two-stage approach:
1. Estimate beta_i using either lm() or glm()
2. Apply the EM algorithm to fit a Gaussian mixture model and determine K

This heuristic performs well in the linear regression setting but tends to overestimate the number of clusters in the logistic regression case.

(3) To address the label switching issue, we use point estimates from lm(), glm(), and the EM algorithm as initial values for the MCMC chain.
