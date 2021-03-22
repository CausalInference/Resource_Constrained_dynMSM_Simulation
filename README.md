This repository contains R code implementing the simulation experiments from “Estimating Optimal Dynamic Treatment Strategies Under Resource Constraints Using Dynamic Marginal Structural Models”. Relevant generated objects and their interpretations are listed below.

•	opt_x_rc2 is the estimated optimal strategy under resource constraints
•	opt_y_rc2 is the estimated expected outcome under the constrained optimal strategy

•	opt_x2 is the estimated unconstrained optimal strategy 
•	opt_y2 is the estimated expected outcome under the unconstrained optimal strategy

•	actual_opt_y_x_hat is the actual (computed via Monte Carlo simulation) counterfactual mean outcome under the estimated unconstrained optimal strategy.

•	actual_opt_y_x_hat_rc is the actual (computed via Monte Carlo simulation) counterfactual mean outcome under the estimated resource constrained optimal strategy

•	actual_naive is the actual (computed via Monte Carlo simulation) counterfactual mean under the strategy: “follow the optimal unconstrained strategy until resources run out, then administer no further treatments”.


