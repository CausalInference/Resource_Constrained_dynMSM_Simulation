times = 2
shrink = .5 #multiplicative effect on outcome
side = .5 #side effect, i.e. countervailing additive effect on outcome
bound = 0 
treat_bound = .25 #resource constraint (kappa)
inv_logit = function(x){
  exp(x)/(exp(x)+1)
}
N=25000
library(nimble)

#make data from observational regime
X1 = rnorm(N)
#treatment more likely when abs(X) is large
A1 = rbinom(N,1,inv_logit(-.5+abs(X1)))
#X2_0 is the counterfactual X2 under no treatment
X2_0 = rnorm(N,X1)
#If A1=1, X2_0 is multiplied by shrink to give X2
X2 = X2_0*shrink^(A1*(abs(X2_0)>bound))
A2 = rbinom(N,1,inv_logit(-.5+abs(X2)))
#Y_0 is the counterfactual Y under A2=0
Y_0 = rnorm(N,X2)
Y = abs(Y_0*shrink^(A2*(abs(Y_0)>bound))) + (A1+A2)*side
D = A1 + A2
mean(Y)
mean(D)

#fit treatment model and compute adjustment weights
A = c(A1,A2)
X = c(X1,X2)
A_mod = glm(A~abs(X),family='binomial')
probs1 = inv_logit(A_mod$coefficients[1] +  A_mod$coefficients[2]*abs(X1))
probs2 = inv_logit(A_mod$coefficients[1] +  A_mod$coefficients[2]*abs(X2))
A_probs1 = probs1*A1 + (1-probs1)*(1-A1)
A_probs2 = probs2*A2 + (1-probs2)*(1-A2)
prob_weights = 1/(A_probs1*A_probs2)
# followed = vector(mode = 'list', length = n+1)

#Fit dyn-MSMs
n=60
x1s = (0:n)*(3/n)
x2s = (0:n)*(3/n)
xs2 = expand.grid(x1s,x2s)
names(xs2) = c('x1','x2')
# followed = vector(mode = 'list', length = n+1)
weights = vector(mode = 'list', length = n+1)
x = vector(mode = 'list', length = n+1)
y = vector(mode = 'list', length = n+1)
d = vector(mode = 'list', length = n+1)
for(i in 1:nrow(xs2)){
  followed1 = ifelse(abs(X1)>xs2[i,1],A1,1-A1)
  followed2 = ifelse(abs(X2)>xs2[i,2],A2,1-A2)
  followed = followed1 & followed2
  weights[[i]] = prob_weights[which(followed)]
  y[[i]] = Y[which(followed)]
  d[[i]] = D[which(followed)]
  x[[i]] = cbind(rep(xs2[i,1],sum(followed)),rep(xs2[i,2],sum(followed)))
}

y = unlist(y)
x = do.call('rbind',x)
d = unlist(d)
weights = unlist(weights)
data = data.frame(cbind(y,x))
names(data) = c('y','x1','x2')
dynmsm_y_2x = lm(y ~ (x1+x2)^2 + I(x1^2) + I(x1^3) + I(x1^4) + I(x1^5) + I(x2^2) + I(x2^3) + I(x2^4) + I(x2^5) + I((x1*x2)^2) + I((x1*x2)^3) + I((x1*x2)^4),weights=unlist(weights),data=data)
dynmsm_d_2x = lm(d ~ (x1+x2)^2 + I(x1^2) + I(x1^3) + I(x1^4) + I(x1^5) + I(x2^2) + I(x2^3) + I(x2^4) + I(x2^5) + I((x1*x2)^2) + I((x1*x2)^3) + I((x1*x2)^4),weights=unlist(weights),data=data)

y_ests2 = predict(dynmsm_y_2x,xs2)
d_ests2 = predict(dynmsm_d_2x,xs2)

#Find best constrained and unconstrained strategies through nonlinear optimization
library(nloptr)
eval_f = function(x){
  list("objective" = dynmsm_y_2x$coefficients %*% c(1,x[1],x[2],x[1]^2,x[1]^3,x[1]^4,x[1]^5,x[2]^2,
                                 x[2]^3,x[2]^4,x[2]^5,(x[1]*x[2])^2,(x[1]*x[2])^3,
                                 (x[1]*x[2])^4,x[1]*x[2]),
       "gradient" = c(dynmsm_y_2x$coefficients[2] + 2*dynmsm_y_2x$coefficients[4]*x[1] + 3*dynmsm_y_2x$coefficients[5]*x[1]^2 +
                    4*dynmsm_y_2x$coefficients[6]*x[1]^3 + 5*dynmsm_y_2x$coefficients[7]*x[1]^4 +
                    2*dynmsm_y_2x$coefficients[12]*(x[2]^2)*x[1] + 3*dynmsm_y_2x$coefficients[13]*(x[2]^3)*x[1]^2 +
                    4*dynmsm_y_2x$coefficients[14]*(x[2]^4)*x[1]^3,
                    dynmsm_y_2x$coefficients[3] + 2*dynmsm_y_2x$coefficients[8]*x[2] + 3*dynmsm_y_2x$coefficients[9]*x[2]^2 +
                      4*dynmsm_y_2x$coefficients[10]*x[2]^3 + 5*dynmsm_y_2x$coefficients[11]*x[2]^4 +
                      2*dynmsm_y_2x$coefficients[12]*(x[1]^2)*x[2] + 3*dynmsm_y_2x$coefficients[13]*(x[1]^3)*x[2]^2 +
                      4*dynmsm_y_2x$coefficients[14]*(x[1]^4)*x[2]^3)
  )
}
eval_g_ineq = function(x){
  constr = dynmsm_d_2x$coefficients %*% c(1,x[1],x[2],x[1]^2,x[1]^3,x[1]^4,x[1]^5,x[2]^2,
                                 x[2]^3,x[2]^4,x[2]^5,(x[1]*x[2])^2,(x[1]*x[2])^3,
                                 (x[1]*x[2])^4,x[1]*x[2]) - .25
  grad = c(dynmsm_d_2x$coefficients[2] + 2*dynmsm_d_2x$coefficients[4]*x[1] + 3*dynmsm_d_2x$coefficients[5]*x[1]^2 +
             4*dynmsm_d_2x$coefficients[6]*x[1]^3 + 5*dynmsm_d_2x$coefficients[7]*x[1]^4 +
             2*dynmsm_d_2x$coefficients[12]*(x[2]^2)*x[1] + 3*dynmsm_d_2x$coefficients[13]*(x[2]^3)*x[1]^2 +
             4*dynmsm_d_2x$coefficients[14]*(x[2]^4)*x[1]^3,
           dynmsm_d_2x$coefficients[3] + 2*dynmsm_d_2x$coefficients[8]*x[2] + 3*dynmsm_d_2x$coefficients[9]*x[2]^2 +
             4*dynmsm_d_2x$coefficients[10]*x[2]^3 + 5*dynmsm_d_2x$coefficients[11]*x[2]^4 +
             2*dynmsm_d_2x$coefficients[12]*(x[1]^2)*x[2] + 3*dynmsm_d_2x$coefficients[13]*(x[1]^3)*x[2]^2 +
             4*dynmsm_d_2x$coefficients[14]*(x[1]^4)*x[2]^3)
  return( list( "constraints"=constr, "jacobian"=grad ) )
}
lb = c(0,0)
ub = c(3,3)
local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                    "xtol_rel" = 1.0e-7 )
opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
              "xtol_rel" = 1.0e-7,
              "maxeval" = 1000,
              "local_opts" = local_opts )
res_con <- nloptr( x0=c(2.95,1.55),
               eval_f=eval_f,
               lb=lb,
               ub=ub,
               eval_g_ineq=eval_g_ineq,
               opts=opts)
#estimated optimal strategy under resource constraints
opt_x_rc2 = res_con$solution
#estimated expected outcome under constrained optimal strategy
opt_y_rc2 = predict(dynmsm_y_2x,data.frame(x1=opt_x_rc2[1],x2=opt_x_rc2[2]))

res = optim(c(2,2),function(x)eval_f(x)[[1]],method="L-BFGS-B",lower=0,upper=3)
#estimated expected outcome under unconstrained optimal strategy
opt_y2 = res$value
#estimated unconstrained optimal strategy
opt_x2 = res$par


#compute actual counterfactual mean outcome under estimated unconstrained optimal strategy
N=1000000
threshes = opt_x2
X1 = rnorm(N)
A1 = abs(X1)>threshes[1]
X2_0 = rnorm(N,X1)
X2 = (shrink^A1)*X2_0
A2 = abs(X2)>threshes[2]
Y_0 = rnorm(N,X2)
Y = abs(Y_0*shrink^(A2*(abs(Y_0)>bound))) + (A1+A2)*side
D = A1 + A2
actual_opt_y_x_hat = mean(Y)

#compute actual counterfactual mean outcome under estimated constrained optimal strategy
N=1000000
threshes = opt_x_rc2
X1 = rnorm(N)
A1 = abs(X1)>threshes[1]
X2_0 = rnorm(N,X1)
X2 = (shrink^A1)*X2_0
A2 = abs(X2)>threshes[2]
Y_0 = rnorm(N,X2)
Y = abs(Y_0*shrink^(A2*(abs(Y_0)>bound))) + (A1+A2)*side
D = A1 + A2
actual_opt_y_x_hat_rc = mean(Y)

#compute actual counterfactual mean outcome under estimated unconstrained optimal strategy 
#such that when resources run out, then no treatment is given
N=1000000
threshes = opt_x2
X1 = rnorm(N)
A1 = abs(X1)>threshes[1]
A1_cumsum = cumsum(A1)
if(sum(A1)>(N*treat_bound)){
  stop_ind = min(which(cumsum(A1)>(N*treat_bound)))
  A1[stop_ind:N] = 0
  X2_0 = rnorm(N,X1)
  X2 = (shrink^A1)*X2_0
  A2 = 0
  Y_0 = rnorm(N,X2)
  Y = abs(Y_0*shrink^(A2*(abs(Y_0)>bound))) + (A1+A2)*side
  D = A1 + A2
}else{
  X2_0 = rnorm(N,X1)
  X2 = (shrink^A1)*X2_0
  A2 = abs(X2)>threshes[2]
  if(sum(A1)+sum(A2) > (N*treat_bound)){
    stop_ind = min(which(sum(A1)+cumsum(A2)>N*treat_bound))
    A2[stop_ind:N]=0
  }
  Y_0 = rnorm(N,X2)
  Y = abs(Y_0*shrink^(A2*(abs(Y_0)>bound))) + (A1+A2)*side
  D = A1 + A2
}
actual_naive = mean(Y)
mean(D)

