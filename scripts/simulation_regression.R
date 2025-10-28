library(broom)
library(dplyr)

set.seed(608)

# define some parameters
n_groups = 4
n_per_group = 40
n_pred = 10
sigma_noise = 1.2
target_preds = c("P3", "P8") # choose significant predictors

# generate dataset
X = matrix(rnorm(n_groups*n_per_group*n_pred), ncol=n_pred)
colnames(X) = paste0("P", 1:n_pred)
  
beta_true = rep(0, n_pred)
beta_true[3] =  0.20
beta_true[8] = -0.20
  
y = as.numeric(X %*% beta_true + rnorm(n_groups*n_per_group, sd=sigma_noise))
  
dat = as.data.frame(X)
dat$y = y
dat$group = rep(paste0("G", 1:n_groups), each=n_per_group)

# get pooled model
formula_str = paste("y ~", paste0("P", 1:n_pred, collapse=" + "))
mod_pooled = lm(as.formula(formula_str), data=dat)
pooled_res = tidy(mod_pooled)
pooled_res = pooled_res[pooled_res$term != "(Intercept)",]

# The model using the entire sample shows P3 and P8 as significant predictors
print(pooled_res)

# term  estimate std.error statistic p.value
# P1     0.16     0.10     1.56      0.12 
# P2     0.05     0.11     0.44      0.66
# P3     0.27     0.11     2.55      0.01 **
# P4    -0.05     0.10    -0.56      0.58
# P5    -0.11     0.10    -1.16      0.25 
# P6     0.09     0.11     0.76      0.45 
# P7     0.18     0.11     1.60      0.11 
# P8    -0.24     0.12    -2.04      0.04 *
# P9    -0.18     0.11    -1.72      0.09
# P10    0.00     0.11     0.01      0.99

# group-level models
group_res_list = list()
for (g in unique(dat$group)) {
  mod = lm(as.formula(formula_str), data=dat[dat$group==g,])
  tmp = tidy(mod)
  tmp = tmp[tmp$term != "(Intercept)",]
  tmp$group = g
  group_res_list[[g]] = tmp
}
group_res = bind_rows(group_res_list)

group_min = group_res %>% group_by(group) %>% summarise(min_p=min(p.value), .groups="drop")

print(group_min)

# Note that no individual group has significant predictors
# group min_p
# G1    0.08
# G2    0.05
# G3    0.10 
# G4    0.19 



