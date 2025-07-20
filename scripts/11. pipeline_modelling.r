library(dplyr)
library(tidyr)
library(tibble)
library(glmnet)
library(caret)
library(emmeans)
library(car)
require(doMC)

# Set-Up and Configuration

cl <- parallel::makeCluster(24)
doParallel::registerDoParallel(cl)

norm = 1

tasK = "AC"

cutoff = 50

root_dir = "E:"
base_dir = file.path(root_dir, "data")
group_dir =  file.path(base_dir, "group")
fig_dir =  file.path(root_dir,"figures")

# load helper functions: normData, scale2, composite_score, etc.
source(file.path(base_dir,"scripts","misc_funs.R"))

# Data Loading and Preprocessing

dat_learn_EDA = read.csv(file.path(base_dir, "datasets", "stat-bf_desc-EDA_df.csv"), sep='\t')

dat_learn_beh = read.csv(file.path(base_dir, "datasets", "stat-logit_desc-beh_df.tsv"), sep='\t')
dat_learn = bind_rows(dat_learn_EDA, dat_learn_beh) %>% filter(!is.na(learning))

covariates = c("age","sex")

demographics = read.csv2(file.path(base_dir, "datasets", "demographics.csv"), sep='\t')

comb_AG = gen_combs(unique(dat_learn$AG))
comb_AG = list(
  # c("A08"), # PL
  # c("A03"), # FLr
  # c("A03","A05"), # FLc
  # c("A03","A05","A09"),# FLs
  # c("A02","A03","A05","A09"), # FLel
  # c("A03","A05","A09","A12"), # FLst
  # c("A02","A03","A05","A09","A12"), # FL
  c("A02","A03","A05","A08","A09","A12")) # All


# dupls = read.csv(file.path(base_dir, "desc-duplicates_df.tsv"), sep='\t')
# excl_subs = dupls[dupls$keep==0,c("participant","AG","study")]
# dat_learn = anti_join(dat_learn, excl_subs, by=c("participant","AG","study"))

dat_learn_corrs = dat_learn %>% group_by(AG,study,task) %>% group_modify(~ normData(.x))

df_learndemo = demographics %>%
  right_join(dat_learn, by=c("participant","AG","study","task"), multiple='all')

# dat_excl_conn = read.csv(file.path(base_dir, "desc-ExcludeSubsConn_table.tsv"), sep='\t')

# Connectivity Data Handling

dat_FC = read.csv(file.path(base_dir, "datasets", "desc-FC_df.tsv"), sep='\t')
# dat_FC = subset(dat_FC, !(participant %in% dat_excl_conn[dat_excl_conn$pipeline=="FC",]$participant_id))

dat_DTI = read.csv(file.path(group_dir, "SC", "desc-SC_df.tsv"), sep='\t')
dat_DTI = subset(dat_DTI, !(participant %in% dat_excl_conn[dat_excl_conn$pipeline=="DTI",]$participant_id))

dat_spDCM = read.csv(file.path(group_dir,"EC","desc-EC_df.tsv"), sep='\t')
dat_spDCM = subset(dat_spDCM, !(participant %in% dat_excl_conn[dat_excl_conn$pipeline=="spDCM",]$participant_id))

df_comb = data.frame(); df_coefs = data.frame(); df_pred = data.frame(); df_counts = data.frame()

# Group-wise and Modality-wise Processing

l = list(FC=dat_FC, SC=dat_DTI, EC=dat_spDCM)
l = list(FC=dat_FC)
for (nAG in 1:length(comb_AG)) {
  AGs = comb_AG[[nAG]]
  
  for (n in 1:length(l)) {
    
    if (names(l[n]) == 'FC') {
      metric = c('corrLW','xcorr','EuclideanDist','ManhattanDist','WassersteinDist','dtw','MI','mscohe','wavcohe')
      inv = c('EuclideanDist','ManhattanDist','WassersteinDist','dtw')
      absl = c('corrLW')
      dat_conn = dat_FC
      pair='pair_und'
    } else if (names(l[n]) == 'SC') {
      metric = 'streamlines'
      dat_conn = dat_DTI
      pair='pair_und'
    } else {
      metric = 'spDCM'
      dat_conn = dat_spDCM
      pair='pair_dir'
    }
    
    dat_conn = subset(dat_conn, hemisphere!="bilateral")
    
    df_join = right_join(df_learndemo, dat_conn, by=c("participant","AG","study"), multiple='all')
    
    ####### CHOOSE WHICH AGs TO INCLUDE IN THE ANALYSIS HERE #######
    df_join = df_join %>% subset(AG %in% AGs)
    if (dim(df_join)[1]==0) {next}
    ################################################################
    
    df_sel = df_join %>% group_by(participant, AG, study) %>%
      filter(if_all(!!metric, ~ all(!is.na(.x))))
    
    # df_sel will contain the final sample - that is, excluding subs with any NAs
    # Calculate the composite score
    if (names(l[n]) == 'FC') {
      df_sel = composite_score(df_sel, cols=metric, inv=inv, absl=absl, keep_ori="corrLW")
      metric = "composite"
    }
    
    selvars = unique(c("participant","AG","study","task",pair,metric,"learning",covariates))
    
    # Standardisation of Learning/connectivity Estimates
    df = df_sel[selvars] %>% group_by(AG,study,task,across(all_of(pair))) %>% group_modify(~ normData(.x)) %>%
      pivot_wider(names_from=all_of(pair), values_from=all_of(metric), values_fn=mean) %>% ungroup()
    
    if (names(l[n]) == 'EC') {
      mVp = df_sel %>% group_by(participant, AG, study) %>% mutate(mVp=mean(spDCM_Var)) %>%
        select(participant, AG, study, mVp) %>% distinct()
      mVp$mVp = 1-range01(mVp$mVp)
      df = left_join(df, select(mVp, mVp))
    }
    
    df$const <- factor(rep(1, each=length(df$participant)))
    
    df_acq = subset(df, task == "acquisition")
    df_ext = subset(df, task == "extinction")
    df_ren = subset(df, task == "renewal")
    
    if (tasK=="AC") {
      mdf = df_acq
    } else if (tasK=="EX") {
      mdf = df_ext
    } else {
      mdf = df_ren
    }
    
    pairs = colnames(mdf %>% dplyr::select(starts_with(
      c('AMY','CEB','HIP','ACC','PFC','lAMY','lCEB','lHIP','lACC','lPFC','rAMY','rCEB','rHIP','rACC','rPFC'))))
    
    # LASSO Regression Model Setup
    
    y <- mdf$learning
    xx <- mdf %>% ungroup() %>% dplyr::select(all_of(c(pairs,covariates)))
    x = data.matrix(makeX(xx, na.impute = TRUE))
    myalpha = 1
    
    if (!is.null(covariates)) {
      force.vars = as.integer(!Reduce('|', lapply(covariates, function(y) startsWith(as.character(colnames(x)), y))))
    } else {
      force.vars = rep(1, ncol(xx))
    }
    
    # Cross-Validation and Lambda Optimization
    
    tymea = "mse"
    lambda_max <- max(abs(colSums(x*y,na.rm=T)))/nrow(x)
    epsilon <- .0001
    K <- 1000
    lambdapath <- round(exp(seq(log(lambda_max), log(lambda_max*epsilon), length.out = K)), digits = 10)
    lbpath = lambdapath
    
    ls = foreach(i = 1:100, .combine='rbind', .packages="glmnet") %dopar% {
      fit <- cv.glmnet(x, y, alpha=myalpha, nfolds=10, standardize=T, penalty.factor=force.vars,
                       type.measure=tymea, parallel=T, lambda=lbpath)
      errors = data.frame(fit$lambda,fit$cvm)
    }
    
    ls <- aggregate(ls[, 2], list(ls$fit.lambda), mean)
    
    bestindex = which(ls[2]==min(ls[2]))[1]
    lbd = ls[bestindex,1]
    
    if (names(l[n]) == 'EC') {
      ws = 1-mdf$mVp
    } else {ws = NULL}
    
	
    # df_res = dataf(x,y,lbd,force.vars,tymea,ws,names(l[n]),df_pred)
    # final_fit = rbind(df_comb,df_res)
    
    # Custom Nested Cross-Validation
    fit_lasso_cus = function (x, y, lambda, penalty.factor, type.measure, weights) {
      
      results = data.frame()
      final_df = data.frame()
      pred_df = data.frame()
      
      for (i in 1:10) {
        cv_folds <- vfold_cv(data.frame(x), v = 10)
        
        for (j in seq_along(cv_folds$splits)) {
          split <- cv_folds$splits[[j]]
          
          train_data <- analysis(split)
          test_data  <- assessment(split)
          
          x_train <- data.matrix(train_data)
          y_train <- y[as.integer(rownames(train_data))]
          
          x_test <- data.matrix(test_data)
          test_idx <- as.integer(rownames(test_data))
          y_test <- y[test_idx]
          
          # Fit model
          fit <- glmnet(x_train, y_train, alpha = 1, lambda = lambda, penalty.factor = penalty.factor,
                        type.measure=type.measure, weights=weights)
          
          # Predict on test set
          y_pred <- predict(fit, newx = x_test)
          
          # Save predictions
          pred_df = bind_rows(pred_df,
                              data.frame(
                                iteration = i,
                                fold = j,
                                subject = test_idx,
                                y_true = y_test,
                                y_pred = as.vector(y_pred)
                              )
          )
          
          # Save coefficients
          coefs = coef(fit)
          coef_df = as.data.frame(as.matrix(coefs)) %>%
            rownames_to_column("connection") %>%
            rename(coefficient = 2) %>%
            filter(coefficient != 0) %>%
            mutate(iteration = i, fold = j)
          
          final_df <- bind_rows(final_df, coef_df)
        }
      }
      
      group_df = final_df %>%
        group_by(connection) %>%
        summarise(
          times_selected = n(),
          avg_coefficient = mean(coefficient)
        ) %>%
        arrange(desc(times_selected))
      
      return(group_df)
      
    }
    
    ffit = fit_lasso_cus(x,y,lbpath,force.vars,tymea,ws) %>% filter(times_selected > cutoff)
    ffit$mod=names(l[n])
    ffit$group=paste(AGs,collapse='_')
    
    # Poisson regression
    
    # Set-Up and Parallelized Simulation
    
    tp = foreach(niter = 1:10, .combine='rbind', .packages=c("dplyr","doMC","tibble","glmnet")) %dopar% {
      
      mdf2 = mdf[sample(nrow(mdf), 80),]
      
      lambda_max <- max(abs(colSums(x*y,na.rm=T)))/nrow(x)
      lbpath <- round(exp(seq(log(lambda_max), log(lambda_max*epsilon), length.out = K)), digits = 10)
      
      # LASSO Regression and Feature Selection
      
      ls = foreach(i = 1:100, .combine='rbind', .packages="glmnet") %dopar% {
        fit <- cv.glmnet(x, y, alpha=myalpha, nfolds=10, standardize=T, penalty.factor=force.vars,
                         type.measure=tymea, parallel=T, lambda=lbpath)
        errors = data.frame(fit$lambda,fit$cvm)
      }
      
      ls <- aggregate(ls[, 2], list(ls$fit.lambda), mean)
      
      bestindex = which(ls[2]==min(ls[2]))[1]
      lbd = ls[bestindex,1]
      
      if (names(l[n]) == 'EC') {
        ws = 1-mdf2$mVp
      } else {ws = NULL}
      
      fit.lasso = glmnet(x=x,y=y,lambda=lbd, alpha=myalpha, penalty.factor=force.vars, type.measure=tymea,
                         weights=ws)
      results = coef(fit.lasso)
      
      df_res = as.data.frame(as.matrix(results)) %>% rownames_to_column(var="connection")
      
      preds = df_res$connection
      
      # Feature Region Analysis
      df_tempc = data.frame(G = niter, AMY=length(grep("AMY", preds)), ACC=length(grep("ACC", preds)), HIP=length(grep("HIP", preds)), PFC=length(grep("PFC", preds)),
                            CEB=length(grep("CEB", preds)), mod=names(l[n]))
      
    } 
    df_counts = rbind(df_counts,tp)
  }
}
# Statistical Comparison

m = df_counts %>% pivot_longer(cols = -c(G,mod)) %>% group_by(G,name)

model <- glmer(value ~ name + (1|G), family=poisson, data=subset(m, mod==mg))
summary(model)
em = emmeans(model, pairwise ~ name, adjust="fdr")
