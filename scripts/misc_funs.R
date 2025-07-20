normData <- function(x, cols=NULL, inv=NULL, absl=NULL, na.rm=TRUE) {
  x_std = x
  if (is.null(cols)) { cols=colnames(x)[sapply(x,is.numeric)] }
  x2 = x[cols]
  if (!is.null(inv)) { x2[inv] = max(x2[inv]) - x2[inv] }
  if (!is.null(absl)) { x2[absl] = abs(x2[absl]) }
  temp = (x2 - colMeans(x2,na.rm=na.rm)[col(x2)]) / sapply(x2, sd, na.rm)[col(x2)]
  x_std[colnames(temp)] = temp
  return(x_std)
}

scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)

composite_score <- function(x, cols=NULL, inv=NULL, absl=NULL, keep_ori=NULL) {
  x2 = x
  if (is.null(cols)) { cols=colnames(x2)[sapply(x2,is.numeric)] }
  if (!is.null(inv)) { x2[inv] = max(x2[inv]) - x2[inv] }
  if (!is.null(absl)) { x2[absl] = abs(x2[absl]) }
  x2$composite = rowSums(x2[cols])
  x2[keep_ori] = x[keep_ori]
  return(x2)
}


gen_combs = function(x, n='all') {
  if (n=='all') { n = length(x) }
  if (length(n)==1) {
    l = do.call(c, lapply(n, combn, x=x, simplify=FALSE))
  } else {
    l = do.call(c, lapply(1:n, combn, x=x, simplify=FALSE))
  }
  return(l)
}


get_demogs = function(files, demo) {
  df_demo = data.frame()
  for (f in files) {
    temp = read.csv(f, sep="\t")
    demo = unique(demo[demo %in% colnames(temp)])
    df_demo = rbind(df_demo, temp[,demo])
  }
  return(df_demo)
}


em_con_vec = function(pairs_em, con) {
  b = str_split(con,"-")
  v1 = as.integer(pairs_em==b[[1]][1])
  v2 = as.integer(pairs_em==b[[1]][2])*-1
  con_vec = v1 + v2
  return(con_vec)
}

# Define the function to calculate the P values
ridge.pvalues <- function(data, indices, lambda) {
  x <- data$x[indices, ]
  y <- data$y[indices]
  ridge.mod <- glmnet(x, y, alpha = 0, lambda = lambda)
  beta <- as.numeric(ridge.mod$beta)
  se <- sqrt(sum((y - x %*% beta)^2) / (nrow(x) - ncol(x))) * sqrt(diag(solve(t(x) %*% x + lambda * diag(ncol(x)))))
  pvalues <- 2 * (1 - pnorm(abs(beta / se)))
  return(pvalues)
}

## inspired by `lme4:::coef.merMod`
coef_ordinal <- function (object, ...) 
{
  if (length(list(...))) 
    warning("arguments named \"", paste(names(list(...)), 
                                        collapse = ", "), "\" ignored")
  fef <- data.frame(rbind(object$beta), check.names = FALSE)  ## adapted
  ref <- ordinal::ranef(object, condVar = FALSE)  ## adapted
  refnames <- unlist(lapply(ref, colnames))
  nmiss <- length(missnames <- setdiff(refnames, names(fef)))
  if (nmiss > 0) {
    fillvars <- setNames(data.frame(rbind(rep(0, nmiss))), 
                         missnames)
    fef <- cbind(fillvars, fef)
  }
  val <- lapply(ref, function(x) fef[rep.int(1L, nrow(x)), 
                                     , drop = FALSE])
  for (i in seq(a = val)) {
    refi <- ref[[i]]
    row.names(val[[i]]) <- row.names(refi)
    nmsi <- colnames(refi)
    if (!all(nmsi %in% names(fef))) 
      stop("unable to align random and fixed effects")
    for (nm in nmsi) val[[i]][[nm]] <- val[[i]][[nm]] + refi[, 
                                                             nm]
  }
  val
}

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

dataf <- function(x,y,lbpath,force.vars,tymea,ws,ns,df_pred) {

	fit.lasso = glmnet(x,y,lambda=lbd, alpha=1, penalty.factor=force.vars, type.measure=tymea, weights=ws)
    
	results = coef(fit.lasso)

    yhat = predict(fit.lasso, newx=x)

    temp = data.frame(modality=ns), yhat=yhat, y=y)
    df_pred = rbind(df_pred, temp)

    preds = rownames(results)
    exclpreds = preds[grepl(paste(c('(Intercept)',covariates), collapse="|"),preds)]

    df_res = as.data.frame(as.matrix(results)) %>%
      filter_all(any_vars(.!=0)) %>% rownames_to_column(var="connection") %>%
      filter(!(connection %in% exclpreds))
    if (dim(df_res)[1]==0) {df_res[1,]=NA}
    df_res['modality'] = names(l[n])
    df_res['groups'] = paste(AGs,collapse='_')
    df_res['N'] = nrow(mdf)
    df_res['devrat'] = fit.lasso$dev.ratio

    df_comb = rbind(df_comb,df_res)
}