
make_rhs <- function(y, p , type, season = NULL, exogen = NULL) {
  y <- as.matrix(y) 
  obs <- dim(y)[1]
  K <- dim(y)[2]
  sample <- obs - p
  # making a matrix of lags of y for rhs with the embed function.  
  # ylags becomes a matrix of lag regressors.
  # notice that embed cuts down the number of rows by "dimension minus 1", 
  # e.g. one lag p=1 cuts down by 1 observations row
  # [, -(1:K)] eliminates first K rows that contain lag zero y values (original y values). 
  ylags <- embed(y, dimension = p + 1)[, -(1:K)]
  
  # making names for lagged variable columns e.g. "variable.l1" 
  temp1 <- NULL
  for (i in 1:p) {
    temp <- paste(colnames(y), ".l", i, sep = "")
    temp1 <- c(temp1, temp)
  }
  colnames(ylags) <- temp1
  
  # making addtional columns for const, trend, seasonality
  if (type == "const") {
    rhs <- cbind(ylags, rep(1, sample))
    colnames(rhs) <- c(colnames(ylags), "const")
  }
  else if (type == "trend") {
    rhs <- cbind(ylags, seq(p + 1, length = sample))
    colnames(rhs) <- c(colnames(ylags), "trend")
  }
  else if (type == "both") {
    rhs <- cbind(ylags, rep(1, sample), seq(p + 1, length = sample))
    colnames(rhs) <- c(colnames(ylags), "const", "trend")
  }
  else if (type == "none") {
    rhs <- ylags
    colnames(rhs) <- colnames(ylags)
  }
  
  # binding seasonality 
  if (!(is.null(season))) {
    season <- abs(as.integer(season))
    dum <- (diag(season) - 1/season)[, -season]
    dums <- dum
    while (nrow(dums) < obs) {
      dums <- rbind(dums, dum)
    }
    dums <- dums[1:obs, ]
    colnames(dums) <- paste("sd", 1:ncol(dums), sep = "")
    rhs <- cbind(rhs, dums[-c(1:p), ])
  }
  
  # binding exogenous variables 
  if (!(is.null(exogen))) {
    exogen <- as.matrix(exogen)
    if (!identical(nrow(exogen), nrow(y))) {
      stop("\nDifferent row size of y and exogen.\n")
    }
    if (is.null(colnames(exogen))) {
      colnames(exogen) <- paste("exo", 1:ncol(exogen), 
                                sep = "")
      warning(paste("No column names supplied in exogen, using:", 
                    paste(colnames(exogen), collapse = ", "), ", instead.\n"))
    }
    colnames(exogen) <- make.names(colnames(exogen))
    tmp <- colnames(rhs)
    rhs <- cbind(rhs, exogen[-c(1:p), ])
    colnames(rhs) <- c(tmp, colnames(exogen))
  }
  return(rhs)
}

make_datamat <-  function(y, p , type, season = NULL, exogen = NULL) {
  rhs_matrix <-  make_rhs(y, p = p, type = type, season = season, exogen = exogen )
  rhs_df <- as.data.frame(rhs_matrix)
  colnames(rhs_df ) <- colnames(rhs_matrix)  
  return(rhs_df )
}

VAR_alternative <- function (y, p = 1, type = c("const", "trend", "both", "none"), 
                     season = NULL, exogen = NULL, lag.max = NULL, ic = c("AIC", 
                                                                          "HQ", "SC", "FPE")) 
{
  # embed needs a matrix
  y <- as.matrix(y)
  # check: NAs should be eliminated before running VARS
  if (any(is.na(y))) 
    stop("\nNAs in y.\n")
  # check: should be multivariate
  if (ncol(y) < 2) 
    stop("The matrix 'y' should contain at least two variables. For univariate analysis consider ar() and arima() in package stats.\n")
  # check: if has no  colnames then make some default named.
  if (is.null(colnames(y))) {
    colnames(y) <- paste("y", 1:ncol(y), sep = "")
    warning(paste("No column names supplied in y, using:", 
                  paste(colnames(y), collapse = ", "), ", instead.\n"))
  }
  
  # Make syntactically valid names out of character vectors.
  colnames(y) <- make.names(colnames(y))
  
  #backup y for export at end
  y.orig <- y
  
  #match.arg matches a character arg against a table of candidate values as specified by choices.
  type <- match.arg(type)
  
  # obs is number of rows, K is number of columns 
  obs <- dim(y)[1]
  K <- dim(y)[2]
  
  # lag.max overrides p 
  if (!is.null(lag.max)) {
    lag.max <- abs(as.integer(lag.max))
    ic <- paste(match.arg(ic), "(n)", sep = "")
    p <- VARselect(y, lag.max = lag.max, type = type, season = season, 
                   exogen = exogen)$selection[ic]
  }
  
  # the number of data rows for linear regression is number of rows minus number of lags
  sample <- obs - p
  
  # cutting y for lhs by p observations
  # was this necessary? when we had first K columns of ylags?
  yend <- y[-c(1:p), ]
  
  rhs <-   make_rhs(y=y, p=p, type = type, season = season, exogen = exogen)
  
  # Coerce to a Data Frame 
  #datamat <- as.data.frame(rhs)
  #colnames(datamat) <- colnames(rhs)
  
  datamat  <- make_datamat  (y=y, p=p, type = type, season = season, exogen = exogen)
  
    # empty list to gather the result of lm regressions 
  equation <- list()
  for (i in 1:K) {
    # this loop runs 1 to K linear regressions with lm , for 1 to K and a common regressor dataframe
    # yend is renamed to y to allow the formula in lm() to have "y" on the left hand side
    y <- yend[, i]
    y_colname <- colnames(yend)[i]
    equation[[y_colname]] <- lm(y ~ -1 + ., data = datamat)
    # attach attribute to the lm object that there was an intercept=1, 
    # despite the "~ -1" eliminating the default lm() intercept
    if (any(c("const", "both") %in% type)) {
      attr(equation[[y_colname]]$terms, "intercept") <- 1
    }
  }
  
  # block: returning the result object 
  
  # prepare call description for export
  call <- match.call()
  if ("season" %in% names(call)) 
    call$season <- eval(season)
  
  # collecting everything for the result 
  result <- list(varresult = equation, 
                 datamat = data.frame(cbind(yend, rhs)), 
                 y = y.orig, 
                 type = type, 
                 p = p, 
                 K = K, 
                 obs = sample, 
                 totobs = sample + p, 
                 restrictions = NULL, 
                 call = call)
  class(result) <- "varest"
  return(result)
}

make_y <- function (df, p, i )  {
  as.matrix(df)[-c(1:p), i]
}
