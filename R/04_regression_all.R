utils::globalVariables(c("yvar", "chans"))

#' Wrapper to SVM training. Defined separetely to avoid passing too many objects in parLapplyLB
#' @param x passed from fit_regressions
#' @param params passed from fit_regressions
#' @importFrom stats predict
#' @export
#' @return A list with two elements: predictions and a fitted model
#' @examples
#' fitter_svm()
fitter_svm <- function(x = NULL, params = NULL){
  if(!requireNamespace("e1071", quietly = TRUE)){
    stop("Please run install.packages('e1071')")
  }
  if(!is.null(x) & !is.null(params)){
    w <- x[,"train_set"]==1
    model <- do.call(function(...){e1071::svm(...,x=x[w,chans],y=x[w,yvar])},params)
    pred <- predict(model,x[,chans])
    rm(list=setdiff(ls(),c("pred","model")))
    return(list(pred=pred,model=model))
  }
}

#' Wrapper to XGBoost training. Defined separetely to avoid passing too many objects in parLapplyLB
#' @param x passed from fit_regressions
#' @param params passed from fit_regressions
#' @importFrom stats predict
#' @export
#' @return A list with two elements: predictions and a fitted model
#' @examples
#' fitter_xgboost()

fitter_xgboost <- function(x = NULL, params = NULL){
  if(!requireNamespace("xgboost", quietly = TRUE)){
    stop("Please run install.packages(\"xgboost\")")
  }
  
  if(!is.null(x) & !is.null(params)){
    w <- x[,"train_set"]==1
    args <- c(list(data = x[w, chans], label = x[w, yvar], nthread = 1L, verbose = 0), params)
    model <- do.call(xgboost::xgboost, args)
    pred <- predict(model, x[, chans])
    rm(list=setdiff(ls(),c("pred","model")))
    return(list(pred=pred,model=model))
  }
}


#' Wrapper to linear model training. Defined separetely to avoid passing too many objects in parLapplyLB
#' @param x passed from fit_regressions
#' @param params passed from fit_regressions
#' @importFrom stats lm predict
#' @export
#' @return A list with two elements: predictions and a fitted model
#' @examples
#' fitter_linear()
fitter_linear <- function(x = NULL, params = NULL){
  if(!is.null(x) & !is.null(params)){
    w <- x[,"train_set"]==1
    fmla <- paste0(make.names(yvar),"~",polynomial_formula(variables=chans,degree=params$degree))
    model <- lm(formula=fmla,data=as.data.frame(x[w,c(chans,yvar)]))
    pred <- predict(model,as.data.frame(x[,chans]))
    rm(list=setdiff(ls(),c("pred","model")))
    model$model=NULL ## Trim down for slimmer objects
    model$qr$qr=NULL ## Trim down for slimmer objects
    return(list(pred=pred,model=model))
  }
}

#' Wrapper to glmnet. Defined separetely to avoid passing too many objects in parLapplyLB
#' @param x passed from fit_regressions
#' @param params passed from fit_regressions 
#' @importFrom stats as.formula predict
#' @importFrom utils getS3method
#' @export
#' @return A list with two elements: predictions and a fitted model
#' @examples
#' fitter_glmnet()
fitter_glmnet <- function(x = NULL, params = NULL){
  if(!requireNamespace("glmnetUtils", quietly = TRUE)){
    stop("Please run install.packages(\"glmnetUtils\")")
  }
  if(!is.null(x) & !is.null(params)){
    w <- x[,"train_set"] == 1
    fmla <- paste0(make.names(yvar), "~", polynomial_formula(variables = chans, degree = params$degree))
    flma <- as.formula(fmla)
    params <- params[setdiff(names(params), "degree")]
    params <- c(
      params,
      list(
        formula = fmla,
        data = as.data.frame(x[w, c(chans, yvar)]),
        use.model.frame = TRUE
      )
    )
    
    fun <- getS3method("cv.glmnet", "formula", envir = asNamespace("glmnetUtils"))
    model <- do.call(fun, params)
    model$call <- NULL ## Slimming down object
    model$glmnet.fit$call <- NULL ## Slimming down object
    attributes(model$terms)[[".Environment"]] <- NULL ## Slimming down object
    pred <- predict(model, as.data.frame(x[, chans]), s = model$lambda.min)
    
    rm(list = setdiff(ls(), c("pred", "model")))
    return(list(pred = pred, model = model))
  }
}

polynomial_formula <- function(variables,degree){
  n <- length(variables)
  polys <- lapply(
    seq_len(degree),
    function(deg){
      res <- apply(gtools::combinations(n,deg,variables,repeats.allowed=TRUE),1,paste,collapse="*")
      res <- paste0("I(",res,")")
    }
  )
  paste(do.call(c,lapply(polys,paste,sep="+")),collapse="+")
}

#' Wrapper to predict. Defined separetely to avoid passing too many objects in parLapplyLB
#' @param x passed from predict_from_models
#' @noRd
predict_wrapper <- function(x){
  if(is(x, "lm")){
    xp <- as.data.frame(xp)
  }
  if(is(x, "cv.glmnet")){
    requireNamespace("glmnetUtils")
    xp <- as.data.frame(xp)
    return(predict(x,xp,s=x$lambda.min)[,1])
  }
  if(is(x, "raw")){
    requireNamespace("keras")
    x = keras::unserialize_model(x)
  }
  if(is(x, "xgb.Booster")){
    requireNamespace("xgboost")
    x <- xgboost::xgb.Booster.complete(x)
    xgboost::xgb.parameters(x) <- list(nthread = 1)
  }
  if(is(x, "svm")){
    requireNamespace("e1071")
  }
  res <- predict(x, xp)
}

#' Wrapper to Neural Network training. Defined separetely to avoid passing too many objects in parLapplyLB
#' @param x passed from fit_regressions. Defines model architecture
#' @param params passed from fit_regressions
#' @importFrom generics fit
#' @importFrom stats predict
#' @export
#' @return A list with two elements: predictions and a fitted model
#' @examples
#' fitter_xgboost()
fitter_nn <- function(x,params){
  if(!requireNamespace("tensorflow", quietly = TRUE) & !requireNamespace("keras", quietly = TRUE)){
    stop("Please run install.packages(c(\"tensorflow\", \"keras\")) and make sure that install_tensorflow() and install_keras() have been run")
  }
  
  if(!is.null(x) & !is.null(params)){
    model <- keras::unserialize_model(params$object)
    params <- params[setdiff(names(params),"object")]
    
    early_stop <- keras::callback_early_stopping(monitor = "val_loss", patience = 20)
    callbacks <- list(early_stop)
    
    w <- x[,"train_set"]==1
    
    fit_history <- do.call(
      function(...){
        keras::fit(
          ...,
          object=model,
          x=x[w,chans],
          y=x[w,yvar],
          callbacks=callbacks
        )
      },
      params
    )
    
    pred <- predict(model, x[, chans])
    rm(list=setdiff(ls(),c("pred","model","fit_history")))
    return(list(pred = pred, model = keras::serialize_model(model), fit_history = fit_history))
  }
}

#' Train SVM regressions
#' @param yvar name of the exploratory measurement
#' @param params Named list of arguments passed to regression fitting functions
#' @param paths Character vector of paths to store intput, intermediary results, outputs...
#' @param cores Number of cores to use for parallel computing
#' @param xp Logicle-transformed backbone expression matrix
#' @param chans vector of backbone channels' names
#' @param events.code vector of length nrow(xp) specifying from which well each event originates
#' @param transforms_chan named list of logicle-transformations for backbone channels
#' @param transforms_pe named list of logicle-transformations for Infinity channels
#' @param regression_functions named list of fitter_* functions, passed from infinity_flow()
#' @param neural_networks_seed Seed for computational reproducibility when using neural networks. Passed from infinity_flow()
#' @param verbose Verbosity
#' @importFrom parallel makeCluster clusterExport clusterEvalQ stopCluster
#' @importFrom pbapply pblapply
#' @importFrom rhdf5 h5read
#' @noRd

fit_regressions2 <- function(
    yvar,
    params,
    paths,
    cores,
    verbose=TRUE,
    xp=h5read(file.path(paths["rds"],"xp_transformed_scaled.h5"), "/xp"),
    chans=readRDS(file.path(paths["rds"],"chans.Rds")),
    events.code=readRDS(file.path(paths["rds"],"pe.Rds")),
    transforms_chan=readRDS(file.path(paths["rds"],"transforms_chan.Rds")),
    transforms_pe=readRDS(file.path(paths["rds"],"transforms_pe.Rds")),
    regression_functions,
    neural_networks_seed
)

{   
  if(verbose){
    message("Fitting regression models")
  }
  library(parallel) #NEW
  chans <- make.names(chans)
  yvar <- make.names(yvar)
  # colnames(xp) <- make.names(colnames(xp))
  
  requireNamespace("parallel")
  
  cl <- makeCluster(min(cores,length(unique(events.code))))
  
  RNGkind("L'Ecuyer-CMRG")
  
  if(.Platform$OS.type == "windows") {
    mc.reset.stream <- function() return(invisible(NULL))
  } else {
    mc.reset.stream <- parallel::mc.reset.stream
  }
  
  mc.reset.stream()
  
  env <- environment()
  
  
  if(verbose){
    message("\tRandomly selecting 50% of the subsetted input files to fit models")
  }

  col_names <- c(
    "FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W",
    "CD69-CD301b", "Zombie", "MHCII", "CD4", "CD44", "CD8",
    "CD11c", "CD11b", "F480", "Ly6C", "Lineage", "CD45a488",
    "FJComp-PE(yg)-A", "CD24", "CD103", "Time"
  )
  
  colnames(xp) <- make.names(col_names)
  

  fcs_files <- unique(events.code) 
  
  
  # Fonction de regression (dans le fichier args)
  # sourcer le fichier args
  
  regression_functions=c(
    XGBoost = fitter_xgboost
  );
  extra_args_regression_params = list(
    list(nrounds = 500, eta=0.05)
  )
  

  ## Test avec fonction FUN 
  models <- list()  # Initialisation d'1 liste pour stocker les modèles
  timings <- numeric()
  for (i in seq_along(regression_functions)) {
    cat("\t\t", names(regression_functions)[i], "\n\n", sep="")
    t0 <- Sys.time()
    
    fun <- regression_functions[[i]]
    environment(fun) <- environment()  ## Fixing issue with scoping when cores = 1L
    
    ## MODIF
    models[[i]] <- list()  # Initialisation de la liste pour stocker les modèles
    
    for (evenement in fcs_files) {
      booleen <- events.code == evenement
      # booleen <- which(events.code == evenement)
      event_data <- xp[booleen, ]  # Filtrer les données pour l'événement actuel
      w <- sample(rep(c(TRUE, FALSE), times = c(floor(nrow(event_data)/2), nrow(event_data)-floor(nrow(event_data)/2))))
      event_data <- cbind(event_data, train_set = w)
      
      model <- fun(event_data, params = extra_args_regression_params[[1]])  # fonction de régression
      
      models[[i]][[evenement]] <- model  #stocker le modèle pour cet événement
      
    }
    
    # Enregistrement du temps 
    t1 <- Sys.time()
    dt <- difftime(t1, t0, units = "secs")
    cat("\t", dt, " seconds","\n", sep="")
    timings <- c(timings, dt)
    
  }
  
  # saveRDS(models,file=file.path(paths["rds"],"models_sha.Rds"))
  # model_sha = readRDS(file.path("/home/students","models_sha.Rds"))
  
  
  
  
  clusterExport(
    cl,
    c("yvar","chans","neural_networks_seed", "regression_functions", "fitter_nn"),
    envir=env
  )
  
  clusterEvalQ(
    cl,
    {
      chans <- make.names(chans)
      yvar <- make.names(yvar)

    }
  )
  
  if(verbose){
    message("\tFitting...")
  }
  
  
  names(models) <- names(regression_functions)
  stopCluster(cl)
  saveRDS(models,file=file.path(paths["rds"],"regression_models.Rds"))
  saveRDS(train_set,file=file.path(paths["rds"],"train_set.Rds"))
  list(timings=timings)
}

#' Predict missing measurements using fitted SVM regressions
#' @param paths Character vector of paths to store intput, intermediary results, outputs...
#' @param prediction_events_downsampling Number of events to predict data for per file
#' @param cores Number of cores to use for parallel computing
#' @param chans vector of backbone channels' names
#' @param events.code vector of length nrow(xp) specifying from which well each event originates
#' @param xp Logicle-transformed backbone expression matrix
#' @param models list of list of machine learning models, created by infinityFlow:::fit_regressions
#' @param train_set data.frame with nrow(xp) rows, specifying which events were used to train the models. Generated by infinityFlow:::fit_regressions
#' @param regression_functions named list of fitter_* functions, passed from infinity_flow()
#' @param neural_networks_seed Seed for computational reproducibility when using neural networks. Passed from infinity_flow()
#' @param verbose Verbosity
#' @noRd
predict_from_models <- function(
    paths,
    prediction_events_downsampling,
    cores,
    chans=readRDS(file.path(paths["rds"],"chans.Rds")),
    events.code=readRDS(file.path(paths["rds"],"pe.Rds")),
    xp=h5read(file.path(paths["rds"],"xp_transformed_scaled.h5"), "/xp"),
    models=readRDS(file.path(paths["rds"],"regression_models.Rds")),
    train_set=readRDS(file.path(paths["rds"],"train_set.Rds")),
    regression_functions = NULL,
    verbose=TRUE,
    neural_networks_seed
)
{   
  
  
  if(verbose){
    message("Imputing missing measurements")
  }
  
  xp_raw <- xp
  chans <- make.names(chans)
  
  col_names <- c(
    "FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W",
    "CD69-CD301b", "Zombie", "MHCII", "CD4", "CD44", "CD8",
    "CD11c", "CD11b", "F480", "Ly6C", "Lineage", "CD45a488",
    "FJComp-PE(yg)-A", "CD24", "CD103", "Time"
  )
  
  colnames(xp) <- make.names(col_names)
  
  # colnames(xp) <- make.names(colnames(xp))
  xp <- xp[, chans]
  
  requireNamespace("parallel")
  cl <- makeCluster(cores)
  
  if(.Platform$OS.type == "windows") {
    mc.reset.stream <- function() return(invisible(NULL))
  } else {
    mc.reset.stream <- parallel::mc.reset.stream
  }
  mc.reset.stream()
  
  env <- environment()
  
  if(verbose){
    message("\tRandomly drawing events to predict from the test set")
  }
  
  spl <- split(train_set[,1],events.code)
  spl <- lapply(
    spl,
    function(x){
      res <- rep(FALSE,length(x))
      w <- x==0
      res[w][sample(seq_len(sum(w)),min(prediction_events_downsampling,sum(w)))] <- TRUE
      res
    }
  )
  
  pred_set <- which(do.call(c,spl))
  
  xp <- xp[pred_set,]
  xp_raw <- xp_raw[pred_set, ]
  
  clusterExport(
    cl,
    c("xp","chans","neural_networks_seed", "regression_functions", "fitter_nn"),
    envir=env
  )
  invisible(clusterEvalQ(
    cl,
    {
      colnames(xp) <- make.names(colnames(xp))
      xp <- xp[,make.names(chans)]
    }
  ))
  
  for(i in seq_along(models)){
    models[[i]] <- lapply(models[[i]],"[[",2)
  }
  
  if(verbose){
    message("\tImputing...")
  }
  preds <- list()
  timings <- numeric()
  fun <- predict_wrapper
  environment(fun) <- environment() ## Fixing issue with scoping when cores = 1L
  for(i in seq_along(models)){
    cat("\t\t",names(models)[i],"\n\n",sep="")
    t0 <- Sys.time()
    preds[[i]] <- do.call(
      cbind,
      pblapply(
        cl=cl,
        X=models[[i]],
        FUN=fun
      )
    )
    t1 <- Sys.time()
    dt <- difftime(t1,t0,units="secs")
    cat("\t",dt," seconds","\n",sep="")
    timings <- c(timings,dt)
    colnames(preds[[i]]) <- names(models[[i]])
  }
  
  stopCluster(cl)
  
  if(verbose){
    message("\tConcatenating predictions")
  }
  preds <- lapply(preds,function(x){cbind(xp_raw,x)})
  preds <- lapply(preds,as.matrix)
  names(preds) <- names(models)
  
  if(verbose){
    message("\tWriting to disk")
  }
  saveRDS(preds,file=file.path(paths["rds"],"predictions.Rds"))
  saveRDS(pred_set,file=file.path(paths["rds"],"sampling_preds.Rds"))   
  
  
  list(timings=timings)
  
}