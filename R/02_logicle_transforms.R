#' Transform FCS files
#' @description Backbone measurements use a common transformation across the whole dataset. Exploratory measurements are transformed well/file-wise.
#' @param yvar name of the exploratory measurement
#' @param paths Character vector of paths to store intput, intermediary results, outputs...
#' @param xp Logicle-transformed backbone expression matrix
#' @param chans vector of backbone channels' names
#' @param events.code vector of length nrow(xp) specifying from which well each event originates
#' @param annot annotation table generated by infinityFlow:::intialize()
#' @param verbose Verbosity
#' @importFrom utils read.table
#' @importFrom stats quantile
#' @importFrom flowCore logicleTransform inverseLogicleTransform
#' @noRd

library(rhdf5)


logicle_transform_input <- function(
    yvar,
    paths,
    # xp=readRDS(file.path(paths["rds"],"xp.Rds")),
    xp = h5read(file.path(paths["rds"],"xp.h5"),"/matrix"),
    chans=readRDS(file.path(paths["rds"],"chans.Rds")),
    events.code=readRDS(file.path(paths["rds"],"pe.Rds")),
    annot=read.table(paths["annotation"],sep=",",header=TRUE,stringsAsFactors=FALSE),
    verbose=TRUE
){
  
  ## ##################
  ## Computing parameters for each channel for each project using the code from flowCore's estimateLogicle
  ## ##################
  if(verbose){
    message("Logicle-transforming the data")
    message("\tBackbone data")
  }
  
  
  colnames(xp) <- c(
    "FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W",
    "CD69-CD301b", "Zombie", "MHCII", "CD4", "CD44", "CD8",
    "CD11c", "CD11b", "F480", "Ly6C", "Lineage", "CD45a488",
    "FJComp-PE(yg)-A", "CD24", "CD103", "Time"
  )
  
  transforms_chan <- vector("list", length(chans))
  names(transforms_chan) <- chans
  
  for (i in seq_along(chans)) {
    x <- chans[i]
    data <- xp[, x]
    t <- max(data)
    m <- 4.5
    q <- 0.05
    r <- .Machine$double.eps + quantile(data, q)
    w <- max((m - log10(t / abs(r))) / 2, 0.1)
    w <- min(w, m / 2)
    a <- 0
    transforms_chan[[x]] <- logicleTransform(w = w, t = t, m = m, a = a)
  }
  
  saveRDS(transforms_chan,file=file.path(paths["rds"],"transforms_chan.Rds"))
  
  
  # transforms_chan <- setNames(
  #   lapply(
  #     chans,
  #     function(x){
  #       data <- xp[,x]
  #       t <- max(data)
  #       m <- 4.5
  #       q <- 0.05
  #       r <- .Machine$double.eps + quantile(data, q)
  #       w <- max((m-log10(t/abs(r)))/2,0.1)
  #       w <- min(w,m/2)
  #       a <- 0
  #       logicleTransform(w=w,t=t,m=m,a=a) ##Just use summary() to retrive the parameters
  #     }
  #   ),
  #   chans
  # )
  
  if(verbose){
    message("\tExploratory data")
  }
  
  unique_events <- unique(events.code)
  transforms_pe <- setNames(vector("list", length(unique_events)), unique_events)
  
  for (event in unique_events) {
    data <- xp[events.code == event, yvar]
    t <- max(data)
    m <- 4.5
    q <- 0.05
    r <- .Machine$double.eps + quantile(data, q)
    w <- max((m - log10(t / abs(r))) / 2, 0.1)
    w <- min(w, m / 2)
    a <- 0
    transforms_pe[[event]] <- logicleTransform(w = w, t = t, m = m, a = a)
  }
  
  saveRDS(transforms_pe,file=file.path(paths["rds"],"transforms_pe.Rds"))
  
  
  if(verbose){
    message("\tWriting to disk")
  }

  ## ##################
  ## Exporting transformed expression matrices
  ## ##################
  if(verbose){
    message("\tTransforming expression matrix")
  }
  
  xp_h5_path <- file.path(paths["rds"], "xp_transformed_scaled.h5")
  h5createFile(xp_h5_path)
  h5createDataset(xp_h5_path, "xp", dims=dim(xp), storage.mode="double")
  
  for (chan in chans) {
    xp[, chan] <- transforms_chan[[chan]](xp[, chan])
  }
  
  unique_events <- unique(events.code)
  
  for (event in unique_events) {
    subset_xp <- xp[events.code == event, ]
    subset_xp_matrix <- as.matrix(subset_xp)
    
    subset_xp_matrix[, yvar] <- transforms_pe[[event]](subset_xp_matrix[, yvar])
    
    start_idx <- which(events.code == event)[1]
    end_idx <- start_idx + nrow(subset_xp_matrix) - 1
    
    h5write(subset_xp_matrix, xp_h5_path, "xp", index=list(start_idx:end_idx, 1:ncol(xp)))
  }
  
  
  
  
  if(verbose){
    message("\tWriting to disk")
  }
  
  # saveRDS(xp,file=file.path(paths["rds"],"xp_transformed.Rds"))
  # xp_h5_path <- file.path(paths["rds"], "xp_transformed.h5")
  # h5write(xp, xp_h5_path, "/xp")
  
  invisible()
}


