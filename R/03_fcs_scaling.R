#' Scaling of the backbone measurements
#' @param yvar name of the exploratory measurement
#' @param paths Character vector of paths to store intput, intermediary results, outputs...
#' @param scale_function Scaling function to apply. Should apply to a matrix (events x backbone measurements) and return a matrix of similar size. Defaults to scale(,center=TRUE,scale=FALSE)
#' @param xp Logicle-transformed backbone expression matrix
#' @param chans vector of backbone channels' names
#' @param events.code vector of length nrow(xp) specifying from which well each event originates
#' @param verbose Verbosity
#' @noRd
standardize_backbone_data_across_wells <- function(
    yvar,
    paths,
    scale_function=function(x){scale(x,center=TRUE,scale=TRUE)},
    # xp=readRDS(file.path(paths["rds"],"xp_transformed.Rds")),
    xp=h5read(file.path(paths["rds"],"xp_transformed_6aout.h5"), "/xp"),
    chans=readRDS(file.path(paths["rds"],"chans.Rds")),
    events.code=readRDS(file.path(paths["rds"],"pe.Rds")),
    verbose=TRUE
){
  if(verbose){
    message("Harmonizing backbone data")
    message("\tScaling expression matrices")
  }
  
  # xp <- split(as.data.frame(xp),events.code)
  # xp <- lapply(xp,as.matrix)
  # xp <- lapply(xp,function(x){
  #   x[,chans] <- scale_function(x[,chans])
  #   x
  # })
  # xp <- do.call(rbind,xp)
  # 
  
  # colnames(xp) <- c(
  #   "FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W",
  #   "CD69-CD301b", "Zombie", "MHCII", "CD4", "CD44", "CD8",
  #   "CD11c", "CD11b", "F480", "Ly6C", "Lineage", "CD45a488",
  #   "FJComp-PE(yg)-A", "CD24", "CD103", "Time"
  # )
  # 
  # #stocke matrices transformées
  # result_list <- list()
  # 
  # unique_events <- unique(events.code)
  # 
  # for (event in unique_events) {
  #   subset_xp <- xp[events.code == event, ]
  # 
  #   subset_xp_matrix <- as.matrix(subset_xp)
  # 
  #   subset_xp_matrix[, chans] <- scale_function(subset_xp_matrix[, chans])
  #   result_list[[event]] <- subset_xp_matrix
  # }
  # 
  # xp <- do.call(rbind, result_list)
  # 
  # 
  # 
  # if(verbose){
  #   message("\tWriting to disk")
  # }
  # 
  # # saveRDS(xp,file=file.path(paths["rds"],"xp_transformed_scaled.Rds"))
  # xp_h5_path <- file.path(paths["rds"], "xp_transformed_scaled_5aout.h5")
  # h5write(xp, xp_h5_path, "/xp")
  # 
  # 
  # 
  # invisible()
  
  colnames(xp) <- c(
    "FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W",
    "CD69-CD301b", "Zombie", "MHCII", "CD4", "CD44", "CD8",
    "CD11c", "CD11b", "F480", "Ly6C", "Lineage", "CD45a488",
    "FJComp-PE(yg)-A", "CD24", "CD103", "Time"
  )
  
  xp_h5_path <- file.path(paths["rds"], "xp_transformed_scaled.h5")
  h5createFile(xp_h5_path)
  h5createDataset(xp_h5_path, "xp", dims=dim(xp), storage.mode="double")
  
  unique_events <- unique(events.code)
  
  for (event in unique_events) {
    subset_xp <- xp[events.code == event, ] ## subset of the matriix
    subset_xp_matrix <- as.matrix(subset_xp)
    subset_xp_matrix[, chans] <- scale_function(subset_xp_matrix[, chans])
    
    ## Sequential
    start_idx <- which(events.code == event)[1]
    end_idx <- start_idx + nrow(subset_xp_matrix) - 1
    h5write(subset_xp_matrix, xp_h5_path, "xp", index=list(start_idx:end_idx, 1:ncol(xp)))
  }
  
  if (verbose) {
    message("\t Finish")
  }
  
  invisible()
}



