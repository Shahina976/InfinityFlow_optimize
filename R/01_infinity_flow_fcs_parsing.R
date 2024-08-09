#' Parsing FCS files
#' @param input_events_downsampling See \code{\link{infinity_flow}}
#' @param paths Character vector of paths to store input, intermediary results, outputs...
#' @param extra_args_read_FCS passed to flowCore:read.FCS
#' @param name_of_PE_parameter Name of the exploratory measurement
#' @description This function reads the input FCS files, downsamples to a user-defined number of events if necessary, harmonizes the input data names, saves a concatenated expression matrix and the corresponding vector mapping events to their file of origin.
#' @param verbose Verbosity
#' @noRd
#' @importFrom flowCore read.FCS
#' @importFrom Biobase pData exprs
#' @importFrom rhdf5 h5createFile h5createDataset h5write h5writeAttribute


subsample_data <- function(
    input_events_downsampling,
    paths,
    extra_args_read_FCS,
    name_of_PE_parameter,
    verbose=TRUE
){
  ## Subsampling
  if(verbose){
    message("Parsing and subsampling input data")
    message("\tDownsampling to ", input_events_downsampling, " events per input file")
  }
  files <- list.files(paths["input"], full.names=TRUE, recursive=TRUE, pattern=".fcs")
  invisible(
    lapply(
      files,
      function(file){
        res <- do.call(read.FCS, c(list(filename=file), extra_args_read_FCS))
        w <- sort(sample(seq_len(nrow(res)), min(input_events_downsampling, nrow(res))))
        res <- res[w,]
        write.FCS(res, file.path(paths["subset"], basename(file)))
      }
    )
  )
  
  ## Convert 
  if(verbose){
    message("\tConcatenating expression matrices")
  }
  files <- list.files(paths["subset"], full.names=TRUE, recursive=FALSE, include.dirs=FALSE, pattern=".fcs")
  ns <- setNames(integer(length(files)), files)
  xp <- lapply(
    files,
    function(file){
      xp <- do.call(read.FCS, c(list(filename=file), extra_args_read_FCS))
      annot <- pData(xp@parameters)
      xp <- exprs(xp)
      targets <- annot$desc
      targets[is.na(targets)] <- annot$name[is.na(targets)]
      colnames(xp)[colnames(xp) != name_of_PE_parameter] <- targets[colnames(xp) != name_of_PE_parameter]
      colnames(xp)[colnames(xp) == name_of_PE_parameter] <- name_of_PE_parameter
      xp
    }
  )
  names(xp) <- files
  ns <- vapply(xp, nrow, 1L)
  xp <- do.call(rbind, xp)
  
  ## Map which events originate from which file.
  if(verbose){
    message("\tWriting to disk")
  }
  events.code <- unlist(
    lapply(
      names(ns),
      function(x){
        rep(tail(strsplit(x, "/")[[1]], 1), ns[x])
      }
    )
  )
  
  ## Save to HDF5 format
  library(rhdf5)
  h5file <- file.path(paths["rds"], "xp.h5")
  if (file.exists(h5file)) {
    file.remove(h5file)
  }
  h5createFile(h5file)
  h5createDataset(h5file, "data", dims=dim(xp), storage.mode="numeric")
  h5write(xp, h5file, "data")
  h5writeAttribute("Expression matrix", h5file, "data", "description")
  
  saveRDS(events.code, file=file.path(paths["rds"], "pe.Rds"))
  invisible()
}







## Script bis pour optimiser 
# library(flowCore)
# library(Biobase)
# library(rhdf5)
# 
# subsample_data <- function(
    #     input_events_downsampling,
#     paths,
#     extra_args_read_FCS,
#     name_of_PE_parameter,
#     verbose = TRUE
# ){
#   if(verbose){
#     message("Parsing input data")
#   }
#   
#   files <- list.files(paths["input"], full.names = TRUE, recursive = TRUE, pattern = "\\.fcs$")
#   xp_h5_path <- file.path(paths["rds"], "xp.h5")
#   
#   # Remove existing file if it exists
#   if (file.exists(xp_h5_path)) {
#     file.remove(xp_h5_path)
#   }
#   
#   # Initialize HDF5 file
#   h5createFile(xp_h5_path)
#   first_file <- TRUE
#   total_rows <- 0
#   col_names <- NULL
#   
#   for (file in files) {
#     res <- do.call(read.FCS, c(list(filename = file), extra_args_read_FCS))
#     xp <- exprs(res)
#     annot <- pData(res@parameters)
#     targets <- annot$desc
#     targets[is.na(targets)] <- annot$name[is.na(targets)]
#     
#     # Renommer les colonnes
#     colnames(xp) <- targets
#     colnames(xp)[colnames(xp) == name_of_PE_parameter] <- name_of_PE_parameter
#     
#     if (first_file) {
#       col_names <- colnames(xp)
#       h5createDataset(xp_h5_path, "matrix", dims = dim(xp), maxdims = c(Inf, length(col_names)), storage.mode = "double")
#       h5write(xp, xp_h5_path, "matrix")
#       first_file <- FALSE
#     } else {
#       current_dims <- dim(h5read(xp_h5_path, "matrix"))
#       new_dims <- c(current_dims[1] + nrow(xp), current_dims[2])
#       h5set_extent(xp_h5_path, "matrix", new_dims)
#       h5write(xp, xp_h5_path, "matrix", start = c(current_dims[1] + 1, 1), count = dim(xp))
#     }
#     
#     total_rows <- total_rows + nrow(xp)
#   }
#   
#   if (verbose) {
#     message("\tWriting completed. Total rows: ", total_rows)
#   }
#   saveRDS(events.code,file=file.path(paths["rds"],"pe.Rds"))

#   
#   invisible()
# }


