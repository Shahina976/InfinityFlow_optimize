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
                                                xp=h5read(file.path(paths["rds"],"xp_transformed.h5"), "/xp"),
                                                chans=readRDS(file.path(paths["rds"],"chans.Rds")),
                                                events.code=readRDS(file.path(paths["rds"],"pe.Rds")),
                                                verbose=TRUE
                                                ){
    if(verbose){
        message("Harmonizing backbone data")
        message("\tScaling expression matrices")
    }
    
    colnames(xp) <- c(
      "FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W",
      "CD69-CD301b", "Zombie", "MHCII", "CD4", "CD44", "CD8",
      "CD11c", "CD11b", "F480", "Ly6C", "Lineage", "CD45a488",
      "FJComp-PE(yg)-A", "CD24", "CD103", "Time"
    )
    
    xp <- split(as.data.frame(xp),events.code)
    xp <- lapply(xp,as.matrix)
    xp <- lapply(xp,function(x){
        x[,chans] <- scale_function(x[,chans])
        x
    })
    xp <- do.call(rbind,xp)

    if(verbose){
        message("\tWriting to disk")
    }    
    
    # saveRDS(xp,file=file.path(paths["rds"],"xp_transformed_scaled.Rds"))
    xp_h5_path <- file.path(paths["rds"], "xp_transformed_scaled.h5")
    h5write(xp, xp_h5_path, "/xp")
    
    invisible()
}
