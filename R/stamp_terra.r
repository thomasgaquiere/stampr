# ---- roxygen documentation ----
#
#' @title Spatial temporal analysis of moving polygons (terra version)
#'@importFrom tidyterra rename
#'@importFrom sf st_as_sf st_distance st_area
#'@importFrom dplyr select
#'@import terra
#'@import parallel
#'@import doParallel
#'@import foreach
#'@import spdep
#' @description
#' This function generates a \code{SpatVector} that can be used for spatial temporal analysis of moving polygons
#' as described in the paper Robertson et al. (2007).
#'
#' @details
#'  The \code{stamp_terra} function can be used to perform spatial temporal analysis of moving polygons (STAMP)
#'  as outlined in the paper by Robertson et al., (2007). Polygon movement "groups" are delineated based on
#'  polygon connectedness defined by the distance threshold \code{dc}. That is, if polygon
#'  boundaries (in T1 or T2) are within distance \code{dc} of one another they will be designated
#'  to the same group. STAMP events are reported at four levels of increasing complexity: \cr
#'  LEV1 -- disappearance (\code{DISA}), stable (\code{STBL}), and generation (\code{GENA}); \cr
#'  LEV2 -- disappearance (\code{DISA}), contraction (\code{CONT}), stable (\code{STBL}),
#'          expansion (\code{EXPN}), and generation (\code{GENR}); \cr
#'  LEV3 -- disappearance (\code{DISA}), T1 displacement (\code{DISP1}), convergence (\code{CONV}),
#'          concentration (\code{CONC}), contraction (\code{CONT}), stable (\code{STBL}),
#'          expansion (\code{EXP}), fragmentation (\code{FRAG}), divergence (\code{DIV}),
#'          T2 displacement (\code{DISP2}), and generation (\code{GENR}); \cr
#'  LEV4 -- LEV4 is different from other levels. It is used to identify those groups where
#'          union (\code{UNION}), division (\code{DIVISION}), and both union and division
#'          (\code{BOTH}) events occur. These events occur when there are more than one
#'          stable event in a group. Groups with one or no stable events recieve an \code{NA}
#'          value for LEV4. \cr
#'  See Robertson et al. (2007; especially Figure 1) for complete descriptions of all STAMP movement
#'  event types.
#'  
#'  Note also that there must be a globally unique \code{ID} column in each data frame passed to the stamp function
#'
#' @param T1 a \code{SpatVector} object of polygons from time 1.
#' @param T2 a \code{SpatVector} object of polygons from time 2.
#' @param dc spatial distance threshold for determining groupings (see \bold{Details}) in appropriate units.
#' @param direction logical, whether or not to perform directional analysis. See documentation for
#'    \code{stamp.direction} for further details.
#' @param distance logical, whether or not to perform distance analysis. See documentation for 
#'    \code{stamp.distance} for further details.
#' @param ... additional parameters to be passed to functions if \code{direction}, or \code{distance} are
#'    set to \code{TRUE}. 
#'
#'
#' @return
#'  This function returns a \code{SpatVector} with the following data columns:
#'  \item{ID1}{Polygon ID from T1 polygons; \code{NA} if it did not exist,}
#'  \item{ID2}{Polygon ID from T2 polygons; \code{NA} if it did not exist,}
#'  \item{LEV1}{Level 1 STAMP designation,}
#'  \item{LEV2}{Level 2 STAMP designation,}
#'  \item{LEV3}{Level 3 STAMP designation,}
#'  \item{LEV4}{Level 4 STAMP designation,}
#'  \item{GROUP}{Group ID signifying group membership,}
#'  \item{AREA}{Polygon area in appropriate areal units,}
#'  \item{--}{(optional) Additional columns from directional analysis if \code{direction = TRUE},}
#'  \item{--}{(optional) Additional columns from distance analysis if \code{distance = TRUE},}
#  \item{--}{(optional) Additional columns from shape analysis if \code{shape = TRUE}.}
#'
#' @references
#'  Robertson, C., Nelson, T., Boots, B., and Wulder, M. (2007) STAMP: Spatial-temporal analysis of moving polygons.
#'  \emph{Journal of Geographical Systems}, 9:207-227.
#'
#' @keywords stamp_terra
#' @seealso stamp.direction stamp.distance stamp.shape stamp.map stamp.group.summary
#' @export
#
# ---- End of Documentation ----

stamp_terra <- function(T1, T2, dc=0, direction=FALSE, distance=FALSE,cores=1, ...){ 
  # intersection b/w T1 and T2
  if (!'ID' %in% names(T1))
    stop("Need a unique 'ID' column.")
  
  if (!'ID' %in% names(T2))
    stop("Need a unique 'ID' column.")
  
  T1 <- rename(T1, 'ID1' = 'ID')
  T2 <- rename(T2, 'ID2'= 'ID')
  
  pI <- terra::intersect(T1, T2)
  
  
  if (!is.null(pI)) {
    pI$LEV1 <- "STBL"
  } 
  
  
  #T1 difference
  #res <- list()
  dfD1 <- data.frame(ID1 = rep(NA,length(T1)),ID2 = rep(NA,length(T1)))
  #This is slow, can we improve?
  
  res <- terra::erase(T1,T2)
  
  
  
  for (i in seq(along=res)) {
    if (!is.null(res[i])){                                          
      dfD1[i,1] <- res[i]$ID1
    }}
  
  
  #Get rid of problem scenarios
  ind <- which(is.na(dfD1$ID1) & is.na(dfD1$ID2))
  if (length(ind) > 0){
    dfD1 <- dfD1[-ind,]
  }
  
  res1 <- res[!sapply(res, is.null)]
  pD1 <- NULL
  
  if (length(res1)>0){
    
    pD1 <- res1
    pD1$LEV1 <- "DISA"
    
  }
  
  #T2 difference
  
  dfD2 <- data.frame(ID1 = rep(NA,length(T2)),ID2 = rep(NA,length(T2)))
  #This is slow, can we improve?
  
  res <- terra::erase(T2, T1)
  
  for (i in seq(along=res)) {
    if (!is.null(res[i])){                                          
      
      dfD2[i,2] <- res[i,]$ID2
    }
    
  }
  
  #Get rid of problem scenarios
  ind <- which(is.na(dfD2$ID1) & is.na(dfD2$ID2))
  if (length(ind) > 0){
    dfD2 <- dfD2[-ind,]
  }
  
  
  res1 <- res[!sapply(res, is.null)]
  pD2 <- NULL
  #if (!is.null(res1)){
  if (length(res1) > 0){
    #ouT1 <- do.call("rbind", res1)
    pD2 <- res1
    pD2$LEV1 <- "GENR"
    #row.names(pD2) <- paste("GENR",seq(1:length(pD2)),sep="")
  }
  
  #Piece them together
  stmp <- rbind(pD1,pI,pD2)
  
  stmp <- st_as_sf(stmp)
  #stmp <- as.data.frame(stmp)
  
  #assign event types ---
  stmp$LEV2 <- stmp$LEV1
  
  #get contraction events
  id.stab1 <- unique(stmp$ID1[which(stmp$LEV1 == "STBL")])
  stmp$LEV2[which(stmp$LEV1 == "DISA" & stmp$ID1 %in% id.stab1)] <- "CONT"
  
  #get expansion events
  id.stab2 <- unique(stmp$ID2[which(stmp$LEV1 == "STBL")])
  stmp$LEV2[which(stmp$LEV1 == "GENR" & stmp$ID2 %in% id.stab2)] <- "EXPN"
  
  
  #Delineate contiguous bases for groups
  
  #stmp <- as_Spatial(stmp)
  stmp$TMP <- 1
  if(nrow(stmp) > 1) {
    nbl <- poly2nb(stmp)
    for(i in 1:nrow(stmp)) {
      nbl[[i]] <- c(unlist(nbl[i]), i)
    }
    stmp$TMP <- n.comp.nb(nbl)$comp.id
  }
  #stmp <- st_as_sf(stmp)
  #stmp <- terra::vect(stmp)
  #terra::crs(stmp) <- terra::crs(T1)
  #Label all other LEV2 movement types...
  gdInd <- which(stmp$LEV2 == "GENR" | stmp$LEV2 == "DISA")
  tempLev <- stmp$LEV2
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = length(gdInd), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  #stmpL <- nrow(stmp)
  
  for(i in gdInd) {
    
    
    
    evenT1 <- stmp$LEV2[i]
    #find D of all appropriate polys
    #dists <- vector(length=length(stmp), mode="numeric")
    #dists[] <- NA
    #stmp <- sf::st_as_sf(stmp)
    #stmp <- wrap(stmp)
    
    dists <-  foreach(j = 1:nrow(stmp), .packages='sf') %dopar% {
      
      #stmp <- unwrap(stmp)
      
      #Do not include nearest GEN-GEN or DIS-DIS as they do not change names
      if (stmp$LEV2[i] != stmp$LEV2[j]){ st_distance(stmp[j,], stmp[i,])}
    }
    #stmp <- unwrap(stmp)
    #stmp <- sf::as_Spatial(stmp)
    dists[sapply(dists, is.null)] <- NA
    dists <- unlist(dists)
    #sort by D then extract if below dc value
    if (min(dists,na.rm=T) <= dc){
      minInd <- which(dists == min(dists, na.rm=T))[1]
      evenT2 <- stmp$LEV2[minInd]
      if (evenT1 == "DISA"){
        tempLev[i] <- switch(evenT2,
                             GENR = "DISP1",
                             EXPN = "CONV",
                             CONT = "CONC",
                             STBL = "CONC")
      }
      else {
        tempLev[i] <- switch(evenT2,
                             DISA = "DISP2",
                             EXPN = "FRAG",
                             CONT = "DIVR",
                             STBL = "CONC")
      }
      #Group movement event into appropriate contiguous group
      stmp$TMP[i] <- stmp$TMP[minInd]
    }
    setTxtProgressBar(pb, i)
  }
  stopCluster(cl)
  
  
  stmp$LEV3 <- tempLev
  #Rename groups so there are no gaps
  grps <- unique(stmp$TMP)
  for (i in 1:length(grps)){
    stmp$TMP[which(stmp$TMP == grps[i])] <- i
  }
  #Label Groups with Multi-Stable events as union or division
  stmp$LEV4 <- 'N/A'
  for (grp in unique(stmp$TMP)){
    ind <- which(stmp$TMP == grp & stmp$LEV3 == "STBL")
    ind.grp <- which(stmp$TMP == grp)
    if (length(ind) > 1){
      if (length(unique(stmp$ID2[ind])) == 1){stmp$LEV4[ind.grp] <- "UNION"}
      else if (length(unique(stmp$ID1[ind])) == 1){stmp$LEV4[ind.grp] <- "DIVISION"}
      else {stmp$LEV4[ind.grp] <- "BOTH"}
    }
  }
  
  
  #Delete TMP column and make a GROUP column
  stmp$GROUP <- stmp$TMP
  #stmp[stmp$TMP]
  #stmp@data <- stmp@data[,-5]
  #sort by group column
  #stmp <- stmp[order(stmp$GROUP),]
  #rename FID's
  stmp <- select(stmp, 'x', 'y', 'ID1', 'ID2', 'LEV1', 'LEV2', 'LEV3', 'LEV4', 'GROUP')
  #stmp <- spChFIDs(stmp,as.character(seq(0,(length(stmp)-1))))
  #Create a polygon area column
  #stmp$AREA <- gArea(stmp,byid=TRUE)
  stmp$AREA <- st_area(stmp)
  
  
  # others functions not modified yet
  
  #directional analysis
  #if (direction==TRUE){stmp <- stamp.direction(stmp,...)}
  #distance analysis
  #if (distance==TRUE){stmp <- stamp.distance(stmp,...)}
  
  #output
  return(stmp)
}
#-------------- END of stamp ---------------------------------------------------