# ---- roxygen documentation ----
#
#' @title Mapping (plotting) functionality for \code{stamp} output
#'
#' @description
#'  This function maps STAMP output for visual assessment of STAMP events and groupings.
#'  Choice of which aspect of the stamp output to be visualized is controlled by passing
#'  the column name to the \code{stamp.map} function.
#'
#' @details
#'  The \code{stamp.map} function can be used to visualize any of the stamp event designation levels
#'  (e.g., \code{"LEV1"}, \code{"LEV2"}, \code{"LEV3"}, \code{"LEV4"}, or the STAMP groupings
#'  (based off of parameter \code{dc} in the \code{stamp} function).
#'
#' @param stmp output from the \code{stamp} function, i.e., a (\code{SpatialPolygonsDataFrame}).
#' @param by tells the function which attribute to visualize, one of \code{"LEV1"},
#'           \code{"LEV2"}, \code{"LEV3"}, \code{"LEV4"}, or \code{"GROUP"}
#'
#' @return
#'   \code{stamp.map} returns a map of the \code{stamp} output using the \code{spplot} functionality. 
#'   It implements a pre-defined coloring scheme.
#'
#' @keywords stamp plot
#' @seealso stamp
#' @export
#
# ---- End of roxygen documentation ----

stamp.map <- function(stmp,by="LEV1"){
  if (by %in% colnames(stmp@data)){
    if (by == "LEV1"){
      stmp@data[,by] <- factor(stmp@data[,by],levels=c("DISA","STBL","GENR"))
      cols <- c("red","gray25","blue")
      }
    else if (by == "LEV2"){
      stmp@data[,by] <- factor(stmp@data[,by],levels=c("DISA","CONT","STBL","EXPN","GENR"))
      cols <- c("red","darkred","gray25","royalblue4","blue")
      }
    else if (by == "LEV3"){
      stmp@data[,by] <- factor(stmp@data[,by],levels=c("DISA","DISP1","CONV","CONC","CONT","STBL","EXPN","FRAG","DIVR","DISP2","GENR"))
      cols <- c("red","indianred1","darkorange","violetred","darkred","gray25","royalblue4","deepskyblue4","dodgerblue","skyblue","blue")
      }
    else if (by == "LEV4"){
      stmp@data[,by] <- factor(stmp@data[,by],levels=c(NA,"UNION","DIVISION","BOTH"))
      cols <- c("grey","red","blue","purple")
      }
    else {
      stmp@data[,by] <- factor(stmp@data[,by])
      cols <- rainbow(length(unique(stmp@data[,by])))
      }
    spplot(stmp,zcol=by,col.regions=cols)
    }
  else {stop(paste("The column ",by," does not exist.",sep=""))}
  }
#------ end of stamp.map function ---------------------------------------------