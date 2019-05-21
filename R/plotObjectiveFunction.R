#' Plot result of Objective Function
#'
#' This function plots the result of  \code{Lslide::ObjectiveFunction} function using
#' \code{\link[ggplot2]{ggplot}}
#'
#' @param x result of function \code{Lslide::ObjectiveFunction}
#' @param legend.position position of legend in plot. Default: "bottom"
#' @param title plot title. Default: ""
#' @param selected.Scale scale parameter for highlighting. Default: NULL
#' @param selected.Scale.label legend label. Default: "selected scale"
#' @param col.OF line color for objective function. Default: "#DC143C"
#' @param col.PF line color for plateau function. Default: "black"
#' @param col.nIV line color for normalized intrasegment variance. Default: "#6495ED"
#' @param col.nMI line color for normalized Moran's I. Default: "#8B4513"
#' @param col.sSl point color for selected scale parameter. Default: "gold"
#' @param selected.Scale.size scale size of selected point. Default: 4
#' @param selected.Scale.shape shape of selected point. Default: 1
#' @param selected.Scale.stroke width of outline of selected point. Default: 1.2
#' @param ... ggplot-theme options
#'
#' @return
#' \code{\link[ggplot2]{ggplot}} object
#'
#'
#'
#' @keywords ggplot2, Objective Function, plot
#'
#'
#' @export
plotObjectiveFunction <- function(x, legend.position = "bottom", title = "", selected.Scale = NULL, selected.Scale.label = "selected scale",
                                  col.OF = "#DC143C", col.PF = "black", col.nIV = "#6495ED", col.nMI = "#8B4513", col.sSl = "gold",
                                  selected.Scale.size = 4, selected.Scale.shape = 1, selected.Scale.stroke = 1.2, ...)
{
  # get values
  minVal <- min(x$Objective.Function)
  maxVal <- max(x$Objective.Function)
  # platVal <- x$Plateau[1]


  x.melt <- x[, c("Scale.Parameter", "Normalized.Intrasegment.Variance", "Normalized.Morans.I", "Objective.Function", "Plateau")]

  x.melt <- reshape2::melt(x[, c("Scale.Parameter", "Normalized.Intrasegment.Variance", "Normalized.Morans.I", "Objective.Function", "Plateau")],
                           id.vars = "Scale.Parameter", stringAsFactors = FALSE)

  labels <- c("normalized intrasegment variance   ", "normalized Moran's I   ", "objective function", "plateau function")

  # working
  plotOF <- ggplot2::ggplot() +
    geom_line(data = x.melt[grep("Obj|Pla", x.melt$variable),], aes(x = Scale.Parameter, y = value, color = variable, linetype = variable, size = variable)) +
    geom_line(data = x.melt[grep("Norm", x.melt$variable),], aes(x = Scale.Parameter, y = value * (maxVal - minVal) + minVal, color = variable, linetype = variable, size = variable)) +
    scale_linetype_manual(values = c(3, 3, 1, 2), name = "", labels = labels) +
    scale_colour_manual(values = c(col.nIV, col.nMI, col.OF, col.PF), name = "", labels = labels) +
    scale_size_manual(values = c(0.9, 0.9, 1, 1), name = "", labels = labels) +
    scale_y_continuous(sec.axis = sec_axis(trans = y ~ I(. - minVal)/(maxVal-minVal), name = "normalized parameters [indices]")) +
    scale_x_continuous(breaks = unique(x$Scale.Parameter)) +
    labs(y = "objective function [indices]", x = "segmentation scale", title = title) +
    theme_bw() +
    theme(legend.position = legend.position, legend.direction = "horizontal", legend.text = element_text(face = "bold", size = 11),
          axis.title.y = element_text(face = "bold", size = 11, margin = margin(r = 10)), axis.title.x = element_text(face = "bold", size = 11, margin = margin(t = 10)),
          axis.title.y.right = element_text(face = "bold", size = 11, margin = margin(l = 10)),
          legend.margin = margin(t = -5),
          title = element_text(face = "bold", size = 13, margin = margin(b = 10)), legend.key.width = unit(1.5, "line")) +
    guides(linetype = guide_legend(ncol = 2))

  if(!is.null(selected.Scale))
  {
    plotOF <- plotOF + geom_point(data =  x.melt[intersect(grep("Obj", x.melt$variable),   which(x.melt$Scale.Parameter %in% selected.Scale)),],
                        aes(x = Scale.Parameter, y = value, shape = "shape scale"), color = col.sSl, size = selected.Scale.size, stroke = selected.Scale.stroke) +
                        scale_shape_manual(name = "", labels = selected.Scale.label, values = c("shape scale" = selected.Scale.shape))
  }


  return(plotOF)

  # ggplot2::ggplot(data = x, aes(x = Scale.Parameter)) +
  #   geom_line(aes(y = Objective.Function, colour = "objective function", linetype = "objective function"), size = 1) +
  #   geom_line(aes(y = Normalized.Morans.I * (maxVal - minVal) + minVal, colour = "normalized Moran's I", linetype = "normalized Moran's I")) +
  #   geom_line(aes(y =  Normalized.Intrasegment.Variance * (maxVal - minVal) + minVal,  colour = "normalized intrasegment variance", linetype = "normalized intrasegment variance")) +
  #   scale_y_continuous(sec.axis = sec_axis(trans = y ~ I(. - minVal)/(maxVal-minVal), name = "normalized parameters [indices]")) +
  #   scale_x_continuous(breaks = unique(x$Scale.Parameter)) +
  #   geom_line(aes(y = Plateau, colour = "plateau function", linetype = "plateau function"), size = 1) +
  #   scale_linetype_manual(values = c("objective function" = 1, "normalized Moran's I" = 2, "normalized intrasegment variance" = 2, "plateau function" = 2)) +
  #   scale_colour_manual(values = c(col.nIV, col.nMI, col.OF, col.PF)) +
  #   labs(y = "objective function [indices]", x = "segmentation scale", colour = "", title = title) +
  #   theme(legend.position = legend.position)
  #


}


