#' Sky-view factor as a relief visualization technique
#'
#' This function enables the possibility to directly address the
#' stand-alone Relief Visualization Toolbox of ZAKSEK et al. (2011).
#'
#' @param path.software path where software is located
#' @param path.input path to input file
#' @param quiet option to control outputs in console during process. Default: TRUE
#' @param other... see \href{http://iaps.zrc-sazu.si/sites/default/files/rvt_1.3_0.pdf}{Relief Visualization Toolbox} \emph{(last call: 13-04-2017)}
#'
#' @note
#' \itemize{
#'   \item output is located in input directory
#'   \item compatible to Relief Visualization Toolbox version 1.3, \href{http://iaps.zrc-sazu.si/en/rvt#v}{RVT 1.3} \emph{(last call: 13-04-2017)}
#'   \item ZAKSEK, K., K. OSTIR & Z. KOKALJ (2011): Sky-View Factor as a Relief Visualization Technique. - Remote Sensing 3, 12, 398-415.
#' }
#'
#' @keywords sky-view factor
#'
#' @examples
#' skyViewFactor()
#'
#' @export
skyViewFactor <- function(path.software, path.input, options_overwrite = 1, options_exaggeration_factor = 1.0,
                          hillshading = c(0, 315, 35, 0), multiple_hillshading = c(0, 16),
                          pca_hillshading = c(0, 3), slope_gradient = 0, simple_local_relief = c(0, 20),
                          sky_view_factor = c(1, 16, 10, 0, "low"), anisotropic_svf = c(0, "low", 315),
                          positive_openness = 0, negative_openness = 0, sky_illumination = c(0, "overcast", 250, 100),
                          local_dominance = c(0, 10, 20), quiet = TRUE)
{

  # set current working directory
  current_wd <- getwd()

  # # #
  # check add-on files for automatic computation -------------------------------------------------------------

  if(!file.exists(paste0(path.software, "/", "sendKeys.bat")))
  {
    stop("sendKeys.bat not found")
  }


  startBatFile <- paste0(path.software, "/", "RVT_1.3_AutomaticStart.bat")

  if(!file.exists(startBatFile))
  {
    stop("RVT_1.3_AutomaticStart.bat not found")
  }

  # # #
  # create process_files.txt -------------------------------------------------------------
  write(path.input, paste0(path.software, "/settings/process_files.txt"))


  # # #
  # create text for processing -------------------------------------------------------------

  header <- c("#----------------------------------------------------------------------------------------------------------",
              "#----------------------------------------------------------------------------------------------------------",
              "#",
              "#	Relief Visualization Toolbox default settings",
              "#",
              "#	For more information about available settings in this file see RVT manual.",
              "#",
              "#----------------------------------------------------------------------------------------------------------",
              "#----------------------------------------------------------------------------------------------------------")

  txt_options <- c(paste("overwrite", options_overwrite, sep = " = "), "", paste("exaggeration_factor", as.numeric(options_exaggeration_factor), sep = " = "))

  txt_hillshading <- c(paste("hillshading", hillshading[1], sep = " = "),
                       paste("sun_azimuth", hillshading[2], sep = " = "),
                       paste("sun_elevation", hillshading[3], sep = " = "),
                       paste("shadow_modelling", hillshading[4], sep = " = "))


  txt_multiple_hillshading <- c(paste("multiple_hillshading", multiple_hillshading[1], sep = " = "),
                         paste("hillshade_directions", multiple_hillshading[2], sep = " = "))



  txt_pca_hillshading <- c(paste("pca_hillshading", pca_hillshading[1], sep = " = "),
                       paste("number_components", pca_hillshading[2], sep = " = "))



  txt_slope_gradient <- c(paste("slope_gradient", slope_gradient, sep = " = "))



  txt_simple_local_relief <- c(paste("simple_local_relief", simple_local_relief [1], sep = " = "),
                           paste("trend_radius", simple_local_relief [2], sep = " = "))



  txt_sky_view_factor <- c(paste("sky_view_factor", sky_view_factor[1], sep = " = "),
                           paste("svf_directions", sky_view_factor[2], sep = " = "),
                           paste("search_radius", sky_view_factor[3], sep = " = "),
                           paste("remove_noise", sky_view_factor[4], sep = " = "),
                           paste("noise_removal", sky_view_factor[5], sep = " = "))




  txt_anisotropic_svf <- c(paste("anisotropic_svf", anisotropic_svf [1], sep = " = "),
                           paste("anisotropy_level", anisotropic_svf [2], sep = " = "),
                           paste("anisotropy_direction", anisotropic_svf [3], sep = " = "))


  txt_positive_openness <- c(paste("positive_openness", positive_openness, sep = " = "))


  txt_negative_openness <- c(paste("negative_openness", negative_openness, sep = " = "))


  txt_sky_illumination <- c(paste("sky_illumination", sky_illumination[1], sep = " = "),
                       paste("sky_model", sky_illumination[2], sep = " = "),
                       paste("number_points", sky_illumination[3], sep = " = "),
                       paste("max_shadow_dist", sky_illumination[4], sep = " = "))



  txt_local_dominance <- c(paste("local_dominance", local_dominance[1], sep = " = "),
                            paste("min_radius", local_dominance[2], sep = " = "),
                            paste("max_radius", local_dominance[3], sep = " = "))



  # # #
  # write text for processing -------------------------------------------------------------
  fileConn <- file(paste0(path.software, "/settings/default_settings.txt"))

  writeLines(c(header, "", txt_options, "", txt_hillshading, "", txt_multiple_hillshading, "", txt_pca_hillshading, "",
               txt_slope_gradient, "", txt_simple_local_relief, "", txt_sky_view_factor, "", txt_anisotropic_svf, "",
               txt_positive_openness, "", txt_negative_openness, "", txt_sky_illumination, "", txt_local_dominance), fileConn)

  close(fileConn)

  # # #
  # start processing -------------------------------------------------------------
  setwd(path.software)


  if(quiet == FALSE)
  {
    system(command = 'RVT_1.3_AutomaticStart.bat > nul 2>&1')
  } else
  {
    system(command = 'RVT_1.3_AutomaticStart.bat')
  }


  # reset working directory
  setwd(current_wd)

}
