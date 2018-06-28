pacman::p_load(sf)

path.slope <- "F:/Re-Thesis/Data/Input/slope_test.tif"
slp <- raster::raster(path.slope)
env.rsaga <- RSAGA::rsaga.env(path = "C:/OSGeo4W64/apps/saga-6.3.0")

sf.inventory <- sf::st_read(dsn = "F:/Re-Thesis/Data/Input/Landslides/OP14_AOI_landslides_parts.shp")
sf.scarp <- sf.inventory[sf.inventory$LS_PART == "S",]

inventory <- raster::rasterize(x = sf.scarp, y = slp, field = "LS_PART", background = 0) %>%
             raster::calc(x = ., fun = function(x){ifelse(x >= 1, 1, 0)})

hiPThresh <- Lslide::hiPassThresh(x = slp, scale.factor = 10, threshold = 50, env.rsaga = env.rsaga, sieve.thresh = 50,
                             quiet = FALSE, show.output.on.console = T)



# ------ 20: Run of optHPT:  7.404  Minutes ------
# length(c(3:15)) * length(c(2:9)) * 7.5/20 ------ Run of optHiPassThresh:  41.129  Minutes -----
optHPT.50 <- Lslide::optHiPassThresh(x = slp, inventory = inventory, range.scale.factor = c(3:15), range.threshold = c(2:9),
               sieve.thresh = 50, cores = 4, quiet = FALSE, env.rsaga = env.rsaga)


optHPT.50_2 <- Lslide::optHiPassThresh(x = slp, inventory = inventory, range.scale.factor = c(8:20), range.threshold = c(7:15),
                                     sieve.thresh = 50, cores = 4, quiet = FALSE, env.rsaga = env.rsaga)


opt <- dplyr::bind_rows(optHPT.50, optHPT.50_2)

saveRDS(opt, file = "d:/Users/qo23hel/Desktop/opt_test.rds")


# testing functions
seg <- sf::st_read("L:/PBAvsOBIA/Geom/Lososina/10m/Output/Segmentation/L2_seg.shp")
seg.sp <- as(seg, "Spatial")

## Main Direction --------------------
mainDir.c1 <- Lslide::mainDirection(spdf = seg.sp, quiet = FALSE, cores = 1)
mainDir.c4 <- Lslide::mainDirection(spdf = seg.sp, quiet = FALSE, cores = 4)

identical(mainDir.c4$angle, mainDir.c1$angle) # TRUE

seg.sp$MnDir <- mainDir.c4$angle


## Length-Width-Ratio --------------------
# ------ Run of LengthWidthRatio: 0.2318 Minutes ------
LeWiRat.c1 <- Lslide::lengthWidthRatio(spdf = seg.sp, cores = 1, quiet = FALSE)

# ------ Run of LengthWidthRatio: 0.1522 Minutes ------
LeWiRat.c4 <- Lslide::lengthWidthRatio(spdf = seg.sp, cores = 4, quiet = FALSE)

identical(LeWiRat.c1$ratio, LeWiRat.c4$ratio) # TRUE




## Get Bounding Box --------------------
seg.sp$MnDir <- mainDir.c4$angle

gBB <- Lslide::getBoundingBox(spdf = seg.sp, col.name = "MnDir", quiet = FALSE)

# ------ Run of getBoundingBox: 0.278 Minutes ------
gBB.Compl.c1 <- Lslide::getBoundingBox(spdf = seg.sp, scale.factor = c(2, 1.3), col.name = "MnDir",
                                    k = 2, scale.side = "long", centroid = TRUE, cores = 1, quiet = FALSE)

# ------ Run of getBoundingBox: 0.1808 Minutes ------
gBB.Compl.c4 <- Lslide::getBoundingBox(spdf = seg.sp, scale.factor = c(2, 1.3), col.name = "MnDir",
                                       k = 2, scale.side = "long", centroid = TRUE, cores = 4, quiet = FALSE)




## Neighbor Direction --------------------
nb <- spdep::poly2nb(pl = seg.sp)
nbDir.c1 <- Lslide::neighborDirection(spdf = seg.sp, col.name = "MnDir", cores = 1,
                                      tol = 75, modus = "nb", quiet = FALSE)

# ------ Run of ClassNeighborFunction:  0.0667  Minutes ------
nbDir.c4 <- Lslide::neighborDirection(spdf = seg.sp, col.name = "MnDir", cores = 4,
                                   tol = 75, nb =  nb, modus = "nb", quiet = FALSE)

identical(nbDir.c1, nbDir.c4) # TRUE
