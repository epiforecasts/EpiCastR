#' Calculate population centroids
#'
#' Calculate the population centroids of each polygon of a given shape file using a population raster.
#'
#' @param shape The shape file containing the polygons (areas). An sf object. Has a column `id`.
#' @param raster_obj The population raster. A RasterLayer object.
#' @param crs The coordinate reference system to use for both `shape` and `raster_obj`. No default given.
#'
#' @return A data.frame object with three columns: `id` (the ID of the polygon in shape), `lng` (longitude) and `lat` (latitude).
#'
calc_pop_centroids <- function(shape, raster_obj, crs) {

  if(missing(crs)) {
    stop("No CRS given.")
  }

  raster_obj <- raster::projectRaster(raster_obj,
                                      crs = crs)
  shape <- sf::st_transform(shape, crs = crs)

  out <- purrr::map_df(.x = shape$id,
                       .f = ~ {

                         raster_crop <- raster::crop(raster_obj, raster::extent(shape[shape$id == .x, ]))
                         raster_mask <- raster::mask(raster_crop, shape[shape$id == .x, ])
                         raster_points <- raster::rasterToPoints(raster_crop)

                         x_coord <- raster::weighted.mean(x = raster_points[,1],
                                                          w = raster_points[,3],
                                                          na.rm = TRUE)
                         y_coord <- raster::weighted.mean(x = raster_points[,2],
                                                          w = raster_points[,3],
                                                          na.rm = TRUE)

                         out <- tibble::tibble(id = .x,
                                               lng = x_coord,
                                               lat = y_coord)

                       }) %>%
    dplyr::bind_rows()

  return(out)

}
