## data
source("data_wrangling.r")
## fit model
require(stelfi)
## Spatial only model as per .ipynb
fit_spatial <- fit_lgcp(locs = locs, 
                        sf = colombia, 
                        smesh = smesh,
                        covariates = covariate_matrix,
                        parameters = list(beta = c(0, 0),  
                                          log_tau = log(1),
                                          log_kappa = log(1)))

get_coefs(fit_spatial)
## estimatde random field
get_fields(fit_spatial, smesh) |>
    show_field(smesh = smesh, sf = colombia, clip = TRUE)

## now using UTMs rather than longitude and latitude
points <- SV_combined_data_sf_entire_conflict
## downloaded 2020 population tif, for greater resolution
## need to update to include all yearly data
## source https://hub.worldpop.org/geodata/summary?id=40675
pop_raster <- terra::rast("../Data/col_pd_2020_1km.tif")
utm_epsg <- 32618
points_utm <- st_transform(points, crs = utm_epsg)
colombia_utm <- st_transform(colombia, crs = utm_epsg)
## keep points inside
## inside <- sf::st_within(points_utm, colombia_utm, sparse = FALSE)
## points_utm <- points_utm[rowSums(inside) > 0, ]
pop_raster_utm <- terra::project(pop_raster, terra::crs(points_utm))
## very course for illustration only
locs <- sf::st_coordinates(points_utm) |> as.data.frame(); names(locs) <- c("x", "y")
smesh <- fmesher::fm_mesh_2d_inla(boundary = colombia_utm,locs = locs,
                                  max.edge = c(50000, 150000), cutoff = 60000)
smesh$crs <- st_crs(points_utm)
## new population covariate
mesh_sf <- sf::st_as_sf(
   data.frame(x = smesh$loc[,1], y = smesh$loc[,2]),
   coords = c("x","y"),
   crs = st_crs(pop_raster_utm)
 )

## Extract raster population density at mesh nodes, note only 2020
## for now
covariate_matrix <-  terra::extract(pop_raster_utm, mesh_sf)[,2] |> as.matrix()
colnames(covariate_matrix) <- "population_density"
covariate_matrix[is.na(covariate_matrix[,1]), 1] <- 1e-6
## Spatial only model using UTMs and 2020 
fit_spatial_02 <- fit_lgcp(locs = locs, 
                        sf = colombia_utm, 
                        smesh = smesh,
                        covariates = covariate_matrix,
                        parameters = list(beta = c(0, 0),  
                                          log_tau = log(1),
                                          log_kappa = log(1)))
get_coefs(fit_spatial_02)
get_fields(fit_spatial_02, smesh) |>
    show_field(smesh = smesh, sf = colombia_utm, clip = TRUE)

## likelihood
fit_spatial_02$objective

## inlabru
require(INLA)
require(inlabru)

## for hyperpar priors
diameter <- as.numeric(
  st_length(st_cast(st_convex_hull(colombia_utm), "MULTILINESTRING"))
)
## SPDE
spde <- inla.spde2.pcmatern(
    mesh = smesh,
    prior.range = c(diameter/3, 0.5),  
    prior.sigma = c(1, 0.01))

## linear predictor
components <- geometry  ~ Intercept(1) + population_density(pop_raster_utm) +
    spatial(geometry, model = spde)

fit <- lgcp(components, points_utm,
    samplers = colombia_utm,
    domain = list(geometry = smesh),
    options = list(control.inla = list(int.strategy = "eb"))
  )

fit$summary.fixed
fit$summary.hyperpar


## field
show_field(fit$summary.random$spatial$mean, smesh = smesh, sf = colombia_utm, clip = TRUE)


## simulate from GMRFs

field_to_im <- function(x, smesh, sf, dims = c(500,500)){
    nx <- dims[1]
    ny <- dims[2]
    xs <- seq(min(smesh$loc[, 1]), max(smesh$loc[, 1]), length = nx)
    ys <- seq(min(smesh$loc[, 2]), max(smesh$loc[, 2]), length = ny)
    data <- expand.grid(xs = xs, ys = ys)
    pxl <- sf::st_multipoint(as.matrix(data))
    A <- fmesher::fm_basis(smesh, pxl)
    data$colz <-  as.vector(A %*% x)
    xy <- sf::st_as_sf(data, coords = c("xs", "ys"))
    sf::st_crs(xy) <- sf::st_crs(sf)
    idx <- lengths(sf::st_intersects(xy, sf)) > 0
    data[!idx, ] <- 0
    z_grid <- matrix(data$colz, nrow = ny, ncol = nx, byrow = TRUE)
    im_obj <- im(z_grid, xcol = xs, yrow = ys)
    return(im_obj)
}
## inlabru - sim 
image <- field_to_im(x = exp(fit$summary.random$spatial$mean - 3.3), smesh = smesh, sf = colombia_utm)
pixel_area <- (diff(range(image$xcol)) / (length(image$xcol)-1)) *
              (diff(range(image$yrow)) / (length(image$yrow)-1))

## scale intensity 
total_expected <- nrow(points_utm) 
image$v <- image$v / sum(image$v) * total_expected / pixel_area

## simulate from random field, no population effect
pp_sims <- spatstat.random::rpoispp(image, nsim = 100)
## pair correlation & Kest
library(spatstat.explore)
Kobs <- Kest(as.ppp(points_utm))
Ksim <- lapply(pp_sims, function(x) Kest(as.ppp(x)))
plot(Kobs, lwd = 2, main = "inlabru-Kest")
lapply(Ksim, plot, col = "pink", add = TRUE)

## stelfi
## fitted intensity
get_lambda <- function(obj, covariates){
    designmat <- cbind(1, covariates)
    res <- TMB::sdreport(obj)
    field <- res$par.random
    beta <- res$value["beta" == names(res$value)]
    beta <- as.matrix(beta, ncol = 1)
    lambda <- exp(field + designmat%*%beta)
    return(lambda)
}

image <- field_to_im(x = get_lambda(fit_spatial_02, covariate_matrix), smesh = smesh, sf = colombia_utm)

## simulate 9 realisations from the fitted intensity (including field)
pp_sims <- spatstat.random::rpoispp(image, nsim = 9)
lapply(pp_sims, plot, pch = 20, main = "simulated reslisation")
plot(locs, main = "observed data")


