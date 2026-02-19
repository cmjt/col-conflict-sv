source("data_wrangling.r")
## fit model
require(stelfi)
require(spatstat.explore)
require(ggplot2)
## "show" observed data are clustered
## pair correlation & Kest usin UTMs
points <- SV_combined_data_sf_entire_conflict
utm_epsg <- 32618
points_utm <- st_transform(points, crs = utm_epsg)
K <- Kest(as.ppp(points_utm), correction = "iso")
df <- data.frame(r = K$r,K_iso = K$iso, K_theo = K$theo)
ggplot(df, aes(x = r)) +
  geom_line(aes(y = K_iso), linewidth = 1) +
  geom_line(aes(y = K_theo), linetype = "dashed") +
  labs(y = "K(r)", x = "Distance r (m)") +
    theme_minimal()
ggsave("../Output/ripley_plot.png")

## Spatial only model as per .ipynb
## sticking with lat long
t01 <- system.time({fit_spatial <- fit_lgcp(locs = locs, 
                        sf = colombia, 
                        smesh = smesh,
                        covariates = covariate_matrix,
                        parameters = list(beta = c(0, 0),  
                                          log_tau = log(1),
                                          log_kappa = log(1)))})["elapsed"]

coefs_01 <- get_coefs(fit_spatial) ## ensure matches current .ipynb


## inverstgigate mesh resolution
## NOTE, not using covariate data here
## as this is not available at mesh nodes

## coarser mesh
smesh_00 <- fmesher::fm_mesh_2d_inla(loc = locs[, 1:2], max.edge = c(2,3), cutoff = 2.5)
## original mesh
smesh_01  <- smesh
## finer mesh
smesh_02 <- fmesher::fm_mesh_2d_inla(loc = locs[, 1:2], max.edge = 0.6, cutoff = 0.9)
## finer mesh again
smesh_03 <- fmesher::fm_mesh_2d_inla(loc = locs[, 1:2], max.edge = 0.4, cutoff = 0.6)

## plot
png("../Output/meshes.png")
par(mfrow = c(2,2))
## 1
plot(smesh_00, col = "grey")
plot(colombia, col = NA, lwd = 2, add = TRUE)
legend("topleft", legend = "a", bty = "n")
## 2
plot(smesh_01, col = "grey")
plot(colombia, col = NA, lwd = 2, add = TRUE)
legend("topleft", legend = "b", bty = "n")
## 3
plot(smesh_02, col = "grey")
plot(colombia, col = NA, lwd = 2, add = TRUE)
legend("topleft", legend = "c", bty = "n")
## 4
plot(smesh_03, col = "grey")
plot(colombia, col = NA, lwd = 2, add = TRUE)
legend("topleft", legend = "d", bty = "n")
dev.off()
## plot mesh weights
require(patchwork)
get_weights(smesh_00, sf = colombia, plot = TRUE) + get_weights(smesh_01, sf = colombia, plot = TRUE) + get_weights(smesh_02, sf = colombia, plot = TRUE) + get_weights(smesh_03, sf = colombia, plot = TRUE) + plot_annotation(tag_levels = "a")
ggsave("../Output/weights_meshes.png")
## fit and compare models using different mesh
## same starting values
pars <- list(beta = c(0, 0), log_tau = log(1), log_kappa = log(1))
## function to get po_density at new mesh nodes
get_pop <- function(smesh, pop_density){
    mesh_coords <- smesh$loc[, 1:2]  
    mesh_coords_sf <- st_as_sf(as.data.frame(mesh_coords), 
                               coords = c("V1", "V2"), 
                               crs = 4326)
    pop_density_at_nodes <- st_join(mesh_coords_sf, 
                                    municipios %>% dplyr::select(pop_density), 
                                    join = st_intersects)
    covariate_matrix <- matrix(log(pop_density_at_nodes$pop_density + 1e-6), 
                               ncol = 1)
    covariate_matrix[is.na(covariate_matrix)] <- mean(covariate_matrix, na.rm = TRUE)
    return(covariate_matrix)
}
## mesh 00
cov_mat_00 <- get_pop(smesh = smesh_00, pop_density = pop_density)
t00 <- system.time({fit_00 <- fit_lgcp(locs = locs, 
                   sf = colombia, 
                   smesh = smesh_00,
                   covariates = cov_mat_00,
                   parameters = pars)})["elapsed"]
t00

coefs_00 <- get_coefs(fit_00)
## mesh 02
cov_mat_02 <- get_pop(smesh = smesh_02, pop_density = pop_density)
t02 <- system.time({fit_02 <- fit_lgcp(locs = locs, 
                   sf = colombia, 
                   smesh = smesh_02,
                   covariates = cov_mat_02,
                   parameters = pars)})["elapsed"]

coefs_02 <- get_coefs(fit_02)
## mesh 03
cov_mat_03 <- get_pop(smesh = smesh_03, pop_density = pop_density)
t03 <- system.time({fit_03 <- fit_lgcp(locs = locs, 
                   sf = colombia, 
                   smesh = smesh_03,
                   covariates = cov_mat_03,
                   parameters = pars)})["elapsed"]

coefs_03 <- get_coefs(fit_03)
## table of pars and mesh attributes
format_model <- function(coefs) {
  est <- sprintf("%.3f", coefs[, "Estimate"])
  se  <- sprintf("%.3f", coefs[, "Std. Error"])
  paste0(est, " (", se, ")")
}
get_expected <- function(fit,  cov, smesh, sf = colombia){
    w <- get_weights(smesh, sf = sf)$weights
    f <- get_fields(fit, smesh)
    coefs <- get_coefs(fit)
    beta0 <- coefs[1, 1]  
    beta1 <- coefs[2, 1]  
    expected_events_n <- sum(w * exp(beta0 + beta1 * cov[, 1] + f))
    return(round(expected_events_n, 0))
}
tab <- cbind(
  Model_01 = format_model(coefs_00),
  Model_02 = format_model(coefs_01),
  Model_03 = format_model(coefs_02),
  Model_04 = format_model(coefs_03)
)
rownames(tab) <- rownames(coefs_00)
## note observed num events 16258
extra <- rbind(
    'Expected Num Events' = c(get_expected(fit = fit_00, cov = cov_mat_00, smesh = smesh_00),
                            get_expected(fit = fit_spatial, cov = covariate_matrix, smesh = smesh_01),
                            get_expected(fit = fit_02, cov = cov_mat_02, smesh = smesh_02),
                            get_expected(fit = fit_03, cov = cov_mat_03, smesh = smesh_03)),
        'Runtime (s)' = c(
            round(t00["elapsed"], 2),round(t01["elapsed"], 2),
            round(t02["elapsed"], 2),round(t03["elapsed"], 2)),
    'Num Mesh Nodes' = c(smesh_00$n,smesh_01$n, smesh_02$n,smesh_03$n)
)
rbind(tab, extra)
## plot fields
fits <- list(fit_00, fit_spatial, fit_02, fit_03)
meshes  <- list(smesh_00, smesh_01, smesh_02, smesh_03)
plots <- Map(function(fit, smesh) {
  p1 <- show_lambda(fit, smesh, sf = colombia, clip = TRUE)
  p2 <- get_fields(fit, smesh) |>
    show_field(smesh = smesh, sf = colombia, clip = TRUE)
  p3 <- get_fields(fit, smesh, sd = TRUE) |>
      show_field(smesh = smesh, sf = colombia, clip = TRUE)
  list(p1, p2, p3)
  
}, fits, meshes)

plots_flat <- unlist(plots, recursive = FALSE)

p1 <- wrap_elements((plots_flat[[1]] + plots_flat[[2]] + plots_flat[[3]]) +
                    plot_annotation(title = "Model 01", tag_levels = list(c("Intensity Surface", "GMRF", "Spatial Uncertainty"))) & theme_void())
p2 <- wrap_elements((plots_flat[[4]] + plots_flat[[5]] + plots_flat[[6]]) +
                    plot_annotation(title = "Model 02", tag_levels = list(c("Intensity Surface", "GMRF", "Spatial Uncertainty"))) & theme_void())
p3 <- wrap_elements((plots_flat[[7]] + plots_flat[[8]] + plots_flat[[9]]) +
                    plot_annotation(title = "Model 03", tag_levels = list(c("Intensity Surface", "GMRF", "Spatial Uncertainty"))) & theme_void())
p4 <- wrap_elements((plots_flat[[10]] + plots_flat[[11]] + plots_flat[[12]]) +
                    plot_annotation(title = "Model 04", tag_levels = list(c("Intensity Surface", "GMRF", "Spatial Uncertainty"))) & theme_void())

wrap_plots(p1, p2, p3, p4, ncol = 1, axis_titles = "collect")
ggsave("../Output/model_fit.png", width = 14, height = 10)

## simulate GMRF from estimated intensity surface from ms model 

## from image
## ASSUMES UTM, We should change to UTMs overall really..
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
get_lambda <- function(obj, covariates){
    designmat <- cbind(1, covariates)
    res <- TMB::sdreport(obj)
    field <- res$par.random
    beta <- res$value["beta" == names(res$value)]
    beta <- as.matrix(beta, ncol = 1)
    lambda <- exp(field + designmat%*%beta)
    return(lambda)
}
image <- field_to_im(x = get_lambda(fit_03, cov_mat_03), smesh = smesh_03, sf = colombia)
## simulate realisations from the fitted intensity (including field)
pp_sims <- spatstat.random::rpoispp(image, nsim = 1000)
png("../Output/hist_n_est.png")
hist(lapply(pp_sims, function(x) x$n) |> unlist(), main = "", xlab = "Simulated number of events")
abline(v = nrow(locs), col = "red", lwd = 2)
legend("topright", bty = "n", legend = "Observed number of events", text.col = "red")
dev.off()

## plot pp visially
pp_obs <- as.ppp(st_coordinates(points),
                 W = as.owin(image))
sigma_val <- 0.2 ## bandwidth
dens_obs <- density(pp_obs, sigma = sigma_val, dimx = 500)
dens_sims <- lapply(pp_sims, density, sigma = sigma_val, dimx = 500)
sim_array <- simplify2array(lapply(dens_sims, function(x) x$v))
mean_sim <- apply(sim_array, c(1,2), mean)
lo_sim   <- apply(sim_array, c(1,2), quantile, probs = 0.025)
hi_sim   <- apply(sim_array, c(1,2), quantile, probs = 0.975)
mean_im <- dens_obs;mean_im$v <- mean_sim
lo_im <- dens_obs;lo_im$v <- lo_sim
hi_im <- dens_obs;hi_im$v <- hi_sim

png("../Output/simulated_density_surfaces.png")
par(mfrow = c(2,2), mar = c(2,2,3,2))
zlim_range <- range(c(dens_obs$v, mean_sim))
## 1
plot(dens_obs, main = "Observed density, bw = 0.2", zlim = zlim_range)
plot(colombia, add = TRUE, col = NA, lwd = 2)
## 2
plot(mean_im,  main = "Mean simulated density, bw = 0.2", zlim = zlim_range)
plot(colombia, add = TRUE, col = NA, lwd = 2)
## 3
plot(lo_im, main = "2.5%Q, bw = 0.2", zlim = zlim_range)
plot(colombia, add = TRUE, col = NA, lwd = 2)
## 4
plot(hi_im,  main = "97.5% Q, bw = 0.2", zlim = zlim_range)
plot(colombia, add = TRUE, col = NA, lwd = 2)
dev.off()
## Data are likely quite a bit more clustered than model assumes,
## perhaps due to stacking of points in villages/regions

