## data wranglins as per ../Pappa_et_al_EntireConflict.ipynb
library(dplyr)         
library(readxl)
library(janitor)
library(spatstat.geom)
library(sf)
library(rnaturalearth)
library(fmesher)
## devtools::install_github("nebulae-co/colmaps")
library(colmaps)
library(terra)
##  relative path
data_path <- "../Data/"
## setup
required_columns <- c("ID Caso", "Código DANE de Municipio", "Municipio", "Departamento", 
                      "Año", "Mes", "Día", "ID Persona", "Sexo", "Etnia", "Ocupación", 
                      "Calidad de la Víctima o la Baja", "Latitud", "Longitud", "Edad",
                      "Fuerza o Grupo Armado Organizado al que Pertenece el Combatiente", "Situación Actual de la Víctima", 
                      "Días de Cautiverio",  "No. de Veces secuestrado")

file_list <- list.files(path = data_path, pattern = ".xlsx$", full.names = TRUE)

detect_and_convert_types <- function(df) {
  df_converted <- df
  for (col_name in names(df)) {
    col_data <- df[[col_name]]
    if (all(is.na(col_data))) {
      next
    }
    non_na_values <- col_data[!is.na(col_data)]
        numeric_conversion <- suppressWarnings(as.numeric(non_na_values))
    if (sum(!is.na(numeric_conversion)) / length(non_na_values) > 0.8) {
      df_converted[[col_name]] <- suppressWarnings(as.numeric(col_data))
    } else {
      df_converted[[col_name]] <- as.character(col_data)
    }
  }
  return(df_converted)
}

data_list <- lapply(file_list, function(file) {
  df <- read_excel(file, col_types = "text")
  df <- detect_and_convert_types(df)
  required_columns <- as.character(required_columns)
  missing_cols <- setdiff(required_columns, names(df))
  if (length(missing_cols) > 0) {
    df[missing_cols] <- NA
  }
  df <- dplyr::select(df, intersect(names(df), required_columns))
  df$Source_File <- basename(file)
  return(df)
})
     
# Combine into single dataframe
combined_data <- bind_rows(data_list)

# Some data tidying
# Turns to snake_case, removes accents
combined_data <- janitor::clean_names(combined_data)  

# Replace all blank values with NA
combined_data[combined_data == ""] <- NA

# If columns are characters and might include blanks
combined_data <- combined_data %>%
  mutate(across(c(ano, mes, dia), ~ ifelse(. == 0 | . == "" | is.na(.), NA, .)))

# Ensure the columns are numeric
combined_data$ano <- as.numeric(combined_data$ano)
combined_data$mes <- as.numeric(combined_data$mes)
combined_data$dia <- as.numeric(combined_data$dia)

# Remove rows with any missing date parts
combined_data <- combined_data[complete.cases(combined_data[, c("ano", "mes", "dia")]), ]

# Combine year, month, day into a string
combined_data$DateString <- paste(combined_data$ano,
                                  sprintf("%02d", combined_data$mes),
                                  sprintf("%02d", combined_data$dia),
                                  sep = "-")

# Convert to POSIXct
combined_data$Date <- as.POSIXct(combined_data$DateString, format = "%Y-%m-%d", tz = "UTC")

# Recode violence event types
combined_data <- combined_data %>%
  dplyr::mutate(EventType = case_when(
    source_file == "VictimasAB_202409.xlsx" ~ "AB", # Acciones Bélicas
    source_file == "VictimasAP_202409.xlsx" ~ "AP", # Ataques Poblaciones
    source_file == "VictimasAS_202409.xlsx" ~ "AS", # Asesinatos selectivos
    source_file == "VictimasAT_202409.xlsx" ~ "AT", # Atentados Terroristas
    source_file == "VictimasDB_202409.xlsx" ~ "DB", # Daños a bienes civiles
    source_file == "VictimasDF_202409.xlsx" ~ "DF", # Desaparición forzada
    source_file == "VictimasMA_202409.xlsx" ~ "MA", # Masacres
    source_file == "VictimasMI_202409.xlsx" ~ "MI", # Minas
    source_file == "VictimasRU_202409.xlsx" ~ "RU", # Reclutamiento y utilización de niños, niñas y adolescentes
    source_file == "VictimasSE_202409.xlsx" ~ "SE", # Secuestos
    source_file == "VictimasVS_202409.xlsx" ~ "VS", # Violencia sexual
    TRUE ~ NA_character_  # default case if no match
  ))

combined_data$EventType <- as.factor(combined_data$EventType)

# Convert to numeric time (e.g., for some Hawkes models)
combined_data$Timestamp <- as.numeric(combined_data$Date)

# Convert to spatial point pattern (longitude & latitude as spatial features)

# Define observation window (bounding box for spatial analysis)
W <- owin(xrange = c(min(combined_data$longitud), max(combined_data$longitud)), 
          yrange = c(min(combined_data$latitud), max(combined_data$latitud)))

# Convert time dimension to relative time units
combined_data$Time <- as.numeric(difftime(combined_data$Timestamp, min(combined_data$Timestamp), units="days"))

# Remove missing timepoints
combined_data <- combined_data[!is.na(combined_data$Timestamp), ]

# Sort the data by time:
combined_data <- combined_data[order(combined_data$Timestamp), ]

#Convert to numeric relative time (in days or seconds):
combined_data$t <- as.numeric(difftime(combined_data$Timestamp, min(combined_data$Timestamp), units = "days"))
combined_data <- combined_data[order(combined_data$t),]
combined_data$x <- combined_data$longitud
combined_data$y <- combined_data$latitud

# Convert data frame to an sf object
combined_data_sf <- st_as_sf(combined_data, coords = c("longitud", "latitud"), crs = 4326)

# Set the CRS to the desired one
st_crs(combined_data_sf) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Filter out rows where coordinates are 0
combined_data_sf <- combined_data_sf %>%
  filter(st_coordinates(geometry)[,1] != 0 & st_coordinates(geometry)[,2] != 0)

# Filter out rows with non-finite Date values
combined_data_sf <- combined_data_sf %>%
  filter(!is.na(Date))

# Ensure the variable is in Date format
combined_data_sf$DateString <- as.Date(combined_data_sf$DateString)

# Extract the year from the date
combined_data_sf$Year <- format(combined_data_sf$DateString, "%Y")

# Convert year to numeric
combined_data_sf$Year <- as.numeric(combined_data_sf$Year)

# Filter for dates from 1964-05-27 onwards
combined_data_sf_entire_conflict <- combined_data_sf %>%
  filter(DateString >= as.Date("1964-05-27"))

# Only keep sexual violence events
SV_combined_data_sf_entire_conflict <- combined_data_sf_entire_conflict %>%
  filter(EventType == "VS")

colombia <- ne_countries(country = "Colombia", returnclass = "sf") %>%
  sf::st_make_valid()

# dataframe of sighting locations (lat, long)
locs <- sf::st_coordinates(SV_combined_data_sf_entire_conflict) %>%
  as.data.frame() %>%
    rename(., c("x" = "X", "y" = "Y"))

locs <- cbind(locs, year = SV_combined_data_sf_entire_conflict$Year)

# Filter out datapoints in San Andres and Providencia
locs <- locs %>%
  filter(!(x < -80 & y > 11 & y < 14.5))  # Remove points in that box
  
# Delauney triangluation of domain
smesh <- fmesher::fm_mesh_2d_inla(loc = locs[, 1:2], max.edge = 1, cutoff = 1)

pop_file <- paste(data_path, "municipios_pop_1985-2024.csv", sep = "")

# Check if file exists, if not download and save
if (!file.exists(pop_file)) {
  municipios_pop <- download_pop_projections("municipality", 1985, 2024)
  # Fix errors in the fetched data from orginal source!
  # 1985 - 2014 codigo_municipio column is labelled correctly, municipio column is labelled correctly
  # 2005 - 2019 codigo_municipio column and municipio column have been mislabelled. the labels have been swapped
  # 2020 onwards codigo_municipio column is labelled correctly, municipio column is labelled correctly
  municipios_pop <- municipios_pop %>%
    mutate(
      codigo_municipio = ifelse(ano >= 2005 & ano <= 2019, municipio, codigo_municipio),
      municipio = ifelse(ano >= 2005 & ano <= 2019, codigo_municipio, municipio)
    ) 
  write.csv(municipios_pop, pop_file, row.names = FALSE)
}

# Load population data
municipios_pop <- read.csv(pop_file, 
                           colClasses = c(codigo_municipio = "character"))

# Filter for total population (data file has separate entries for rural, urban and total within each municipio)
# and for the correct year range (1985-2024 - the only data range available) for the analysis
municipios_pop <- municipios_pop %>%
  filter(area == "total") %>%
    filter(between(ano, 1985, 2024))
# Calculate average population by municipality
avg_population <- municipios_pop %>%
  group_by(codigo_municipio) %>%
  summarise(avg_population = mean(total, na.rm = TRUE))

# Get municipality spatial data
municipios_sp <- colmaps::municipios
municipios <- st_as_sf(municipios_sp)

# Ensure codigo_municipio is character for proper joining
municipios$codigo_municipio <- as.character(municipios$id)

# Join population data
municipios <- municipios %>%
  dplyr::left_join(avg_population, by = "codigo_municipio")

# Calculate area and population density
municipios$area_km2 <- as.numeric(st_area(municipios)) / 1e6

municipios$pop_density <- municipios$avg_population / municipios$area_km2

# Handle NA and infinite values
municipios$pop_density[is.na(municipios$pop_density) | is.infinite(municipios$pop_density)] <- 1e-6

# Ensure CRS consistency (WGS84 / EPSG:4326)
municipios <- st_transform(municipios, crs = 4326)

# Create raster template with proper CRS
bbox <- st_bbox(municipios)
raster_template <- rast(xmin = bbox["xmin"], xmax = bbox["xmax"],
                        ymin = bbox["ymin"], ymax = bbox["ymax"],
                        resolution = 0.01, crs = "EPSG:4326")

# Convert to terra vector for rasterization
municipios_vect <- vect(municipios)

# Rasterize population density
pop_density_raster <- rasterize(municipios_vect, raster_template, 
                                field = "pop_density", fun = "mean")

# Log-transform population density (recommended for skewed distributions)
# Add small constant to avoid log(0)
pop_density_raster_log <- log(pop_density_raster + 1e-6)

# Convert raster to matrix for fit_lgcp
# The covariate matrix should match the mesh structure
pop_density_matrix <- as.matrix(pop_density_raster_log, wide = TRUE)

# Alternative: Extract values at mesh node locations
mesh_coords <- smesh$loc[, 1:2]  # Get mesh node coordinates
mesh_coords_sf <- st_as_sf(as.data.frame(mesh_coords), 
                           coords = c("V1", "V2"), 
                           crs = 4326)

# Extract population density values at mesh nodes
pop_density_at_nodes <- st_join(mesh_coords_sf, 
                                municipios %>% dplyr::select(pop_density), 
                                join = st_intersects)

# Create covariate matrix (one column for pop density)
covariate_matrix <- matrix(log(pop_density_at_nodes$pop_density + 1e-6), 
                           ncol = 1)

# Handle any remaining NA values
covariate_matrix[is.na(covariate_matrix)] <- mean(covariate_matrix, na.rm = TRUE)
