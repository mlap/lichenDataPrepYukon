## Everything in this file and any files in the R directory are sourced during `simInit()`;
## all functions and objects are put into the `simList`.
## To use objects, use `sim$xxx` (they are globally available to all modules).
## Functions can be used inside any function that was sourced in this module;
## they are namespaced to the module, just like functions in R packages.
## If exact location is required, functions will be: `sim$.mods$<moduleName>$FunctionName`.
defineModule(sim, list(
  name = "lichenDataPrepYukon",
  description = "",
  keywords = "",
  authors = structure(list(list(given = c("Marcus", "Francois"), family = "Lapeyrolerie", role = c("aut", "cre"), email = "mlapeyro@mail.ubc.ca", comment = NULL)), class = "person"),
  childModules = character(0),
  version = list(lichenDataPrepYukon = "0.0.0.9000"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("NEWS.md", "README.md", "lichenDataPrepYukon.Rmd"),
  reqdPkgs = list("PredictiveEcology/SpaDES.core@box (>= 2.1.6.9000)", "reproducible", 
                  "terra", "whitebox"),
  parameters = bindrows(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
    defineParameter(".plots", "character", "screen", NA, NA,
                    "Used by Plots function, which can be optionally used here"),
    defineParameter(".plotInitialTime", "numeric", start(sim), NA, NA,
                    "Describes the simulation time at which the first plot event should occur."),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    "Describes the simulation time interval between plot events."),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "Describes the simulation time at which the first save event should occur."),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between save events."),
    defineParameter(".studyAreaName", "character", NA, NA, NA,
                    "Human-readable name for the study area used - e.g., a hash of the study",
                          "area obtained using `reproducible::studyAreaName()`"),
    ## .seed is optional: `list('init' = 123)` will `set.seed(123)` for the `init` event only.
    defineParameter(".seed", "list", list(), NA, NA,
                    "Named list of seeds to use for each event (names)."),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    "Should caching of events or module be used?")
  ),
  inputObjects = bindrows(
    #expectsInput("objectName", "objectClass", "input object description", sourceURL, ...),
    expectsInput(objectName = NA, objectClass = NA, desc = NA, sourceURL = NA)
  ),
  outputObjects = bindrows(
    #createsOutput("objectName", "objectClass", "output object description", ...),
    createsOutput(objectName = NA, objectClass = NA, desc = NA)
  )
))

doEvent.lichenDataPrepYukon = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      # do stuff for this event
      sim <- Init(sim)
      
      sim <- updateDynamicCovariates(sim)
      
      sim <- scheduleEvent(sim, time(sim) + P(sim)$predictionInterval, "lichenDataPrepYukon", "updateDynamicCovariates")

    },
    updateDynamicCovariates = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event

      sim <- updateDynamicCovariates(sim)

      # e.g.,
      sim <- scheduleEvent(sim, time(sim) + P(sim)$predictionInterval, "lichenDataPrepYukon", "updateDynamicCovariates")

      # ! ----- STOP EDITING ----- ! #
    },
    warning(noEventWarning(sim))
  )
  return(invisible(sim))
}

### template initialization
Init <- function(sim) {
  # # ! ----- EDIT BELOW ----- ! #
  # Loading data collected from Yukon's site, which is hosted on google drive
  targetFile <- "data"
  reproducible::prepInputs(url="https://drive.google.com/file/d/1OZjcd7Ln_SMkrJA-_mpjke2NruIJNe1M/view?usp=drive_link",
                           destinationPath = paths(sim)$inputPath,
                           targetFile = paste0(targetFile, "2.zip"),
                           fun = NA) |>
    reproducible::Cache()
  dataDirPath <- file.path(paths(sim)$inputPath, targetFile)
  whitebox::install_whitebox() |>
    reproducible::Cache()
  
  # Extracting topographic info
  demFile <- "50n150w_20101117_gmted_med075.tif"
  demPath <- file.path(dataDirPath, demFile)
  dem <- terra::rast(demPath)
  aspect <- terra::terrain(dem, v = "aspect", unit = "radians") |>
    reproducible::Cache()
  #wbt_breach_depressions(demPath, file.path(dataDirPath, "filledDEM.tif")) |>
  #  reproducible::Cache()
  #wbt_slope(file.path(dataDirPath, "filledDEM.tif"), file.path(dataDirPath, "slopeDEM.tif"), units = "radians") |>
  #  reproducible::Cache()
  #wbt_d8_flow_accumulation(file.path(dataDirPath, "filledDEM.tif"), file.path(dataDirPath, "d8flowDEM.tif")) |>
  #  reproducible::Cache()
  slope <- rast(file.path(dataDirPath, "../slopeDEM.tif"))
  d8flow <- rast(file.path(dataDirPath, "../d8flowDEM.tif"))
  wbt_twi <- log(d8flow / (tan(slope) + 0.001))
  
  # Projecting rasters to the appropriate CRS
  slope <- terra::project(slope, crs(sim$pixelGroupMap)) |>
    reproducible::Cache()
  elev <- terra::project(dem, crs(sim$pixelGroupMap)) |>
    reproducible::Cache()
  aspect <- terra::project(aspect, crs(sim$pixelGroupMap)) |>
    reproducible::Cache()
  wbt_twi <- terra::project(wbt_twi, crs(sim$pixelGroupMap)) |>
    reproducible::Cache()
  
  # Resampling rasters to match resolution and extent of study area
  elev_rs <- terra::resample(elev, sim$pixelGroupMap, "average") |>
    reproducible::Cache()
  slope_rs <- terra::resample(slope, sim$pixelGroupMap, "average") |>
    reproducible::Cache()
  aspect_rs <- terra::resample(aspect, sim$pixelGroupMap, "average") |>
    reproducible::Cache()
  wbt_twi_rs <- terra::resample(wbt_twi, sim$pixelGroupMap, "average") |>
    reproducible::Cache()
  
  # Removing values outside the study area
  elev_mask <- terra::mask(elev_rs, sim$studyArea) 
  slope_mask <- terra::mask(slope_rs, sim$studyArea)
  aspect_mask <- terra::mask(aspect_rs, sim$studyArea)
  wbt_twi_mask <- terra::mask(wbt_twi_rs, sim$studyArea)
  
  sim$lichenStaticCovariates <- data.table(pixelGroup = as.numeric(terra::values(sim$pixelGroupMap)),
                                           elev = as.numeric(terra::values(elev_mask)),
                                           slope = as.numeric(terra::values(slope_mask)),
                                           aspect = as.numeric(terra::values(aspect_mask)),
                                           wbt_twi = as.numeric(terra::values(wbt_twi_mask))) |>
    reproducible::Cache()
  # ! ----- STOP EDITING ----- ! #

  return(invisible(sim))
}
### template for save events
Save <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  sim <- saveFiles(sim)

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for plot events
plotFun <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  sampleData <- data.frame("TheSample" = sample(1:10, replace = TRUE))
  Plots(sampleData, fn = ggplotFn) # needs ggplot2

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for your event1
updateDynamicCovariates <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  cohortData <- vegTypeGenerator(sim$cohortData)
  cohortData$standType <- ifelse(grepl("Pice|Abie|Pinu", cohortData$leading, ignore.case = TRUE),
                            "C", "D")
  cohortData$standAge <- floor(cohortData$age / 10) * 10
  
  cohortData <- cohortData[, .(pixelGroup, standAge, standType)]
  cohortData <- cohortData[, .(
          standAge = first(standAge),
          standType = first(standType)
      ), by = pixelGroup]
  sim$lichenDynamicCovariates <- data.table(pixelGroup = as.numeric(terra::values(sim$pixelGroupMap)))
  
  sim$lichenDynamicCovariates[
    cohortData,
    on = "pixelGroup",
    `:=`(standAge = i.standAge, standType = i.standType)
  ]
  sim$lichenDynamicCovariates$inStudyArea <- terra::values(sim$rasterToMatch)
  sim$lichenDynamicCovariates[is.na(standType) & inStudyArea == 1, standType := "U"]
  sim$lichenDynamicCovariates[is.na(standAge) & inStudyArea == 1, standAge := 80] # REVISIT THIS
  sim$lichenDynamicCovariates[, inStudyArea := NULL]
  
  
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  # Any code written here will be run during the simInit for the purpose of creating
  # any objects required by this module and identified in the inputObjects element of defineModule.
  # This is useful if there is something required before simulation to produce the module
  # object dependencies, including such things as downloading default datasets, e.g.,
  # downloadData("LCC2005", modulePath(sim)).
  # Nothing should be created here that does not create a named object in inputObjects.
  # Any other initiation procedures should be put in "init" eventType of the doEvent function.
  # Note: the module developer can check if an object is 'suppliedElsewhere' to
  # selectively skip unnecessary steps because the user has provided those inputObjects in the
  # simInit call, or another module will supply or has supplied it. e.g.,
  # if (!suppliedElsewhere('defaultColor', sim)) {
  #   sim$map <- Cache(prepInputs, extractURL('map')) # download, extract, load file from url in sourceURL
  # }

  #cacheTags <- c(currentModule(sim), "function:.inputObjects") ## uncomment this if Cache is being used
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  # ! ----- EDIT BELOW ----- ! #

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

ggplotFn <- function(data, ...) {
  ggplot2::ggplot(data, ggplot2::aes(TheSample)) +
    ggplot2::geom_histogram(...)
}

