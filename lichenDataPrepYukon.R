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
      ### check for more detailed object dependencies:
      ### (use `checkObject` or similar)

      # do stuff for this event
      sim <- Init(sim)
      
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
  # Loading data collected from Yukon's site, which we host on a google drive for
  # efficiency
  browser()
  targetFile <- "data"
  reproducible::prepInputs(url="https://drive.google.com/file/d/1OZjcd7Ln_SMkrJA-_mpjke2NruIJNe1M/view?usp=drive_link",
                           destinationPath = sim$paths$inputPath,
                           targetFile = paste0(targetFile, ".zip")) |>
    reproducible::Cache()
  utils::unzip(zipfile = file.path(sim$paths$inputPath, paste0(targetFile, ".zip")), 
               exdir = file.path(sim$paths$inputPath, targetFile)) |>
    reproducible::Cache()
  dataDirPath <- file.path(sim$paths$inputPath, targetFile, "data")
  
  # Extracting topographic info
  demFile <- "50n150w_20101117_gmted_med075.tif"
  demPath <- file.path(dataDirPath, demFile)
  dem <- terra::rast(demPath)
  aspect <- terra::terrain(dem, v = "aspect", unit = "radians")
  wbt_breach_depressions(demFile, "./filledDEM.tif")
  wbt_slope("./filledDEM.tif", "./slopeDEM.tif", units = "radians")
  wbt_d8_flow_accumulation("./filledDEM.tif", "./d8flowDEM.tif")
  slope <- rast("./slopeDEM.tif")
  d8flow <- rast("./d8flowDEM.tif")
  wbt_twi <- log(d8flow / (tan(slope) + 0.001))
  
  # Linking NTEMS id's to stand type
  speciesKeyPath <- file.path(dataDirPath, "sppEquivalencies_CA.csv")
  speciesKey <- read.csv(speciesKeyPath) |> 
    filter(is.na(NTEMS_Species_Code) == FALSE)
  treeLegend <- speciesKey[, c("NTEMS_Species_Code", "Type")] |> unique()
  treeLegend[nrow(treeLegend) + 1, ] <- c(0, "Unf")
  abbreviationLabels <- c("C", "B", "U")
  abbreviationCodes <- c("Conifer", "Deciduous", "Unf")
  treeLegend$Type <- abbreviationLabels[match(treeLegend$Type, abbreviationCodes)]
  treeLegend$NTEMS_Species_Code <- as.numeric(treeLegend$NTEMS_Species_Code)
  levels(treeInventoryYukon) <- treeLegend
  rm(speciesKey, treeLegend, abbreviationLabels, abbreviationCodes)
  
  # Filling in background of fire year and Stand type
  unforestedVect <- terra::deepcopy(sim$studyArea)
  unforestedVect$zone <- "U"
  unforestedBackground <- rasterize(unforestedVect, sim$rasterToMatch, field = "zone")
  treeInventoryFilled <- terra::merge(treeInventory, unforestedbackground)
  rm(fireBackground, vegBackground, tSinceFireRast, vegInventoryYukon)
  
  # Project to the same CRS
  slope <- terra::project(slope, crs(sim$studyArea))
  elev <- terra::project(dem, crs(sim$studyArea))
  aspect <- terra::project(aspect, crs(sim$studyArea))
  wbt_twi <- terra::project(wbt_twi, crs(sim$studyArea))
  
  elev_rs <- terra::resample(terra::mask(elev, sim$studyArea), sim$rasterToMatch, "average")
  slope_rs <- terra::resample(slope, sim$rasterToMatch, "average")
  aspect_rs <- terra::resample(aspect, sim$rasterToMatch, "average")
  wbt_twi_rs <- terra::resample(wbt_twi, sim$rasterToMatch, "average")
  treeInventory_rs <- terra::resample(treeInventoryFilled, sim$rasterToMatch, "mode")
  tSinceFireRast_rs <- terra::resample(tSinceFireRastFilled, sim$rasterToMatch, "mode")
  rm(slope, elev, aspect, wbt_twi, vegInventoryYukonFilled, tSinceFireRastFilled, lichen_presence_absence)
  
  # RETHINK dynamicInputs HERE
  staticInputs <- c(slope_rs, elev_rs, aspect_rs, wbt_twi_rs) 
  dynamicInputs <- c(tSinceFireRast_rs, vegInventoryYukon_rs)
  staticDT <- as.data.table(staticInputs, na.rm = TRUE, xy = TRUE)
  dynamicDT <- as.data.table(dynamicInputs, na.rm = TRUE, xy = TRUE)
  data.table::fwrite(staticDT, file = file.path(dataDirPath, "staticInputs.csv"), row.names = FALSE)
  data.table::fwrite(dynamicDT, file = file.path(dataDirPath, "dynamicInputs.csv"), row.names = FALSE)
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
  # WRITE CODE TO TAKE STAND AGE AND STAND CLASS FOR COHORT DATA
  
  # OUTPUT A NEW dynamicInputs.csv
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

