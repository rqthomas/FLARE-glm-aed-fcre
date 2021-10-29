#renv::restore()
#diskutil partitionDisk $(hdiutil attach -nomount ram://2048000) 1 GPTFormat APFS 'ramdisk' '100%'
install.packages("remotes")
remotes::install_github("FLARE-forecast/GLM3r")
remotes::install_github("FLARE-forecast/FLAREr")


library(tidyverse)
library(lubridate)

lake_directory <- here::here()
update_run_config <- TRUE #TRUE is used for an iterative workflow

files.sources <- list.files(file.path(lake_directory, "R"), full.names = TRUE)
sapply(files.sources, source)

#start_day <- as_date("2018-12-17")
start_day <- as_date("2018-07-20")
forecast_start<- as_date("2018-09-01")
#forecast_start<- as_date("2018-12-31")
#forecast_start<- as_date("2019-12-31")
#forecast_start<- as_date("2019-12-31")
holder1 <- start_day
holder2 <- forecast_start
while(forecast_start < as_date("2020-01-01")){
  start_day <- forecast_start
  forecast_start <- forecast_start + days(7)
  if(forecast_start < as_date("2020-01-01")){
    holder1 <- c(holder1, start_day)
    holder2 <- c(holder2, forecast_start)
  }
}

forecast_days_vector <- rep(16, length(holder1))
forecast_days_vector[1] <- 0
forecasting_timings <- data.frame(holder1,holder2,forecast_days_vector)

saved_file <- NA
#saved_file <- "/Users/quinn/workfiles/Research/SSC_forecasting/run_flare_package/flarer_aed/forecast_output/test_run_H_2018_09_08_2018_09_15_F_16_20210908T101708.nc"


config_obs <- yaml::read_yaml(file.path(lake_directory,"configuration","observation_processing","observation_processing.yml"))
config <- yaml::read_yaml(file.path(lake_directory,"configuration","FLAREr","configure_flare_glm_aed.yml"))
#Note: lake_directory need to be set prior to running this script
config$file_path$qaqc_data_directory <- file.path(lake_directory, "data_processed")
config$file_path$data_directory <- file.path(lake_directory, "data_raw")
config$file_path$configuration_directory <- file.path(lake_directory, "configuration")
config$file_path$execute_directory <- file.path(lake_directory, "flare_tempdir")
config$file_path$run_config <- file.path(lake_directory, "configuration", "FLAREr","configure_run.yml")
config$file_path$forecast_output_directory <- file.path(lake_directory, "forecasts")
config$file_path$noaa_directory <- file.path(lake_directory, "drivers")


run_config <- yaml::read_yaml(config$file_path$run_config)

setwd(file.path(lake_directory, "data_raw"))
if(!dir.exists(file.path(lake_directory, "data_raw", config_obs$realtime_insitu_location))){
  system(paste0("git clone --depth 1 --single-branch --branch ",config_obs$realtime_insitu_location, " https://github.com/FLARE-forecast/FCRE-data.git ", config_obs$realtime_insitu_location))
}else{
  setwd(file.path(lake_directory, "data_raw", config_obs$realtime_insitu_location))
  system("git pull")
}

setwd(file.path(lake_directory, "data_raw"))
if(!dir.exists(file.path(lake_directory, "data_raw", config_obs$realtime_met_station_location))){
  system(paste0("git clone --depth 1 --single-branch --branch ",config_obs$realtime_met_station_location, " https://github.com/FLARE-forecast/FCRE-data.git ", config_obs$realtime_met_station_location))
}else{
  setwd(file.path(lake_directory, "data_raw", config_obs$realtime_met_station_location))
  system("git pull")
}

setwd(file.path(lake_directory, "data_raw"))
if(!dir.exists(file.path(lake_directory, "data_raw", config_obs$realtime_inflow_data_location))){
  system(paste0("git clone --depth 1 --single-branch --branch ",config_obs$realtime_inflow_data_location, " https://github.com/FLARE-forecast/FCRE-data.git ", config_obs$realtime_inflow_data_location))
}else{
  setwd(file.path(lake_directory, "data_raw", config_obs$realtime_inflow_data_location))
  system("git pull")
}

setwd(file.path(lake_directory, "data_raw"))
if(!dir.exists(file.path(lake_directory, "data_raw", config_obs$manual_data_location))){
  system(paste0("git clone --depth 1 --single-branch --branch ",config_obs$manual_data_location, " https://github.com/FLARE-forecast/FCRE-data.git ", config_obs$manual_data_location))
}else{
  setwd(file.path(lake_directory, "data_raw", config_obs$manual_data_location))
  system("git pull")
}

if(!file.exists(file.path(lake_directory, "data_raw", config_obs$met_raw_obs_fname[2]))){
  download.file("https://pasta.lternet.edu/package/data/eml/edi/389/5/3d1866fecfb8e17dc902c76436239431", destfile = file.path(lake_directory, "data_raw",config_obs$manual_data_location,"/Met_final_2015_2020.csv"), method="curl")
}

if(!file.exists(file.path(lake_directory, "data_raw", config_obs$inflow_raw_file1[2]))){
  download.file("https://pasta.lternet.edu/package/data/eml/edi/202/7/f5fa5de4b49bae8373f6e7c1773b026e", destfile = file.path(lake_directory, "data_raw",config_obs$manual_data_location,"/inflow_for_EDI_2013_10Jan2021.csv"), method="curl")
}

if(!file.exists(file.path(lake_directory, "data_raw", config_obs$insitu_obs_fname[2]))){
  download.file("https://pasta.lternet.edu/package/data/eml/edi/271/5/c1b1f16b8e3edbbff15444824b65fe8f", destfile = file.path(lake_directory, "data_raw",config_obs$manual_data_location,"/Catwalk_cleanedEDI.csv"), method="curl")
}

if(!file.exists(file.path(lake_directory, "data_raw", config_obs$secchi_fname))){
  download.file("https://pasta.lternet.edu/package/data/eml/edi/198/8/336d0a27c4ae396a75f4c07c01652985", destfile = file.path(lake_directory, "data_raw",config_obs$manual_data_location,"/Secchi_depth_2013-2020.csv"), method="curl")
}




met_qaqc(realtime_file = file.path(config$file_path$data_directory, config_obs$met_raw_obs_fname[1]),
         qaqc_file = file.path(config$file_path$data_directory, config_obs$met_raw_obs_fname[2]),
         cleaned_met_file_dir = config$file_path$qaqc_data_directory,
         input_file_tz = "EST",
         nldas = file.path(config$file_path$data_directory, config_obs$nldas))

cleaned_inflow_file <- paste0(config$file_path$qaqc_data_directory, "/inflow_postQAQC.csv")

inflow_qaqc(realtime_file = file.path(config$file_path$data_directory, config_obs$inflow_raw_file1[1]),
            qaqc_file = file.path(config$file_path$data_directory, config_obs$inflow_raw_file1[2]),
            nutrients_file = file.path(config$file_path$data_directory, config_obs$nutrients_fname),
            cleaned_inflow_file ,
            input_file_tz = 'EST')


cleaned_observations_file_long <- paste0(config$file_path$qaqc_data_directory,
                                         "/observations_postQAQC_long.csv")

config_obs$data_location <- config$file_path$data_directory
in_situ_qaqc(insitu_obs_fname = file.path(config$file_path$data_directory,config_obs$insitu_obs_fname),
             data_location = config$file_path$data_directory,
             maintenance_file = file.path(config$file_path$data_directory,config_obs$maintenance_file),
             ctd_fname = file.path(config$file_path$data_directory,config_obs$ctd_fname),
             nutrients_fname =  file.path(config$file_path$data_directory, config_obs$nutrients_fname),
             secchi_fname = file.path(config$file_path$data_directory, config_obs$secchi_fname),
             cleaned_observations_file_long = cleaned_observations_file_long,
             lake_name_code = config$location$lake_name_code,
             config_obs = config_obs)


#for(i in 1:nrow(forecasting_timings)){
  for(i in 1:1){

  #config_obs <- yaml::read_yaml(file.path(lake_directory,"configuration","observation_processing","observation_processing.yml"))
  #config <- yaml::read_yaml(file.path(lake_directory,"configuration","FLAREr","configure_flare_glm_aed.yml"))
  run_config <- yaml::read_yaml(file.path(lake_directory,"configuration","FLAREr","configure_run.yml"))

  run_config$start_datetime <- as.character(forecasting_timings[i,1])
  run_config$forecast_start_datetime <- as.character(forecasting_timings[i, 2])
  run_config$forecast_horizon <- forecasting_timings[i, 3]
  run_config$restart_file <- saved_file
  yaml::write_yaml(run_config, file = file.path(lake_directory,"configuration","FLAREr","configure_run.yml"))

  file.copy(file.path(config$file_path$data_directory,config$management$sss_fname), file.path(config$file_path$qaqc_data_directory,basename(config$management$sss_fname)))

  config$run_config <- run_config

  if(!dir.exists(config$file_path$execute_directory)){
    dir.create(config$file_path$execute_directory)
  }

  pars_config <- readr::read_csv(file.path(config$file_path$configuration_directory, "FLAREr", config$model_settings$par_config_file), col_types = readr::cols())
  obs_config <- readr::read_csv(file.path(config$file_path$configuration_directory, "FLAREr", config$model_settings$obs_config_file), col_types = readr::cols())
  states_config <- readr::read_csv(file.path(config$file_path$configuration_directory, "FLAREr", config$model_settings$states_config_file), col_types = readr::cols())


  #Download and process observations (already done)

  cleaned_observations_file_long <- file.path(config$file_path$qaqc_data_directory,"observations_postQAQC_long.csv")
  cleaned_inflow_file <- file.path(config$file_path$qaqc_data_directory, "/inflow_postQAQC.csv")
  observed_met_file <- file.path(config$file_path$qaqc_data_directory,"observed-met_fcre.nc")


  #Step up Drivers

  start_datetime <- lubridate::as_datetime(config$run_config$start_datetime)
  if(is.na(config$run_config$forecast_start_datetime)){
    end_datetime <- lubridate::as_datetime(config$run_config$end_datetime)
    forecast_start_datetime <- end_datetime
  }else{
    forecast_start_datetime <- lubridate::as_datetime(config$run_config$forecast_start_datetime)
    end_datetime <- forecast_start_datetime + lubridate::days(config$run_config$forecast_horizon)
  }
  forecast_hour <- lubridate::hour(forecast_start_datetime)
  if(forecast_hour < 10){forecast_hour <- paste0("0",forecast_hour)}
  noaa_forecast_path <- file.path(config$file_path$noaa_directory, config$met$forecast_met_model,config$location$site_id,lubridate::as_date(forecast_start_datetime),forecast_hour)


  forecast_files <- list.files(noaa_forecast_path, full.names = TRUE)
  

  if(length(forecast_files) == 0){
    Sys.setenv("AWS_DEFAULT_REGION"="data",
               "AWS_S3_ENDPOINT"="rquinnthomas.com")
    download_s3_objects(lake_directory = lake_directory,
                        bucket = "drivers",
                        prefix = file.path(config$met$forecast_met_model,config$location$site_id,lubridate::as_date(forecast_start_datetime),forecast_hour))
  }


  print("Creating FLARE Met files")
  met_out <- FLAREr::generate_glm_met_files(obs_met_file = observed_met_file,
                                            out_dir = config$file_path$execute_directory,
                                            forecast_dir = noaa_forecast_path,
                                            config = config)

  met_file_names <- met_out$met_file_names
  historical_met_error <- met_out$historical_met_error

  if(config$model_settings$model_name == "glm_aed"){

  file.copy(file.path(config$file_path$data_directory, "fcre-manual-data/FCR_weir_inflow_2013_2019_20200828_allfractions_2poolsDOC.csv"),
            file.path(config$file_path$execute_directory, "FCR_weir_inflow_2013_2019_20200828_allfractions_2poolsDOC.csv"))

  file.copy(file.path(config$file_path$data_directory, "fcre-manual-data/FCR_wetland_inflow_2013_2019_20200828_allfractions_2DOCpools.csv"),
            file.path(config$file_path$execute_directory, "FCR_wetland_inflow_2013_2019_20200828_allfractions_2DOCpools.csv"))

  file.copy(file.path(config$file_path$data_directory, "fcre-manual-data/FCR_SSS_inflow_2013_2019_20200701_allfractions_2DOCpools.csv"),
            file.path(config$file_path$execute_directory, "FCR_SSS_inflow_2013_2019_20200701_allfractions_2DOCpools.csv"))

  file.copy(file.path(config$file_path$data_directory, "fcre-manual-data/FCR_spillway_outflow_SUMMED_WeirWetland_2013_2019_20200615.csv"),
            file.path(config$file_path$execute_directory, "FCR_spillway_outflow_SUMMED_WeirWetland_2013_2019_20200615.csv"))
  }

  print("Creating FLARE Inflow files")

  inflow_forecast_path <- file.path(config$file_path$inflow_directory, config$inflow$forecast_inflow_model,config$location$site_id,lubridate::as_date(forecast_start_datetime),forecast_hour)

  #inflow_outflow_files <- FLAREr::create_glm_inflow_outflow_files(inflow_file_dir = inflow_forecast_path,
  #                                                                inflow_obs = cleaned_inflow_file,
  #                                                                working_directory = config$file_path$execute_directory,
  #                                                                config = config,
  #                                                                state_names = states_config$state_names)
  #inflow_file_names <- inflow_outflow_files$inflow_file_name
  #outflow_file_names <- inflow_outflow_files$outflow_file_name

  if(config$model_settings$model_name == "glm_aed"){

  file1 <- file.path(config$file_path$execute_directory, "FCR_weir_inflow_2013_2019_20200828_allfractions_2poolsDOC.csv")
  file2 <- file.path(config$file_path$execute_directory, "FCR_wetland_inflow_2013_2019_20200828_allfractions_2DOCpools.csv")
  inflow_file_names <- tibble(file1 = file1,
                              file2 = file2,
                              file3 = "sss_inflow.csv")
  outflow_file_names <- tibble(file_1 = file.path(config$file_path$execute_directory, "FCR_spillway_outflow_SUMMED_WeirWetland_2013_2019_20200615.csv"),
                               file_2 = "sss_outflow.csv")

  management <- FLAREr::generate_oxygen_management(config = config)
  }else{
    management <- NULL
  }

  #Create observation matrix
  print("Creating Observation Matrix")
  obs <- FLAREr::create_obs_matrix(cleaned_observations_file_long,
                                   obs_config,
                                   config)

  #full_time_forecast <- seq(start_datetime, end_datetime, by = "1 day")
  #obs[,2:dim(obs)[2], ] <- NA

  if(i > 1){
    obs[3:11,, which(seq(1,28,1) != 6)] <- NA
  }else{
    obs[3:11, 2:dim(obs)[2], which(seq(1,28,1) != 6)] <- NA
  }

  states_config <- FLAREr::generate_states_to_obs_mapping(states_config, obs_config)

  model_sd <- FLAREr::initiate_model_error(config, states_config)

  init <- FLAREr::generate_initial_conditions(states_config,
                                              obs_config,
                                              pars_config,
                                              obs,
                                              config,
                                              restart_file = config$run_config$restart_file,
                                              historical_met_error = met_out$historical_met_error)

  #Run EnKF
  print("Starting EnKF")
  print("-----------------------------------")
  da_forecast_output <- FLAREr::run_da_forecast(states_init = init$states,
                                       pars_init = init$pars,
                                       aux_states_init = init$aux_states_init,
                                       obs = obs,
                                       obs_sd = obs_config$obs_sd,
                                       model_sd = model_sd,
                                       working_directory = config$file_path$execute_directory,
                                       met_file_names = met_out$filenames,
                                       inflow_file_names = inflow_file_names,
                                       outflow_file_names = outflow_file_names,
                                       config = config,
                                       pars_config = pars_config,
                                       states_config = states_config,
                                       obs_config = obs_config,
                                       management,
                                       da_method = config$da_setup$da_method,
                                       par_fit_method = config$da_setup$par_fit_method)

  print("Writing output file")
  saved_file <- FLAREr::write_forecast_netcdf(da_forecast_output = da_forecast_output,
                                              forecast_output_directory = config$file_path$forecast_output_directory)

  #Create EML Metadata
  print("Creating metadata")
  FLAREr::create_flare_metadata(file_name = saved_file,
                                da_forecast_output = da_forecast_output)

  unlist(config$file_path$execute_directory, recursive = TRUE)

  rm(da_output)
  gc()

  #file_name <- saved_file
  print("Generating plot")
  FLAREr::plotting_general(file_name = saved_file,
                           qaqc_data_directory = config$file_path$qaqc_data_directory,
                           ncore = config$model_settings$ncore,
                           plot_profile = FALSE,
                           obs_csv = FALSE)

}


