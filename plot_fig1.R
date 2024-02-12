
## Plotting VS et al. response ##

## Packages ----

library(tidyverse)
library(lubridate)
library(gridExtra)
library(nls.multstart)
library(rTPC)

## VS et al. 2021 ----

## weather station data --

# read weather station data (wsd)
sweden_wsd <- read.delim("data/weather_station_temperatures.txt")

# modify the format of the date.time column and extract hour & day 
sweden_wsd$date.time <- ymd_hms(sweden_wsd$date.time)
sweden_wsd$hour <- hour(sweden_wsd$date.time)
sweden_wsd$day <- day(sweden_wsd$date.time)
sweden_wsd$month <- month(sweden_wsd$date.time)

# determine if an observation was made during daytime 
sweden_wsd$daytime <- ifelse(sweden_wsd$hour %in% c(5:20), 1, 0)

# filter for only daytime observations 
sweden_wsd <- sweden_wsd %>% filter(daytime == 1)

# group by day and average daytime temperatures
sweden_wsd <- sweden_wsd %>% group_by(month, day) %>% summarise(temp = mean(temp)) %>% ungroup()

## thermal performance curve data for eggs --

# read raw egg TPC data
sweden_eggs <- read.delim("data/eggs.txt")

# define function for the Ratkowsky model
ratkowsky <- function(temp, a , b, tmin, tmax){
  
  perf <- ifelse(temp >= tmax|temp <= tmin, 0 ,((a * (temp - tmin)) * (1 - exp(b * (temp - tmax))))^2)
  
  return(perf)
}

# define sequence of temperatures from which to predict Ratkowsky function
temps <- seq(0,47, by = 0.1)

# generate TPC for eggs using Ratkowsky model and parameters estimated by VS et al. 2021
sweden_tpc_eggs <- ratkowsky(temps, 0.023, 0.39, 0.9, 36)

# build egg TPC data
sweden_tpc_eggs <- data.frame(temp = temps, performance = sweden_tpc_eggs)

## thermal performance curve data for larvae --

# read raw larvae tpc data
sweden_larvae <- read.delim("data/larvae.txt")

# generate TPC for larvae using Ratkowsky model and parameters estimated by VS et al. 2021
sweden_tpc_larvae <- ratkowsky(temps, 0.012, 0.31, 0.6, 37)

# build larvae TPC data
sweden_tpc_larvae <- data.frame(temp = temps, performance = sweden_tpc_larvae)

## combine egg and larvae data ---

# add variables to indicate life stage on raw dev.rate data 
sweden_eggs <- sweden_eggs %>% select(temp, dev.rate) %>% mutate(stage = rep("Egg", nrow(.)))
sweden_larvae <- sweden_larvae %>% select(temp, dev.rate) %>% mutate(stage = rep("Larvae", nrow(.)))

# combine raw dev.rate data 
sweden_dev.rate <- rbind(sweden_eggs, sweden_larvae)

# add variables to indicate life stage on TPC data 
sweden_tpc_eggs <- sweden_tpc_eggs %>% mutate(stage = rep("Egg", nrow(.)))
sweden_tpc_larvae <- sweden_tpc_larvae %>% mutate(stage = rep("Larvae", nrow(.)))

# combine TPC data 
sweden_tpc <- rbind(sweden_tpc_eggs, sweden_tpc_larvae)

## plotting --

sweden_plot <- ggplot() +
  geom_line(data = sweden_tpc, aes(x = temp, y = performance, col = stage), lwd = 1) +
  geom_point(data = sweden_dev.rate, aes(x = temp, y = dev.rate, col = stage, fill = stage), shape = 21, size = 3, alpha = 0.5) +
  scale_color_manual(values = c("orange", "forestgreen")) +
  scale_fill_manual(values = c("orange", "forestgreen")) +
  geom_density(data = sweden_wsd, aes(x = temp, y = 1.5 * ..density..), fill = "lightblue", col = "lightblue", lwd = 1, alpha = 0.5) +
  ylab(expression(Developmental~rate~(days^{-1}))) +
  xlab("Temperature (°C)") +
  geom_text(aes(x = 2, y = 0.49, label = "C)"), size = 5) +
  theme_classic() +
  theme(legend.position = c(0.9, 0.9), 
        legend.title = element_blank()) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.52)) +
  scale_x_continuous(expand = c(0,0))

## Logan et al. 2020 & Neel et al. 2020 - Panama ----

## weather station data 

# read weather station data 
panama_wsd <- read.csv("data/monthly_weather_station_data_panama.csv")

# change formatting of the date column
panama_wsd$date <- mdy(panama_wsd$date)

# filter for data from the same months as VS et al. 2021
panama_wsd <- panama_wsd %>% filter(month %in% c("June", "July", "August")) 

# group by date and get mean daily temperature
panama_wsd <- panama_wsd %>% group_by(date) %>% summarise(temp = mean(temp))

## Anolis apletophallus TPC data ---

# read apleto data --
apleto_sprint_speed <- read.csv("data/apletophallus_sprint_speed.csv")

# fit Ratkowsky model to apleto data --

# get start values to fit the model
apleto_start_vals <- get_start_vals(apleto_sprint_speed$Temp, apleto_sprint_speed$Speed, model_name = "ratkowsky_1983")

# fit model
apleto_model <- nls.multstart::nls_multstart(Speed ~ ratkowsky_1983(temp = Temp, tmin, tmax, a, b),
                                             data = apleto_sprint_speed,
                                             iter = c(4,4,4,4),
                                             start_lower = apleto_start_vals - 10,
                                             start_upper = apleto_start_vals + 10,
                                             lower = get_lower_lims(apleto_sprint_speed$Temp, apleto_sprint_speed$Speed, model_name = 'ratkowsky_1983'),
                                             upper = get_upper_lims(apleto_sprint_speed$Temp, apleto_sprint_speed$Speed, model_name = 'ratkowsky_1983'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE)

# generate sprint speed prediction 
apleto_sprint_pred <- predict(apleto_model, newdata = data.frame(Temp = temps))

# generate TPC data 
apleto_tpc <- data.frame(temp = temps, speed = apleto_sprint_pred)

# apply correction as VS et al. does
apleto_tpc$speed <- ifelse(apleto_tpc$temp > 10.71, apleto_tpc$speed, 0)
apleto_tpc$speed <- ifelse(apleto_tpc$temp < 35.93, apleto_tpc$speed, 0)

## plotting --

panama_plot <- ggplot() +
  geom_line(data = apleto_tpc, aes(x = temp, y = speed), col = "gold1", lwd = 1) +
  geom_jitter(data = apleto_sprint_speed, aes(x = Temp, y = Speed), 
              width = 0.25, col = "gold4", fill = "gold1", shape = 21, size = 3, alpha = 0.5) +
  geom_density(data = panama_wsd, aes(x = temp), col = "lightblue", fill = "lightblue", lwd = 1, alpha = 0.5) +
  geom_text(aes(x = 2, y = 0.57, label = "A)"), size = 5) +
  ylab("Sprint speed (m/s)") +
  xlab("Temperature (°C)") +
  theme_classic() +
  theme(axis.title.x = element_blank()) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.6)) +
  scale_x_continuous(expand = c(0,0))

## Logan et al. 2020 & Neel et al. 2020 - Bahamas ----

## weather station data 

# read weather station data 
bahamas_wsd <- read.csv("data/monthly_weather_station_data_bahamas.csv")

# change formatitng of the date column
bahamas_wsd$Date <- mdy_hm(bahamas_wsd$Date)

# filter for data from the same months as VS et al. 
bahamas_wsd <- bahamas_wsd %>% filter(month %in% c("June", "July", "August")) 

# extracting the hour from the complex date column 
bahamas_wsd$hour <- hour(bahamas_wsd$Date)

# extract day 
bahamas_wsd$day <- day(bahamas_wsd$Date)

# determine day-time temperatures 
bahamas_wsd$daytime <- ifelse(bahamas_wsd$hour %in% c(6:20), 1, 0)

# group by date and get mean daily temperature
bahamas_wsd <- bahamas_wsd %>% filter(daytime == 1) %>%
  group_by(month, day) %>% summarise(temp = mean(Temp)) %>% ungroup()

## anolis apletophallus TPC data ---

# read apleto data --
sagrei_sprint_speed <- read.csv("data/sagrei_sprint_speed.csv")

# fit ratkowsky model to apleto data --

# get start values to fit the model
sagrei_start_vals <- get_start_vals(sagrei_sprint_speed$Temp, sagrei_sprint_speed$Speed, model_name = "ratkowsky_1983")

# fit model
sagrei_model <- nls.multstart::nls_multstart(Speed ~ ratkowsky_1983(temp = Temp, tmin, tmax, a, b),
                                             data = sagrei_sprint_speed,
                                             iter = c(4,4,4,4),
                                             start_lower = sagrei_start_vals - 10,
                                             start_upper = sagrei_start_vals + 10,
                                             lower = get_lower_lims(sagrei_sprint_speed$Temp, sagrei_sprint_speed$Speed, model_name = 'ratkowsky_1983'),
                                             upper = get_upper_lims(sagrei_sprint_speed$Temp, sagrei_sprint_speed$Speed, model_name = 'ratkowsky_1983'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE)

# generate sprint speed prediction 
sagrei_sprint_pred <- predict(sagrei_model, newdata = data.frame(Temp = temps))

# generate TPC data 
sagrei_tpc <- data.frame(temp = temps, speed = sagrei_sprint_pred)

# apply correction as VS et al. does
sagrei_tpc$speed <- ifelse(sagrei_tpc$temp > 0.40, sagrei_tpc$speed, 0)
sagrei_tpc$speed <- ifelse(sagrei_tpc$temp < 46.86, sagrei_tpc$speed, 0)

## ploting --

bahamas_plot <- ggplot() +
  geom_line(data = sagrei_tpc, aes(x = temp, y = speed), col = "darkgoldenrod4", lwd = 1) +
  geom_jitter(data = sagrei_sprint_speed, aes(x = Temp, y = Speed), 
              width = 0.25, col =  "darkgoldenrod4",  fill = "darkgoldenrod4", shape = 21, size = 3, alpha = 0.5) +
  geom_density(data = bahamas_wsd, aes(x = temp, y = 1.5 * ..density..), col = "lightblue", fill = "lightblue", lwd = 1, alpha = 0.5) +
  geom_text(aes(x = 2, y = 1.25, label = "B)"), size = 5) +
  ylab("Sprint speed (m/s)") +
  xlab("Temperature (°C)") +
  theme_classic() +
  theme(axis.title.x = element_blank()) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.3), breaks = seq(0,1.2,by = 0.2 )) +
  scale_x_continuous(expand = c(0,0))

## Combine all plots ----

grid.arrange(panama_plot, bahamas_plot, sweden_plot, ncol = 1)

