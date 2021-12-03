# time series project with air pollution data.
library(tidyverse)
library(ggplot2)
library(lubridate)
library(readxl)
library(CADFtest)
library(vars)
library(imputeTS)
library(fGarch)
library(kableExtra)
library(forecast)

# note: will need to go in seasonal differences for multivariate, but i won't be able to 
# do sarima in the multivariate package vars

require(plyr)
#setwd("d:/asus_documents/ku_leuven/courses/sustainability/project/project_resources/data/time_series_brux_8")
# new data scraping:
setwd("d:/linux_documents_11_2021/thesis/code/multi-modal-pollution/data/eea_air/time_series_brux_8")
air_df <- ldply(list.files(), read.csv, header=T)
# this package does not play nice with dplyr, so unload once done.
detach("package:plyr", unload=T)

require(dplyr)
# check memory
disk_usage <- object.size(air_df)
print(disk_usage, units="MB")
dim(air_df)
names(air_df)
unique(air_df$AirQualityStation)

# get highest reading station
air_summary <- air_df %>% group_by(AirQualityStation) %>% 
  summarise(avg_conc=mean(Concentration, na.rm=T), 
            sd_conc=sd(Concentration, na.rm=T), 
            length=n(), 
            min_date=min(DatetimeBegin), 
            max_date=max(DatetimeBegin)) %>% 
  arrange(desc(length), desc(avg_conc))
  #arrange(desc(avg_conc))
air_summary
# sta-betb001 has the most observations, and a relatively high avg, 41.8
# sta-benat39 has the highest average, 44.4
 
# get the stations with the most observations 
highest_station <- air_df %>% dplyr::filter(AirQualityStation %in% air_summary$AirQualityStation[1])
highest_station
nrow(highest_station)

# next fill in missing dates and run imputation
# then select the date range that i want.
highest_station$datetime <- ymd_hms(highest_station$DatetimeBegin)
min(highest_station$datetime)
max(highest_station$datetime)
# need to do le and ge because of how the == is computed.
#station_subset <- highest_station %>% filter(datetime >= "2019-01-01" & datetime <= "2021-09-22")
station_subset <- highest_station %>% dplyr::filter(datetime >= "2020-03-01" & datetime <= "2021-01-01")
#station_subset <- highest_station %>% filter(datetime >= "2020-01-01" & datetime <= "2020-02-01")
# check that the times are corerct
min(station_subset$DatetimeBegin)
max(station_subset$DatetimeBegin)
nrow(station_subset)
# missing some time measurements
station_complete <- station_subset %>% tidyr::complete(datetime = seq(min(datetime), max(datetime), by="hour"))
nrow(station_complete)
# this now contains a complete time series for 2019-01-01 through 2021-09-21
# because of creating new time indices, it is necessary to do imputation.
# using imputets, for the univariate case.
# todo: look at what this does more corefully.
# https://cran.r-project.org/web/packages/imputets/vignettes/imputets-time-series-missing-value-imputation-in-r.pdf
statsna(as.matrix(station_complete$Concentration))
# run imputations
# choice between linear, spline, and stineman interpolations
# need to round for some reason.
station_complete$concentration_inter <- round(na_interpolation(station_complete$Concentration),1)
ggplot_na_imputations(station_complete$Concentration, station_complete$concentration_inter)
statsNA(as.matrix(station_complete$concentration_inter))

# aggregate by day, but need to round up if staying in 2019
air_daily <- station_complete %>% mutate(day=lubridate::ceiling_date(datetime, "day")) %>% group_by(day) %>% 
  summarise(daily_concentration=mean(concentration_inter)) 
air_daily
nrow(air_daily)
min(air_daily$day)
max(air_daily$day)
# create the ts
air_ts <- ts(air_daily$daily_concentration, frequency=7, start=c(2020, 9))

png(filename="d:/asus_documents/ku_leuven/courses/time_series/project/figures/air_ts.png", width=1000)
#ts.plot(air_ts)
p <- ggplot(air_daily, aes(x=day, y=daily_concentration)) +
    geom_line() + 
    ggtitle("daily average no2 in 2020, brussels measuring station sta-betb001") +
    xlab("2020") +
    ylab("daily concentration")
p
dev.off()

# using the xts library instead
# but might not be compatible with at least the cadftest
#time_index <- seq(from = min(station_complete$datetime), 
#                  to = max(station_complete$datetime), by = "hour")
#length(time_index)
#air_ts = xts(station_complete$concentration_inter, order.by = time_index)

# check lags
png(filename="d:/asus_documents/ku_leuven/courses/time_series/project/figures/acf_air.png")
acf(air_ts)
dev.off()
acf(diff(air_ts))
png(filename="d:/asus_documents/ku_leuven/courses/time_series/project/figures/pacf_air.png")
pacf(air_ts)
dev.off()
pacf(diff(air_ts))
acf(diff(air_ts,7))
pacf(diff(air_ts,7))
# notice that there is seasonality.
# take differences and seasonal differences (by day, but has to be in frequency of ts)
air_2_diff <- diff(diff(air_ts), 7)
length(air_2_diff)
air_max_lag <- round(sqrt(length(air_ts)))
air_max_lag
png(filename="d:/asus_documents/ku_leuven/courses/time_series/project/figures/acf_2d_air.png")
acf(air_2_diff)
dev.off()
png(filename="d:/asus_documents/ku_leuven/courses/time_series/project/figures/pacf_2d_air.png")
pacf(air_2_diff)
dev.off()
Box.test(air_ts, lag = air_max_lag, type = "Ljung-Box")
Box.test(diff(air_ts), lag = air_max_lag, type = "Ljung-Box")
Box.test(air_2_diff, lag = air_max_lag, type = "Ljung-Box")
# rejects that is white noise

# test ts for non stationarity
# don't need type of "trend"
# reject non-stationary 
CADFtest(air_ts, type= "drift", criterion= "BIC", max.lag.y=air_max_lag)
CADFtest(diff(air_ts), type= "drift", criterion= "BIC", max.lag.y=air_max_lag)
CADFtest(air_2_diff, type= "drift", criterion= "BIC", max.lag.y=air_max_lag)

# this is for the aggregated daily series
fit1 <- arima(air_ts,order=c(1,1,1), seasonal=c(1,1,1))
fit2 <- arima(air_ts,order=c(2,1,1), seasonal=c(1,1,1)) 
fit3 <- arima(air_ts,order=c(3,1,1), seasonal=c(1,1,1)) 
fit4 <- arima(air_ts,order=c(2,1,1), seasonal=c(0,1,1)) 
fit5 <- arima(air_ts,order=c(2,1,1), seasonal=c(1,1,0)) 
fit6 <- arima(air_ts,order=c(3,1,2), seasonal=c(1,1,0)) 
fit7 <- arima(air_ts,order=c(3,1,1), seasonal=c(2,1,0)) 
fit8 <- arima(air_ts,order=c(2,1,2), seasonal=c(2,1,0))
fit9 <- arima(air_ts,order=c(2,1,1), seasonal=c(2,1,1))
fit10 <- arima(air_ts,order=c(2,1,1), seasonal=c(2,1,2))
fit11 <- arima(air_ts,order=c(3,1,1), seasonal=c(2,1,1))
fit12 <- arima(air_ts,order=c(3,1,0), seasonal=c(2,1,1))
fit13 <- arima(air_ts,order=c(1,0,0)) 
fit14 <- arima(air_ts,order=c(0,0,1)) 
fit15 <- arima(air_ts,order=c(1,0,1)) 
fit16 <- arima(air_ts,order=c(2,0,1)) 
fit16 <- arima(air_ts,order=c(1,1,0)) 
fit17 <- arima(air_ts,order=c(1,1,1)) 
fit18 <- arima(air_ts,order=c(2,1,0)) 
fit19 <- arima(air_ts,order=c(2,1,1)) 
fit20 <- arima(air_ts,order=c(3,1,0)) 
fit21 <- arima(air_ts,order=c(3,1,1)) 
fit22 <- arima(air_ts,order=c(3,1,2)) 
fit23 <- arima(air_ts,order=c(1,0,0), seasonal=c(1,1,1))
fit24 <- arima(air_ts,order=c(1,0,0), seasonal=c(1,1,0))
fit25 <- arima(air_ts,order=c(1,0,0), seasonal=c(0,1,1))
fit26 <- arima(air_ts,order=c(1,0,1), seasonal=c(1,1,1))
fit27 <- arima(air_ts,order=c(1,0,1), seasonal=c(2,1,1))
fit28 <- arima(air_ts,order=c(2,0,1), seasonal=c(1,1,1))
fit29 <- arima(air_ts,order=c(2,0,2), seasonal=c(1,1,1))
fit30 <- arima(air_ts,order=c(3,0,1), seasonal=c(1,1,1))

white_noise <- c()
bp <- c()
for (i in 1:30){
    p_val <- Box.test(eval(as.name(paste("fit", i, sep="")))$residuals, lag = air_max_lag, type = "Ljung-Box")$p.value
    if (p_val >= .05){
        white_noise <- c(white_noise, i)
        bp <- c(bp, p_val)
    }
}
model_comparison <- matrix(nrow=length(white_noise), ncol=4)
model_comparison[,1] <- white_noise
model_comparison[,2] <- bp
for (i in 1:length(white_noise)){
    model <- white_noise[i]
    m_call <- eval(as.name(paste("fit", model, sep="")))$call
    print(paste("model", model, ":", "order:", m_call[3], "seasonal:", m_call[4]))
    model_comparison[i, 3] <- AIC(eval(as.name(paste("fit", model, sep=""))))
    model_comparison[i, 4] <- AIC(eval(as.name(paste("fit", model, sep=""))), k=log(length(air_ts)))
}
model_comparison_df <- as.data.frame(model_comparison)
colnames(model_comparison_df) <- c("Model", "Box", "AIC", "SIC")
model_comparison_df
kbl(model_comparison_df, booktabs = T, format="latex")
# this generates all the in-sample statistics

acf(fit1$residuals)
pacf(fit1$residuals)
acf(fit2$residuals)
pacf(fit2$residuals)
acf(fit3$residuals)
pacf(fit3$residuals)
acf(fit4$residuals)
pacf(fit4$residuals)
acf(fit8$residuals)
pacf(fit8$residuals)
acf(fit9$residuals)
pacf(fit9$residuals)
acf(fit10$residuals)
pacf(fit10$residuals)
acf(fit11$residuals)
pacf(fit11$residuals)
acf(fit12$residuals)
pacf(fit12$residuals)
acf(fit23$residuals)
pacf(fit23$residuals)
acf(fit25$residuals)
pacf(fit25$residuals)
acf(fit29$residuals)
pacf(fit29$residuals)
acf(fit31$residuals)
pacf(fit31$residuals)
acf(fit32$residuals)
pacf(fit32$residuals)
acf(fit34$residuals)
pacf(fit34$residuals)

# testing models for garch effects.
acf(fit1$residuals^2)
pacf(fit1$residuals^2)
garch_eff <- c()
p_vals <- c()
for (i in 1:length(white_noise)){
    model_i <- white_noise[i]
    p_val <- Box.test(eval(as.name(paste("fit", model_i, sep="")))$residuals^2, lag = air_max_lag, type = "Ljung-Box")$p.value
    print(p_val)
    if (p_val >= .05){
        garch_eff <- c(garch_eff, model_i)
        p_vals <- c(pvals, p_val)
    }
}
garch_eff
p_vals
# the squared residuals are all rejected for being white noise.
# this indicates that we would benefit from a garch model.
fit_garch0<-garchFit(~arma(1,1)+garch(1,1),data=air_ts)
summary(fit_garch0)
plot(fit_garch0)
# use ql to deal with non normality.
fit_garch1<-garchFit(~arma(1,1)+garch(1,1),data=air_ts, cond.dist = "QMLE")

# out of sample computations.
# compare error with expanding window forecasts
# Get each model order and season, to fit all of them.
params <- list()
for (i in 1:length(white_noise)){
    model_i <- white_noise[i]
    model <- eval(as.name(paste("fit", model_i, sep="")))
    order_call <- model$call[[3]]
    se_call <- model$call[[4]]
    params[[length(params)+1]] <- list(order=c(order_call[[2]], order_call[[3]], order_call[[4]]), 
                                       seasonal=c(se_call[[2]], se_call[[3]], se_call[[4]]))
}
    
H <- c(1, 10, 30)
s <- round(0.75*length(air_ts))
s
out_of_sample_error <- matrix(nrow=length(white_noise), ncol=7)
out_of_sample_error
# iterate through the models whose residuals were white noise. 
# for each model, fit for 3 h values
for (i in 1:length(params)){
    print(i)
    for (j in 1:length(H)){
        print(j)
        if (j==1){
            error_index <- 2
        }
        else if (j==2){
            error_index <- 4
        }
        else{
            error_index <- 6
        }
        error.h<-c()
        h <- H[j]
        # Need to leave h off the end of the time series to be able to predict h ahead.
        for (k in s:(length(air_ts)-h)){
          # train on expanding window
          mymodel.sub<-arima(air_ts[1:k], order=params[[i]]$order, seasonal=params[[i]]$seasonal)
          # predict h step ahead
          predict.h<-predict(mymodel.sub,n.ahead=h)$pred[h]
          # subtract the predicted value predict.h from the observed value, y[i+h]
          error.h<-c(error.h,air_ts[k+h]-predict.h)
        }
        mae <- mean(abs(error.h))
        # Divide by all y_t
        mape <- mean(abs(error.h)/air_ts[(s+h):length(air_ts)])
        rmse <- sqrt(mean(error.h^2))
        out_of_sample_error[i, 1] <- white_noise[i]
        #out_of_sample_error[i, error_index] <- mae
        out_of_sample_error[i, error_index] <- mape
        out_of_sample_error[i, error_index+1] <- rmse
    }
}

oos_df <- as.data.frame(out_of_sample_error) 
oos_df
colnames(oos_df) <- c("Model", "MAPE, h=1", "RMSE, h=1", 
                      "MAPE, h=10", "RMSE, h=10", 
                      "MAPE, h=30", "RMSE, h=30")
oos_df
kbl(oos_df, booktabs = T, format="latex")

# test difference model 23 and 25
# absolute and squared error
error_23<-c()
error_25<-c()
h <- 10
for (k in s:(length(air_ts)-h)){
          fit23_sub <- arima(air_ts[1:k],order=c(1,0,0), seasonal=c(1,1,1))
          fit25_sub <- arima(air_ts[1:k],order=c(1,0,0), seasonal=c(0,1,1))
          # predict h step ahead
          predict_23<-predict(fit23_sub,n.ahead=h)$pred[h]
          predict_25<-predict(fit25_sub,n.ahead=h)$pred[h]
          # subtract the predicted value predict.h from the observed value, y[i+h]
          error_23<-c(error_23,air_ts[k+h]-predict_23)
          error_25<-c(error_25,air_ts[k+h]-predict_25)
}
dm.test(error_23, error_25,h=h,power=1)
dm.test(error_23, error_25,h=h,power=2)

summary(fit23)
summary(fit25)

# Generate plot of forecast
# using the best model, fit23.
# I can get the actual data observed during this period (not included in the analysis above.)
station_observed <- highest_station %>% dplyr::filter(datetime >= "2021-01-02" & datetime <= "2021-02-01")
# missing some time measurements
station_obs_complete <- station_observed %>% tidyr::complete(datetime = seq(min(datetime), max(datetime), by="hour"))
station_obs_complete$concentration_inter <- round(na_interpolation(station_obs_complete$Concentration),1)
statsNA(as.matrix(station_obs_complete$concentration_inter))

# aggregate by day, but need to round up if staying in 2019
obs_daily <- station_obs_complete %>% mutate(day=lubridate::ceiling_date(datetime, "day")) %>% group_by(day) %>% 
  summarise(daily_concentration=mean(concentration_inter)) 
obs_daily$daily_concentration
min(obs_daily$day)
max(obs_daily$day)
nrow(obs_daily)
# create the ts
obs_ts <- ts(obs_daily$daily_concentration, frequency=7, start=c(2021, 1))
length(obs_ts)
ts.plot(obs_ts)

fcast <- predict(fit23,n.ahead=31)
expected <- fcast$pred
lower<-fcast$pred-qnorm(0.975)*fcast$se
upper<-fcast$pred+qnorm(0.975)*fcast$se
forecast_df <- data.frame(pred=expected, day=obs_daily$day, lower=lower, upper=upper)

png(filename="d:/asus_documents/ku_leuven/courses/time_series/project/figures/univariate_forecast.png")
colors <- c("Observed"="green", "Predicted"="red", "lower"="black", "upper"="black")
p <- ggplot() +
    geom_line(data=obs_daily, aes(x=day, y=daily_concentration, color="Observed")) + 
    geom_line(data=forecast_df, aes(x=day, y=pred, color="Predicted")) + 
    geom_line(data=forecast_df, aes(x=day, y=lower, color="lower")) + 
    geom_line(data=forecast_df, aes(x=day, y=upper, color="upper")) + 
    ggtitle("Predicted and observed concentration for January-February, 2021") +
    labs(y = "Daily. conc", x = "Day",
         color = "Legend") +
    scale_color_manual(values = colors)
p
dev.off()


########################
# multivariate analysis
########################


# https://www.ifo.de/docdl/cesifo1_wp1939.pdf
# https://aip.scitation.org/doi/pdf/10.1063/1.5016666
# load covid data.
covid_df <- read_excel("d:/asus_documents/ku_leuven/courses/time_series/project/covid19be.xlsx")
head(covid_df)
names(covid_df) <- tolower(names(covid_df))
head(covid_df)
# aggregate the classes.
covid_df %>% group_by()
covid_flanders <- covid_df %>% dplyr::filter(region=="Flanders") %>% 
    mutate(day=lubridate::ceiling_date(date, "day")) %>% 
    group_by(day) %>% summarise(total_cases=sum(cases))
covid_flanders
covid_subset <- covid_flanders %>% dplyr::filter(day > "2020-03-01" & day <= "2021-01-02")
covid_subset
nrow(covid_subset)
min(covid_subset$day)
max(covid_subset$day)
# there are no missing dates or nas, so the following code is not necessary.
#covid_complete <- covid_subset %>% complete(day = seq(min(day), max(day), by="day"))
#nrow(covid_complete)
## finish with imputation
#statsna(as.matrix(covid_complete$total_cases))
## run imputations
## choice between linear, spline, and stineman interpolations
## need to round for some reason.
#covid_complete$cases_inter <- round(na_interpolation(covid_complete$total_cases),1)
#ggplot_na_imputations(covid_complete$total_cases, covid_complete$cases_inter)
#statsna(as.matrix(covid_complete$cases_inter))

# create covid ts
cov_ts <- ts(covid_subset$total_cases, frequency=7, start=c(2020, 9))
#plot.ts(cov_ts)
png(filename="d:/asus_documents/ku_leuven/courses/time_series/project/figures/cov_ts.png", width=1000)
p <- ggplot(covid_subset, aes(x=day, y=total_cases)) +
    geom_line() + 
    ggtitle("Daily COVID-19 cases, Flanders") +
    xlab("2020") +
    ylab("Daily cases")
p
dev.off()

plot.ts(cov_ts)

# test for non stationarity
cov_max_lag=round(sqrt(length(cov_ts)))
CADFtest(cov_ts, type= "drift", criterion= "BIC", max.lag.y=cov_max_lag)
CADFtest(diff(cov_ts), type= "drift", criterion= "BIC", max.lag.y=cov_max_lag)
# both are stationary

png(filename="d:/asus_documents/ku_leuven/courses/time_series/project/figures/cov_acf.png")
acf(cov_ts)
dev.off()
png(filename="d:/asus_documents/ku_leuven/courses/time_series/project/figures/cov_pacf.png")
pacf(cov_ts)
dev.off()
d_cov <- diff(cov_ts) 
acf(d_cov)
pacf(d_cov)
acf(diff(d_cov, 7))
pacf(diff(d_cov, 7))


# then load the traffic data.
traffic_df <- read_excel(
        path="d:/linux_documents_11_2021/thesis/code/multi-modal-pollution/data/traffic/20210921-tunnels2019-sept2021.xlsx"
        ) 
names(traffic_df) <- tolower(names(traffic_df))
head(traffic_df)
nrow(traffic_df)
# first need to select the detector to use.
traffic_summary <- traffic_df %>% dplyr::group_by(detector) %>% 
  summarise(avg_tf=mean(volume, na.rm=T), sd_tf=sd(volume, na.rm=T), len=n()) %>% 
  arrange(desc(avg_tf))
traffic_summary
# using tun bel in because it has data for longest period with a high average
detector_df <- traffic_df %>% dplyr::filter(detector=="Tun Bel IN")
nrow(detector_df)
detector_df$from <- ymd_hms(detector_df$from)
min(detector_df$from)
max(detector_df$from)
head(detector_df)
# this detector starts at 2019-01-01 00:00:00
# day duration counter https://www.timeanddate.com/date/duration.html?d1=01&m1=01&y1=2019&d2=&m2=&y2=&ti=on&
detector_subset <- detector_df %>% dplyr::filter(from > "2020-03-01 00:00:00" & from <= "2021-01-01")
min(detector_subset$from)
max(detector_subset$from)
nrow(detector_subset)

# fill in missing dates
# https://blog.exploratory.io/populating-missing-dates-with-complete-and-fill-functions-in-r-and-exploratory-79f2a321e6b5
detector_complete <- detector_subset %>% complete(from = seq(min(from), max(from), by="hour"))
nrow(detector_complete)
detector_complete[1:3, "from"]

# next do imputation. variations options are available.
# https://cran.r-project.org/web/packages/amelia/index.html
# http://schd.ws/hosted_files/user2017/e5/user2017_steffen_moritz.pdf
# https://cran.r-project.org/web/packages/imputets/
# https://stats.stackexchange.com/questions/261271/imputation-methods-for-time-series-data
# https://stats.stackexchange.com/questions/261271/imputation-methods-for-time-series-data
# imputets is for univariate, it is probably ok to use for the sake of this
# project
statsNA(as.matrix(detector_complete$volume))
# run imputations
# choice between linear, spline, and stineman interpolations
# need to round for some reason.
detector_complete$volume_inter <- round(na_interpolation(detector_complete$volume),1)
ggplot_na_imputations(detector_complete$volume, detector_complete$volume_inter)
statsNA(as.matrix(detector_complete$volume_inter))

# aggregate by day, but need to round up if staying in 2019
traffic_daily <- detector_complete %>% mutate(day=lubridate::ceiling_date(from, "day")) %>% group_by(day) %>% 
  summarise(daily_volume=mean(volume_inter)) 
traffic_daily
nrow(traffic_daily)
min(traffic_daily$day)
max(traffic_daily$day)
# create the ts
traffic_ts <- ts(traffic_daily$daily_volume, frequency=7, start=c(2020, 9))
#plot.ts(traffic_ts)
png(filename="d:/asus_documents/ku_leuven/courses/time_series/project/figures/traffic_ts.png", width=1000)
p <- ggplot(traffic_daily, aes(x=day, y=daily_volume)) +
    geom_line() + 
    ggtitle("Daily traffic volume, Brussels Belliard Tunnel") +
    xlab("2020") +
    ylab("Daily volume")
p
dev.off()

# test for non stationarity
tf_max_lag=round(sqrt(length(traffic_ts)))
CADFtest(traffic_ts, type= "drift", criterion= "BIC", max.lag.y=tf_max_lag)
CADFtest(diff(traffic_ts), type= "drift", criterion= "BIC", max.lag.y=tf_max_lag)
# in differences we have stationarity. in levels we do not
# for the air time series, it was stationary for both.

acf(traffic_ts)
pacf(traffic_ts)
d_traffic <- diff(traffic_ts) 
acf(d_traffic)
pacf(d_traffic)
acf(diff(d_traffic, 7))
pacf(diff(d_traffic))
# remember to correct for residuals when testing for cointegration.

# Testing for cointegration
length(cov_ts)
length(air_ts)
length(diff(traffic_ts))
# Set up dfs
var_df_cov <- data.frame(air_ts, traffic_ts, cov_ts)
var_df_no_cov <- data.frame(air_ts, traffic_ts)
var_df_no_air <- data.frame(traffic_ts, cov_ts)
# lag of 9
names(var_df_cov)<-c("air","traffic", "covid")
names(var_df_no_cov)<-c("air","traffic")
names(var_df_no_air)<-c("traffic","covid")
VARselect(var_df_cov,lag.max=10,type="const")
# order 8
VARselect(var_df_no_cov,lag.max=10,type="const")
# by sic we select order 8
VARselect(var_df_no_air,lag.max=10,type="const")
# order 9

# Testing cointegration
# This first test is the most important, because we have to check if the two 
# non-stationary series are cointegrated.
trace1 <- ca.jo(var_df_no_air,type="trace",K=9,ecdet="const",spec="transitory")
summary(trace1)
# The result supposts cointegration with 1 test equation.
trace2 <- ca.jo(var_df_cov,type="trace",K=8,ecdet="const",spec="transitory")
summary(trace2)
# When including that stationary air series, we get 2 test equations, though the second
# is rather on the boundary b/w signficant and not
# This test below is not necessary:
trace3 <- ca.jo(var_df_no_cov,type="trace",K=8,ecdet="const",spec="transitory")
summary(trace3)

#repeat the same procedure using johansen’s maximum eigenvalue test statistic.
# const is constant term in the long run relationshp
# transitory means to not include a deterministic trend in the test equation.
# we interpret this test as we did for the prvious test
maxeigen_test <- ca.jo(var_df_no_air,type="eigen",K=9,ecdet="const",spec="transitory")
summary(maxeigen_test)
# we see that we reject no cointegration at the 1% level

#############
# VAR model
############
# Going in differences for the non-stationary series
# Need to remove first measurement from stationary series
var_df_cov <- data.frame(air_ts[2:length(air_ts)], diff(traffic_ts), diff(cov_ts))
var_df_no_cov <- data.frame(air_ts[2:length(air_ts)], diff(traffic_ts))
names(var_df_cov)<-c("air","traffic", "covid")
names(var_df_no_cov)<-c("air","traffic")
# there may be cases where using nonstationary series in var is ok.
fit_var1 <- VAR(var_df_no_cov,type="const",p=8)
fit_var2 <- VAR(var_df_cov,type="const",p=8)
summary(fit_var1)
summary(fit_var2)

# checking the residuals
var2_residuals<-resid(fit_var2)
acf(var2_residuals[,1])
acf(var2_residuals[,2])
acf(var2_residuals[,3])
ccf(var2_residuals[,1],var2_residuals[,2])
ccf(var2_residuals[,1],var2_residuals[,3])
ccf(var2_residuals[,2],var2_residuals[,3])

# we then examine the impulse function
irf_var<-irf(fit_var2,ortho=F,boot=T)
plot(irf_var)

# Generating predictions
H <- c(1, 5, 10)
s <- round(0.75*length(air_ts))
dm_tests <- matrix(nrow=length(H), ncol=3)
oos_errors <- matrix(nrow=3, ncol=5)
for (i in 1:length(H)){
    h <- H[i]
    print(h)
    error_cov <- c()
    error_nocov <- c()
    for (k in s:(length(air_ts)-h)){
        var_df_cov_sub <- data.frame(air_ts[1:k], diff(traffic_ts)[1:k], diff(cov_ts)[1:k])
        var_df_no_cov_sub <- data.frame(air_ts[1:k], diff(traffic_ts)[1:k])
        names(var_df_cov_sub)<-c("air","traffic", "covid")
        names(var_df_no_cov_sub)<-c("air","traffic")
        fit_var_cov <- VAR(var_df_cov_sub,type="const",p=8)
        fit_var_nocov <- VAR(var_df_no_cov_sub,type="const",p=8)
        # predict h step ahead
        predict_cov_out <- predict(fit_var_cov,n.ahead=h)
        predict_cov <- predict_cov_out[[1]][1]$air[h]
        predict_nocov_out <- predict(fit_var_nocov,n.ahead=h)
        predict_nocov <- predict_nocov_out[[1]][1]$air[h]
        # subtract the predicted value predict.h from the observed value, y[i+h]
        error_cov <- c(error_cov,air_ts[k+h]-predict_cov)
        error_nocov <- c(error_nocov,air_ts[k+h]-predict_nocov)
    }
    cov_mape <- mean(abs(error_cov)/air_ts[(s+h):length(air_ts)])
    nocov_mape <- mean(abs(error_nocov)/air_ts[(s+h):length(air_ts)])
    cov_rmse <- sqrt(mean(error_cov^2))
    nocov_rmse <- sqrt(mean(error_nocov^2))
    oos_errors[i, 1] <- h
    oos_errors[i, 2] <- cov_mape
    oos_errors[i, 3] <- cov_rmse
    oos_errors[i, 4] <- nocov_mape
    oos_errors[i, 5] <- nocov_rmse
    dm_tests[i, 1] <- h
    dm_tests[i, 2] <- dm.test(error_cov, error_nocov,h=h,power=1)$p.value[[1]]
    dm_tests[i, 3] <- dm.test(error_cov, error_nocov,h=h,power=2)$p.value[[1]]
}
oos_errors
dm_tests
oos_df <- as.data.frame(oos_errors)
dm_tests
dm_df <- as.data.frame(dm_tests)
colnames(oos_df) <- c("Horizon", "W/ COVID MAPE", "W/COVID RMSE", "W/out COVID MAPE", "W/out COVID RMSE")
colnames(dm_df) <- c("Horizon", "MAE", "MSE")
kbl(oos_df, booktabs = T, format="latex")
kbl(dm_df, booktabs = T, format="latex")

#########
# The analysis is then repeated for October to the end of the year, to account for the period
# of the most intense covid
# First traffic
detector_subset <- detector_df %>% dplyr::filter(from > "2020-10-01 00:00:00" & from <= "2021-01-01")
detector_complete <- detector_subset %>% complete(from = seq(min(from), max(from), by="hour"))
nrow(detector_complete)
detector_complete$volume_inter <- round(na_interpolation(detector_complete$volume),1)
# aggregate by day, but need to round up if staying in 2019
traffic_daily <- detector_complete %>% mutate(day=lubridate::ceiling_date(from, "day")) %>% group_by(day) %>% 
  summarise(daily_volume=mean(volume_inter)) 
traffic_daily
nrow(traffic_daily)
min(traffic_daily$day)
max(traffic_daily$day)
traffic_ts <- ts(traffic_daily$daily_volume, frequency=7, start=c(2020, 39))
# Covid
covid_subset <- covid_flanders %>% dplyr::filter(day > "2020-10-01" & day <= "2021-01-02")
covid_subset
nrow(covid_subset)
min(covid_subset$day)
max(covid_subset$day)
cov_ts <- ts(covid_subset$total_cases, frequency=7, start=c(2020, 39))
# And NO2
station_subset <- highest_station %>% dplyr::filter(datetime >= "2020-10-01" & datetime <= "2021-01-01")
station_complete <- station_subset %>% tidyr::complete(datetime = seq(min(datetime), max(datetime), by="hour"))
nrow(station_complete)
station_complete$concentration_inter <- round(na_interpolation(station_complete$Concentration),1)
air_daily <- station_complete %>% mutate(day=lubridate::ceiling_date(datetime, "day")) %>% group_by(day) %>% 
  summarise(daily_concentration=mean(concentration_inter)) 
air_daily
# create the ts
air_ts <- ts(air_daily$daily_concentration, frequency=7, start=c(2020, 39))

ts.plot(air_ts)
ts.plot(cov_ts)
ts.plot(traffic_ts)

# Then reselect lags and remodel
length(cov_ts)
length(air_ts)
length(traffic_ts)
# Set up dfs
var_df_cov <- data.frame(air_ts, traffic_ts, cov_ts)
var_df_no_cov <- data.frame(air_ts, traffic_ts)
names(var_df_cov)<-c("air","traffic", "covid")
names(var_df_no_cov)<-c("air","traffic")
VARselect(var_df_cov,lag.max=10,type="const")
VARselect(var_df_no_cov,lag.max=10,type="const")
# both still order 8

# REgenerating predictions
H <- c(1, 5, 10)
s <- round(0.75*length(air_ts))
dm_tests <- matrix(nrow=length(H), ncol=3)
oos_errors <- matrix(nrow=3, ncol=5)
for (i in 1:length(H)){
    h <- H[i]
    print(h)
    error_cov <- c()
    error_nocov <- c()
    for (k in s:(length(air_ts)-h)){
        var_df_cov_sub <- data.frame(air_ts[1:k], diff(traffic_ts)[1:k], diff(cov_ts)[1:k])
        var_df_no_cov_sub <- data.frame(air_ts[1:k], diff(traffic_ts)[1:k])
        names(var_df_cov_sub)<-c("air","traffic", "covid")
        names(var_df_no_cov_sub)<-c("air","traffic")
        fit_var_cov <- VAR(var_df_cov_sub,type="const",p=8)
        fit_var_nocov <- VAR(var_df_no_cov_sub,type="const",p=8)
        # predict h step ahead
        predict_cov_out <- predict(fit_var_cov,n.ahead=h)
        predict_cov <- predict_cov_out[[1]][1]$air[h]
        predict_nocov_out <- predict(fit_var_nocov,n.ahead=h)
        predict_nocov <- predict_nocov_out[[1]][1]$air[h]
        # subtract the predicted value predict.h from the observed value, y[i+h]
        error_cov <- c(error_cov,air_ts[k+h]-predict_cov)
        error_nocov <- c(error_nocov,air_ts[k+h]-predict_nocov)
    }
    cov_mape <- mean(abs(error_cov)/air_ts[(s+h):length(air_ts)])
    nocov_mape <- mean(abs(error_nocov)/air_ts[(s+h):length(air_ts)])
    cov_rmse <- sqrt(mean(error_cov^2))
    nocov_rmse <- sqrt(mean(error_nocov^2))
    oos_errors[i, 1] <- h
    oos_errors[i, 2] <- cov_mape
    oos_errors[i, 3] <- cov_rmse
    oos_errors[i, 4] <- nocov_mape
    oos_errors[i, 5] <- nocov_rmse
    dm_tests[i, 1] <- h
    dm_tests[i, 2] <- dm.test(error_cov, error_nocov,h=h,power=1)$p.value[[1]]
    dm_tests[i, 3] <- dm.test(error_cov, error_nocov,h=h,power=2)$p.value[[1]]
}
oos_errors
dm_tests
oos_df <- as.data.frame(oos_errors)
dm_tests
dm_df <- as.data.frame(dm_tests)
colnames(oos_df) <- c("Horizon", "W/ COVID MAPE", "W/COVID RMSE", "W/out COVID MAPE", "W/out COVID RMSE")
colnames(dm_df) <- c("Horizon", "MAE", "MSE")
kbl(oos_df, booktabs = T, format="latex")
kbl(dm_df, booktabs = T, format="latex")


####### ###
# NOT USING FOR REPORT
###########
# VECM
# see exercise 6
# k is from code above

# Test for cointegration with engel grangers
fit_ci <- lm(cov_ts~traffic_ts)
CADFtest(fit_ci$residuals, max.lag.y = cov_max_lag, criterion="BIC", type="drift")
# The test is rejected (indicating stationarity), therefore we have cointegration

# estimate a vector error correcting model, using the 3 time series (trace2 from above)
# we set r=2 because that is the result from the previous tests.
fit_vecm<-cajorls(trace2,r=2)
fit_vecm
# the cointegration vector will be:
# ect1          ect2
# air_ts.l1      1.000000e+00    0.00000000
# traffic_ts.l1 -1.734723e-18    1.00000000
# cov_ts.l1     -1.072716e-03   -0.09431364
# constant      -3.288128e+01 -959.71766568
# see about 1:50 in final lecture for explanation of interpretation of some coefficients.
# i'm still not sure exactly what is what.
# ect is the delta on slide 249.

# then we can do predition with the ecm
# Using observed air values from univariate analysis
vec_forecast <- predict(vec2var(trace_test,r=2),n.ahead=32)
air_forecast<-ts(myforecast$fcst$air[,1],frequency=1,start=c(2021,1))
air_lower<-ts(myforecast$fcst$air[,2],frequency=1,start=c(2021,1))
air_upper<-ts(myforecast$fcst$air[,3],frequency=1,start=c(2021,1))
length(air_forecast)
ts.plot(air_forecast,air_lower,air_upper, obs_ts, col=c("black","red","red", "blue"))

myforecast <- predict(vec2var(trace_test,r=2),n.ahead=30)
cov_forecast<-ts(myforecast$fcst$cov_ts[,1],frequency=1,start=c(2021,1))
cov_lower<-ts(myforecast$fcst$cov_ts[,2],frequency=1,start=c(2021,1))
cov_upper<-ts(myforecast$fcst$cov_ts[,3],frequency=1,start=c(2021,1))
ts.plot(cov_forecast,cov_lower,cov_upper,col=c("black","red","red"))





