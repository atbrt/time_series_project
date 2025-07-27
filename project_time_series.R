# Loading the data
datafile <- "valeurs_mensuelles.csv"

#delete the headers
data <- read.csv(datafile, sep=";", skip = 4, header = FALSE)
data <- data[, 1:2]

# Rename the columns for simpler use
colnames(data) <- c("date", "IPI")

#ordering the data using date
data <- data[order(data$date), ]



require(zoo)
dates <- as.yearmon(data$date)
IPI <- zoo(data$IPI, order.by=dates)

#Plot the series before cutting
plot(index(IPI), IPI, type="l", col="blue", lwd=2,
     xlab="Date", ylab="IPI", main="Industrial Production Index over time")

# Filter to select only data starting in January 2010 after 2008 financial crisis
IPI <- window(IPI, start = as.yearmon("Jan 2010"), end=as.yearmon("Mar 2025"))

filtered_dates <- index(IPI)

#Plotting the IPI and its differenciated series
dIPI <- diff(IPI, 1)

y_range <- range(c(as.numeric(IPI), as.numeric(dIPI)), na.rm = TRUE)
plot(index(IPI), IPI, type="l", col="blue", lwd=2,
     xlab="Date", ylab="Value", main="IPI and dIPI over time", 
     ylim = y_range)
lines(index(dIPI), dIPI, col="red", lwd=2)
legend("topright", legend=c("IPI", "dIPI"), col=c("blue", "red"), lwd=2, cex=0.6)


# Convert to numeric for regression to detect trends
dates_numeric <- as.numeric(filtered_dates)

#Regression to detect trends
summary(lm(IPI ~ dates_numeric))


require(fUnitRoots)
adf <- adfTest(IPI, lag=0, type="ct") #Using type=ct to include constant and linear trend

#First we need to test the residuals autocorrelation for lags until 24
Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}

Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))

#Find the minimum lag order without residual autocorrelation in ADF
series <- IPI; kmax <- 24; adftype="ct"
adfTest_valid <- function(series, kmax, adftype){
  k <- 0
  noautocorr <- 0
  while (noautocorr==0){
    cat(paste0("ADF with ",k," lags: residuals OK? "))
    adf <- adfTest(series, lags=k, type=adftype)
    pvals <- Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))[,2]
    if (sum(pvals<0.05,na.rm=T)==0) {
      noautocorr <- 1; cat("OK \n")
    } else cat("nope \n")
    k <- k+1
  }
  return(adf)
}
adf <- adfTest_valid(IPI,24,adftype="ct")

#Run the ADF test ==> not stationary, we try the differenciated series dIPI
adf


filtered_dates_diff <- index(dIPI)
dates_numeric_diff <- as.numeric(filtered_dates_diff)

#Run the LR to test trend ==> no significancy ==> trend corrected by differenciation
summary(lm(dIPI ~ dates_numeric_diff))

#This time, we run with type="nc", as the trend is corrected with the differentiation
adf <- adfTest_valid(dIPI,24,"nc")
Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))

#ADF test: stationary
adf


#Plot ACF and PACF on the stationary series
par(mfrow=c(1,2))
pacf(dIPI,24);acf(dIPI,24) 
library(lmtest)
pmax=5; qmax=2

####### WITHOUT COVID CORRECTION #############

#Function to estimate ARMA(p, q) models, and check adjustment and validity
modelchoice <- function(p,q,data=dIPI, k=24){
  estim <- try(arima(data, c(p, 0,q),optim.control=list(maxit=20000)))
  if (class(estim)=="try-error") return(c("p"=p,"q"=q,"arsignif"=NA,"masignif"=NA,"resnocorr"=NA, "ok"=NA))
  arsignif <- if (p==0) NA else coeftest(estim)[1:p, 4][p] <= 0.05
  masignif <- if (q==0) NA else coeftest(estim)[(p+1):(p+q), 4][q] <= 0.05
  resnocorr <- sum(Qtests(estim$residuals,24,length(estim$coef)-1)[,2]<=0.05,na.rm=T)==0
  checks <- c(arsignif,masignif,resnocorr)
  ok <- as.numeric(sum(checks,na.rm=T)==(3-sum(is.na(checks))))
  return(c("p"=p,"q"=q,"arsignif"=arsignif,"masignif"=masignif,"resnocorr"=resnocorr,"ok"=ok))
}


#Function to evaluate and validate all the arma(p, 0, q) for p<=pmax and q<=qmax on dIPI
armamodelchoice <- function(pmax,qmax){
  pqs <- expand.grid(0:pmax,0:qmax)
  t(apply(matrix(1:dim(pqs)[1]),1,function(row) {
    p <- pqs[row,1]; q <- pqs[row,2]
    cat(paste0("Computing ARMA(",p,",",q,") \n"))
    modelchoice(p,q)
  }))
}

armamodels <- armamodelchoice(pmax,qmax) 
colnames(armamodels) <- c("p", "q", "arsignif", "masignif", "resnocorr", "ok")
selec <- armamodels[armamodels[,"ok"]==1&!is.na(armamodels[,"ok"]),] 
#Selec contains all the valid models
selec

#Among the valid models, choose the one minimizing the AIC and BIC criteria
pqs <- apply(selec,1,function(row) list("p"=as.numeric(row[1]),"q"=as.numeric(row[2]))) 
names(pqs) <- paste0("arima(",selec[,1],", 0,",selec[,2],")") 
models <- lapply(pqs, function(pq) arima(dIPI, c(pq[["p"]], 0, pq[["q"]])))
vapply(models, FUN.VALUE=numeric(2), function(m) c("AIC"=AIC(m),"BIC"=BIC(m))) 

#Estimate the RMSE out of sample to evaluate the model
n_obs <- length(dIPI)
dIPI_train <- window(dIPI, end = index(dIPI)[n_obs-4])
dIPI_test <- tail(dIPI, 4)

models_train <- lapply(pqs, function(pq) arima(dIPI_train, c(pq[["p"]], 0, pq[["q"]])))

forecasts <- lapply(models_train, function(m) as.zoo(predict(m,4)$pred))

rmse <- vapply(forecasts, FUN.VALUE=numeric(1), function(forecast) {
  sqrt(mean((forecast - dIPI_test)^2))
})

rmse

#Best model according to the criteria: ARMA(1, 2)
selected_model_no_covid <- arima(dIPI, c(1, 0, 2))

library(tseries)

#Test of the normality of the residuals using Jarque Bera test
jb_test_no_covid <- jarque.bera.test(selected_model_no_covid$residuals)
print(jb_test_no_covid)


#Extract parameters
phi1 <- selected_model_no_covid$coef["ar1"]
theta1 <- selected_model_no_covid$coef["ma1"]
sigma_sq_no_covid <- selected_model_no_covid$sigma2

#Calculate a=φ₁ + θ₁, as in the report 
a <- phi1 + theta1

#Determine the covariance matrix
Sigma_no_covid <- sigma_sq_no_covid * matrix(c(1, a, a, 1 + a^2), nrow = 2)

#Estimate forecasts for dIPI without COVID correction
forecasts_no_covid <- predict(selected_model_no_covid, n.ahead = 2)
X_hat_no_covid <- forecasts_no_covid$pred  

par(mfrow = c(1,1))
library(ellipse)

#Generate ellipse points using the correct covariance matrix
ellipse_points_no_covid <- ellipse(Sigma_no_covid, centre = as.numeric(X_hat_no_covid), level = 0.95)

#Plot the 95% confidence region
plot(ellipse_points_no_covid, type = 'l', lwd = 2, col = 'red',
     xlab = expression(dIPI[T+1]), ylab = expression(dIPI[T+2]),
     main = "95% Confidence Region for dIPI",
     xlim = range(ellipse_points_no_covid[,1]) + c(-0.1, 0.1)*diff(range(ellipse_points_no_covid[,1])),
     ylim = range(ellipse_points_no_covid[,2]) + c(-0.1, 0.1)*diff(range(ellipse_points_no_covid[,2])))

#Add the forecast
points(X_hat_no_covid[1], X_hat_no_covid[2], pch = 19, col = 'blue', cex = 1.5)
text(X_hat_no_covid[1], X_hat_no_covid[2], "Point Forecast", pos = 3, col = 'blue', cex = 0.8)

grid()

#Add individual confidence intervals using correct standard errors
se_1_corrected <- sqrt(Sigma_no_covid[1,1])
se_2_corrected <- sqrt(Sigma_no_covid[2,2])

ci_1_no_covid <- X_hat_no_covid[1] + c(-1, 1) * qnorm(0.975) * se_1_corrected
ci_2_no_covid <- X_hat_no_covid[2] + c(-1, 1) * qnorm(0.975) * se_2_corrected

abline(v = ci_1_no_covid, lty = 2, col = 'gray')
abline(h = ci_2_no_covid, lty = 2, col = 'gray')

legend("topright", 
       legend = c("95% Joint CI", "Point Forecast", "Individual 95% CIs"),
       col = c("red", "blue", "gray"),
       lty = c(1, NA, 2),
       pch = c(NA, 19, NA),
       lwd = c(2, NA, 1),
       cex = 0.6)



######### WITH COVID (using dummy variables) ########

#Creation of dummy variables for March and April 2020, where covid effect is strong
covid_dates <- c(as.yearmon("Mar 2020"), as.yearmon("Apr 2020"))
covid_indicator_dIPI <- as.numeric(index(dIPI) %in% covid_dates)
covid_xreg_dIPI <- matrix(covid_indicator_dIPI, ncol=1) #We need to take into account the fewer amount of observation
colnames(covid_xreg_dIPI) <- "COVID"

#Modification to use the dummy variables as external regressors in our ARMA models for dIPI
modelchoice_covid <- function(p,q,data=dIPI, xreg=covid_xreg_dIPI, k=24){
  estim <- try(arima(data, c(p, 0,q), xreg=xreg, optim.control=list(maxit=20000)))
  if (class(estim)=="try-error") return(c("p"=p,"q"=q,"arsignif"=NA,"masignif"=NA,"resnocorr"=NA, "ok"=NA))
  n_xreg <- if(is.null(xreg)) 0 else ncol(xreg)
  coef_names <- names(estim$coef)
  arsignif <- if (p==0) NA else {
    ar_indices <- grep("^ar", coef_names)
    if(length(ar_indices) > 0) coeftest(estim)[ar_indices[p], 4] <= 0.05 else NA
  }
  masignif <- if (q==0) NA else {
    ma_indices <- grep("^ma", coef_names)
    if(length(ma_indices) > 0) coeftest(estim)[ma_indices[q], 4] <= 0.05 else NA
  }
  resnocorr <- sum(Qtests(estim$residuals,24,length(estim$coef)-1)[,2]<=0.05,na.rm=T)==0
  checks <- c(arsignif,masignif,resnocorr)
  ok <- as.numeric(sum(checks,na.rm=T)==(3-sum(is.na(checks))))
  return(c("p"=p,"q"=q,"arsignif"=arsignif,"masignif"=masignif,"resnocorr"=resnocorr,"ok"=ok))
}

#Modification to use the new modelchoice function for dIPI
armamodelchoice_covid <- function(pmax,qmax, xreg=covid_xreg_dIPI){
  pqs <- expand.grid(0:pmax,0:qmax)
  t(apply(matrix(1:dim(pqs)[1]),1,function(row) {
    p <- pqs[row,1]; q <- pqs[row,2]
    cat(paste0("Computing ARMA(",p,",0,",q,") with COVID indicators for dIPI \n"))
    modelchoice_covid(p,q,data=dIPI,xreg=xreg)
  }))
}

armamodels_covid <- armamodelchoice_covid(pmax,qmax, xreg=covid_xreg_dIPI) 
colnames(armamodels_covid) <- c("p", "q", "arsignif", "masignif", "resnocorr", "ok")
selec_covid <- armamodels_covid[armamodels_covid[,"ok"]==1&!is.na(armamodels_covid[,"ok"]),] 
selec_covid

pqs_covid <- apply(selec_covid,1,function(row) list("p"=as.numeric(row[1]),"q"=as.numeric(row[2]))) 
names(pqs_covid) <- paste0("arima(",selec_covid[,1],",0,",selec_covid[,2],")") 

#Estimate the models with indicators for dIPI
models_covid <- lapply(pqs_covid, function(pq) arima(dIPI, c(pq[["p"]], 0, pq[["q"]]), xreg=covid_xreg_dIPI))
vapply(models_covid, FUN.VALUE=numeric(2), function(m) c("AIC"=AIC(m),"BIC"=BIC(m))) 

#RMSE out of sample evaluation for dIPI with COVID
n_obs_dIPI <- length(dIPI)
dIPI_train <- window(dIPI, end = index(dIPI)[n_obs_dIPI-4])
dIPI_test <- tail(dIPI, 4)

covid_xreg_dIPI_train <- covid_xreg_dIPI[1:(n_obs_dIPI-4), , drop=FALSE]
covid_xreg_dIPI_test <- covid_xreg_dIPI[(n_obs_dIPI-3):n_obs_dIPI, , drop=FALSE]

models_train_covid <- lapply(pqs_covid, function(pq) arima(dIPI_train, c(pq[["p"]], 0, pq[["q"]]), xreg=covid_xreg_dIPI_train))

forecasts_covid <- lapply(models_train_covid, function(m) as.zoo(predict(m, 4, newxreg=covid_xreg_dIPI_test)$pred))

rmse_covid <- vapply(forecasts_covid, FUN.VALUE=numeric(1), function(forecast) {
  sqrt(mean((forecast - dIPI_test)^2))
})

rmse_covid

#Fit the AR(5) model with COVID dummy
selected_model_covid <- arima(dIPI, order = c(5, 0, 0), xreg = covid_xreg_dIPI)

#Test of the normality of the residuals using Jarque Bera Test
jb_test_covid <- jarque.bera.test(selected_model_covid$residuals)
print(jb_test_covid)

#Extract AR coefficients
phi_coeffs <- selected_model_covid$coef[grep("^ar", names(selected_model_covid$coef))]
sigma_sq_covid <- selected_model_covid$sigma2

phi1 <- phi_coeffs[1] 

#Covariance matrix for AR(5) 2-step forecasts (using same calculation as in the report but for AR(5))
Sigma_covid <- sigma_sq_covid * matrix(c(1, phi1, phi1, 1 + phi1^2), nrow = 2)

#Same code than without covid correction
future_covid_xreg_dIPI <- matrix(0, nrow = 2, ncol = 1)
colnames(future_covid_xreg_dIPI) <- "COVID"

forecasts_covid_final <- predict(selected_model_covid, n.ahead = 2, newxreg = future_covid_xreg_dIPI)
X_hat_covid <- forecasts_covid_final$pred  

ellipse_points_covid <- ellipse(Sigma_covid, centre = as.numeric(X_hat_covid), level = 0.95)

plot(ellipse_points_covid, type = 'l', lwd = 2, col = 'red',
     xlab = expression(dIPI[T+1]), ylab = expression(dIPI[T+2]),
     main = "95% Confidence Region for dIPI\n(With COVID correction)",
     xlim = range(ellipse_points_covid[,1]) + c(-0.1, 0.1)*diff(range(ellipse_points_covid[,1])),
     ylim = range(ellipse_points_covid[,2]) + c(-0.1, 0.1)*diff(range(ellipse_points_covid[,2])))

points(X_hat_covid[1], X_hat_covid[2], pch = 19, col = 'blue', cex = 1.5)
text(X_hat_covid[1], X_hat_covid[2], "Point Forecast", pos = 3, col = 'blue', cex = 0.8)

grid()

se_1_corrected <- sqrt(Sigma_covid[1,1])
se_2_corrected <- sqrt(Sigma_covid[2,2])

ci_1_covid <- X_hat_covid[1] + c(-1, 1) * qnorm(0.975) * se_1_corrected
ci_2_covid <- X_hat_covid[2] + c(-1, 1) * qnorm(0.975) * se_2_corrected

abline(v = ci_1_covid, lty = 2, col = 'gray')
abline(h = ci_2_covid, lty = 2, col = 'gray')

legend("topright", 
       legend = c("95% Joint CI", "Point Forecast", "Individual 95% CIs"),
       col = c("red", "blue", "gray"),
       lty = c(1, NA, 2),
       pch = c(NA, 19, NA),
       lwd = c(2, NA, 1),
       cex = 0.6)
