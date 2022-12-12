#######################################
# Select the n-th element in list lst
  # input: lst is a R list
  # input: n is an integer
#######################################
select_n <- function(lst, n){
  sapply(lst, `[`, n)
}

###################################################
# calculate RMSE, use 2 norm distance
###################################################
#' Title
#'
#' @param a A numeric vector
#' @param b A numeric vector that has same length as $a$
#'
#' @return The rooted mean squared distance between $a$ and $b$
#' @export
RMSE <- function(a,b){
  return(sqrt(mean((a-b)^2)))
}

#######################################
# Return a vector of b_n = a_{n+1}-a_n
  # for input vector a_n
#######################################
diff_vec <- function(a_n){
  n <- length(a_n)
  b_n <- rep(NA, n-1)
  for (k in 1:(n-1)){
    b_n[k] = a_n[k+1] - a_n[k]
  }
  return(b_n)
}

#######################################
# Check whether there is negative value
  # in input vector a_n
#######################################
#' Check whether there is negative values in a vector
#'
#' @param a_n A number array
#'
#' @return Boolean
#' @export
#'
#' @examples
#' check_negative(c(1,2,-1))
#' check_negative(c(1,2))
check_negative <- function(a_n){
  #initial
  i=1; K=FALSE
  #loop stop when detect a negative value
  while ((K==FALSE)&(i<=length(a_n))){
    if(a_n[i]<0){
      K <- TRUE
    }
    i=i+1
  }
  return(K)
}

##########################################################
# generate simulate data
  # Use this function to generate the data matrix y
  # The Observation y(j,i) represents the function
  #     value observed at time j and location i.
  # y(j,i) = Y_j(i) + epsilon_{j,i}
  # Y_j is cdf of N(noise, change_speed*j)
  # The true value of function at time j
  #     is cdf of N(0,change_speed*j)
  # input:
    # t: time index, used to record the time,
    #    must be uniform design
    # locations: location, the sample locations on function
    #    must (much) less than t to aviod sigularity
    # func_test_locations: used to test the function
    #                      estimation accurancy
    # sigma_noise: the variance for noise term,
    #              must smaller than change_speed
    # sigma_error: the variance for measurement error,
    #              default as 0.
    # change_speed: the variance change speed over time for
    #               variance of true cdf function
###########################################################
#' Generate Simulate Data Sets
#'
#' @param t time length, the default value is 40 unit time
#' @param locations The observed locations (x's) on functions
#' @param func_test_locations The locations needs to evaluate in simulation
#' @param sigma_noise The variance of noise for mean shift on normal distributions, the default value is 0.05
#' @param sigma_error The variance of measurement error, the default value is 0
#' @param change_speed The speed of normal cdf rotates
#'
#' @return The function returns 2 components:
#' "data" a data matrix contains the values of normal cdf at observed locations;
#' "test" a data matrix contains the values of normal cdf at evaluate locations
#' @export
#'
#' @examples
#' Sim_data(t=40,locations=c(1,2,3),func_Test_locations=c(4,5,6),sigma_noise=0.05,sigma_error=0.02,change_speed=0.2)$data
#' Sim_data(t=200,locations=c(1,2,3),func_Test_locations=c(4,5,6),sigma_noise=0.05,sigma_error=0,change_speed=0.235)$test
Sim_data <- function(t=40, locations, func_test_locations, sigma_noise = 0.05, sigma_error=0, change_speed=0.235) {

  l = length(locations)
  l2 = length(func_test_locations)

  data.sim <- matrix(0,t,l)
  data.sim_test <- matrix(0,t,l2)

  for(j in 1:t){
    #noise term for model at time j
    noise <- rnorm(1,0,sigma_noise)

    for(i in 1:l){
      data.sim[j,i] <- pnorm(locations[i], noise, (0+change_speed*j))
      #measurement error for model at location l
      measure.error <- rnorm(1, 0, sigma_error)
      data.sim[j,i] <- data.sim[j,i] + measure.error

      #truncated the value so that it belongs to [0,1]
      data.sim[j,i] <- min(1,data.sim[j,i]); data.sim[j,i] <- max(0,data.sim[j,i])

    }
    for(i in 1:l2){
      data.sim_test[j,i] <- pnorm(func_test_locations[i], noise, (0+change_speed*j))
      #measurement error for model at location l
      measure.error <- rnorm(1, 0, sigma_error)
      data.sim_test[j,i] <- data.sim_test[j,i] + measure.error

      #truncated the value so that it belongs to [0,1]
      data.sim_test[j,i] <- min(1,data.sim_test[j,i]); data.sim_test[j,i] <- max(0,data.sim_test[j,i])
    }
  }
  my_list <- list("data" = data.sim, "test" = data.sim_test)
  return(my_list)
}

#####################################################
# Vector Auto-Regression Predictions
#####################################################
#' Vector Auto-Regression Predictions
#'
#' @param y A data matrix with each column is values of function on observed locations
#' @param p Integer for the lag order (default is p=1)
#' @param ahead An integer specifying the number of forecast steps.
#' @param ci The forecast confidence interval (default is 0.95)
#'
#' @importFrom vars VAR
#'
#' @return A list with 3 items: "forecast" forecast value of the function for each step at locations; "lower" lower confidence band; "upper" upper confidence band
#' @export
#'
#' @examples
#' y_ram = matrix(runif(48),nrow=3)
#' Var_Pred(y_ram, p=1, ahead=5, ci=0.95)
Var_Pred <- function(y, p = 1, ahead, ci = 0.95){
  Var_est <- vars::VAR(y, p)
  forecast <- vars::predict(VAR_est, n.ahead = ahead, ci = ci)

  forecast.sim <- sapply(forecast$fcst, `[`, 1)
  for(k in 2:ahead){
    forecast.sim <- rbind(forecast.sim, sapply(forecast$fcst, `[`, k))
  }

  lower.sim <- sapply(forecast$fcst, `[`, (ahead+1))
  for(k in (ahead+1):(2*ahead)){
    lower.sim <- rbind(lower.sim, sapply(forecast$fcst, `[`, k))
  }

  upper.sim <- sapply(forecast$fcst, `[`, ((2*ahead)+1))
  for(k in ((2*ahead)+1):(3*ahead)){
    upper.sim <- rbind(upper.sim, sapply(forecast$fcst, `[`, k))
  }

  my_list <- list("forecast" = forecast.sim, "lower" = lower.sim, "upper" = upper.sim)
  return(my_list)
}


######################################################
# Truncated between [0,1]
######################################################
#' Truncate vector between 0 and 1
#'
#' @param an A numeric array
#'
#' @return The array with every entry bigger than 1 shrinked into 1 and every entry smaller than 0 shrinked into 0
#' @export
#'
#' @examples
#' tru_func(c(0,2,1,-2))
tru_func <- function(an){
  n <- length(an)
  for(i in 1:n){
    if (an[i]<0){an[i]=0}
    if (an[i]>1){an[i]=1}
  }
  return(an)
}


######################################################
# isotonic regression (monotone non-decreasing)
# Input: vector a, return an monotone non-decreasing
#        vector by implementing PAVA.
######################################################
#' Isotonic (Increasing) Regression
#'
#' @param a A numeric vector
#'
#' @return A vector with increasing terms using Pool-Adjacent-Violators Algorithm
#' @export
#'
#' @examples iso_reg(c(1,5,2,3,1))
#'
#'
iso_reg <- function(a){
   C <- TRUE
   n <- length(a)
   while(C==TRUE){
     block = c(0); bk_n = 1
     for(i in 1:(n-1)){
       if (a[i]<a[i+1]) {
         block <- append(block, i)
         bk_n <- bk_n+1
       } else if (a[i]>a[i+1]){
         block <- append(block, (i+1))
         bk_n <- bk_n+1
         s=0
         for(j in (block[bk_n-1]+1):block[bk_n]){s=s+a[j]}
         for(j in (block[bk_n-1]+1):block[bk_n]){a[j] = s/(block[bk_n]-block[bk_n-1])}
       }
     }
     C = check_negative(diff_vec(a))
   }
   return(a)
}


###################################################
# Monotone Cubic Hermite Interpolation
# input: x, same as location defined previously
#        y, estimated value at each locations
#        xnew, new evaluation point.
###################################################

#main function
#' Monotone Cubic Hermite Interpolation
#'
#' @param x Observed locations on a function
#' @param y Estimated Value of the function at each observed locations
#' @param xnew Evaluate locations
#'
#' @return Estimated Value of the function at each evaluate locations
#' @export
#'
#' @examples
#' mchi(c(1,2,3),c(0.6,0.7,0.77),c(4,5))
mchi <- function(x,y,xnew){
  #Compute the linear slopes between successive points
  n <- length(y)
  delta <- rep(NA,(n-1))
  delta <- diff_vec(y)/diff_vec(x)

  #Initialize the tangents at every interior data point
  # as the average of the secants
  m <- rep(NA,(n-1))
  m[1] <- delta[1]
  m[n] <- delta[n-1]
  for (i in 2:(n-1)){
    m[i] = (delta[i-1]+delta[i])/2
    if (delta[i]==0){
      m[i] <- 0
      m[i+1] <- 0
    }
  }

  #update m if it is not equal to 0
  for(i in 1:(n-1)){
    if(delta[i]!=0){
      alpha <- m[i]/delta[i]
      beta <- m[i+1]/delta[i]
      if ((alpha>3)|(beta>3)){
        m[i] = 3*delta[i]
      }
    }
  }

  #cubic interpolation
  #record the location of xnew
  for(i in 1:(n-1)){
    if ((x[i]<=xnew)&(x[i+1]>=xnew)){k<-i}
  }
  dif <- x[k+1] - x[k]
  t <- (xnew-x[k])/dif

  #Basis functions for cubic Hermite spline
  h00 <- 2*(t^3)-3*(t^2)+1
  h10 <- t^3-2*(t^2)+t
  h01 <- -2*(t^3)+3*(t^2)
  h11 <- (t^3)-(t^2)

  ynew <- y[k]*h00+dif*m[k]*h10+y[k+1]*h01+dif*m[k+1]*h11

  return(unname(ynew))
}

#####################################################
# Function Auto-Regression Predictions
#####################################################
#' Functional Auto-Regression Predictions
#'
#' @param x Observed locations of function
#' @param y Function values at observed locations through time, row is location, column is time.
#' @param xnew Evaluate locations of function. When this is same as 'x', the return value is same as vector auto regression.
#' @param p Integer for the lag order (default is p=1)
#' @param ahead An integer specifying the number of forecast steps
#' @param ci The forecast confidence interval (default is 0.95)
#'
#' @return forecast value of the function for each step at evaluate locations
#' @export
#'
#' @examples
#' x = c(-1,0,1)
#' y = cbind(pnorm(x,0,0.05),pnorm(x,0,0.1),pnorm(x,0,0.15),pnorm(x,0,0.2))
#' Fun_Pred(x, y, xnew = c(0.5), p=1, ahead=1, ci=0.95)
Fun_Pred <- function(x, y, xnew,p = 1, ahead, ci = 0.95){
  Dis_pred <- Var_Pred(y, p, ahead, ci)
  Var_Pred_Value <- Dis_pred$forecast
  Fun_Pred_Value <- mchi(x, Pred_Value, xnew)
  return(Fun_Pred_Value)
}

#####################################################
# Function Auto-Regression for CDF Predictions
#####################################################
#' Monotone Functional Auto-Regression Predictions
#'
#' @param x Observed locations of function
#' @param y Function values at observed locations through time, row is location, column is time.
#' @param xnew Evaluate locations of function. When this is same as 'x', the return value is same as vector auto regression.
#' @param p Integer for the lag order (default is p=1)
#' @param ahead An integer specifying the number of forecast steps
#' @param ci The forecast confidence interval (default is 0.95)
#'
#' @return forecast value of the function for each step at evaluate locations
#' @export
#'
#' @examples
#' x = c(-1,0,1)
#' y = cbind(pnorm(x,0,0.05),pnorm(x,0,0.1),pnorm(x,0,0.15),pnorm(x,0,0.2))
#' Fun_cdf_Pred(x, y, xnew = c(0.5), p=1, ahead=1, ci=0.95)
Fun_cdf_Pred <- function(x, y, xnew,p = 1, ahead, ci = 0.95){
  Dis_pred <- Var_Pred(y, p, ahead, ci)
  Var_Pred_Value <- Dis_pred$forecast
  for(i in 1:ahead){
    Var_Pred_Value[i,] <- iso_reg(tru_func(Var_Pred_Value[i,]))
  }
  Fun_Pred_Value <- mchi(x, Pred_Value, xnew)
  return(Fun_Pred_Value)
}





