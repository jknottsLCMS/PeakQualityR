#' Calculate the first moment (mean retention time) using LC trace data
#'
#' @param RT A numeric vector - retention time from your LC trace
#' @param intensity A numeric vector - intensity from your LC trace
#' @return Mean retention time (numeric)
#' @export
muFun = function(RT,intensity){

  #variables are separated for ease of use

  # create a weighted sum consisting of E(RT * intensity)
  wt = sum(RT * intensity)

  # create a population by adding up intensity
  pop = sum(intensity)

  # calculate a weighted average
  return(wt/pop)

}

#' Calculate the second moment and return standard deviation
#'
#' @param RT A numeric vector - retention time from your LC trace
#' @param intensity A numeric vector - intensity from your LC trace
#' @return A numeric result for standard deviation
#' @export
sigmaFun = function(RT, intensity){

  # variables are separated for ease of use

  #calculate mu using muFun
  mu = muFun(RT,intensity)

  # create weighted sum consisting of E((RT-mu)^2 * intensity)
  wt = sum((RT-mu)^2 * intensity)

  # tally population
  pop = sum(intensity)

  # calculate result (stdev)
  return(sqrt(wt/pop))

}

#' Calculate the third moment and return a numeric value for skewness
#'
#' @param RT A numeric vector - retention time from LC trace
#' @param intensity A numeric vector - intensity from LC trace
#' @return A numeric value - calculated skewness
#' @export
skewFun = function(RT,intensity){

  # calculate mu using muFun
  mu = muFun(RT,intensity)

  # calculate sigma using sigmaFun
  sigma = sigmaFun(RT,intensity)

  # calculate a weighted sum consisting of E((RT-mu)^3 * intensity)
  wt = sum((RT-mu)^3 * intensity)

  # tally population
  pop = sum(intensity)

  return((wt/pop)/sigma^3)
}

#' Calculate the third moment and return a numeric value for kurtosis
#'
#' @param RT A numeric vector - retention time from LC trace
#' @param intensity A numeric vector - intensity from LC trace
#' @return A numeric value - calculated kurtosis
#' @export
kurtFun = function(RT,intensity){

  # calculate mu using muFun
  mu = muFun(RT,intensity)

  # calculate sigma using sigmaFun
  sigma = sigmaFun(RT,intensity)

  # create a weighted sum consisting of E((RT-mu)^4 * intensity)
  wt = sum((RT-mu)^4*intensity)

  # tally population
  pop = sum(intensity)

  return((wt/pop) / sigma^4)

}

#'Calculate and display all 4 moments for LC trace data
#'
#' @param RT A numeric vector - retention time from LC trace
#' @param intensity A numeric vector - intensity from LC trace
#' @return A concatenated string of results
#' @export
moments = function(RT,intensity){

  df = tibble(RT = {{RT}},
              intensity = {{intensity}})

  # mew i'm a cat nya =(^.^)=
  mew = muFun(RT,intensity)
  sig = sigmaFun(RT,intensity)
  skew = skewFun(RT,intensity)
  kurt = kurtFun(RT,intensity)

  tblOut = tibble("Mean Retention Time" = mew,
              "Standard deviation" = sig,
              "Skewness" = skew,
              "Kurtosis" = kurt)

  # out = paste("Mean retention time: ",mew,
  #             "\nStandard deviation: ",sig,
  #             "\nSkewness: ",skew,
  #             "\nKurtosis: ",kurt)

  # cat(out)
  return(tblOut)
}

