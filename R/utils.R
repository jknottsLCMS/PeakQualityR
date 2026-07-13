#' Performs a baseline correction on LC trace data
#'
#' Applies a default baseline of 5%, but can be modified by the user.
#' Important to avoid error from baseline noise.
#'
#' @param RT A numeric vector - retention time from LC trace
#' @param intensity A numeric vector - intensity from LC trace
#' @param baseline A numeric value between 0 and 1; default is 0.05 (5%)
#' @return A filtered tibble of LC trace data keeping the data above the set baseline
#' @export
.thresholdR = function(RT,intensity,baseline=0.05){

  df = tibble("RT" = {{RT}},
              "intensity" = {{intensity}})

  df = df |>
    filter(intensity>max(intensity)*baseline)

  return(df)
}
