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

#' Skyline data import utility
#'
#' @param skydata A character string consisting of the filepath to a Skyline file
#' @param output A character string consisting of a name for the exported Skyline chromatograms
#' @param .skyRunLoc A character string to specify filepath to the SkylineRunner shim
#' if not located in the default C:/User/Downloads location
#' @return A data frame containing data nested by peptide, isotope label (H/L), and transition
#' @export
.skyline_auto_import = function(skydata, output=NULL,.skyRunLoc=NULL){

  # save wd for easy modifications
  curDir = getwd()

  # if an output isn't provided, supply a default
  if (is.null(output) == TRUE){
    output = paste0(curDir,'/skyline_chromatograms.tsv')
  }
  else{
    output =  paste0(curDir,output)
  }

  # point to location of skyline cmd tool and skyline data

  if (is.null(.skyRunLoc)){
    skyline_exe = file.path(Sys.getenv("USERPROFILE"),"Downloads","SkylineRunner.exe") |>
      normalizePath(winslash = "/")
  }

  else {
    skyline_exe = .skyRunLoc
  }
  skyline_doc = skydata

  # run the command to export data automatically

  df <- system2(skyline_exe, args = c(
    paste0('--in="',skyline_doc,'"'),
    paste0('--chromatogram-file="',output,'"')
  ))

  # import file into R

  df = read_tsv(output,name_repair = "universal",col_types = c("c","c","f","f","f","f","f","d","d","d"))

  # fix time and intensity columns so that they are interpreted as a list of numbers

  df = df |>
    mutate(Times = str_split(Times,",") |> map(as.numeric)) |>
    mutate(Intensities = str_split(Intensities,",") |> map(as.numeric))

  # create a nested list based on peptide, isotope label, and transition

  df = df |>
    unnest(cols = c(Times,Intensities)) |>
    nest(.by = c(PeptideModifiedSequence,IsotopeLabelType,FragmentIon,ProductMz))

  df = df |>
    mutate(Threshold=map(data,~.x |>
                           group_by(FileName) |>
                           reframe(.thresholdR(RT = Times, intensity = Intensities))
    ))

  return(df)
}

#' Create a benchmark group for calculating MD
#'
#' @param data A data frame containing properly formatted data
#' @return A data frame with computed moments, mu, and sigma for benchmark samples
#' @export
.setBenchmark = function(data){
  # default name because easier
  df = data

  # create a list of all available injection names
  injNames = levels(as.factor(df["data"][[1]][[1]][[1]]))

  # allow user to select their preferred injections
  injSelect = dlg_list(
    injNames,
    multiple = TRUE,
    title = "Select benchmark injections"
  )$res

  # filter data with user selected injection names
  Bench = df |>
    mutate(
      Benchmark = map(
        Threshold,~.x |>
          filter(FileName %in% injSelect)
      )
    )

  # compute moments for filtered results
  Bench = Bench |>
    mutate(
      BenchmarkMoments = map(
        Benchmark,~.x |>
          reframe(
            moments(RT = RT,intensity = intensity),
            .by = FileName
          )
      )
    )

  # create new columns for mu and sigma values

  Bench = Bench |>
    mutate(
      Mu = map(
        BenchmarkMoments,~.x |>
          select(!FileName) |>
          colMeans()
      )) |>
    mutate(
      Sigma = map(
        BenchmarkMoments,~.x |>
          select(!FileName) |>
          cov()
      ))

  df = left_join(df,Bench)

  return(df)
}

#' Environemnt to store path to SkylineRunner for automated imports
#'
#' @export
.skyRunLocCache <- new.env(parent = emptyenv())
