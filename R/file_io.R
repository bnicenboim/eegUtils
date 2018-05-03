#' Function for reading raw data.
#'
#' Currently BDF/EDF and 32-bit .CNT files are supported. Filetype is determined
#' by the file extension.The \code{edfReader} package is used to load BDF files.
#' The function creates an eeg_data structure for subsequent use.
#'
#' @author Matt Craddock, \email{matt@mattcraddock.com}
#' @param file_name File to import. Should include file extension.
#' @param file_path Path to file name, if not included in filename.
#' @param chan_nos Channels to import. All channels are included by default.
#' @import edfReader
#' @import tools
#' @importFrom purrr map_df
#' @importFrom tibble tibble as_tibble
#' @export

import_raw <- function(file_name, file_path = NULL, chan_nos = NULL, ...) {
  file_type <- tools::file_ext(file_name)

  if (file_type == "bdf" | file_type == "edf") {
    data <- edfReader::readEdfSignals(edfReader::readEdfHeader(file_name))
    sigs <- purrr::map_df(data, "signal")
    srate <- data[[1]]$sRate
    events <- sigs$Status %% (256)
    timings <- tibble::tibble(sample = 1:dim(sigs)[[1]])
    timings$time <- (timings$sample - 1) / srate
    if (is.null(chan_nos)) {
      chan_nos <- 1:(dim(sigs)[[2]] - 1)
    }
    sigs <- tibble::as_tibble(sigs[, chan_nos])
    events_diff <- diff(events)
    event_table <- tibble::tibble(event_onset = which(events_diff > 0) + 1,
                                  event_time = which(events_diff > 0) / srate,
                              event_type = events[which(events_diff > 0) + 1])
    data <- eeg_data(data = sigs, srate = srate,
                     events = event_table, timings = timings,
                     continuous = TRUE)
  } else if (file_type == "cnt") {
    data <- import_cnt(file_name)
    sigs <- tibble::as_tibble(t(data$chan_data))
    names(sigs) <- data$chan_info$chan_name
    srate <- data$head_info$samp_rate
    timings <- tibble::tibble(sample = 1:dim(sigs)[[1]])
    timings$time <- (timings$sample - 1) / srate
    event_table <- tibble::tibble(event_onset = data$event_list$offset + 1,
                                  event_time = (data$event_list$offset + 1) / srate,
                                  event_type = data$event_list$event_type)
    data <- eeg_data(data = sigs, srate = srate,
                     chan_info = data$chan_info[1:4],
                     events = event_table, timings = timings,
                     continuous = TRUE)
  } else if (file_type == "vhdr") {
    # Takes the files from the header:
    header_info <- load_vhdr(paste0(file_path, file_name))
    data_file <- header_info$common_info$data_file
    data_ext <- tools::file_ext(data_file)
    # It only accepts .dat files (for now)
    if(data_ext == "dat"){
      vmrk_file <- header_info$common_info$vmrk_file
      srate  <- header_info$common_info$srate
      event_table <- load_vmrk(paste0(file_path, vmrk_file), srate)
      data <- import_dat(paste0(file_path, data_file), 
                        header_info$common_info$data_points,
                        header_info$chan_info, 
                        srate, 
                        header_info$common_info$orientation, 
                        event_table)


    } else {
      warning(paste0(".",data_ext, " files are unsupported."))
    }

  } else{
    warning("Unsupported filetype")
    return()
  }
  return(data)
}


#' Import Neuroscan .CNT file
#'
#' Beta version of function to import Neuroscan .CNT files. Only intended for
#' import of 32-bit files.
#'
#' @param file_name Name of .CNT file to be loaded.
#' @importFrom tibble tibble

import_cnt <- function(file_name) {

  cnt_file <- file(file_name, "rb")

  # Read in meta-info - number of channels, location of event table, sampling
  # rate...

  pos <- seek(cnt_file, 12)
  next_file <- readBin(cnt_file, integer(), size = 4, n = 1, endian = "little")
  pos <- seek(cnt_file, 353)
  n_events <- readBin(cnt_file, integer(), n = 1, endian = "little")
  pos <- seek(cnt_file, 370)
  n_channels <- readBin(cnt_file, integer(), n = 1, size = 2,
                        signed = FALSE, endian = "little")
  pos <- seek(cnt_file, 376)
  samp_rate <-  readBin(cnt_file, integer(), n = 1, size = 2,
                        signed = FALSE, endian = "little")
  pos <- seek(cnt_file, 864)
  n_samples <- readBin(cnt_file, integer(), size = 4, n = 1, endian = "little")
  pos <- seek(cnt_file, 886)
  event_table_pos <- readBin(cnt_file, integer(), size = 4,
                             n = 1, endian = "little") # event table
  pos <- seek(cnt_file, 900)

  data_info <- tibble::tibble(n_events,
                              n_channels,
                              samp_rate,
                              event_table_pos)

  chan_df <- tibble::tibble(chan_name = character(n_channels),
                            chan_no = numeric(n_channels),
                            x = numeric(n_channels),
                            y = numeric(n_channels)
  )

  # Read channel names and locations

  for (i in 1:n_channels) {
    chan_start <- seek(cnt_file)
    chan_df$chan_name[i] <- readBin(cnt_file, character(),
                                    n = 1, endian = "little")
    chan_df$chan_no[i] <- i
    pos <- seek(cnt_file, chan_start + 19)
    chan_df$x[i] <- readBin(cnt_file, double(), size = 4,
                            n = 1, endian = "little") # x coord
    chan_df$y[i] <- readBin(cnt_file, double(), size = 4,
                            n = 1, endian = "little") # y coord

    pos <- seek(cnt_file, chan_start + 47)
    chan_df$baseline[i] <- readBin(cnt_file, integer(), size = 1,
                                   n = 1, endian = "little")
    pos <- seek(cnt_file, chan_start + 59)
    chan_df$sens[i] <- readBin(cnt_file, double(), size = 4, n = 1,
                               endian = "little")
    pos <- seek(cnt_file, chan_start + 71)
    chan_df$cal[i] <- readBin(cnt_file, double(), size = 4, n = 1,
                              endian = "little")
    pos <- seek(cnt_file, (900 + i * 75))
  }

  beg_data <- seek(cnt_file) # beginning of actual data
  real_n_samples <- event_table_pos - (900 + 75 * n_channels) / (2 * n_channels)

  frames <- floor((event_table_pos - beg_data) / n_channels / 4)

  chan_data <- matrix(readBin(cnt_file,
                              integer(),
                              size = 4,
                              n = n_channels * frames,
                              endian = "little"),
                      nrow = n_channels, ncol = frames)

  # rescale chan_data to microvolts
  mf <- chan_df$sens * (chan_df$cal / 204.8)
  chan_data <- (chan_data - chan_df$baseline) * mf

  # Read event table

  pos <- seek(cnt_file, event_table_pos)
  teeg <- readBin(cnt_file, integer(), size = 1, n = 1, endian = "little")
  tsize <- readBin(cnt_file, integer(), n = 1, endian = "little")
  toffset <- readBin(cnt_file, integer(), n = 1, endian = "little")
  ev_table_start <- seek(cnt_file)

  ev_list <- tibble::tibble(event_type = integer(n_events),
                            keyboard = character(n_events),
                            keypad_accept = integer(n_events),
                            accept_evl = integer(n_events),
                            offset = integer(n_events),
                            type = integer(n_events),
                            code = integer(n_events),
                            latency = numeric(n_events),
                            epochevent = integer(n_events),
                            accept = integer(n_events),
                            accuracy = integer(n_events)
                            )

  for (i in 1:n_events) {
    ev_list$event_type[i] <- readBin(cnt_file, integer(), size = 2,
                                     n = 1, endian = "little")
    ev_list$keyboard[i] <- readBin(cnt_file, integer(), size = 1,
                                   n = 1, endian = "little")
    temp <- readBin(cnt_file, integer(), size = 1, n = 1, signed = FALSE,
                    endian = "little")
    ev_list$keypad_accept[i] <- bitwAnd(15, temp)
    ev_list$accept_evl[i] <- bitwShiftR(temp, 4)
    ev_list$offset[i] <- readBin(cnt_file, integer(), size = 4,
                                 n = 1, endian = "little")
    ev_list$type[i] <- readBin(cnt_file, integer(), size = 2,
                               n = 1, endian = "little")
    ev_list$code[i] <- readBin(cnt_file, integer(), size = 2,
                               n = 1, endian = "little")
    ev_list$latency[i] <- readBin(cnt_file, double(), size = 4,
                                  n = 1, endian = "little")
    ev_list$epochevent[i] <- readBin(cnt_file, integer(), size = 1,
                                     n = 1, endian = "little")
    ev_list$accept[i] <- readBin(cnt_file, integer(), size = 1,
                                 n = 1, endian = "little")
    ev_list$accuracy[i] <- readBin(cnt_file, integer(), size = 1,
                                   n = 1, endian = "little")
  }

  ev_list$offset <- (ev_list$offset - beg_data) / (4 * n_channels) + 1

  close(cnt_file)
  out <- list(chan_info = chan_df, head_info = data_info,
              chan_data = chan_data, event_list = ev_list)
}


#' Load EEGLAB .set files
#'
#' EEGLAB .set files are standard Matlab .mat files, but EEGLAB can be set to
#' export either v6.5 or v7.3 format files. Only v6.5 files can be read with
#' this function. v7.3 files (which use HDF5 format) are not currently
#' supported, as they cannot be fully read with existing tools.
#'
#' @param file_name Filename (and path if not in present working directory)
#' @param df_out Defaults to FALSE - outputs an object of class eeg_data. Set to
#'   TRUE for a normal data frame.
#' @author Matt Craddock \email{matt@mattcraddock.com}
#' @importFrom R.matlab readMat
#' @importFrom dplyr group_by mutate select rename
#' @importFrom tibble tibble as_tibble
#' @importFrom purrr is_empty
#' @export


load_set <- function(file_name, df_out = FALSE) {

  temp_dat <- R.matlab::readMat(file_name)
  var_names <- dimnames(temp_dat$EEG)[[1]]

  n_chans <- temp_dat$EEG[[which(var_names == "nbchan")]]
  n_trials <- temp_dat$EEG[[which(var_names == "trials")]]
  times <- temp_dat$EEG[[which(var_names == "times")]]

  chan_info <- temp_dat$EEG[[which(var_names == "chanlocs")]]
  col_names <- dimnames(chan_info)[1]
  size_chans <- dim(chan_info)
  chan_info <- lapply(chan_info, function(x) ifelse(purrr::is_empty(x), NA, x))
  dim(chan_info) <- size_chans
  dimnames(chan_info) <- col_names
  chan_info <- tibble::as_tibble(t(as.data.frame(chan_info)))
  chan_info <- tibble::as_tibble(data.frame(lapply(chan_info, unlist),
                                            stringsAsFactors = FALSE))

  # check if the data is stored in the set or in a separate .fdt
  if (is.character(temp_dat$EEG[[which(var_names == "data")]])) {
    message("loading from .fdt")
    fdt_file <- paste0(tools::file_path_sans_ext(file_name), ".fdt")
    fdt_file <- file(fdt_file, "rb")

    # read in data from .fdt
    # do this in chunks to avoid memory errors for large files...?
    signals <- readBin(fdt_file,
                       "double",
                       n = n_chans * n_trials * length(times),
                       size = 4,
                       endian = "little")
    close(fdt_file)

    dim(signals) <- c(n_chans, length(times) * max(n_trials, 1))
    times <- rep(times, max(n_trials, 1))

    if (n_trials == 1) {
      continuous <- TRUE
    } else {
      continuous <- FALSE
    }

  } else {

    # if the data is in the .set file, load it here instead of above
    signals <- temp_dat$EEG[[which(dimnames(temp_dat$EEG)[[1]] == "data")]]
    dim_signals <- dim(signals)

    if (length(dim_signals) == 3) {
      dim(signals) <- c(dim_signals[1], dim_signals[2] * dim_signals[3])
      times <- rep(times, n_trials)
      continuous <- FALSE
    } else {
      continuous <- TRUE
    }
  }

  signals <- data.frame(cbind(t(signals), times))
  srate <- c(temp_dat$EEG[[which(var_names == "srate")]])
  names(signals) <- c(unique(chan_info$labels), "time")
  signals <- dplyr::group_by(signals, time)
  signals <- dplyr::mutate(signals, epoch = 1:n())
  signals <- dplyr::ungroup(signals)

  event_info <- temp_dat$EEG[[which(var_names == "event")]]

  event_table <- tibble::as_tibble(t(matrix(as.integer(event_info),
                               nrow = dim(event_info)[1],
                               ncol = dim(event_info)[3])))

  names(event_table) <- unlist(dimnames(event_info)[1])
  event_table$event_time <- (event_table$latency - 1) / srate
  event_table <- dplyr::select(event_table, latency, event_time, type, epoch)
  event_table <- dplyr::rename(event_table, event_type = "type",
                               event_onset = "latency")
  event_table$time <- NA

  if (df_out) {
    return(signals)
  } else {
    signals$time <- signals$time / 1000 # convert to seconds - eeglab uses milliseconds
    timings <- tibble::tibble(time = signals$time,
                              epoch = signals$epoch,
                              sample = 1:length(signals$time))
    event_table$time <- timings[which(timings$sample %in% event_table$event_onset, arr.ind = TRUE), ]$time
    out_data <- eeg_data(signals[, 1:n_chans],
                         srate = srate,
                         timings = timings,
                         continuous = continuous,
                         chan_info = chan_info,
                         events = event_table)
    if (!continuous) {
      class(out_data) <- c("eeg_epochs", "eeg_data")
    }
    out_data
  }
}

#' Import binary .dat file
#'
#' Beta version of a function to import binary .dat files. Only intended for
#' importing of vectorized files in time domain.
#'
#' @param file_name Name of .dat file to be loaded.
#' @param data_points Number of data points.
#' @param srate Sampling rate.
#' @param orientation Data orientation: VECTORIZED=ch1,pt1, ch1,pt2..., MULTIPLEXED=ch1,pt1, ch2,pt1.
#' @param domain needs to be time domain.
#' @param event_table Table with markers.
#' @author Bruno Nicenboim \email{bruno.nicenboim@uni-potsdam.de}
#' @import tidyverse
#' @export


import_dat <- function(file_name, data_points, chan_info, 
                   srate, orientation = "vectorized", event_table = NULL) {
  
  n_chan <- nrow(chan_info)
  dat_bin <- readBin(file_name, "double", data_points * n_chan, size = 4)

  if(stringr::str_sub(orientation,1,nchar("vector")) %>% 
      stringr::str_to_lower() == "vector"){
    sigs <- matrix(dat_bin, ncol = n_chan)  %>% 
                tibble::as.tibble() 
  } else {
    stop("orientiation needs to be 'vectorized'")
  }

  colnames(sigs)  <- chan_info$labels
  n_chan <- nrow(chan_info)
  dat_bin <- readBin(file_name, "double", data_points * n_chan, size = 4)
  timings <- tibble::tibble(sample = 1:dim(sigs)[[1]]) %>%
              dplyr::mutate(time = (sample - 1) / srate)
  data <- eeg_data(data = sigs, srate = srate,
                     chan_info = chan_info,
                     events = event_table, timings = timings,
                     continuous = TRUE)

  return(data)
}


#' Load header file (.vhdr) produced by BrainVision Analyzer v 2.0
#'
#' Beta version of a function to load header file (.vhdr) produced by BrainVision Analyzer v 2.0.
#'
#' @param file_name Name of .vhdr file to be loaded.
#' @author Bruno Nicenboim \email{bruno.nicenboim@uni-potsdam.de}
#' @import tidyverse
#' @export
#' 
load_vhdr <- function(file_name) {

  content_vhdr <- readr::read_file(file_name) %>% 
            stringr::str_match_all(stringr::regex("\\[(.*?)\\]\n(.*?\n)(\n|\\Z)", 
              dotall = TRUE, multiline=TRUE) ) %>% .[[1]]

  read_metadata <- function(tag) {
   readr::read_delim(content_vhdr[content_vhdr[,2] == tag,3],
                 delim = "=", comment = ";", col_names=c("type", "value"))
 }           

  out <- list()
  binary_info <- read_metadata("Binary Infos")
  channel_info <- read_metadata("Channel Infos") %>% 
              separate(value, c("labels","ref","res","unit"), sep=",")
  coordinates <- read_metadata("Coordinates") %>% 
                separate(value, c("radius","theta","phi"), sep=",")
  
  
  common_info <- read_metadata("Common Infos") %>% 
                 tidyr::spread(type, value) %>% 
                 dplyr::transmute(data_points = as.numeric(DataPoints),
                            seg_data_points = as.numeric(SegmentDataPoints),
                            orientation = DataOrientation,
                            domain = DataType,
                            srate =  1000000 / as.numeric(SamplingInterval),
                            data_file = DataFile,
                            vmrk_file = MarkerFile)

  if(stringr::str_sub(common_info$domain,1,nchar("time")) %>% 
      stringr::str_to_lower() != "time") {
    stop("DataType needs to be 'time'")
  } else 
  chan_info <- dplyr::full_join(channel_info, coordinates, by = "type") %>% 
                dplyr::mutate(type = "EEG", sph.theta = NA, sph.phi = NA, 
                        sph.radius = NA, urchan = NA,X=NA,Y=NA,Z=NA) %>%
                dplyr::select(labels, type, theta, radius, X, Y, Z, sph.theta,
                 sph.phi, sph.radius, urchan, ref)

  out <- list()
  out$chan_info <- chan_info
  out$common_info <- common_info
  return(out)
}

#' Load marker file (.vmrk) produced by BrainVision Analyzer v 2.0
#'
#' Beta version of a function to load marker files (.vhdr) produced by BrainVision Analyzer v 2.0.
#'
#' @param file_name Name of .vhdr file to be loaded.
#' @param srate Sampling rate.
#' @author Bruno Nicenboim \email{bruno.nicenboim@uni-potsdam.de}
#' @import tidyverse
#' @export
#' 


load_vmrk <- function(file_name, srate) {
  # Each entry looks like this in the vmrk file:
  # Mk<Marker number>=<Type>,<Description>,<Position in data
  # points>, <Size in data points>, <Channel number (0 = marker is related to
  # all channels)>, <Date (YYYYMMDDhhmmssuuuuuu)>
  # More information can be found in
  # http://pressrelease.brainproducts.com/markers/
  markers_info_lines <- readr::read_lines(file_name) %>% 
                        stringr::str_detect("Mk[0-9]*?=") %>% which 
  start <- markers_info_lines %>% min - 1
  end <- markers_info_lines %>% max - start 

  col_names = c("Mk_number=Type","description","pos_dpoints", 
              "size_dp", "channel","date")

  tibble_vmrk <- suppressWarnings(readr::read_csv(file_name, col_names = col_names,
                    col_types = readr::cols(
                                `Mk_number=Type` = readr::col_character(),
                                description = readr::col_character(),
                                pos_dpoints = readr::col_integer(),
                                size_dp = readr::col_integer(),
                                channel = readr::col_integer(),
                                date = readr::col_double()
                                ),
                    skip = start, n_max= end,
                    trim_ws = TRUE)) %>%
                  tidyr::separate(`Mk_number=Type`, 
                            c("mk","type"), sep ="=") %>%
                  dplyr::mutate(mk = as.numeric(stringr::str_remove(mk,"Mk")))

  event_table <- tibble_vmrk %>% 
                  dplyr::transmute(event_type = 
                        stringr::str_replace(description,"s","") %>% 
                                          as.integer, 
                          event_time = (pos_dpoints - 1 ) / srate,
                          event_onset = pos_dpoints - 1) %>% 
                  dplyr::filter(!is.na(event_type))
           
  return(event_table)              
}
