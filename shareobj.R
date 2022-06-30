log_msg <- function(level = "INFO", ...) {
  y <- paste("\n[", level, "] ", ..., sep = "")
  writeLines(y)
  return(y)
}

check_path <- function(x, quit = FALSE, mkdir = FALSE) {
  if (file.exists(x)) {
    return(TRUE)
  } else {
    if (quit) {
      log_msg("ERROR", "path does not exist:", x)
      quit(status = 1)
    } else {
      if (mkdir) {
        log_msg("INFO", "mkdir:", x)
        dir.create(x)
        return(TRUE)
      } else {
        return(FALSE)
      }
   }
  }
}

make_path <- function(dest_dir = ".", suffix = "txt", ...) {
  y <- paste(..., suffix, sep = ".")
  z <- file.path(dest_dir, y)
  return(z)
}

