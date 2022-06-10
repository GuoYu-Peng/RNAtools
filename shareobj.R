check_path <- function(x, quit = FALSE, mkdir = FALSE) {
  if (!file.exists(x)) {
    writeLines("\n[ERROR] path does not exist:")
    print(x)
    if (quit) {
      quit(status = 1)
    } else {
      if (mkdir) {
        writeLines("\n[INFO] mkdir:")
        print(x)
        dir.create(x)
      }
    }
  }
}

make_path <- function(dest_dir = ".", suffix = "txt", ...) {
  y <- paste(..., suffix, sep = ".")
  z <- file.path(dest_dir, y)
  return(z)
}

log_msg <- function(level = "INFO", ...) {
  y <- paste("\n[", level, "] ", ..., sep = "")
  return(y)
}