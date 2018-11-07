.onUnload <- function (libpath) {
  library.dynam.unload("GUSLD", libpath)
}
