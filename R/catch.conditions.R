catch.conditions <- function (expr) {
  W <- NULL
  w.handler <- function(w) {
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                   warning = w.handler), warning = W)
}
