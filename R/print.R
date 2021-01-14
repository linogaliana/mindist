print.mindist <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  nsteps <- ifelse(attr(x, "approach") == "two_step",
                   "2 steps",
                   "1 step")


  p <- length(x$estimates$theta_hat)
  m <- nrow(x$moments$moment_first_step)

  cat(sprintf("Minimum distance estimation\n%s estimation:    \n", nsteps))
  cat("---------------------------\n")

  cat(sprintf("Number of parameters:    %s\n", p))
  cat(sprintf("Number of moments:       %s\n", m))

  if (p == m){
    cat("(Exact identification)\n")
  } else if (m > p){
    cat("(Over-identification)\n")
  } else{
    cat("(Under-identification)\n")
  }


  cat("\nE(\\epsilon'\\epsilon) :\n")
  cat(sprintf("   Step 1: %s\n",
              format(mean(x$moments$moment_first_step$epsilon),
                     digits = digits)))
  if (attr(x, "approach") == "two_step"){
    cat(sprintf("   Step 2: %s\n",
                format(mean(x$moments$moment_optimum$epsilon),
                       digits = digits)))
  }

  cat("Number of iterations:\n ")
  cat(sprintf("  Step 1: %s\n",
              format(as.numeric(
                x$NelderMead$NM_step1$counts['function']
              ),
              digits = digits)))
  if (attr(x, "approach") == "two_step"){
    cat(sprintf("   Step 2: %s\n",
                format(as.numeric(
                  x$NelderMead$NM_step2$counts['function']
                ),
                digits = digits)))
  }


  cat("---------------------------\n")

  cat("   Estimated parameters:\n")

  std <- x$estimates$se_theta_hat

  formated_coef <- paste0(
    format(x$estimates$theta_hat, digits = digits),
    signif_stars_vectorized(2*pnorm( -abs(x$estimates$theta_hat/std)),
                            type = "none"),
    sprintf(" (%s)",format(std, digits = digits))
  )
  names(formated_coef) <- names(x$estimates$theta_hat)

  print.default(
    formated_coef,
    print.gap = 2,
    quote = FALSE
  )


  if (attr(x, "approach") == "two_step"){
    cat("\n First step estimated parameters:\n")
    print.default(
      x$estimates$theta_1,
      digits = digits,
      print.gap = 2,
      quote = FALSE
    )
  }

  invisible(x)
}


