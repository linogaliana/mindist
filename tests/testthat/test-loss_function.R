context("test-loss_function")


# true_theta <- c(4,2)
#
# loss <- function(est_theta, true_theta = c(4,2), weights = NULL){
#   x <- rnorm(1000, mean = est_theta[1], sd = est_theta[2])
#   epsilon <- c(mean(x) - true_theta[1], sd(x) - true_theta[2])
#   # UNWEIGHTED GMM ESTIMATES =============
#   if (is.null(weights)){
#     W <- diag(length(epsilon))
#   } else{
#     W <- optimal_weight_matrix(epsilon)
#   }
#   return(t(epsilon) %*% W %*% epsilon)
# }
#
# NM_step1 <- optim(loss, par = c(0,1), weights = NULL)
#
# W <- optimal_weight_matrix(NM_step1$par - true_theta)
#
#
# NM_step2 <- optim(loss, par = NM_step1$par, weights = W)
#
#
#
# epsilon_f <- function(est_theta, true_theta = c(4,2)){
#   x <- rnorm(1000, mean = est_theta[1], sd = est_theta[2])
#   epsilon <- c(mean(x) - true_theta[1], sd(x) - true_theta[2])
#   return(epsilon)
# }
#
# approx_jacobian_epsilon <- function(theta,
#                                     step = 1e-6,
#                                     weights = 1L, N_moments = 2L){
#
#
#
#
#   all_theta <- lapply(seq_along(theta), function(i) list(
#     c(theta[i]+step, theta[-i]),
#     c(theta[i]-step, theta[-i]))
#   )
#
#   all_theta <- unlist(all_theta, recursive = FALSE)
#
#   grad <- lapply(all_theta, function(param) ({
#     eps <- epsilon_f(est_theta = param, true_theta = c(4,2))
#     return(
#       eps
#     )
#   }))
#
#   names(theta) <- c("mu","sigma")
#
#   names(grad) <- do.call(
#     c,
#     lapply(names(theta), function(g) paste0(g, 1:2))
#   )
#
#
#   H <- lapply(names(theta), function(param) ({
#     (grad[[paste0(param,"1")]] - grad[[paste0(param,"2")]])/(2*step)
#   })
#   )
#
#   names(H) <- names(theta)
#
#   data.table::setDT(H)
#
#   return(H/N_moments)
# }
#
#
#
# Gamma <-  approx_jacobian_epsilon(NM_step2$par)
#
#
# # We must use Gamma = Jacobian/N_moments
# Gamma <- Gamma/2
#
#
#
# df_moment_optim <- epsilon_f(est_theta = NM_step2$par)
#
#
# Delta <- optimal_weight_matrix(df_moment_optim)
#
# matinv <- solve(t(as.matrix(Gamma)) %*% W %*% as.matrix(Gamma))
#
# Vtheta <- matinv %*% t(as.matrix(Gamma)) %*% W %*% Delta %*% W %*% as.matrix(Gamma) %*% matinv
# Vtheta <- Vtheta/2
#
# se_estimator <- sqrt(diag(Vtheta))
#
#
# x <- rnorm(1000, mean = true_theta[1], sd = true_theta[2])
#
# g1 <- function(tet,x){
#   m1 <- (tet[1]-x)
#   m2 <- (tet[2]^2 - (x - tet[1])^2)
#   f <- cbind(m1,m2)
#   return(f)
# }
#
# gmm::gmm(g1, x1, c(mu = 0, sig = 0), wmatrix = "ident")
# NM_step1$par
# summary(
#   gmm::gmm(g1, x1, c(mu = 0, sig = 0), wmatrix = "ident")
# )
