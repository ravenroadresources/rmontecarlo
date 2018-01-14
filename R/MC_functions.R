utils::globalVariables(c("nn", "x" , "y", "y_", "point", "flag"))


#' Generate Random Distribution
#'
#' @param n number of realizations
#' @param mean mean
#' @param sd standard deviation
#' @param lower lower limit of the truncated distribution
#' @param upper upper limit of the truncated distribution
#' @param method sampling method, either "RND" or "LHS". Default = "RND"
#' @param lhs vector of quatiles (obtained from LH sampling functions). Not required
#'     if \code{method = "RND"}. Default = NULL
#' @param distro distribution
#' @param seed integer serving as seed number, work only for the RND option
#'
#' @importFrom truncdist rtrunc qtrunc
#' @importFrom stats rlnorm rnorm runif qnorm qlnorm qunif
#' @importFrom mc2d rtriang qtriang
#'
#' @return an array with n random values for the given parameters
#'
#' @examples
#' \dontrun{
#' distribution_1 <- var_distro(100, 10, 4)
#' hist(distribution_1)
#'
#' distribution_2 <- var_distro(100, 10, 4, "truncated normal", 7, 13, seed = 123)
#' hist(distribution_2)
#' }
#'
#' @export
var_distro <- function(n, mean, sd = 0, lower = -Inf, upper = Inf,
                         method = "RND", lhs = NULL,
                         distro = "normal",
                         seed = 1234) {

  X <- seq(mean, mean, length.out = n) # constant value
  if (method == "RND") {
    set.seed <- seed
    if (distro == "normal") X <- rnorm(n, mean, sd)
    else if (distro == "truncated normal") X <- truncdist::rtrunc(spec = "norm", n, linf = lower, lsup = upper, mean = mean, sd = sd)
    else if (distro == "lognormal") X <- rlnorm(n, log(mean), sd)
    else if (distro == "truncated lognormal") X <- truncdist::rtrunc(spec = "lnorm", n, linf = lower, lsup = upper, mean = log(mean), sd = sd)
    else if (distro == "uniform") X <- runif(n, min = lower, max = upper)
    else if (distro == "triangular") X <- mc2d::rtriang(n, min = lower, mode = mean, max = upper)
  }
  else if (method == "LHS") {
    if (distro == "normal") X <- qnorm(lhs, mean, sd)
    else if (distro == "truncated normal") X <- truncdist::qtrunc(lhs, spec = "norm", mean = mean, sd = sd, a = lower, b = upper)
    else if (distro == "lognormal") X <- qlnorm(lhs, log(mean), sd)
    else if (distro == "truncated lognormal") X <- truncdist::qtrunc(lhs, spec = "lnorm", mean = mean, sd = sd, a = lower, b = upper)
    else if (distro == "uniform") X <- qunif(lhs, min = lower, max = upper)
    else if (distro == "triangular") X <- mc2d::qtriang(lhs, min = lower, mode = mean, max = upper)
  }

  return(X)
}



#' Aproximate Integral Solution by Monte Carlo Method
#'
#' Aproximate the solution of the integral of a function using Monte Carlo Method.
#'
#' @param func either a regular function of x, eithe a character string representing
#'     the function to be integrated. Examples: \code{function(x) x + 2} OR \code{"x + 2"}
#' @param interval length 2 array, giving the x limits of the inegral.
#' @param y_interval length 2 array, giving the y limits. Default is NULL: in this case the function
#'     calculate the minimum and maximum of the function whitin the \code{interval}.
#' @param n integer representing the number of repetition of the Monte Carlo method. Default = 1e4.
#' @param plot logical defining if the output should be plotted or not.
#'
#' @importFrom dplyr mutate %>%
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual
#'
#' @return a list fo 2 items: \code{$plot} which is the plot of the function and
#'     \code{$solution} with the value aproximating the exact solution of the integral.
#'
#' @examples
#' \dontrun{
#'   ff <- "x * sin(x^2)"
#'   xx <- c(0, 2 * sqrt(pi))
#'   yy <- c(-10, 10)
#'   nn <- 1e5
#'
#'   integrate_mc(ff, xx, n = nn)
#'
#'   sol <- integrate_mc(function(x) x * sin(x^2), xx, y_interval = yy, n = nn)
#'   sol[[2]]
#'   sol$solution
#' }
#'
#' @export
integrate_mc <- function(func, interval, y_interval = NULL, n = 1e4, plot = TRUE) {

  if(is.character(func)) ff <- function(x) eval(parse(text = func))
  else if(is.function(func)) ff <- func
  else stop("func nedds to be either a regular function of x or a character string.")

  if(is.null(y_interval)) {
    temp_x <- seq(interval[1], interval[2], length.out = 1001)
    yy <- c(min(ff(temp_x)), max(ff(temp_x)))

    if(is.infinite(yy[1])) yy[1] <- -100
    if(is.infinite(yy[2])) yy[2] <- 100
  } else {
    yy <- y_interval
  }

  df <- data.frame(nn = 1:n) %>%
    dplyr::mutate(
      x = runif(nn, interval[1], interval[2]),
      y = runif(nn, yy[1], yy[2]),
      y_ = ff(x),
      point = ifelse(y < y_ & y >= 0, 1,
                    ifelse(y > y_ & y <= 0, -1, 0)),
      flag = factor(point, levels = sort(unique(point), decreasing = TRUE)))

  sol <- sum(df$point) / n * diff(interval) * diff(yy)

  if(plot) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point(aes(color = flag), alpha = 0.8) +
      ggplot2::scale_color_manual(values = c("red", "grey", "blue")) +
      ggplot2::geom_line(data = data.frame(x = temp_x, y = ff(temp_x)))

    return(list(plot = p, solution = sol))
  }
  else return(list(plot = NULL, solution = sol))
}

