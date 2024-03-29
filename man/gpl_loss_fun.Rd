% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loss.R
\name{gpl_loss_fun}
\alias{gpl_loss_fun}
\title{Create generalized piecewise linear (gpl) scoring/loss function,
which in general need not be piecewise linear}
\usage{
gpl_loss_fun(g = "x", kappa = 1, alpha, O, U = NA, offset = 0)
}
\arguments{
\item{g}{a non-decreasing increment function}

\item{kappa}{scale factor}

\item{alpha}{normalized loss when outcome y exceeds forecast x}

\item{O}{cost when forecast x exceeds outcome y; equals kappa*(1-alpha)}

\item{U}{loss when outcome y exceeds forecast x; equals kappa*alpha}

\item{offset}{an added function of y, defaults to 0 for which gpl loss function \code{L} will have \code{L(x,x) = 0}}
}
\value{
function with arguments \code{x} and \code{y} giving loss
}
\description{
Create generalized piecewise linear (gpl) scoring/loss function,
which in general need not be piecewise linear
}
