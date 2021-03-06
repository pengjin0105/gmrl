% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/newton.R
\name{newton}
\alias{newton}
\title{newton method}
\usage{
newton(data_order, z_mat, beta_ini, model = "exp", iter_max = 500,
  tol = 1e-04)
}
\arguments{
\item{data_order}{a ordered dataset by observed time}

\item{z_mat}{covariate matrix}

\item{beta_ini}{given initial parameters}

\item{model}{specify additive or proportional MRL model}

\item{iter_max}{iteration maximum}

\item{tol}{tolerance}
}
\value{
estimated value
}
\description{
newton method
}
\author{
Peng Jin and Mengling Liu
}
