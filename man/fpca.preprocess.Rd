% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fpca.preprocess.R
\name{fpca.preprocess}
\alias{fpca.preprocess}
\title{fpca.preprocess}
\usage{
fpca.preprocess(
  data,
  workinggrid = NULL,
  observationgrid = NULL,
  method = c("smoothsplines", "pace", "naturalsplines")
)
}
\arguments{
\item{data}{A multivariate time series of class 'ts' or 'data.frame'.}

\item{workinggrid}{An optional list, sequence or numeric vector which contains a predetermined workinggrid for the analysis.
If NULL, the sequence with gridsize 1 is used, ranging from the smallest to the highest observed domain point. NULL is default.}

\item{observationgrid}{An optional list, sequence or numeric vector which contains a predetermined observationgrid for the
observations. If NULL, the column names will be used. NULL is default.}

\item{method}{Preprocessing method to be used. If 'smoothsplines', smoothing splines from the package 'fda' are used (see Kokoszka
and Reimherr, 2017). If 'pace', the PACE method from the package 'fdapace' is applied (see Yao et al., 2005).
If 'natralsplines', natural spline interpolation is used.
In case of irregularly observed functional data please set to 'pace' or 'naturalsplines'.
'smoothsplines' is default.}
}
\value{
fpca.preprocess returns an object of class 'FPCAobj' with:
\item{eigenvalues}{A vector containing the eigenvalues of the sample covariance operator.}
\item{fpcascores}{A matrix/ts-matrix of the computed FPCA-scores/projection coefficients.}
\item{eigenfunctions.obsgrid}{A matrix containing the orthonormal eigenfunctions of the sample covariance operator, which are evaluated on the workinggrid.}
\item{meanfunction.obsgrid}{A vector containing the sample mean function, which is evalauted on the observationgrid.}
\item{eigenfunctions.workgrid}{A matrix containing the orthonormal eigenfunctions of the sample covariance operator, which are evaluated on the observationgrid.}
\item{meanfunction.workgrid}{A vector containing the sample mean function, which is evalauted on the workinggrid.}
\item{raw.data}{The raw data used for the analysis.}
\item{observationgrid}{Observationgrid of the raw data used for the analysis.}
\item{workinggrid}{The dense workinggrid of the analysis.}
\item{dense.fts}{The dense representation of the functional time series using the empirical Karhunen-Loeve expansion on the workinggrid.}
\item{method}{Preprocessing method to obtain the dense functional time series.}
}
\description{
fpca.preprocess converts discrete point observations of a functional time series into functional time series evaluated on a dense
workinggrid. It also computes its eigendecomposition and empirical Karhunen-Loeve expansion. The data is preprocessed using
either natural splines, smoothing splines, or the PACE-method.
}
\examples{
# standard workinggrid
fed = load.fed()
fpca.preprocess(data = fed)
# A tighter workinggrid
wg = (2:720)/2
fpca.preprocess(fed, workinggrid = wg)
}
\references{
\itemize{
\item Kokoszka P. and Reimherr M. (2017), "Introduction to Functional Data Analysis". CRC Press.
\item Otto S. and Salish N. (2022). "Dynamic Factor Model for Functional Time Series: Identification, Estimation, and Prediction". ...
\item Yao, F., Mueller, H.-G., and Wang, J.-L. (2005). Functional data analysis for sparse longitudinal
data. Journal of the American Statistical Association, 100:577-590.
}
}
