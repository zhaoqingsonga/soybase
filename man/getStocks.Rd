% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{getStocks}
\alias{getStocks}
\title{This is a function that can get stock data from yahoo web}
\usage{
getStocks(
  stock = "000001.SS",
  start = "1979-01-01",
  end = Sys.Date(),
  freq = "daily",
  offset = 28800
)
}
\arguments{
\item{stock}{is the name of stock.}

\item{start}{is the start day, for example "2022-01-01".}

\item{end}{is the end day, for example "2022-05-01".}

\item{offset}{set default 28800}

\item{filename}{is t, for example "2022-05-01".}
}
\description{
This is a function that can get stock data from yahoo web
}
