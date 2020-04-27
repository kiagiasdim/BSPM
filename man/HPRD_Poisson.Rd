\name{HPRD_Poisson}
\alias{HPRD_Poisson}

\title{HPRD for Poisson process
}

\description{
HPRD_Poisson is used to calculate the Highest Predictive Density (HPrD) for the PCC process for Poisson data, asssuming the rate parameter unknown.
}

\usage{
HPRD_Poisson(far, R, P)
}

\arguments{
   \item{far}{ scalar (between 0 and 1); the false alarm rate, which reflects to an 100*(1-far)\% coverage for the HPrD. The default value is 1/370.4.
}
  \item{R}{scalar (positive); posterior predictive parameter.
}
  \item{P}{scalar (between 0 and 1); posterior predictive parameter.
}
}
\details{
HPRD_Poisson provides 100*(1-far)\% the HPrD region of the PCC process for Poisson data, assumming the rate parameter unknown.

More specifically the function provides the HPrD of a negative Binomial distribution with number of failures R and failure probability P.}

\examples{
HPRD_Poisson(far = 0.0017, R = 10, P = 0.9)
}

