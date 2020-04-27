\name{HPRD_Normal_BothUnknown}
\alias{HPRD_Normal_BothUnknown}

\title{HPRD for Normal process with both parameters unknown
}
\description{
HPRD_Normal_BothUnknown is used to calculate the Highest Predictive Density (HPrD) for the PCC process for Normal data, asssuming both the mean and the variance unknown.
}
\usage{
HPRD_Normal_BothUnknown(far, Mpr, Lp, Ap, Bp)
}

\arguments{
  \item{far}{ scalar (between 0 and 1); the false alarm rate, which reflects to an 100*(1-far)\% coverage for the HPrD. The default value is 1/370.4.
}
  \item{Mpr}{scalar; posterior predictive parameter.
}
  \item{Lp}{scalar (positive); posterior predictive parameter.
}
  \item{Ap}{scalar (positive); posterior predictive parameter.
}
  \item{Bp}{scalar (positive); posterior predictive parameter.
}
}
\details{
HPRD_Normal_BothUnknown provides (1-far)\% the HPrD region of the PCC process for Normal data, assumming both the mean and the variance unknown.

More specifically the function provides the HPrD of a non-standardized Student's t-distribution with location parameter Mpr and scale parameter  sqrt( Bp * (Lp+1) / (Ap * Lp ) ).

}

\examples{
HPRD_Normal_BothUnknown(far = 0.0017, Mpr = 30, Lp = 12, Ap = 6, Bp = 5)
}

