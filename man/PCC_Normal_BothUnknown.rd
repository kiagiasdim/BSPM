\name{PCC_Normal_BothUnknown}
\alias{PCC_Normal_BothUnknown}

\title{PCC for Normal with both parameters unknown}

\description{
PCC_Normal_BothUnknown is used to fit the PCC process for Normal data, assuming both the mean and the variance unknown.
}

\usage{
PCC_Normal_BothUnknown(data = NULL, historical_data = NULL,
                       mu0 = 0, l0 = 0, a0 = -1/2, b0 = 0, alpha_0 = NULL,
                       ARL_0 = 370.4, FAP = NULL, FIR = FALSE, fFIR = .99,
                       aFIR=1/8, summary_list = TRUE, PCC_PLOT = TRUE,
                       historical_data_PLOT = FALSE, pdf_report = FALSE,
                       xlab = "Observations", ylab = "Values", main="PCC Normal")
}

\arguments{
  \item{data}{ vector; a dataset for PCC implementation. Data needs to be in a vector form.
}
  \item{historical_data}{vector; an optional dataset of historical data. Historical data needs to be in a vector form.
}
  \item{mu0}{ scalar; prior hyperparameter. It needs to be a number. The default is 0 and it refers to the initial reference prior.
}
  \item{l0}{ scalar (positive); prior hyperparameter. It needs to be a number. The default is 0 and it refers to the initial reference prior.
}
  \item{a0}{ scalar (positive); prior hyperparameter. It needs to be a number. The default is -1/2 and it refers to the initial reference prior.
}
  \item{b0}{ scalar (positive); prior hyperparameter. It needs to be a number. The default is 0 and it refers to the initial reference prior.
}
  \item{alpha_0}{ scalar (non negative); It is a power prior parameter controls the influence of the historical data on the posterior distribution. The default is 1/n_0, where n_0 is the size of the historical data.
}
  \item{ARL_0}{scalar (positive); In Control (IC) Average Run Length (ARL). It is average number of IC data points that we will plot in the PCC before a false alarm occurs. The default value is 370.4
}
  \item{FAP}{scalar (between 0 and 1); False Alarm Probability (FAP). It is the probability of raising at least one false alarm out of a pre-determined number of N hypothesis tests and it can be used instead of ARL_0. It is based on the Sidak's correction.
}
  \item{FIR}{logical; If TRUE, then the Fast Initial Response (FIR) PCC is applied, which is an adjustment for the fisrt tests by narrowing the PCC limits.
}
  \item{fFIR}{a number between 0 and 1; It is used if FIR=TRUE. The default value is 0.99 and represents the proportion of the adjusted PCC region over the initial one for the first test.
}
  \item{aFIR}{non-negative number; It is used if FIR=TRUE. The default value is 0.125 and it is a smoothing parameter for the FIR adjustment.
}
  \item{summary_list}{logical; If it is TRUE, then a data frame is provided, containing the data sequence, the PCC limits and the occurence of an alarm. It is TRUE by default.
}
  \item{PCC_PLOT}{logical; if TRUE, the PCC plot is displayed. It is TRUE by default.
}
 \item{historical_data_PLOT}{logical; if TRUE, the historical data are plotted along with the PCC plot is displayed. It is FALSE by default.
}
 \item{pdf_report}{logical; if TRUE then the summary list and the PCC plot reported in a pdf file.
}
  \item{xlab, ylab, main}{The titles of the x-axis, y-axis and the overall title for the PCC plot.
}
}

\details{

PCC_Normal_BothUnknown provides the PCC process for Normal data, assumming both the mean and the variance unknown. The PCC process is based on the sequential testing of the future observable against the Highest Predictive Density (HPrD), which is obtained by the posterior predictive distribution. The first test is available for the second observation.

The joint initial prior for the unknown parameters is a NIG(mu0,l0,a0,b0). Furthermore, the direct use of historical data is possible via the power prior, if they are available. In this case, the default value for the power prior parameter alpha_0 is the reciprocal of the length of the historical data, which conveys the weight of a single observation to the prior information. The default prior is the non-informative reference prior NIG(0,0,-1/2,0), without the use of historical data. In this case, the first test is available for the third observation.

A FIR option is available by narrowing the first few control limits. The metrics that can be used for the false alarms tolerance are ARL_0 and FAP.

}

\examples{
# 30 Normal observations introducing an outlier of 3*sd at the 15th observation
set.seed(1234)
out <- rnorm(30)
out[15] <- out[15] + 3
PCC_Normal_BothUnknown(out)

# Real data application
attach(aPTT)
PCC_Normal_BothUnknown(data = aPTT_current, historical_data = aPTT_historical)
}


