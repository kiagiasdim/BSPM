\name{aPTT}
\alias{aPTT}
\docType{data}
\title{Dataset for PCC process for Normal with both parameters unknown
}

\description{
The activated Partial Thromboplastin Time (aPTT) dataset is a data frame that contains two columns. The first column (aPTT_current) repsesents 30 sequentially gathered current observations from the daily routine of a medical laboratory. The values refer to the activated Partial Thromboplastin Time (aPTT), measured in seconds. The second column (aPTT_historical) refers to 30 historical In Control measurements of the same variable. APTT is a blood test that characterizes coagulation of the blood. More specifically, it is the time that until the coagulation when reagents are added to plasma.
}

\usage{data("aPTT")}

\examples{
# Loading data
attach(aPTT)

# Plotting data
plot(1:length(aPTT_current), aPTT_current, type = "b", lty = 1, xlab = "", ylab = "",
     ylim = c(27.3, 33.4), xlim = c(-30, 30), lwd = 1, pch = 16, axes = FALSE, yaxs = "i",
     main = "aPTT dataset")

# x - axis for historical and current data
pastx <- c(-30, -20, -10, 0) ; currentx <- c(0, 10, 20, 30)
# Adding points
points(-length(aPTT_historical):(-1), aPTT_historical,
       type = "b", lty = 2, xlab = "", ylab = "", lwd = 1, pch = 21, col = "gray55")
# Adding axis with names
axis(2)
mtext("Current Data", side = 1, at = 15, line = 2.2, cex = 1.1)
mtext("Historical Data", side = 1, at = -15, line = 2.2, cex = 1.1, col = "gray55")
mtext("aPTT [sec]", side = 2, line = 2.2, cex = 1.1)
axis(1, at = currentx, labels = currentx)
axis(1, at = pastx, labels = pastx, col.axis = "gray55", col = "gray55", lty = 2)
segments(0, 27.5, 0, 33.5, lwd = 1)

}

