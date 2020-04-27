\name{ECE}
\alias{ECE}
\docType{data}
\title{ECE dataset for the PCC process for Poisson with rate parameter unknown
}

\description{
The Electrical Complex Equipment (ECE) dataset repsesents 25 paired observations in a data frame. The first column (defect_counts) refers the number of defects per inspected number of units, encountered in a complex electrical equipment of an assembly line. The second column (inspected_units) contains the corresponding inspected number of units for the counting of defects. The data were also analyzed by Hansen and Ghare (1987) and Bayarri and Garcia-Donato (2005).}

\usage{data("ECE")}

\references{
[1] Bayarri, M. J., and Garcia-Donato, G. (2005), "A Bayesian Sequential Look at u-Control Charts", Technometrics, 47, 2, pp. 142-151

[2] Hansen, B., and Ghare, P. (1987), "Quality Control and Application, Prentice-Hall", Englewood Cliffs, NJ
}

\examples{
# Loading data
attach(ECE)

# Plotting data
par( oma = c(1,3,2,3) )
plot(1:length(defect_counts), defect_counts/inspected_units, type = "b", lty = 1,
     xlab = "Observations", ylab = "", xlim = c(0, 25), ylim = c(1.5, 10.5),
     lwd = 1, pch = 16, axes = FALSE, yaxs = "i", main = "ECE dataset")

# Adding points
points(1:length(defect_counts), inspected_units, type = "b", lty = 2, lwd = 1, pch = 21, col = "gray55")
# Adding legend
legend("topleft", legend=c(expression(paste(s[i])), expression(paste(x[i]/s[i])) ), bty = "n",
        cex = 0.8, lty = c(2, 1), lwd = 1, col = c ("gray55", "black") , pch = c(21, 16))
# Adding axis with names
axis(1) ; axis(2) ; axis(4, col.axis = "gray55", col = "gray55")
mtext("Number of Defects \n per unit", side = 2,  line = 2.2, cex = 0.9)
mtext("Inspected units", side = 4,  line = 2, cex = 0.9, col = "gray55")

}

