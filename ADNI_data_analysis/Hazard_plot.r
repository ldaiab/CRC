##########################################################
#Protein
##########################################################
#srun -t 8:00:00 -p interact -N 1 -n 1 --mem=32gb --pty /bin/bash
cd /overflow/projects/KernelLogrankTest-master/data
module load r/4.4.0;R
rm(list = ls())
library(MASS)
library(reticulate)
library(RANN)
library(FOCI)
library(cubature)
library(matrixStats)
library(dcortools)
library(RcppArmadillo)
library(pracma)
library(survival)
library(ggplot2)
setwd('./data') #or 'PseudoRealData/'
A2=read.csv('surv.csv');M1=read.csv('protein.csv')
##########################################################
A2$B7H1=M1[,5186]
A2$VEGF=M1[,4344]
mfit_B7H1 <- coxph(Surv(time, censor) ~ gender + age + I(age^2) + pspline(B7H1, df=4), data = A2)
mfit_VEGF <- coxph(Surv(time, censor) ~ gender + age + I(age^2) + pspline(VEGF, df=4), data = A2)
tp_B7H1 <- termplot(mfit_B7H1, term = 4, se = TRUE, plot = FALSE)
tp_VEGF <- termplot(mfit_VEGF, term = 4, se = TRUE, plot = FALSE)
normalize <- function(x) (x - min(x)) / (max(x) - min(x))
prepare_data <- function(tp, var_name) {
  data.frame(
    x = normalize(tp[[var_name]]$x),
    y = tp[[var_name]]$y,
    se = tp[[var_name]]$se
  )
}

B7H1_data <- prepare_data(tp_B7H1, "B7H1")
VEGF_data <- prepare_data(tp_VEGF, "VEGF")

# Determine y-axis limits for the plot
y_ranges <- range(c(
                    B7H1_data$y - 1.96*B7H1_data$se,
                    B7H1_data$y + 1.96*B7H1_data$se,
                    VEGF_data$y - 1.96*VEGF_data$se,
                    VEGF_data$y + 1.96*VEGF_data$se))

# Create the base plot
pdf('effect_VEGF.pdf')
plot(NA, xlim = c(0, 1), ylim = y_ranges, 
     xlab = "Normalized VEGF", ylab = "Log Hazard Ratio",cex=)

# Function to add each term with confidence intervals
add_term <- function(data, color) {
  lines(data$x, data$y, col = color, lwd = 2)
  polygon(c(data$x, rev(data$x)), 
          c(data$y + 1.96*data$se, rev(data$y - 1.96*data$se)), 
          col = adjustcolor(color, alpha.f = 0.2), border = NA)
}

# Add each spline to the plot
#add_term(B7H1_data, "red")
add_term(VEGF_data, "blue")

# Add a legend
#legend("topright", legend = c( "B7H1", "VEGF"),
#       col = c("red", "blue"), lty = 1, lwd = 2, bty = "n")
dev.off()
