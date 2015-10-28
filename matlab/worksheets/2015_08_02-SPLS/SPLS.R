#Prepare  the response matrices in MATLAB, then save as mat file, then
#import into here to run SPLS

#Load data
library("R.matlab")
library("spls")

a = readMat('./20130117SpankyUtah001.mat')
nU = 35;
data = a$data
spikes <- t(matrix(unlist(data[3]), ncol = 33367, byrow = TRUE))
cursor <- t(matrix(unlist(data[1]), ncol = 33367, byrow = TRUE))
torque <- t(matrix(unlist(data[4]), ncol = 33367, byrow = TRUE))

# SPLS with eta=0.9 & 8 hidden components
# Should set eta through cross validation
f <- spls( spikes, cursor, K=8, eta=0.9, kappa = 0.5, scale.y = TRUE, trace = TRUE, maxstep = 500 )
print(f)
# Print out coefficients
coef.f <- coef(f)
coef.f[,1]
# Coefficient path plot
pdf(file="./coeffspathplot.pdf")
plot( f, yvar=1 )
dev.off()
pdf(file="./coeffs.pdf")
# Coefficient plot of selected variables
coefplot.spls( f, xvar=c(1:4) )
dev.off()

#Test code
data(yeast)
# SPLS with eta=0.7 & 8 hidden components
fy <- spls( yeast$x, yeast$y, K=8, eta=0.7, kappa = 0.3, scale.y = TRUE, trace = TRUE )
#spls( spikes, cursor, K=8, eta=0.7, kappa = 0.3, scale.y = TRUE, trace = TRUE )
print(fy)
# Print out coefficients
coef.fy <- coef(fy)
coef.fy[,1]
# Coefficient path plot
pdf(file="./yeastpathplot.pdf")
plot( fy, yvar=1 )
dev.off()
pdf(file="./yeastcoeffs.pdf")
# Coefficient plot of selected variables
coefplot.spls( fy, xvar=c(1:4) )
dev.off()