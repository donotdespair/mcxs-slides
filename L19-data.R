############################################################
# Macroeconometrics: ECOM90007, ECOM40003
# prepared by Tomasz Woźniak
# R file for Lecture 19: Modeling trend inflation
############################################################

# Define colors
############################################################
mcxs1  = "#05386B"
mcxs2  = "#379683"
mcxs3  = "#5CDB95"
mcxs4  = "#8EE4AF"
mcxs5  = "#EDF5E1"
purple = "#b02442"

mcxs1.rgb       = col2rgb(mcxs1)
mcxs1.shade1    = rgb(mcxs1.rgb[1],mcxs1.rgb[2],mcxs1.rgb[3], alpha = 120, maxColorValue = 255)
mcxs2.rgb       = col2rgb(mcxs2)
mcxs2.shade1    = rgb(mcxs2.rgb[1],mcxs2.rgb[2],mcxs2.rgb[3], alpha = 120, maxColorValue = 255)


# Download and create series
############################################################
CPI_downloaded  = readabs::read_abs(series_id = "A2325846C")
CPI_tmp         = as.matrix(CPI_downloaded[,6]) 
cpi             = ts( log(CPI_tmp[111:295,]), start = c(1976, 1), frequency = 4)
pi              = ts( 100 * diff(log(CPI_tmp[107:295, ]), lag = 4), 
                      start = c(1976, 1), frequency = 4)
save(cpi, pi, file = "cpi_au.rda")

pdf("results/data-CPI-downloaded.pdf", height = 5, width = 8)
plot(ts(CPI_tmp, start = c(1948, 3), frequency = 4), bty = "n", lwd = 2, col = mcxs2, xlab = "", ylab = expression(CPI[t]))
dev.off()

# autocorrelation plots
############################################################
cpi_acf         = FinTS::Acf(cpi, lag.max = 20, plot = FALSE)
pi_acf          = FinTS::Acf(pi, lag.max = 20, plot = FALSE)

pdf("results/data-cpi.pdf", height = 5, width = 8)
par(mfrow = c(2,2),mar = c(4.2,4.2,1,1))
plot(cpi, bty = "n", lwd = 2, col = mcxs2, xlab = "", ylab = expression(cpi[t]))
plot(x = 1:20, y = cpi_acf$acf[2:21], ylim = c(0,1), type = "h", lwd = 5, col = mcxs2, bty = "n", xlab = "", ylab = expression(acf(cpi[t])))
plot(pi, bty = "n", lwd = 2, col = mcxs3, xlab = "time", ylab = expression(pi[t]))
abline(h = 0)
plot(x = 1:20, y =  pi_acf$acf[2:21], ylim = c(0,1), type = "h", lwd = 5, col = mcxs3, bty = "n", xlab = "lags", ylab = expression(acf(pi[t])))
dev.off()

# integration order
############################################################
cpi_ar   = ar(cpi, order.max = 24, aic = FALSE)
round(cpi_ar$ar, 3)
round(sqrt(diag(cpi_ar$asy.var.coef)), 3)

FinTS::AutocorTest(cpi_ar$resid, lag = 10)
FinTS::ArchTest(cpi_ar$resid, lags = 10)

fUnitRoots::adfTest(cpi, lags = 24, type = "ct")
fUnitRoots::adfTest(cpi, lags = 24, type = "c")
fUnitRoots::adfTest(pi, lags = 23, type = "c")
fUnitRoots::adfTest(pi, lags = 23, type = "nc")
fUnitRoots::adfTest(diff(pi), lags = 22, type = "nc")

