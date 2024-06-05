

mcxs1  = "#05386B"
mcxs2  = "#379683"
mcxs3  = "#5CDB95"
mcxs4  = "#8EE4AF"
mcxs5  = "#EDF5E1"
purple = "#b02442"

mcxs1.rgb   = col2rgb(mcxs1)
mcxs1.shade1= rgb(mcxs1.rgb[1],mcxs1.rgb[2],mcxs1.rgb[3], alpha=80, maxColorValue=255)
mcxs2.rgb   = col2rgb(mcxs2)
mcxs2.shade1= rgb(mcxs2.rgb[1],mcxs2.rgb[2],mcxs2.rgb[3], alpha=120, maxColorValue=255)
mcxs3.rgb   = col2rgb(mcxs3)
mcxs3.shade1= rgb(mcxs3.rgb[1],mcxs3.rgb[2],mcxs3.rgb[3], alpha=120, maxColorValue=255)
mcxs4.rgb   = col2rgb(mcxs4)
mcxs4.shade1= rgb(mcxs4.rgb[1],mcxs4.rgb[2],mcxs4.rgb[3], alpha=100, maxColorValue=255)


S         = 100000
s         = 10
ig2_prior = s/rchisq(S, 5)
g_prior   = rgamma(S, 0.5, scale = 2*s)

pdf(file = "prior_sigma2eta.pdf", width = 7, height = 5)
hist(ig2_prior, breaks = 500, freq = FALSE, ylim = c(0, 0.4), xlim = c(0, 40), col = mcxs1.shade1, border = mcxs1.shade1, main = expression(p(sigma[eta]^2)), xlab = expression(sigma[eta]^2), ylab = "")
hist(g_prior, breaks = 500, freq = FALSE, add = TRUE, col = mcxs4.shade1, border = mcxs4.shade1)
legend(25, 0.35, legend = c("IG2(10, 5)", "G(20, 0.5"), col = c(mcxs1.shade1, mcxs4.shade1), lwd = 8, bty = "n")
dev.off()
