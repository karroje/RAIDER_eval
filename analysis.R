ce1 = "ce10.chrI.sim"
ce2 = "ce10.chrII.sim"
ce3 = "ce10.chrIII.sim"
ce4 = "ce10.chrIV.sim"
ce5 = "ce10.chrV.sim"
tair5 = "tair10.chr5.sim"
hg22 = "hg38.chr22.sim"
mm19 = "mm10.chr19.sim"

best_file = "SEEDS1/stats.best1.txt"

#  A plot by length with density of 1
plot.by.length <- function(sorg, file = "SEEDS1/stats.txt") {
  par(mfrow=c(1,2))
  T <- subset(read.table(file, header=TRUE), tool == "phRAIDER" & org==sorg & w.l == 1)  
  R <- subset(read.table(file, header=TRUE), tool == "RepeatScout" & org==sorg)
  
  xlim = range(T$l, na.rm = TRUE)
  ylim = with(T, range(c(tpr, QuCoverage, R$tpr, R$QuCoverage), na.rm = TRUE))
  plot(c(), c(), xlim = xlim, ylim = ylim, xlab = "seed length", ylab = "Sensitivity", main = "density = 1")
  points(T$l, T$tpr, col = 'blue', pch = 1)
  points(T$l, T$QuCoverage, col = 'red', pch = 2)
  abline(h=R$tpr, col = 'blue')
  abline(h=R$QuCoverage, col = 'red')

  
  xlim = range(T$l, na.rm = TRUE)
  ylim = with(T, range(c(tnr, ConCoverage, R$tnr, R$ConCoverage)), na.rm = TRUE)
  plot(c(), c(), xlim = xlim, ylim = ylim, xlab = "seed length", ylab = "Speficitiy", main = "density = 1")
  points(T$l, T$tnr, col = 'blue', pch = 1)
  points(T$l, T$ConCoverage, col = 'red', pch = 2)
  abline(h=R$tnr, col = 'blue')
  abline(h=R$ConCoverage, col = 'red')
  
  abline(v = 32)
  abline(v = 40)
}

plot.by.density <- function(sorg, file = "SEEDS1/stats.txt") {
  par(mfrow=c(1,2))
  T <- subset(read.table(file, header=TRUE), tool == "phRAIDER" & org==sorg & w.l != 1)  
  #T <- T[is.na(T$tpr),]
  R <- subset(read.table(file, header=TRUE), tool == "RepeatScout" & org==sorg)
  
  xlim = range(T$w.l, na.rm = TRUE, na.rm = TRUE)
  ylim = with(T, range(c(tpr, QuCoverage, R$tpr, R$QuCoverage), na.rm = TRUE))
  plot(c(), c(), xlim = xlim, ylim = ylim, xlab = "density", ylab = "Sensitivity", main = "length = 40")
  points(T$w.l, T$tpr, col = 'blue', pch = 1)
  points(T$w.l, T$QuCoverage, col = 'red', pch = 2)
  abline(h=R$tpr, col = 'blue')
  abline(h=R$QuCoverage, col = 'red')
  
  xlim = range(T$w.l, na.rm = TRUE)
  ylim = with(T, range(c(tnr, ConCoverage, R$tnr, R$ConCoverage), na.rm = TRUE))
  plot(c(), c(), xlim = xlim, ylim = ylim, xlab = "density", ylab = "Speficitiy", main = "length = 40")
  points(T$w.l, T$tnr, col = 'blue', pch = 1)
  points(T$w.l, T$ConCoverage, col = 'red', pch = 2)
  abline(h=R$tnr, col = 'blue')
  abline(h=R$ConCoverage, col = 'red')
  
}


plotFig2 <- function(file = "SEEDS1/stats.txt", sorg = ce5) {
  #setEPS()
  #postscript("fig2.eps")
  pdf(file="manuscript/JCB/PreformanceFig.pdf", height=7, width=14)
  par(mfrow=c(1,2))
  
  T <- subset(read.table(file, header=TRUE), tool == "phRAIDER" & org==sorg & w.l == 1)  
  R <- subset(read.table(file, header=TRUE), tool == "RepeatScout" & org==sorg)
  
  xlim = range(T$l, na.rm = TRUE)
  #ylim = with(T, range(c(tpr, QuCoverage, R$tpr, R$QuCoverage), na.rm = TRUE))
  ylim = range(0.6, 0.8)
  plot(c(), c(), xlim = xlim, ylim = ylim, xlab = "seed length", ylab = "Sensitivity", main = "density = 1")
  points(T$l, T$tpr, col = 'blue', pch = 1)
  points(T$l, T$QuCoverage, col = 'red', pch = 2)
  abline(h=R$tpr, col = 'blue')
  abline(h=R$QuCoverage, col = 'red')
  
#   xlim = range(T$l, na.rm = TRUE)
#   ylim = with(T, range(c(tnr, ConCoverage, R$tnr, R$ConCoverage)), na.rm = TRUE)
#   plot(c(), c(), xlim = xlim, ylim = ylim, xlab = "seed length", ylab = "Speficitiy", main = "density = 1")
#   points(T$l, T$tnr, col = 'blue', pch = 1)
#   points(T$l, T$ConCoverage, col = 'red', pch = 2)
#   abline(h=R$tnr, col = 'blue')
#   abline(h=R$ConCoverage, col = 'red')
  
  T <- subset(read.table(file, header=TRUE), tool == "phRAIDER" & org==sorg & w.l != 1)  
  #T <- T[is.na(T$tpr),]
  R <- subset(read.table(file, header=TRUE), tool == "RepeatScout" & org==sorg)
  
  xlim = range(T$w.l, na.rm = TRUE, na.rm = TRUE)
  #ylim = with(T, range(c(tpr, QuCoverage, R$tpr, R$QuCoverage), na.rm = TRUE))
  ylim = range(0.6, 0.8)
  plot(c(), c(), xlim = xlim, ylim = ylim, xlab = "density", ylab = "Sensitivity", main = "length = 40")
  points(T$w.l, T$tpr, col = 'blue', pch = 1)
  points(T$w.l, T$QuCoverage, col = 'red', pch = 2)
  abline(h=R$tpr, col = 'blue')
  abline(h=R$QuCoverage, col = 'red')
  
#   xlim = range(T$l, na.rm = TRUE)
#   ylim = with(T, range(c(tnr, ConCoverage, R$tnr, R$ConCoverage)), na.rm = TRUE)
#   plot(c(), c(), xlim = xlim, ylim = ylim, xlab = "seed length", ylab = "Speficitiy", main = "density = 1")
#   points(T$l, T$tnr, col = 'blue', pch = 1)
#   points(T$l, T$ConCoverage, col = 'red', pch = 2)
#   abline(h=R$tnr, col = 'blue')
#   abline(h=R$ConCoverage, col = 'red')
  
  
  dev.off()
}
