resourcePlot <- function() {
  # Produce pdf plotting memory and runtime as a function of genome size.
  RT <- read.table("runtime.txt", header=TRUE)
  RT <- subset(RT, size > 20000)
  RT$size <- RT$size / 1000000
  RT$mem <- RT$mem
  
  #setEPS()
  #postscript("RuntimeFig.eps")
  #pdf(file="RuntimeFig.pdf", height=7, width=14)
  par(mfrow=c(1,2))
  
  # Make RS table.  (Need to add lmer and RS times/mem.)
  RS = subset(RT, tool == "RepeatScout")
  LM = subset(RT, tool == "lmer")
  RSLM = merge(RS, LM, by = c("label","size","seed","f"), suffixes = c(".RS", ".LM"))  
  
  PH0 <- subset(RT, tool == "phRAIDER" & seed == 0)
  PH1 <- subset(RT, tool == "phRAIDER" & seed == 0)
  PPH1 <- subset(RT, tool == 'pre-phRAIDER' & seed == 0)
  return(list(RSLM,PH1, PPH1))
  
  # First: runtimes for phRAIDER and RptScout
  plot(c(), c(), xlim = range(c(RS$size, PH1$size)), ylim = range(c(RSLM$time.RS, RSLM$timeLM, PH1$time)), xlab = "Genome size (Mb)", ylab = "Runtime (seconds)", main = "Runtime: phRAIDER v. RptScout")
  with(RSLM, points(size, time.RS + time.LM, col = 'red'))
  with(PH1, points(size, time, col = 'blue', pch = 2))
  with(PPH1, points(size, time, col = 'green', pch = 3))
  legend("topleft", c("RepeatScout", "phRAIDER", "phRAIDER.LM"), col = c('red', 'blue', 'green'), pch = c(1,2,3))
  
  # Second: memory usage for phRAIDER and RptScout
  plot(c(), c(), xlim = range(c(RS$size, PH1$size)), ylim = range(c(RSLM$mem.RS, RSLM$mem.LM, PH1$mem)), xlab = "Genome size (Mb)", ylab = "Memory required (Gb)", main = "Memory: phRAIDER v. RptScout")
  with(RSLM, points(size, mem.RS + mem.LM, col = 'red', pch = 1))
  with(PH1, points(size, mem, col = 'blue', pch = 2))
  with(PPH1, points(size, mem, col = 'green', pch = 3))
  #legend("topleft", c("RepeatScout", "phRAIDER", "phRAIDER.LM"), col = c('red', 'blue', 'green'), pch = c(1,2,3))
  dev.off()
}

seedPerf.plot <- function(sorg = "ce10.chrV") {
  # Plot seed performance.
  # First column: sensitivity / specificity as a function of seed length
  # Second column: sensitivity / speficifity as a function of seed weight/length
  # Third plot: runtime / memory as a function of seed weight/length
  
  #pdf(file="PreformanceFig.pdf", height=4, width=8)
  par(mfcol = c(1,2))
 
  
  
  # First: w/l = (sensitivity)
  T <- subset(read.table("seeds1.stats.txt", header=TRUE), tool == "phRAIDER" & org==sorg)   # TODO: NEED TO CHAGE THIS!!!
  R <- subset(read.table("seeds1.stats.txt", header=TRUE), tool == "RepeatScout" & org==sorg)

  xlim = range(T$l, na.rm = TRUE)
  ylim = with(T, range(c(tpr, QuCoverage, R$tpr, R$QuCoverage)), na.rm = TRUE)
  plot(c(), c(), xlim = xlim, ylim = ylim, xlab = "seed length", ylab = "Sensitivity", main = "density = 1")
  points(T$l, T$tpr, col = 'blue', pch = 1)
  points(T$l, T$QuCoverage, col = 'red', pch = 2)
  abline(h=R$tpr, col = 'blue')
  abline(h=R$QuCoverage, col = 'red')
  
  # Second: w/l = 1 (specificity)
  #xlim = range(T$l, na.rm = TRUE)
  #ylim = with(T, range(c(tnr, ConCoverage, R$tnr, R$ConCoverage)), na.rm = TRUE)
  #plot(c(), c(), xlim = xlim, ylim = ylim, xlab = "seed length", ylab = "Specificity", main = "density = 1")
  #points(T$l, T$tnr, col = 'blue')
  #points(T$l, T$ConCoverage, col = 'red', pch = 2) 
  #abline(h=R$tnr, col = 'blue', pch = 1)
  #abline(h=R$ConCoverage, col = 'red', pch = 2)

  #legend("topright", c("Masking", "RepeatMasker"), col = c('red', 'blue'), pch = c(2,1))
  
  # Third: Plot sensitivity on w/l for w==32
  T <- subset(read.table("seeds32.stats.txt", header=TRUE), tool == "phRAIDER" & org==sorg)
  R <- subset(read.table("RptSct.stats.txt", header=TRUE), tool == "RepeatScout" & org==sorg)
  
  xlim = range(T$w.l, na.rm = TRUE)
  ylim = with(T, range(c(tpr, QuCoverage, R$tpr, R$QuCoverage)), na.rm = TRUE)
  plot(c(), c(), xlim = xlim, ylim = ylim, xlab = "seed density", ylab = "Sensitivity", main = "weight = 32")
  points(T$w.l, T$tpr, col = 'blue', pch = 1)
  points(T$w.l, T$QuCoverage, col = 'red', pch = 2)
  abline(h=R$tpr, col = 'blue')
  abline(h=R$QuCoverage, col = 'red')
  
  # Fourth: Plot specificity on w/l for w==32
  #T <- subset(read.table("seeds32.stats.txt", header=TRUE), tool == "phRAIDER" & org==sorg)
  
  #xlim = range(T$w.l, na.rm = TRUE)
  #ylim = with(T, range(c(tnr, ConCoverage, R$tnr, R$ConCoverage)), na.rm = TRUE)
  #plot(c(), c(), xlim = xlim, ylim = ylim, xlab = "seed density", ylab = "Specificity", main = "weight = 32")
  #points(T$w.l, T$tnr, col = 'blue', pch = 1)
  #points(T$w.l, T$ConCoverage, col = 'red', pch = 2)
  #abline(h=R$tnr, col = 'blue')
  #abline(h=R$ConCoverage, col = 'red')
  
  
  # Sixth: Plot mem usage on w/l
  #RT <- read.table("runtime2.txt", header=TRUE)
  #RT$mem <- RT$mem / 1000000
  #T1 <- subset(RT, w==40)
  #with(T1, plot(w/l, time, xlim = range(w/l), ylim = range(time), xlab = "w/l (w = 40)", ylab = "Runtime (seconds)", main = "Runtime v. seed density"))
  #with(T1, plot(w/l, mem, xlim = range(w/l), ylim = range(mem), xlab = "w/l (w = 40)", ylab = "Memory (Gb)", main = "Memory usage v. seed sensity"))  
  
  title("\nphRAIDER applied to c. elegans (chr. V)", outer = TRUE)
  #dev.off()
}

seed.runtime <- function() {
  pdf(file="Runtime2.pdf", height=8, width=8)
  T <- read.table("runtime.long.txt", header = TRUE)
  T$mem <- T$mem / 1000000
  T1 <- subset(T, d == 1)
  T2 <- subset(T, len == 40)
  
  par(mfcol = c(2,2))
  l <- lm(T1$time ~ T1$len)
  r = cor(T1$len, T1$time)**2
  s = sprintf("%s = %5.3f", expression("r^2"), r)
  plot(T1$len, T1$time, xlab = "Seed weight (density = 1)", ylab = "Runtime (s)", main = "runtime v. seed length", 
       sub = s)
  abline(l, col = 'red')
  
  l <- lm(T2$time ~ T2$d)
  r = cor(T2$d, T2$time)**2
  s = sprintf("%s = %5.3f", expression("r^2"), r)
  plot(T2$d, T2$time, xlab = "Seed density (weight = 40)", ylab = "Runtime (s)", main = "runtime v. seed density", 
       sub = s)
  abline(l, col = 'red')
  
  l <- lm(T1$mem ~ T1$len)
  r = cor(T1$len, T1$mem)**2
  s = sprintf("%s = %5.3f", expression("r^2"), r)
  plot(T1$len, T1$mem, xlab = "Seed weight (density = 1)", ylab = "Memory (Gb)", main = "memory v. seed length",
        sub = s)
  abline(l, col = 'red')

  l <- lm(T2$mem ~ T2$d)
  r = cor(T2$d, T2$mem)
  s = sprintf("%s = %5.3f", expression("r^2"), r)
  plot(T2$d, T2$mem, xlab = "Seed density (width = 40)", ylab = "Runtime (Gb)", main = "memory v. seed density",
        sub = s)
  abline(l, col = 'red')
  dev.off()
  #return(list(T1,T2))
}
  


