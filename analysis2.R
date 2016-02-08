read.file <- function(file, pal, test_num) {
  T <- read.table(file, header = TRUE)
  T$pal = pal
  T$test_num = test_num
  return(T)
}

read.stats <- function() {
  stats1 <<- read.file("seeds1.stats.txt", FALSE, 1) # Vary length
  stats2 <<- read.file("seeds2.stats.txt", FALSE, 2)
  stats3 <<- read.file("seeds3.stats.txt", FALSE, 3)
  stats4 <<- read.file("seeds4.stats.txt", FALSE, 4)
  stats5 <<- read.file("seeds5.stats.txt", TRUE,  5)
  stats6 <<- read.file("seeds6.stats.txt", FALSE, 6)
  stats7 <<- read.file("seeds7.stats.txt", TRUE,  7)
  stats8 <<- read.file("seeds8.stats.txt", FALSE, 8)
  stats9 <<- read.file("seeds9.stats.txt", TRUE,  9)
  stats10 <<- read.file("seeds10.stats.txt", FALSE, 10)
  stats11 <<- read.file("seeds11.stats.txt", TRUE,  11)  
  stats12 <<- read.file("seeds12.stats.txt", FALSE, 12)
  stats13 <<- read.file("seeds13.stats.txt", FALSE, 13)
  stats15 <<- read.file("seeds15.stats.txt", FALSE, 15)
  stats.all <<-rbind(stats1, stats2, stats3, stats4, stats5, 
                     stats6, stats7, stats8, stats9, stats10, 
                     stats11, stats12, stats13)
  statsM1 <<- read.file("seeds1m.stats.txt", FALSE, 1)
  best <<- read.file("best.stats.txt", FALSE, -1)
}

add.legend <- function(position) {
  legend(position, c("RM", "PRA"), col = c('blue', 'green'), 
                     lty = c('solid', 'solid', 'solid', 'solid'))
}

plot.gen <- function(T, org, x_col, y_col, color, xlab = NA, ylab = NA, main = NA, 
                     xlim = NA, ylim = NA, legend_position = "topleft", new_plot = TRUE,
                     reverse_x = FALSE, reverse_y = FALSE) {
  T <- droplevels(T[T$org==org,])

  if (reverse_x) {
    xlim = rev(xlim)
  }
  if (reverse_y) {
    ylim = rev(ylim)
  }
  
  if (new_plot) {
    plot(c(), c(), xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main)
  }
  abline(h = subset(T, tool == "RepeatScout", select = y_col)[[1]], col = color)
  points(subset(T, tool=phRAIDER, select=c(x_col, y_col)), col = color)
}
  
  
plot.ROC <- function(T, org, tpr, tnr, pivot, color) {
  T <- T[T$org == org & T$tool == 'phRAIDER' & !is.na(T[tpr]) & !is.na(T[tnr]), ]
  X <- c()
  Y <- c()
  for (v in unique(T[[pivot]])) {
    T2 = T[T[pivot] == v,]
    Y <- c(Y, mean(T2[[tpr]], na.rm=TRUE))
    X <- c(X, 1 - mean(T2[[tnr]]))
  }
  plot(X, Y, xlim = c(0,1), ylim = c(0,1), col=color, main = "ROC curve", xlab = "FNR", ylab = "TPR")
}


plot.test <- function(T, org, x_axis = "w.l", x_label = "w/l") {
  par(mfrow=c(1,2))
  # Sensitity v. weight/length
  xlim = range(T[T$org==org, x_axis], na.rm=TRUE)
  ylim = range(T[T$org==org, c('tpr', 'QuCoverage')], na.rm=TRUE)
  plot.gen(T, org, x_axis, "tpr", "blue", xlab = x_label, ylab = "", xlim =xlim,
           ylim = ylim, main = "Sensitivity", new_plot = TRUE, reverse_x = FALSE)
  plot.gen(T, org, x_axis, "QuCoverage", "red", new_plot = FALSE, reverse_x = FALSE)
  #add.legend("bottomleft")
  
  # Specificity v. weight/length
  xlim = range(T[T$org==org, x_axis], na.rm=TRUE)
  ylim = range(T[T$org==org, c('tnr', 'ConCoverage')], na.rm=TRUE)
  plot.gen(T, org, x_axis, "tnr", "blue", xlab =  x_label, ylab = "", xlim =xlim,
           ylim = ylim, main = "Specificity", new_plot = TRUE, reverse_x = FALSE)
  plot.gen(T, org, x_axis, "ConCoverage", "red", new_plot = FALSE, reverse_x = FALSE) 
  
#   # Runtime v. weight/length
#   xlim = range(T[T$org==org, "w.l"], na.rm=TRUE)
#   ylim = range(T[T$org==org, c('ToolCpuTime', 'ConCoverage')], na.rm=TRUE)
#   plot.gen(T, org, "w.l", "ToolCpuTime", "blue", xlab = "w/l", ylab = "runtime", xlim = xlim,
#            ylim = ylim, main = "Runtime", new_plot = TRUE, reverse_x = TRUE)
#   
#   # ROC
#   plot.ROC(T, org, "tpr", "tnr", "w.l", "blue")
#   #plot.ROC(T, org, "QuCoverage", "ConCoverage", "w.l", "green")

  
}
  


fixed_plots <- function() {
  par(mfrow=c(3,3))
  plot.sen.by.len("ce10.chrV", stats1)
  plot.sen.by.ratio("ce10.chrV", stats2)
  plot.sen.by.ratio("ce10.chrV", stats3)
  plot.sen.by.ratio("ce10.chrV", stats4)
  plot.sen.by.ratio("ce10.chrV", stats5)
  plot.sen.by.ratio("ce10.chrV", stats6)
  plot.sen.by.ratio("ce10.chrV", stats7)
}

############################################
plot.runtimes <- function() {
  RT <- read.table("runtime.txt", header=TRUE)
  RT$size = RT$size / 1000000
  RT <- RT[RT$size > 0,]
  plot(c(), c(), xlim = range(RT$size), ylim = range(RT$time, na.rm = TRUE), main="Runtime",
       xlab = "Genome size (Mb)", ylab = "Tool runtime (seconds)")

  #LMER = RT[RT$tool=="LMER",]
  #points(LMER$size, LMER$time, col = 'green')  
  
  RS = RT[RT$tool=="RepeatScout",]
  points(RS$size, RS$time, col='red')
    
  PHRA = RT[RT$tool=="phRAIDER" & RT$seed==1,]
  points(PHRA$size, PHRA$time, col='blue')
  
  RA = RT[RT$tool=="RAIDER" & RT$seed==1,]
  points(RA$size, RA$time, col='purple')
  
  PRE = RT[RT$tool=="pre-phRAIDER" & RT$seed==1,]
  points(PRE$size, PRE$time, col='black')  
  
}


plot.mem <- function() {
  RT <- read.table("runtime.txt", header=TRUE)
  RT$size = RT$size / 1000000
  RT$mem = RT$mem / 1000000
  RT <- RT[RT$size > 0,]
  plot(c(), c(), xlim = range(RT$size), ylim = range(RT$mem, na.rm = TRUE), main="Memory",
       xlab = "Genome size (Mb)", ylab = "Tool memory (GB)")
  
  #LMER = RT[RT$tool=="LMER",]
  #points(LMER$size, LMER$mem, col = 'green')  
  
  RS = RT[RT$tool=="RepeatScout",]
  points(RS$size, RS$mem, col='red')

  PHRA = RT[RT$tool=="phRAIDER" & RT$seed==1,]
  points(PHRA$size, PHRA$mem, col='blue')
  
  RA = RT[RT$tool=="RAIDER" & RT$seed==1,]
  points(RA$size, RA$mem, col='purple')
  
  PRE = RT[RT$tool=="pre-phRAIDER" & RT$seed==1,]
  points(PRE$size, PRE$mem, col='black')  
  
}

plot.by_seed <- function(tool, column = "time") {
  RT <- read.table("runtime.txt", header=TRUE)
  RT$size = RT$size / 1000000
  RT$mem = RT$mem / 1000000
  RT <- RT[RT$size > 0,]
  
  RT <- RT[RT$tool == tool,]
  plot(c(), c(), xlim = range(RT$size), ylim = range(RT[[column]]))
  R1 <- RT[RT$seed == 0,]
  points(R1$size, R1[[column]], col = 'blue')
  R2 <- RT[RT$seed == 1,]
  points(R2$size, R2[[column]], col = 'red')
  return(list(RT, R1, R2))
}



fig1Plot <- function() {
  RT <- read.table("runtime.txt", header=TRUE)
  RT <- subset(RT, size > 20000)
  RT$size <- RT$size / 1000000
  RT$mem <- RT$mem / 1000000
  
  #setEPS()
  #postscript("RuntimeFig.eps")
  pdf(file="RuntimeFig.pdf", height=7, width=14)
  par(mfrow=c(1,2))
  
  # Make RS table.  (Need to add lmer and RS times/mem.)
  RS = subset(RT, tool == "RepeatScout")
  LM = subset(RT, tool == "LMER")
  RSLM = merge(RS, LM, by = c("label","size","seed","f"), suffixes = c(".RS", ".LM"))  
  
  PH0 <- subset(RT, tool == "phRAIDER" & seed == 0)
  PH1 <- subset(RT, tool == "phRAIDER" & seed == 1)
  PPH1 <- subset(RT, tool == 'pre-phRAIDER' & seed == 1)
  
  # First: runtimes for phRAIDER and RptScout
  plot(c(), c(), xlim = range(c(RS$size, PH1$size)), ylim = range(c(RSLM$time.RS, RSLM$timeLM, PH1$time)), xlab = "Genome size (Mb)", ylab = "Runtime (seconds)", main = "Runtime: phRAIDER v. RptScout")
  with(RSLM, points(size, time.RS + time.LM, col = 'red'))
  with(PH1, points(size, time, col = 'blue', pch = 2))
  with(PPH1, points(size, time, col = 'green', pch = 3))
  legend("topleft", c("RepeatScout", "phRAIDER", "phRAIDER.LM"), col = c('red', 'blue', 'green'), pch = c(1,2,3))
  
  # Second: memory usage for phRAIDER and RptScout
  plot(c(), c(), xlim = range(c(RS$size, PH1$size)), ylim = range(c(RSLM$mem.RS, RSLM$mem.LM, PH1$mem)), xlab = "Genome size (Mb)", ylab = "Memory required (Gb)", main = "Memory: phRAIDER v. RptScout")
  with(RSLM, points(size, mem.RS + mem.LM, col = 'red'))
  with(PH1, points(size, mem, col = 'blue', pch = 2))
  with(PPH1, points(size, mem, col = 'green', pch = 3))
  legend("topleft", c("RepeatScout", "phRAIDER", "phRAIDER.LM"), col = c('red', 'blue', 'green'), pch = c(1,2,3))
  dev.off()
  return(NA)
  
  # Third: runtime for phRAIDER on different seeds
  plot(c(), c(), xlim = range(RT$size), ylim = range(c(PH0$time, PH1$time), na.rm = TRUE), xlab = "Genome size (Mb)", ylab = "Runtime (seconds)", main = "Runtime for phRAIDER on different seeds")
  colors = c('blue', 'black')
  for (s in c(0,1)) {
    with(subset(RT, tool=='phRAIDER' & seed==s), points(size, time, col = colors[[s+1]], pch = s))
  }
  legend("topleft", c("phRAIDER seed 1", "phRAIDER seed 2"), col = c('blue', 'black'), pch = c(1,2))
  
  plot(c(), c(), xlim = range(RT$size), ylim = range(c(PH0$mem, PH1$mem), na.rm = TRUE), xlab = "Genome size Mb)", ylab = "Runtime (seconds)", main = "Memory usage for phRAIDER on different seeds")
  colors = c('blue', 'black')
  for (s in c(0,1)) {
    with(subset(RT, tool=='phRAIDER' & seed==s), points(size, mem, col = colors[[s+1]], pch = s+1))
  }
  legend("topleft", c("phRAIDER seed 1", "phRAIDER seed 2"), col = c('blue', 'black'), pch = c(1,2))
  
  RT <- read.table("runtime2.txt", header=TRUE)
  RT$size <- RT$size / 1000000
  RT$mem <- RT$mem / 1000000
  T1 <- subset(RT, w==l)
  with(T1, plot(l, time, xlim = range(l), ylim = range(time), xlab = "Seed length (w/l = 1)", ylab = "Runtime (seconds)", main = "Runtime for phRAIDER as a function of seed length (ungapped)"))
  with(T1, plot(l, mem, xlim = range(l), ylim = range(mem), xlab = "Seed length (w/l = 1)", ylab = "Memory (Gb)", main = "Memory usage for phRAIDER as a function of seed length (ungapped)"))

  RT <- read.table("runtime2.txt", header=TRUE)
  RT$size <- RT$size / 1000000
  RT$mem <- RT$mem / 1000000
  T1 <- subset(RT, w==40)
  with(T1, plot(w/l, time, xlim = range(w/l), ylim = range(time), xlab = "w/l (w = 40)", ylab = "Runtime (seconds)", main = "Runtime for phRAIDER as a function of seed width/length ratio"))
  with(T1, plot(w/l, mem, xlim = range(w/l), ylim = range(mem), xlab = "w/l (w = 40)", ylab = "Memory (Gb)", main = "Memory usage for phRAIDER as a function of seed width/length ratio"))
  
  return(T1);
}
    

