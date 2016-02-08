orgList <- c("ce10.chrV", "hg38.chr22", "ce10.chrV", "hg38.chr18.sim")

read.file <- function(file) {
  T <- read.table(file, header = TRUE)
  return(T)
}

read.stats <- function() {
  seeds32 <<- read.table("seeds32.stats.txt", header = TRUE)
  best32 <<- read.table("best.stats.txt", header = TRUE)
  RptSct <<- read.table("RptSct.stats.txt", header = TRUE)
  
  TList <<- lapply(orgList, function(o) subset(seeds32, org == o))
  names(TList) <<- orgList;
}

add.legend <- function(position) {
  legend(position, c("RM", "PRA"), col = c('blue', 'green'), 
                     lty = c('solid', 'solid', 'solid', 'solid'))
}

plot.gen <- function(org_name, x_col, y_col, xlab = NA, ylab = NA, main = NA, 
                     xlim = NA, ylim = NA, legend_position = "topleft", new_plot = TRUE,
                     reverse_x = FALSE, reverse_y = FALSE) {

  PH = subset(seeds32, org == org_name, select = c(x_col,y_col))
  RS = subset(RptSct, org == org_name, select = c(x_col,y_col))
  
  if (is.na(xlim)) {
    xlim = range(PH[[x_col]], rm.na = TRUE)
  }
  if (is.na(ylim)) {
    ylim = range(c(PH[[y_col]], RS[[y_col]]));
  }
  
  if (reverse_x) {
    xlim = rev(xlim)
  }
  if (reverse_y) {
    ylim = rev(ylim)
  }
  
  if (new_plot) {
    plot(c(), c(), xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main)
  }
  abline(h = RS[1,y_col], col = 'red')
  points(PH[[x_col]], PH[[y_col]], col = 'blue')
}
  
  


plot.test <- function(PH, RS, org_name, x_axis = "w.l", x_label = "w/l") {
  par(mfcol=c(1,2))
  
  PH = subset(PH, org == org_name & tool == "phRAIDER")
  RS = subset(RS, org == org_name)


  # First the RM-metric
  xlim = range(PH[[x_axis]])
  print(xlim)
  ylim = range(c(PH$tpr, PH$QuCoverage, RS$tpr, RS$QuCoverage))
  print(ylim)
  plot(c(), c(), xlim = xlim, ylim = ylim, xlab = x_label, ylab = "", main = "Sensitivity")
  points(PH[[x_axis]], PH[["tpr"]], col = 'blue')
  abline(h = RS[[1,"tpr"]], col = 'blue')
  points(PH[[x_axis]], PH[["QuCoverage"]], col = 'red')
  abline(h = RS[[1,"QuCoverage"]], col = 'red')
  
  xlim = range(PH[[x_axis]])
  ylim = range(c(PH$tnr, PH$ConCoverage, RS$tnr, RS$ConCoverage))
  plot(c(), c(), xlim = xlim, ylim = ylim, xlab = x_label, ylab = "", main = "Sensitivity")
  points(PH[[x_axis]], PH[["tnr"]], col = 'blue')
  abline(h = RS[[1,"tnr"]], col = 'blue')
  points(PH[[x_axis]], PH[["ConCoverage"]], col = 'red')
  abline(h = RS[[1,"ConCoverage"]], col = 'red') 
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
  #quartz()
  par(mfrow=c(4,2))
  
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
  legend("topleft", c("RepeatScout", "phRAIDER", "pre-phRAIDER"), col = c('red', 'blue', 'green'), pch = c(1,2,3))
  
  # Second: memory usage for phRAIDER and RptScout
  plot(c(), c(), xlim = range(c(RS$size, PH1$size)), ylim = range(c(RSLM$mem.RS, RSLM$mem.LM, PH1$mem)), xlab = "Genome size (Mb)", ylab = "Memory required (Gb)", main = "Memory: phRAIDER v. RptScout")
  with(RSLM, points(size, mem.RS + mem.LM, col = 'red'))
  with(PH1, points(size, mem, col = 'blue', pch = 2))
  with(PPH1, points(size, mem, col = 'green', pch = 3))
  legend("topleft", c("RepeatScout", "phRAIDER", "pre-phRAIDER"), col = c('red', 'blue', 'green'), pch = c(1,2,3))
  
  
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
  T1 <- subset(RT, w==l)
  with(T1, plot(l, time, xlim = range(l), ylim = range(time), xlab = "Seed length (w/l = 1)", ylab = "Runtime (seconds)", main = "Runtime for phRAIDER as a function of seed length (ungapped)"))
  with(T1, plot(l, mem, xlim = range(l), ylim = range(mem), xlab = "Seed length (w/l = 1)", ylab = "Memory (Gb)", main = "Memory usage for phRAIDER as a function of seed length (ungapped)"))

  RT <- read.table("runtime2.txt", header=TRUE)
  T1 <- subset(RT, w==40)
  with(T1, plot(w/l, time, xlim = range(w/l), ylim = range(time), xlab = "w/l (w = 40)", ylab = "Runtime (seconds)", main = "Runtime for phRAIDER as a function of seed width/length ratio"))
  with(T1, plot(w/l, mem, xlim = range(w/l), ylim = range(mem), xlab = "w/l (w = 40)", ylab = "Memory (Gb)", main = "Memory usage for phRAIDER as a function of seed width/length ratio"))
  
  return(T1);
}
    

  
  