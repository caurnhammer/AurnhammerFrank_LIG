#### Join plots other scripts
library(ggplot2)
library(gridExtra)
library(plyr)

setwd("/Users/chr/R/prediction/prediction_NWE/Code/")

# Set x axis limits
#xl <- c(-7.646470, -4.830918)

################### NWE

# EEG Plot
results_lmem = read.csv("../Output_data/eeg_snapshots_LMEM.csv", sep = " ")
dt <- NULL
#dt$average_surprisal <- rep(c(7.525636, 7.177759, 7.035030, 6.647705, 6.233933, 5.762444, 5.274674, 4.948322, 4.850785), 2) * -1
dt$model <- revalue(results_lmem$model, c("1k"="1K", "3k"="3K", "10k"="10K", "30k"="30K",
                                      "100k"="100K", "300k"="300K",
                                      "1M"="1M", "3M"="3M", "6-46M"="6.47M"))[seq(1, 35, 4)]
dt$metric_value      <- c(as.numeric(results_lmem[seq(1, length(results_lmem[,1]), 4),4]), 
                          as.numeric(results_lmem[seq(2, length(results_lmem[,1]), 4),4])) * -1
dt$metric_name       <- c(rep("Surprisal", 9), rep("Next-word Entropy", 9))
dt <- as.data.frame(dt)
dt <- dt[-c(1, 2, 10, 11),]
dt$model <- ordered(dt$model, levels=c("10K", "30K", "100K", "300K", "1M", "3M", "6.47M"))
dt = dt[8:14,]

eeg_plt <- ggplot(dt, aes(x=model, y=metric_value, color=metric_name, shape = metric_name)) + geom_point(size=3)
eeg_plt <- eeg_plt + labs(title = "N400", x = "Training Sentences", y = "Goodness-of-fit")
eeg_plt <- eeg_plt + theme(legend.position = "none") 
#eeg_plt <- eeg_plt + geom_hline(yintercept = 5.991, size = 0.35,  linetype = "dashed") + geom_hline(yintercept = -5.991, size = 0.35, linetype = "dashed")
eeg_plt <- eeg_plt + geom_hline(yintercept = 0, color = "black")
eeg_plt <- eeg_plt + theme(panel.background = element_rect(fill = "#FFFFFF", color="#000000", size = , linetype = "solid"), 
                           panel.grid.major = element_line(size = 0.05, linetype = "solid", color = "#d3d3d3"), 
                           panel.grid.minor = element_line(size = 0.025, linetype = "solid", color = "#d3d3d3"))
eeg_plt <- eeg_plt + scale_colour_manual(values=c("#E69F00", "#0072B2"))
eeg_plt <- eeg_plt + scale_color_grey(start=0.2, end=0.6)

# ET Plto
results_lmem = read.csv("../Output_data/PREV_et_snapshots_LMEM.csv", sep = " ")
dt <- NULL
#dt$average_surprisal <- rep(c(8.145235, 7.816667, 7.646470, 7.281974, 6.880306, 6.407189, 5.878788, 5.531965, 5.423403), 2) * -1
dt$model <- revalue(results_lmem$model, c("1k"="1K", "3k"="3K", "10k"="10K", "30k"="30K",
                                          "100k"="100K", "300k"="300K",
                                          "1M"="1M", "3M"="3M", "6-46M"="6.47M"))[seq(1, 35, 4)]
dt$metric_value      <- c(as.numeric(results_lmem[seq(1, length(results_lmem[,1]), 4),4]),
                          as.numeric(results_lmem[seq(2, length(results_lmem[,1]), 4),4]))
dt$metric_name       <- c(rep("Surprisal", 9), rep("Next-word Entropy", 9))
dt <- as.data.frame(dt)
dt <- dt[-c(1, 2, 10, 11),]
dt$model <- ordered(dt$model, levels=c("10K", "30K", "100K", "300K", "1M", "3M", "6.47M"))
dt = dt[8:14,]

et_plt <- ggplot(dt, aes(x=model, y=metric_value, color=metric_name, shape = metric_name)) + geom_point(size=3)
et_plt <- et_plt + labs(title = "First-pass RT", x = "Training Sentences", y = "", color = "", shape = "")
et_plt <- et_plt + theme(legend.position = "bottom")
#et_plt <- et_plt + geom_hline(yintercept = 5.991, size = 0.35, color = "#F8766D",  linetype = "dashed") + geom_hline(yintercept = -5.991, size = 0.35, color = "#F8766D",  linetype = "dashed")
#et_plt <- et_plt + geom_hline(yintercept = 9.488, size = 0.35, color = "#00BFC4",  linetype = "dashed") + geom_hline(yintercept = -9.488, size = 0.35, color = "#00BFC4",  linetype = "dashed")
et_plt <- et_plt + geom_hline(yintercept = 0, color = "black")
et_plt <- et_plt + theme(panel.background = element_rect(fill = "#FFFFFF", color="#000000", size = , linetype = "solid"), 
                         panel.grid.major = element_line(size = 0.05, linetype = "solid", color = "#d3d3d3"), 
                         panel.grid.minor = element_line(size = 0.025, linetype = "solid", color = "#d3d3d3"))
et_plt <- et_plt + scale_colour_manual(values=c("#E69F00", "#0072B2"))
et_plt <- et_plt + scale_color_grey(start=0.2, end=0.6)

# SPR Plot
# Plot LMEM results
results_lmem = read.csv("../Output_data/PREV_spr_snapshots_LMEM.csv", sep=" ")
dt <- NULL
#dt$average_surprisal <- rep(c(7.465115, 7.110349, 6.952130, 6.556464, 6.160724, 5.717861, 5.265369, 4.936375, 4.830918), 2) * -1
dt$model <- revalue(results_lmem$model, c("1k"="1K", "3k"="3K", "10k"="10K", "30k"="30K",
                                          "100k"="100K", "300k"="300K",
                                          "1M"="1M", "3M"="3M", "6-46M"="6.47M"))[seq(1, 35, 4)]
dt$metric_value      <- c(as.numeric(results_lmem[seq(1, length(results_lmem[,1]), 4),4]), 
                          as.numeric(results_lmem[seq(2, length(results_lmem[,1]), 4),4]))
dt$metric_name       <- c(rep("Surprisal", 9), rep("Next-word Entropy", 9))
dt <- as.data.frame(dt)
dt <- dt[-c(1, 2, 10, 11),]
dt$model <- ordered(dt$model, levels=c("10K", "30K", "100K", "300K", "1M", "3M", "6.47M"))
dt = dt[8:14,]

spr_plt <- ggplot(dt, aes(x=model, y=metric_value, color=metric_name, shape = metric_name)) + geom_point(size=3)
spr_plt <- spr_plt + labs(title = "Self-paced RT", x = "Training Sentences", y = "")
spr_plt <- spr_plt + theme(legend.position = "none")
#spr_plt <- spr_plt + geom_hline(yintercept = 5.991, size = 0.35, color = "#F8766D",  linetype = "dashed") + geom_hline(yintercept = -5.991, size = 0.35, color = "#F8766D",  linetype = "dashed")
#spr_plt <- spr_plt + geom_hline(yintercept = 9.488, size = 0.35, color = "#00BFC4",  linetype = "dashed") + geom_hline(yintercept = -9.488, size = 0.35, color = "#00BFC4",  linetype = "dashed")
spr_plt <- spr_plt + geom_hline(yintercept = 0, color = "black")
spr_plt <- spr_plt + theme(panel.background = element_rect(fill = "#FFFFFF", color="#000000", size = , linetype = "solid"), 
                           panel.grid.major = element_line(size = 0.05, linetype = "solid", color = "#d3d3d3"), 
                           panel.grid.minor = element_line(size = 0.025, linetype = "solid", color = "#d3d3d3"))
spr_plt <- spr_plt + scale_colour_manual(values=c("#E69F00", "#0072B2"))
spr_plt <- spr_plt + scale_color_grey(start=0.2, end=0.6)

# Get legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

lgnd <- g_legend(et_plt)

# Arrange plots
grid.arrange(arrangeGrob(eeg_plt, et_plt + theme(legend.position = "none"), spr_plt, nrow=1), 
             lgnd, heights=c(10, 1))







########## SURPRISAL
# EEG Plot
results_lmem = read.csv("../Output_data/eeg_snapshots_LMEM.csv", sep = " ")
dt <- NULL
#dt$average_surprisal <- rep(c(7.525636, 7.177759, 7.035030, 6.647705, 6.233933, 5.762444, 5.274674, 4.948322, 4.850785), 2) * -1
dt$model <- revalue(results_lmem$model, c("1k"="1K", "3k"="3K", "10k"="10K", "30k"="30K",
                                          "100k"="100K", "300k"="300K",
                                          "1M"="1M", "3M"="3M", "6-46M"="6.47M"))[seq(1, 35, 4)]
dt$metric_value      <- c(as.numeric(results_lmem[seq(1, length(results_lmem[,1]), 4),4]), 
                          as.numeric(results_lmem[seq(2, length(results_lmem[,1]), 4),4])) * -1
dt$metric_name       <- c(rep("Surprisal", 9), rep("Next-word Entropy", 9))
dt <- as.data.frame(dt)
dt <- dt[-c(1, 2, 10, 11),]
dt$model <- ordered(dt$model, levels=c("10K", "30K", "100K", "300K", "1M", "3M", "6.47M"))
dt = dt[1:7,]

eeg_plt <- ggplot(dt, aes(x=model, y=metric_value, color=metric_name, shape = metric_name)) + geom_point(size=3)
eeg_plt <- eeg_plt + labs(title = "N400", x = "Training Sentences", y = "Goodness-of-fit")
eeg_plt <- eeg_plt + theme(legend.position = "none") 
#eeg_plt <- eeg_plt + geom_hline(yintercept = 5.991, size = 0.35,  linetype = "dashed") + geom_hline(yintercept = -5.991, size = 0.35, linetype = "dashed")
eeg_plt <- eeg_plt + geom_hline(yintercept = 0, color = "black")
eeg_plt <- eeg_plt + theme(panel.background = element_rect(fill = "#FFFFFF", color="#000000", size = , linetype = "solid"), 
                           panel.grid.major = element_line(size = 0.05, linetype = "solid", color = "#d3d3d3"), 
                           panel.grid.minor = element_line(size = 0.025, linetype = "solid", color = "#d3d3d3"))
eeg_plt <- eeg_plt + scale_colour_manual(values=c("#E69F00", "#0072B2"))
eeg_plt <- eeg_plt + scale_color_grey(start=0.2, end=0.6)

# ET Plto
results_lmem = read.csv("../Output_data/PREV_et_snapshots_LMEM.csv", sep = " ")
dt <- NULL
#dt$average_surprisal <- rep(c(8.145235, 7.816667, 7.646470, 7.281974, 6.880306, 6.407189, 5.878788, 5.531965, 5.423403), 2) * -1
dt$model <- revalue(results_lmem$model, c("1k"="1K", "3k"="3K", "10k"="10K", "30k"="30K",
                                          "100k"="100K", "300k"="300K",
                                          "1M"="1M", "3M"="3M", "6-46M"="6.47M"))[seq(1, 35, 4)]
dt$metric_value      <- c(as.numeric(results_lmem[seq(1, length(results_lmem[,1]), 4),4]),
                          as.numeric(results_lmem[seq(2, length(results_lmem[,1]), 4),4]))
dt$metric_name       <- c(rep("Surprisal", 9), rep("Next-word Entropy", 9))
dt <- as.data.frame(dt)
dt <- dt[-c(1, 2, 10, 11),]
dt$model <- ordered(dt$model, levels=c("10K", "30K", "100K", "300K", "1M", "3M", "6.47M"))
dt = dt[1:7,]

et_plt <- ggplot(dt, aes(x=model, y=metric_value, color=metric_name, shape = metric_name)) + geom_point(size=3)
et_plt <- et_plt + labs(title = "First-pass RT", x = "Training Sentences", y = "", color = "", shape = "")
et_plt <- et_plt + theme(legend.position = "bottom")
#et_plt <- et_plt + geom_hline(yintercept = 5.991, size = 0.35, color = "#F8766D",  linetype = "dashed") + geom_hline(yintercept = -5.991, size = 0.35, color = "#F8766D",  linetype = "dashed")
#et_plt <- et_plt + geom_hline(yintercept = 9.488, size = 0.35, color = "#00BFC4",  linetype = "dashed") + geom_hline(yintercept = -9.488, size = 0.35, color = "#00BFC4",  linetype = "dashed")
et_plt <- et_plt + geom_hline(yintercept = 0, color = "black")
et_plt <- et_plt + theme(panel.background = element_rect(fill = "#FFFFFF", color="#000000", size = , linetype = "solid"), 
                         panel.grid.major = element_line(size = 0.05, linetype = "solid", color = "#d3d3d3"), 
                         panel.grid.minor = element_line(size = 0.025, linetype = "solid", color = "#d3d3d3"))
et_plt <- et_plt + scale_colour_manual(values=c("#E69F00", "#0072B2"))
et_plt <- et_plt + scale_color_grey(start=0.2, end=0.6)

# SPR Plot
# Plot LMEM results
results_lmem = read.csv("../Output_data/PREV_spr_snapshots_LMEM.csv", sep=" ")
dt <- NULL
#dt$average_surprisal <- rep(c(7.465115, 7.110349, 6.952130, 6.556464, 6.160724, 5.717861, 5.265369, 4.936375, 4.830918), 2) * -1
dt$model <- revalue(results_lmem$model, c("1k"="1K", "3k"="3K", "10k"="10K", "30k"="30K",
                                          "100k"="100K", "300k"="300K",
                                          "1M"="1M", "3M"="3M", "6-46M"="6.47M"))[seq(1, 35, 4)]
dt$metric_value      <- c(as.numeric(results_lmem[seq(1, length(results_lmem[,1]), 4),4]), 
                          as.numeric(results_lmem[seq(2, length(results_lmem[,1]), 4),4]))
dt$metric_name       <- c(rep("Surprisal", 9), rep("Next-word Entropy", 9))
dt <- as.data.frame(dt)
dt <- dt[-c(1, 2, 10, 11),]
dt$model <- ordered(dt$model, levels=c("10K", "30K", "100K", "300K", "1M", "3M", "6.47M"))
dt = dt[1:7,]

spr_plt <- ggplot(dt, aes(x=model, y=metric_value, color=metric_name, shape = metric_name)) + geom_point(size=3)
spr_plt <- spr_plt + labs(title = "Self-paced RT", x = "Training Sentences", y = "")
spr_plt <- spr_plt + theme(legend.position = "none")
#spr_plt <- spr_plt + geom_hline(yintercept = 5.991, size = 0.35, color = "#F8766D",  linetype = "dashed") + geom_hline(yintercept = -5.991, size = 0.35, color = "#F8766D",  linetype = "dashed")
#spr_plt <- spr_plt + geom_hline(yintercept = 9.488, size = 0.35, color = "#00BFC4",  linetype = "dashed") + geom_hline(yintercept = -9.488, size = 0.35, color = "#00BFC4",  linetype = "dashed")
spr_plt <- spr_plt + geom_hline(yintercept = 0, color = "black")
spr_plt <- spr_plt + theme(panel.background = element_rect(fill = "#FFFFFF", color="#000000", size = , linetype = "solid"), 
                           panel.grid.major = element_line(size = 0.05, linetype = "solid", color = "#d3d3d3"), 
                           panel.grid.minor = element_line(size = 0.025, linetype = "solid", color = "#d3d3d3"))
spr_plt <- spr_plt + scale_colour_manual(values=c("#E69F00", "#0072B2"))
spr_plt <- spr_plt + scale_color_grey(start=0.2, end=0.6)

# Get legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

lgnd <- g_legend(et_plt)

# Arrange plots
grid.arrange(arrangeGrob(eeg_plt, et_plt + theme(legend.position = "none"), spr_plt, nrow=1), 
             lgnd, heights=c(10, 1))