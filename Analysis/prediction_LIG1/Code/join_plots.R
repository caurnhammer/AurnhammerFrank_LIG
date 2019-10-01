#### Join plots other scripts
library(ggplot2)
library(gridExtra)
library(plyr)

setwd("/Users/chr/R/prediction/prediction_UNI/Code/")

# Set metric labels
metrics = c(expression("Lookahead information gain"[1]), "Surprisal")

# EEG Plot
results_lmem = read.csv("../Output_data/eeg_snapshots_LMEM.csv", sep = " ")
dt <- NULL
#dt$average_surprisal <- rep(c(7.4309, 7.1081, 6.9922, 6.6013, 6.1707, 5.6739, 5.1746, 4.8661, 4.7802), 2) * -1
dt$model <- revalue(results_lmem$model, c("1k"="1K", "3k"="3K", "10k"="10K", "30k"="30K",
                                          "100k"="100K", "300k"="300K",
                                          "1M"="1M", "3M"="3M", "6-46M"="6.47M"))[seq(1, 35, 4)]
dt$metric_value      <- c(as.numeric(results_lmem[seq(3, length(results_lmem[,1]), 4),4]), 
                          as.numeric(results_lmem[seq(4, length(results_lmem[,1]), 4),4])) * -1    
dt$metric_name       <- c(rep("Surprisal", 9), rep("Lookahead information gain₁", 9))
dt <- as.data.frame(dt)
dt <- dt[-c(1, 2, 10, 11),]
dt$model <- ordered(dt$model, levels=c("10K", "30K", "100K", "300K", "1M", "3M", "6.47M"))
        
eeg_plt <- ggplot(dt, aes(x=model, y=metric_value, color=metric_name, shape = metric_name)) + geom_point(size=3)
eeg_plt <- eeg_plt + labs(title = "N400", x = "Training Sentences", y = "Goodness-of-fit")
eeg_plt <- eeg_plt + theme(legend.position = "none") 
#eeg_plt <- eeg_plt + geom_hline(yintercept = 5.991, size = 0.35,  linetype = "dashed") + geom_hline(yintercept = -5.991, size = 0.35, linetype = "dashed")
eeg_plt <- eeg_plt + geom_hline(yintercept = 0, color = "black")
eeg_plt <- eeg_plt + theme(panel.background = element_rect(fill = "#FFFFFF", color="#000000", size = , linetype = "solid"), 
                           panel.grid.major = element_line(size = 0.05, linetype = "solid", color = "#d3d3d3"), 
                           panel.grid.minor = element_line(size = 0.025, linetype = "solid", color = "#d3d3d3"))
eeg_plt <- eeg_plt + scale_colour_manual(values=c("#E69F00", "#0072B2"))
#eeg_plt <- eeg_plt + scale_color_grey(start=0.2, end=0.6)

# ET Plto
results_lmem = read.csv("../Output_data/PREV_et_snapshots_LMEM.csv", sep = " ")
dt <- NULL
#dt$average_surprisal <- rep(c(8.0437, 7.7381, 7.6013, 7.2359, 6.8171, 6.3182, 5.7768, 5.4480, 5.3508), 2) * -1
dt$model <- revalue(results_lmem$model, c("1k"="1K", "3k"="3K", "10k"="10K", "30k"="30K",
                                          "100k"="100K", "300k"="300K",
                                          "1M"="1M", "3M"="3M", "6-46M"="6.47M"))[seq(1, 35, 4)]
dt$metric_value      <- c(as.numeric(results_lmem[seq(3, length(results_lmem[,1]), 4),4]),
                          as.numeric(results_lmem[seq(4, length(results_lmem[,1]), 4),4]))
dt$metric_name       <- c(rep("Surprisal", 9), rep("Lookahead information gain₁", 9))
dt <- as.data.frame(dt)
dt <- dt[-c(1, 2, 10, 11),]
dt$model <- ordered(dt$model, levels=c("10K", "30K", "100K", "300K", "1M", "3M", "6.47M"))

et_plt <- ggplot(dt, aes(x=model, y=metric_value, color=metric_name, shape = metric_name)) + geom_point(size=3)
et_plt <- et_plt + labs(title = "First-pass RT", x = "Training Sentences", y = "", color = "", shape = "")
et_plt <- et_plt + theme(legend.position = "bottom") 
#et_plt <- et_plt + geom_hline(yintercept = 5.991, size = 0.35, color = "#F8766D",  linetype = "dashed") + geom_hline(yintercept = -5.991, size = 0.35, color = "#F8766D",  linetype = "dashed")
#et_plt <- et_plt + geom_hline(yintercept = 9.488, size = 0.35, color = "#00BFC4",  linetype = "dashed") + geom_hline(yintercept = -9.488, size = 0.35, color = "#00BFC4",  linetype = "dashed")
et_plt <- et_plt + geom_hline(yintercept = 0, color = "black")
et_plt <- et_plt + theme(panel.background = element_rect(fill = "#FFFFFF", color="#000000", size = , linetype = "solid"), 
                           panel.grid.major = element_line(size = 0.05, linetype = "solid", color = "#d3d3d3"), 
                           panel.grid.minor = element_line(size = 0.025, linetype = "solid", color = "#d3d3d3"))
et_plt <- et_plt + scale_colour_manual(values = c("#E69F00", "#0072B2"), labels = metrics)
et_plt <- et_plt + scale_shape_manual(values = c(16, 17), labels = metrics)
#et_plt <- et_plt + scale_color_grey(start=0.2, end=0.6)

# SPR Plot
# Plot LMEM results
results_lmem = read.csv("../Output_data/PREV_spr_snapshots_LMEM.csv", sep=" ")
dt <- NULL
#dt$average_surprisal <- rep(c(7.3679, 7.0417, 6.9108, 6.5104, 6.0985, 5.6292, 5.1624, 4.8508, 4.7545), 2) * -1
dt$model <- revalue(results_lmem$model, c("1k"="1K", "3k"="3K", "10k"="10K", "30k"="30K",
                                          "100k"="100K", "300k"="300K",
                                          "1M"="1M", "3M"="3M", "6-46M"="6.47M"))[seq(1, 35, 4)]
dt$metric_value      <- c(as.numeric(results_lmem[seq(3, length(results_lmem[,1]), 4),4]), 
                          as.numeric(results_lmem[seq(4, length(results_lmem[,1]), 4),4]))
dt$metric_name       <- c(rep("Surprisal", 9), rep("Lookahead information gain₁", 9))
dt <- as.data.frame(dt)
dt <- dt[-c(1, 2, 10, 11),]
dt$model <- ordered(dt$model, levels=c("10K", "30K", "100K", "300K", "1M", "3M", "6.47M"))

spr_plt <- ggplot(dt, aes(x=model, y=metric_value, color=metric_name, shape = metric_name)) + geom_point(size=3)
spr_plt <- spr_plt + labs(title = "Self-paced RT", x = "Training Sentences", y = "")
spr_plt <- spr_plt + theme(legend.position = "none")
#spr_plt <- spr_plt + geom_hline(yintercept = 5.991, size = 0.35, color = "#F8766D",  linetype = "dashed") + geom_hline(yintercept = -5.991, size = 0.35, color = "#F8766D",  linetype = "dashed")
#spr_plt <- spr_plt + geom_hline(yintercept = 9.488, size = 0.35, color = "#00BFC4",  linetype = "dashed") + geom_hline(yintercept = -9.488, size = 0.35, color = "#00BFC4",  linetype = "dashed")
spr_plt <- spr_plt + geom_hline(yintercept = 0, color = "black")
spr_plt <- spr_plt + theme(panel.background = element_rect(fill = "#FFFFFF", color="#000000", size = , linetype = "solid"), 
                         panel.grid.major = element_line(size = 0.05, linetype = "solid", color = "#d3d3d3"), 
                         panel.grid.minor = element_line(size = 0.025, linetype = "solid", color = "#d3d3d3"))
spr_plt <- spr_plt + scale_colour_manual(values=c( "#E69F00", "#0072B2"))
#spr_plt <- spr_plt + scale_color_grey(start=0.2, end=0.6)


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
