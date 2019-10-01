# Christoph Aurnhammer, October 29 2018
# Analyse snapshots on SPR data (Frank et al., 2013) in a nonlinear and linear model
rm(list=ls())

data = read.csv('../Input_data/PREV_spr_entropy_avg_snapshots.csv', sep = '\t')
data$subj_nr  = as.factor(data$subj_nr)
data$sent_nr  = as.factor(data$sent_nr)
data$item     = as.factor(data$item)
data$word     = as.character(data$word)
data$reject_data = as.factor(data$reject_data)
data          = na.omit(data)

# Prepare output file
results = data.frame(model=character(), baseline=character(), interest=character(), chi2=numeric(), df=numeric(), pval=numeric(),
                     stringsAsFactors = FALSE)

library(lme4)

# Prepare output file
results_lmem = data.frame(model=character(), baseline=character(), interest=character(), chi2=numeric(), df=numeric(), pval=numeric(),
                          stringsAsFactors = FALSE)

# Loop: for snapshot in snapshots
mod_num = seq(1, 9)
snapshots = seq(0, 9)
prevs = seq(32, 49)
models = c("1k", "3k", "10k", "30k", "100k", "300k", "1M", "3M", "6-46M")

# build baseline, surprisal, entropy, complete 
for (mod in mod_num){
  print(models[mod])
  print(snapshots[mod])
  print(colnames(data)[14+snapshots[mod]])
  print(colnames(data)[prevs[mod]])
  print(colnames(data)[14+snapshots[mod]+9])
  print(colnames(data)[prevs[mod]+9])
  
  # Build models
  m_baseline = lmer(RT ~ log_freq
                    +log_freq_prev
                    +nr_char
                    +nr_char_prev
                    +word_pos
                    +RT_prev
                    +log_freq:log_freq_prev
                    +log_freq:nr_char
                    +log_freq:nr_char_prev
                    +log_freq:word_pos
                    +log_freq:RT_prev
                    +log_freq_prev:nr_char
                    +log_freq_prev:nr_char_prev
                    +log_freq_prev:word_pos
                    +log_freq_prev:RT_prev
                    +nr_char:nr_char_prev
                    +nr_char:word_pos
                    +nr_char:RT_prev
                    +nr_char_prev:word_pos
                    +nr_char_prev:RT_prev
                    +word_pos:RT_prev
                    +(1|subj_nr)
                    +(1|item)
                    +(0+log_freq|subj_nr)
                    +(0+log_freq_prev|subj_nr)
                    +(0+nr_char|subj_nr)
                    +(0+nr_char_prev|subj_nr)
                    +(0+word_pos|subj_nr)
                    +(0+RT_prev|word),
                    data = data, REML = FALSE, verbose = 1,
                    control=lmerControl(optCtrl=list(75000)))
  
  m_surprisal = lmer(RT ~ data[,14+snapshots[mod]]
                     +data[,prevs[mod]]
                     +log_freq
                     +log_freq_prev
                     +nr_char
                     +nr_char_prev
                     +word_pos
                     +RT_prev
                     +log_freq:log_freq_prev
                     +log_freq:nr_char
                     +log_freq:nr_char_prev
                     +log_freq:word_pos
                     +log_freq:RT_prev
                     +log_freq_prev:nr_char
                     +log_freq_prev:nr_char_prev
                     +log_freq_prev:word_pos
                     +log_freq_prev:RT_prev
                     +nr_char:nr_char_prev
                     +nr_char:word_pos
                     +nr_char:RT_prev
                     +nr_char_prev:word_pos
                     +nr_char_prev:RT_prev
                     +word_pos:RT_prev
                     +(1|subj_nr)
                     +(1|item)
                     +(0+log_freq|subj_nr)
                     +(0+log_freq_prev|subj_nr)
                     +(0+nr_char|subj_nr)
                     +(0+nr_char_prev|subj_nr)
                     +(0+word_pos|subj_nr)
                     +(0+RT_prev|word)
                     +(0+data[,14+snapshots[mod]]|subj_nr)
                     +(0+data[,prevs[mod]]|subj_nr),
                     data = data, REML = FALSE, verbose = 1,
                     control=lmerControl(optCtrl=list(75000)))
  
  m_entropy = lmer(RT ~ data[,14+snapshots[mod]+9]
                   +data[,prevs[mod]+9]
                   +log_freq
                   +log_freq_prev
                   +nr_char
                   +nr_char_prev
                   +word_pos
                   +RT_prev
                   +log_freq:log_freq_prev
                   +log_freq:nr_char
                   +log_freq:nr_char_prev
                   +log_freq:word_pos
                   +log_freq:RT_prev
                   +log_freq_prev:nr_char
                   +log_freq_prev:nr_char_prev
                   +log_freq_prev:word_pos
                   +log_freq_prev:RT_prev
                   +nr_char:nr_char_prev
                   +nr_char:word_pos
                   +nr_char:RT_prev
                   +nr_char_prev:word_pos
                   +nr_char_prev:RT_prev
                   +word_pos:RT_prev
                   +(1|subj_nr)
                   +(1|item)
                   +(0+log_freq|subj_nr)
                   +(0+log_freq_prev|subj_nr)
                   +(0+nr_char|subj_nr)
                   +(0+nr_char_prev|subj_nr)
                   +(0+word_pos|subj_nr)
                   +(0+RT_prev|word)
                   +(0+data[,14+snapshots[mod]+9]|subj_nr)
                   +(0+data[,prevs[mod]+9]|subj_nr),
                   data = data, REML = FALSE, verbose = 1,
                   control=lmerControl(optCtrl=list(75000)))
  
  m_complete = lmer(RT ~ data[,14+snapshots[mod]]
                    +data[,prevs[mod]]
                    +data[,14+snapshots[mod]+9]
                    +data[,prevs[mod]+9]
                    +log_freq
                    +log_freq_prev
                    +nr_char
                    +nr_char_prev
                    +word_pos
                    +RT_prev
                    +log_freq:log_freq_prev
                    +log_freq:nr_char
                    +log_freq:nr_char_prev
                    +log_freq:word_pos
                    +log_freq:RT_prev
                    +log_freq_prev:nr_char
                    +log_freq_prev:nr_char_prev
                    +log_freq_prev:word_pos
                    +log_freq_prev:RT_prev
                    +nr_char:nr_char_prev
                    +nr_char:word_pos
                    +nr_char:RT_prev
                    +nr_char_prev:word_pos
                    +nr_char_prev:RT_prev
                    +word_pos:RT_prev
                    +(1|subj_nr)
                    +(1|item)
                    +(0+log_freq|subj_nr)
                    +(0+log_freq_prev|subj_nr)
                    +(0+nr_char|subj_nr)
                    +(0+nr_char_prev|subj_nr)
                    +(0+word_pos|subj_nr)
                    +(0+RT_prev|word)
                    +(0+data[,14+snapshots[mod]]|subj_nr)
                    +(0+data[,prevs[mod]]|subj_nr)
                    +(0+data[,14+snapshots[mod]+9]|subj_nr)
                    +(0+data[,prevs[mod]+9]|subj_nr),
                    data = data, REML = FALSE, verbose = 1,
                    control=lmerControl(optCtrl=list(75000)))
  
  # compare models & append to output file & write to disc
  # surprisal over baseline
  comp = anova(m_baseline, m_surprisal)
  chisq = round(comp$Chisq[2],4)
  if (coef(summary(m_surprisal))[2] < 0){
    chisq = chisq * -1
  }
  results_lmem[nrow(results_lmem) + 1,] = c(models[mod],
                                            "m_baseline",
                                            "m_surprisal",
                                            chisq,
                                            comp$`Chi Df`[2],
                                            round(comp$`Pr(>Chisq)`[2],4))
  
  # entropy over baseline
  comp = anova(m_baseline, m_entropy)
  chisq = round(comp$Chisq[2],4)
  if (coef(summary(m_entropy))[2] < 0){
    chisq = chisq * -1
  }
  results_lmem[nrow(results_lmem) + 1,] = c(models[mod],
                                            "m_baseline",
                                            "m_entropy",
                                            chisq,
                                            comp$`Chi Df`[2],
                                            round(comp$`Pr(>Chisq)`[2],4))
  
  # surprisal over entropy
  comp = anova(m_entropy, m_complete)
  chisq = round(comp$Chisq[2],4)
  if (coef(summary(m_complete))[2] < 0){
    chisq = chisq * -1
  }
  results_lmem[nrow(results_lmem) + 1,] = c(models[mod],
                                            "m_entropy",
                                            "m_complete",
                                            chisq,
                                            comp$`Chi Df`[2],
                                            round(comp$`Pr(>Chisq)`[2],4))
  
  # entropy over surprisal
  comp = anova(m_surprisal, m_complete)
  chisq = round(comp$Chisq[2],4)
  if (coef(summary(m_complete))[4] < 0){
    chisq = chisq * -1
  }
  results_lmem[nrow(results_lmem) + 1,] = c(models[mod],
                                            "m_surprisal",
                                            "m_complete",
                                            chisq,
                                            comp$`Chi Df`[2],
                                            round(comp$`Pr(>Chisq)`[2],4))
  
  # Write results to disc
  write.table(results_lmem, file="../Output_data/PREV_spr_snapshots_LMEM.csv", row.names = FALSE)
  
  print("-------------------------------------")
  # End loop
}

# Plot LMEM results
results_lmem = read.csv("../Output_data/PREV_spr_snapshots_LMEM.csv", sep=" ")
dev.off()
plot.new()
avg_surprisal = c(7.465115, 7.110349, 6.952130, 6.556464, 6.160724, 5.717861, 5.265369, 4.936375, 4.830918)
avg_entropy   = c(6.917207, 6.484064, 6.656342, 6.408130, 5.966319, 5.352943, 4.990034, 4.710546, 4.558530)
chi_surprisal_lmem = as.numeric(results_lmem[seq(3, length(results_lmem[,1]), 4),4])
chi_entropy_lmem   = as.numeric(results_lmem[seq(4, length(results_lmem[,1]), 4),4])

plot(avg_surprisal*-1,chi_entropy_lmem, col="red", ylab="Chi-squared", xlab="Average surprisal", 
     ylim=c(min(chi_entropy_lmem,chi_surprisal_lmem),max(chi_entropy_lmem,chi_surprisal_lmem)))
points(avg_surprisal*-1, chi_surprisal_lmem, col="blue")
grid (NULL,NULL, lty = 6, col = "grey") 
abline(a=	9.488,b=0, lty=2)
abline(a=	-9.488,b=0, lty=2)
title('Self-paced RT (LMEM)')
legend(-median(avg_surprisal)-1, max(chi_entropy_lmem, chi_surprisal_lmem), legend=c("NW-Entropy", "Surprisal"),
       col=c("red", "blue"), lty=1:2, cex=0.7)
