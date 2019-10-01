# Christoph Aurnhammer, October 30 2018
# Analyse snapshots on SPR data (Frank et al., 2013) in a nonlinear and linear model
rm(list=ls())

data = read.csv('../Input_data/PREV_et_entropy_avg_snapshots.csv', sep = '\t')
data$word = as.character(data$word)
data$item = as.factor(data$item)
data$sent_nr = as.factor(data$sent_nr)
data$subj_nr = as.factor(data$subj_nr)
data$prevfix = as.factor(data$prevfix)
data = na.omit(data)

################ as lmem #######################
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
  m_baseline = lmer(RTfirstpass ~ 
                      # fixed effects
                      log_freq +
                      log_freq_prev +
                      nr_char +
                      nr_char_prev +
                      word_pos +
                      prevfix +
                      # interactions 
                      log_freq:log_freq_prev + log_freq:nr_char + log_freq:nr_char_prev +
                      log_freq:word_pos + log_freq:prevfix +
                      log_freq_prev:nr_char + log_freq_prev:nr_char_prev + log_freq_prev:word_pos +
                      log_freq_prev:prevfix + 
                      nr_char:nr_char_prev + nr_char:word_pos + nr_char:prevfix +
                      nr_char_prev:word_pos + nr_char_prev:prevfix +
                      word_pos:prevfix +
                      # random effects
                      (1|subj_nr) +
                      (1|item) +
                      (0+log_freq|subj_nr) +
                      (0+log_freq_prev|subj_nr) +
                      (0+nr_char|subj_nr) +
                      (0+nr_char_prev|subj_nr) +
                      (0+word_pos|subj_nr) +
                      (0+prevfix|subj_nr),
                    data = data, REML = FALSE, verbose = 1,
                    control=lmerControl(optCtrl=list(75000)))
  
  m_surprisal = lmer(RTfirstpass ~ 
                       # fixed effects
                       log_freq +
                       log_freq_prev +
                       nr_char +
                       nr_char_prev +
                       word_pos +
                       prevfix +
                       data[,14+snapshots[mod]] + 
                       data[,prevs[mod]] +
                       # interactions 
                       log_freq:log_freq_prev + log_freq:nr_char + log_freq:nr_char_prev +
                       log_freq:word_pos + log_freq:prevfix +
                       log_freq_prev:nr_char + log_freq_prev:nr_char_prev + log_freq_prev:word_pos +
                       log_freq_prev:prevfix + 
                       nr_char:nr_char_prev + nr_char:word_pos + nr_char:prevfix +
                       nr_char_prev:word_pos + nr_char_prev:prevfix +
                       word_pos:prevfix +
                       # random effects
                       (1|subj_nr) +
                       (1|item) +
                       (0+log_freq|subj_nr) +
                       (0+log_freq_prev|subj_nr) +
                       (0+nr_char|subj_nr) +
                       (0+nr_char_prev|subj_nr) +
                       (0+word_pos|subj_nr) +
                       (0+prevfix|subj_nr) +
                       (0+data[,14+snapshots[mod]]|subj_nr) +
                       (0+data[,prevs[mod]]|subj_nr),
                     data = data, REML = FALSE, verbose = 1,
                     control=lmerControl(optCtrl=list(75000)))
  
  m_entropy = lmer(RTfirstpass ~ 
                     # fixed effects
                     log_freq +
                     log_freq_prev +
                     nr_char +
                     nr_char_prev +
                     word_pos +
                     prevfix +
                     data[,14+snapshots[mod]+9] +
                     data[,prevs[mod]+9] +
                     # interactions 
                     log_freq:log_freq_prev + log_freq:nr_char + log_freq:nr_char_prev +
                     log_freq:word_pos + log_freq:prevfix +
                     log_freq_prev:nr_char + log_freq_prev:nr_char_prev + log_freq_prev:word_pos +
                     log_freq_prev:prevfix + 
                     nr_char:nr_char_prev + nr_char:word_pos + nr_char:prevfix +
                     nr_char_prev:word_pos + nr_char_prev:prevfix +
                     word_pos:prevfix +
                     # random effects
                     (1|subj_nr) +
                     (1|item) +
                     (0+log_freq|subj_nr) +
                     (0+log_freq_prev|subj_nr) +
                     (0+nr_char|subj_nr) +
                     (0+nr_char_prev|subj_nr) +
                     (0+word_pos|subj_nr) +
                     (0+prevfix|subj_nr) +
                     (0+data[,14+snapshots[mod]+1]|subj_nr) +
                     (0+data[,prevs[mod]+9]|subj_nr),
                   data = data, REML = FALSE, verbose = 1,
                   control=lmerControl(optCtrl=list(75000)))
  
  m_complete = lmer(RTfirstpass ~ 
                      # fixed effects
                      log_freq +
                      log_freq_prev +
                      nr_char +
                      nr_char_prev +
                      word_pos +
                      prevfix +
                      data[,14+snapshots[mod]] + 
                      data[,prevs[mod]] +
                      data[,14+snapshots[mod]+9] +
                      data[,prevs[mod]+9] +
                      # interactions 
                      log_freq:log_freq_prev + log_freq:nr_char + log_freq:nr_char_prev +
                      log_freq:word_pos + log_freq:prevfix +
                      log_freq_prev:nr_char + log_freq_prev:nr_char_prev + log_freq_prev:word_pos +
                      log_freq_prev:prevfix + 
                      nr_char:nr_char_prev + nr_char:word_pos + nr_char:prevfix +
                      nr_char_prev:word_pos + nr_char_prev:prevfix +
                      word_pos:prevfix +
                      # random effects
                      (1|subj_nr) +
                      (1|item) +
                      (0+log_freq|subj_nr) +
                      (0+log_freq_prev|subj_nr) +
                      (0+nr_char|subj_nr) +
                      (0+nr_char_prev|subj_nr) +
                      (0+word_pos|subj_nr) +
                      (0+prevfix|subj_nr) +
                      (0+data[,14+snapshots[mod]]|subj_nr) + 
                      (0+data[,prevs[mod]]|subj_nr) +
                      (0+data[,14+snapshots[mod]+1]|subj_nr) +
                      (0+data[,prevs[mod]+9]|subj_nr),
                    data = data, REML = FALSE, verbose = 1,
                    control=lmerControl(optCtrl=list(75000)))
  
  # compare models & append to output file & write to disc
  # surprisal over baseline
  comp = anova(m_baseline, m_surprisal)
  chisq = round(comp$Chisq[2],4)
  if (coef(summary(m_surprisal))[8] < 0){
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
  if (coef(summary(m_entropy))[8] < 0){
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
  if (coef(summary(m_complete))[8] < 0){
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
  if (coef(summary(m_complete))[10] < 0){
    chisq = chisq * -1
  }
  results_lmem[nrow(results_lmem) + 1,] = c(models[mod],
                                            "m_surprisal",
                                            "m_complete",
                                            chisq,
                                            comp$`Chi Df`[2],
                                            round(comp$`Pr(>Chisq)`[2],4))
  
  # Write results to disc
  write.table(results_lmem, file="../Output_data/PREV_et_snapshots_LMEM.csv", row.names = FALSE)
  
  print("-------------------------------------")
  # End loop
}

# Plot LMEM results
results_lmem = read.csv("../Output_data/PREV_et_snapshots_LMEM.csv", sep = " ")
dev.off()
plot.new()
avg_surprisal = c(8.145235, 7.816667, 7.646470, 7.281974, 6.880306, 6.407189, 5.878788, 5.531965, 5.423403)
avg_entropy   = c(7.114768, 6.509349, 6.658838, 6.312363, 5.885305, 5.258519, 4.861018, 4.608529, 4.445464)
chi_surprisal_lmem = as.numeric(results_lmem[seq(3, length(results_lmem[,1]), 4),4])
chi_entropy_lmem   = as.numeric(results_lmem[seq(4, length(results_lmem[,1]), 4),4])

plot(avg_surprisal*-1,chi_entropy_lmem, col="red", ylab="Chi-squared", xlab="Average surprisal", 
     ylim=c(min(chi_entropy_lmem,chi_surprisal_lmem),max(chi_entropy_lmem,chi_surprisal_lmem)))
points(avg_surprisal*-1, chi_surprisal_lmem, col="blue")
grid (NULL,NULL, lty = 6, col = "grey") 
abline(a=	9.488,b=0, lty=2)
abline(a=	-9.488,b=0, lty=2)
title('First-pass RT (LMEM)')
legend(-median(avg_surprisal), max(chi_entropy_lmem, chi_surprisal_lmem), legend=c("NW-Entropy", "Surprisal"),
       col=c("red", "blue"), lty=1:2, cex=0.7)
