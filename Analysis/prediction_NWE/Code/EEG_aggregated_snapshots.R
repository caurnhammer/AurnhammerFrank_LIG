# Christoph Aurnhammer, September 18 2018
# Analyse snapshots on aggregated N400 data (Frank et al., 2015) in a nonlinear and linear model
rm(list=ls())
library(mgcv)
library(itsadug)

data = read.csv('../Input_data/eeg_entropy_avg_snapshots.csv', sep = '\t')
data$subj_nr  = as.factor(data$subj_nr)
data$sent_nr  = as.factor(data$sent_nr)
data$item     = as.factor(data$item)
data$reject_data = as.factor(data$reject_data)
data$reject_word = as.factor(data$reject_word)

# Prepare output file
results = data.frame(model=character(), baseline=character(), interest=character(), chi2=numeric(), df=numeric(), pval=numeric(),
                     stringsAsFactors = FALSE)

# Loop: for snapshot in snapshots
mod_num = c(1, 2, 3, 4, 5, 6, 7, 8, 9)
snapshots = c(0, 2, 4, 6, 8, 10, 12, 14, 16)
models = c("1k", "3k", "10k", "30k", "100k", "300k", "1M", "3M", "6-46M")

# build baseline, surprisal, entropy, complete 
for (mod in mod_num){
  print(models[mod])
  print(snapshots[mod])
  print(colnames(data)[12+snapshots[mod]])
  print(colnames(data)[12+snapshots[mod]+1])
  # Build models
  m_baseline = bam(N400 ~ 
                     # fixed effects
                     s(word_pos, k = length(unique(data$wordpos))) + 
                     s(nr_char, k = length(unique(data$nrchar))) +
                     s(log_freq) +
                     s(baseline) +
                     # interactions
                     ti(word_pos, nr_char, k = min(length(unique(data$word_pos)), length(unique(data$nr_char)))) +
                     ti(word_pos, log_freq, k = length(unique(data$word_pos))) +
                     ti(word_pos, baseline) +
                     ti(nr_char, log_freq, k = length(unique(data$nr_char))) +
                     ti(nr_char, baseline) +
                     ti(log_freq, baseline) +
                     # random effects
                     s(subj_nr, bs = 're') +
                     s(item, bs = 're') +
                     s(word_pos, subj_nr, bs = 'fs', m=1) +
                     s(nr_char, subj_nr, bs = 'fs', m=1) +
                     s(log_freq, subj_nr, bs = 'fs', m=1) +
                     s(baseline, subj_nr, bs = 'fs', m=1),
                     data = data, discrete=TRUE, nthreads = 4)
  
  m_surprisal = bam(N400 ~ 
                     # fixed effects
                     s(word_pos, k = length(unique(data$wordpos))) + 
                     s(nr_char, k = length(unique(data$nrchar))) +
                     s(log_freq) +
                     s(baseline) +
                     s(data[,12+snapshots[mod]]) +
                     # interactions
                     ti(word_pos, nr_char, k = min(length(unique(data$word_pos)), length(unique(data$nr_char)))) +
                     ti(word_pos, log_freq, k = length(unique(data$word_pos))) +
                     ti(word_pos, baseline) +
                     ti(nr_char, log_freq, k = length(unique(data$nr_char))) +
                     ti(nr_char, baseline) +
                     ti(log_freq, baseline) +
                     # random effects
                     s(subj_nr, bs = 're') +
                     s(item, bs = 're') +
                     s(word_pos, subj_nr, bs = 'fs', m=1) +
                     s(nr_char, subj_nr, bs = 'fs', m=1) +
                     s(log_freq, subj_nr, bs = 'fs', m=1) +
                     s(baseline, subj_nr, bs = 'fs', m=1) +
                     s(data[,12+snapshots[mod]], subj_nr, bs = 'fs', m=1),
                   data = data, discrete=TRUE, nthreads = 4)
  
  m_entropy = bam(N400 ~ 
                    # fixed effects
                    s(word_pos, k = length(unique(data$wordpos))) + 
                    s(nr_char, k = length(unique(data$nrchar))) +
                    s(log_freq) +
                    s(baseline) +
                    s(data[,12+snapshots[mod]+1]) +
                    # interactions
                    ti(word_pos, nr_char, k = min(length(unique(data$word_pos)), length(unique(data$nr_char)))) +
                    ti(word_pos, log_freq, k = length(unique(data$word_pos))) +
                    ti(word_pos, baseline) +
                    ti(nr_char, log_freq, k = length(unique(data$nr_char))) +
                    ti(nr_char, baseline) +
                    ti(log_freq, baseline) +
                    # random effects
                    s(subj_nr, bs = 're') +
                    s(item, bs = 're') +
                    s(word_pos, subj_nr, bs = 'fs', m=1) +
                    s(nr_char, subj_nr, bs = 'fs', m=1) +
                    s(log_freq, subj_nr, bs = 'fs', m=1) +
                    s(baseline, subj_nr, bs = 'fs', m=1) +
                    s(data[,12+snapshots[mod]+1], subj_nr, bs = 'fs', m=1),
                  data = data, discrete=TRUE, nthreads = 4)
  
  m_complete = bam(N400 ~ 
                     # fixed effects
                     s(word_pos, k = length(unique(data$wordpos))) + 
                     s(nr_char, k = length(unique(data$nrchar))) +
                     s(log_freq) +
                     s(baseline) +
                     s(data[,12+snapshots[mod]]) + 
                     s(data[,12+snapshots[mod]+1]) +
                     # interactions
                     ti(word_pos, nr_char, k = min(length(unique(data$word_pos)), length(unique(data$nr_char)))) +
                     ti(word_pos, log_freq, k = length(unique(data$word_pos))) +
                     ti(word_pos, baseline) +
                     ti(nr_char, log_freq, k = length(unique(data$nr_char))) +
                     ti(nr_char, baseline) +
                     ti(log_freq, baseline) +
                     # random effects
                     s(subj_nr, bs = 're') +
                     s(item, bs = 're') +
                     s(word_pos, subj_nr, bs = 'fs', m=1) +
                     s(nr_char, subj_nr, bs = 'fs', m=1) +
                     s(log_freq, subj_nr, bs = 'fs', m=1) +
                     s(baseline, subj_nr, bs = 'fs', m=1) +
                     s(data[,12+snapshots[mod]], subj_nr, bs = 'fs', m=1) +
                     s(data[,12+snapshots[mod]+1], subj_nr, bs = 'fs', m=1),
                   data = data, discrete=TRUE, nthreads = 4)
  
  # compare models & append to output file & write to disc
  # surprisal over baseline
  comp = compareML(m_baseline, m_surprisal)
  results[nrow(results) + 1,] = c(models[mod],
                                  as.character(comp$table[1,1][1]),
                                  as.character(comp$table[2,1][1]),
                                  levels(comp$table[2,4])[2],
                                  levels(comp$table[2,5])[2],
                                  levels(comp$table[2,6])[2])
  
  # entropy over baseline
  comp = compareML(m_baseline, m_entropy)
  results[nrow(results) + 1,] = c(models[mod],
                                  as.character(comp$table[1,1][1]),
                                  as.character(comp$table[2,1][1]),
                                  levels(comp$table[2,4])[2],
                                  levels(comp$table[2,5])[2],
                                  levels(comp$table[2,6])[2])
  
  # surprisal over entropy
  comp = compareML(m_entropy, m_complete)
  results[nrow(results) + 1,] = c(models[mod],
                                  as.character(comp$table[1,1][1]),
                                  as.character(comp$table[2,1][1]),
                                  levels(comp$table[2,4])[2],
                                  levels(comp$table[2,5])[2],
                                  levels(comp$table[2,6])[2])

  # entropy over surprisal
  comp = compareML(m_surprisal, m_complete)
  results[nrow(results) + 1,] = c(models[mod],
                                  as.character(comp$table[1,1][1]),
                                  as.character(comp$table[2,1][1]),
                                  levels(comp$table[2,4])[2],
                                  levels(comp$table[2,5])[2],
                                  levels(comp$table[2,6])[2])
  
  # Write results to disc
  write.table(results, file="../Output_data/eeg_snapshots_GAMM.csv", row.names = FALSE)
  
  print("-------------------------------------")
  # End loop
}

# Plot avg surprisal vs chi-square (for surprisal and entropy; over and above each other)
results = read.csv("../Output_data/eeg_snapshots_GAMM.csv", sep=" ")
avg_surprisal = c(7.525636, 7.177759, 7.035030, 6.647705, 6.233933, 5.762444, 5.274674, 4.948322, 4.850785)
avg_entropy   = c(7.091436, 6.507045, 6.689682, 6.445902, 6.021535, 5.415326, 5.034600, 4.763739, 4.608751)
chi_surprisal = as.numeric(results[seq(3, length(results[,1]), 4),4])
chi_entropy   = as.numeric(results[seq(4, length(results[,1]), 4),4])

plot(avg_surprisal*-1,chi_entropy, col="red", ylab="Chi-squared", xlab="Average surprisal", 
     ylim=c(min(chi_entropy,chi_surprisal),max(chi_entropy,chi_surprisal)))
points(avg_surprisal*-1, chi_surprisal, col="blue")
grid (NULL,NULL, lty = 6, col = "grey") 
abline(a=9.488, b=0, lty=2)
abline(a=-9.488, b=0, lty=2)
title('N400 (GAMM)')
legend(-median(avg_surprisal), max(chi_entropy, chi_surprisal), legend=c("NW-Entropy", "Surprisal"),
       col=c("red", "blue"), lty=1:2, cex=0.7)   

################ as lmem #######################
library(lme4)

# Prepare output file
results_lmem = data.frame(model=character(), baseline=character(), interest=character(), chi2=numeric(), df=numeric(), pval=numeric(),
                     stringsAsFactors = FALSE)

# Loop: for snapshot in snapshots
mod_num = c(1, 2, 3, 4, 5, 6, 7, 8, 9)
snapshots = c(0, 2, 4, 6, 8, 10, 12, 14, 16)
models = c("1k", "3k", "10k", "30k", "100k", "300k", "1M", "3M", "6-46M")

# build baseline, surprisal, entropy, complete 
for (mod in mod_num){
  print(models[mod])
  print(snapshots[mod])
  print(colnames(data)[12+snapshots[mod]])
  print(colnames(data)[12+snapshots[mod]+1])
  # Build models
  m_baseline = lmer(N400 ~ 
                     # fixed effects
                     word_pos + 
                     nr_char +
                     log_freq +
                     baseline +
                     # interactions
                     word_pos * nr_char +
                     word_pos * log_freq +
                     word_pos * baseline +
                     nr_char * log_freq +
                     nr_char * baseline +
                     log_freq * baseline +
                     # random effects
                     (1|subj_nr) + 
                     (1|item) +
                     (0+word_pos|subj_nr) +
                     (0+nr_char|subj_nr) +  
                     (0+log_freq|subj_nr) + 
                     (0+baseline|subj_nr),
                   data = data, REML = FALSE)
  
  m_surprisal = lmer(N400 ~ 
                      # fixed effects
                      word_pos + 
                      nr_char +
                      log_freq +
                      baseline +
                      data[,12+snapshots[mod]] +
                      # interactions
                      word_pos * nr_char +
                      word_pos * log_freq +
                      word_pos * baseline +
                      nr_char * log_freq +
                      nr_char * baseline +
                      log_freq * baseline +
                      # random effects
                      (1|subj_nr) + 
                      (1|item) +
                      (0+word_pos|subj_nr) +
                      (0+nr_char|subj_nr) +  
                      (0+log_freq|subj_nr) + 
                      (0+baseline|subj_nr) +
                      (0+data[,12+snapshots[mod]]|subj_nr),
                      data = data, REML = FALSE)
  
  m_entropy = lmer(N400 ~ 
                    # fixed effects
                    word_pos + 
                    nr_char +
                    log_freq +
                    baseline +
                    data[,12+snapshots[mod]+1] +
                    # interactions
                    word_pos * nr_char +
                    word_pos * log_freq +
                    word_pos * baseline +
                    nr_char * log_freq +
                    nr_char * baseline +
                    log_freq * baseline +
                    # random effects
                    (1|subj_nr) + 
                    (1|item) +
                    (0+word_pos|subj_nr) +
                    (0+nr_char|subj_nr) +  
                    (0+log_freq|subj_nr) + 
                    (0+baseline|subj_nr) +
                    (0+data[,12+snapshots[mod]+1]|subj_nr),
                  data = data, REML = FALSE)
  
  m_complete = lmer(N400 ~ 
                     # fixed effects
                     word_pos + 
                     nr_char +
                     log_freq +
                     baseline +
                     data[,12+snapshots[mod]] +
                     data[,12+snapshots[mod]+1] +
                     # interactions
                     word_pos * nr_char +
                     word_pos * log_freq +
                     word_pos * baseline +
                     nr_char * log_freq +
                     nr_char * baseline +
                     log_freq * baseline +
                     # random effects
                     (1|subj_nr) + 
                     (1|item) +
                     (0+word_pos|subj_nr) +
                     (0+nr_char|subj_nr) +  
                     (0+log_freq|subj_nr) + 
                     (0+baseline|subj_nr) +
                     (0+data[,12+snapshots[mod]]|subj_nr) +
                     (0+data[,12+snapshots[mod]+1]|subj_nr),
                   data = data, REML = FALSE)
  
  # compare models & append to output file & write to disc
  # surprisal over baseline
  comp = anova(m_baseline, m_surprisal)
  chisq = round(comp$Chisq[2],4)
  if (coef(summary(m_surprisal))[6] < 0){
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
  if (coef(summary(m_entropy))[6] < 0){
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
  if (coef(summary(m_complete))[6] < 0){
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
  if (coef(summary(m_complete))[7] < 0){
    chisq = chisq * -1
  }
  results_lmem[nrow(results_lmem) + 1,] = c(models[mod],
                                            "m_surprisal",
                                            "m_complete",
                                            chisq,
                                            comp$`Chi Df`[2],
                                            round(comp$`Pr(>Chisq)`[2],4))
  
  # Write results to disc
  write.table(results_lmem, file="../Output_data/eeg_snapshots_LMEM.csv", row.names = FALSE)
  
  print("-------------------------------------")
  # End loop
}

# Plot
results_lmem = read.csv("../Output_data/eeg_snapshots_LMEM.csv", sep = " ")
dev.off()
plot.new()
avg_surprisal = c(7.525636, 7.177759, 7.035030, 6.647705, 6.233933, 5.762444, 5.274674, 4.948322, 4.850785)
avg_entropy   = c(7.091436, 6.507045, 6.689682, 6.445902, 6.021535, 5.415326, 5.034600, 4.763739, 4.608751)
chi_surprisal_lmem = as.numeric(results_lmem[seq(3, length(results_lmem[,1]), 4),4])
chi_entropy_lmem   = as.numeric(results_lmem[seq(4, length(results_lmem[,1]), 4),4])

plot(avg_surprisal*-1,chi_entropy_lmem, col="red", ylab="Chi-squared", xlab="Average surprisal", 
     ylim=c(min(chi_entropy_lmem,chi_surprisal_lmem),max(chi_entropy_lmem,chi_surprisal_lmem)))
points(avg_surprisal*-1, chi_surprisal_lmem, col="blue")
grid (NULL,NULL, lty = 6, col = "grey") 
abline(a=	5.991,b=0, lty=2)
abline(a=	-5.991,b=0, lty=2)
title('N400 (LMEM)')
legend(-median(avg_surprisal), max(chi_surprisal_lmem, chi_entropy_lmem)+2, legend=c("NW-Entropy", "Surprisal"),
       col=c("red", "blue"), lty=1:2, cex=0.7)
