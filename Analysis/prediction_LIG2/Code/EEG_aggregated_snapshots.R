# Christoph Aurnhammer, Jan 17 2019
# Analyse snapshots on aggregated N400 data (Frank et al., 2015) in a nonlinear and linear model

rm(list=ls())
library(mgcv)
library(itsadug)

data = read.csv('../Input_data/eeg_srp_kld_avg_snapshots.csv', sep = '\t')
data$subj_nr  = as.factor(data$subj_nr)
data$sent_nr  = as.factor(data$sent_nr)
data$item     = as.factor(data$item)
data$reject_data = as.factor(data$reject_data)
data$reject_word = as.factor(data$reject_word)
data = na.omit(data)

################ as lmem #######################
library(lme4)

# Prepare output file
results_lmem = data.frame(model=character(), baseline=character(), interest=character(), chi2=numeric(), df=numeric(), pval=numeric(),
                     stringsAsFactors = FALSE)

# Loop: for snapshot in snapshots
mod_num = c(1, 2, 3, 4, 5, 6, 7, 8, 9)
snapshots = c(0, 2, 4, 6, 8, 10, 12, 14, 16)
models = c("1k", "3k", "10k", "30k", "100k", "300k", "1M", "3M", "6-46M")

# build baseline, surprisal, kld, complete 
for (mod in mod_num){
  print(models[mod])
  print(snapshots[mod])
  print(colnames(data)[12+snapshots[mod]])
  print(colnames(data)[12+snapshots[mod]+1])
  # Build models
  print("baseline")
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
  
  print("surprisal")
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
  
  print("kld")
  m_kld = lmer(N400 ~ 
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
  
  print("complete")
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
  
  # kld over baseline
  comp = anova(m_baseline, m_kld)
  chisq = round(comp$Chisq[2],4)
  if (coef(summary(m_kld))[6] < 0){
    chisq = chisq * -1
  }
  results_lmem[nrow(results_lmem) + 1,] = c(models[mod],
                                            "m_baseline",
                                            "m_kld",
                                            chisq,
                                            comp$`Chi Df`[2],
                                            round(comp$`Pr(>Chisq)`[2],4))
  
  # surprisal over kld
  comp = anova(m_kld, m_complete)
  chisq = round(comp$Chisq[2],4)
  if (coef(summary(m_complete))[6] < 0){
    chisq = chisq * -1
  }
  results_lmem[nrow(results_lmem) + 1,] = c(models[mod],
                                            "m_kld",
                                            "m_complete",
                                            chisq,
                                            comp$`Chi Df`[2],
                                            round(comp$`Pr(>Chisq)`[2],4))
  
  # kld over surprisal
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
avg_surprisal = c(7.4309, 7.1081, 6.9922, 6.6013, 6.1707, 5.6739, 5.1746, 4.8661, 4.7802)
chi_surprisal_lmem = as.numeric(results_lmem[seq(3, length(results_lmem[,1]), 4),4])
chi_kld_lmem   = as.numeric(results_lmem[seq(4, length(results_lmem[,1]), 4),4])

plot(avg_surprisal*-1,chi_kld_lmem, col="red", ylab="Chi-squared", xlab="Average surprisal", 
     ylim=c(min(chi_kld_lmem,chi_surprisal_lmem),max(chi_kld_lmem,chi_surprisal_lmem)))
points(avg_surprisal*-1, chi_surprisal_lmem, col="blue")
grid (NULL,NULL, lty = 6, col = "grey") 
abline(a=	5.991,b=0, lty=2)
abline(a=	-5.991,b=0, lty=2)
title('N400 (LMEM)')
legend(-median(avg_surprisal), max(chi_surprisal_lmem, chi_kld_lmem)+2, legend=c("KLD", "Surprisal"),
       col=c("red", "blue"), lty=1:2, cex=0.5)

