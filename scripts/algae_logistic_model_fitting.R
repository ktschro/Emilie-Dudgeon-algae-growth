#load in ----
setwd("/Users/katieschroeder/Documents/GitHub/Emilie-Dudgeon-algae-growth")
algae <- read.csv("raw-data/algae_counts.csv")

library(tidyverse)
library(magrittr)

# data cleaning ---- 
# Convert raw counts to actual counts
colnames(algae)[5] <- "raw_count"
algae$count <- ifelse(algae$dilution=="none",algae$raw_count,algae$raw_count*10)

# maybe do this at some point? I'm not sure how much it changes the
# initial parameters and I don't want to mess with it too much
# convert count on haemocytometer to count in flask
# algae$count <- algae$haemo_count / 8 * 10000 * 250

algae_all <- algae %>% filter(!is.na(count))

# logistic model 1: each replicate individually ---- 
fit_logistic_model_all <- function(algae_all) {
  time <- algae_all$exp_day
  population <- algae_all$count
  
  logistic_growth <- function(time, r, K) {
    N0 <- population[1]
    N <- N0 * K / (N0 + (K - N0) * exp(-r * time))
    return(N)
  }
  
  model <- optim(par = c(r = 1.5, K = 120),
                 fn = function(par) sum((population - logistic_growth(time, par[1], par[2]))^2),
                 method = "BFGS",
                 control = list(maxit = 200))
  
  r_hat <- model$par[1]
  K_hat <- model$par[2]
  
  # Predict population size for each day using the fitted model
  predicted_population <- logistic_growth(time, r_hat, K_hat)
  
  return(data.frame(temp = unique(algae_all$temp),
                    exp_day = time,
                    observed_population = population,
                    predicted_population = predicted_population,
                    r = r_hat, 
                    K = K_hat,
                    N0 = population[1]))
}

results_all <- algae_all %>%
  group_by(temp,rep) %>%
  do(fit_logistic_model_all(.))

# get predictions and plot each replicaate
temps <- unique(algae$temp)
sequence_values <- seq(1,50,by=0.25)

results_all_pred <- data.frame(
  exp_day = rep(sequence_values,length(temps)),
  temp = rep(temps,each=length(sequence_values)))

randK_all <- results_all %>% 
  select(r,K,N0,temp,rep) %>% 
  unique()

results_all_pred <- merge(results_all_pred,randK_all,by="temp")

# use r, K, and N0 to get model values
results_all_pred$N <- results_all_pred$N0 * results_all_pred$K / 
  (results_all_pred$N0 + 
     (results_all_pred$K - results_all_pred$N0) * 
     exp(-results_all_pred$r * results_all_pred$exp_day))

#plot individual replicates and model fit
#constant temps
results_all_pred %>% filter(!str_detect(temp,"V")) %>%
  ggplot(aes(x=exp_day,y=N)) +
  geom_line() +
  facet_grid(rows=vars(temp),cols=vars(rep),scales = "free") +
  geom_hline(aes(yintercept=K),color="darksalmon",
             linetype="dashed")+
  geom_point(data = filter(algae_all,!str_detect(temp,"V")),
             mapping = aes(x=exp_day,y=count),
             size=2,alpha=0.2)+
  labs(y="Population size",x="Time (days)") +
  theme_bw()


#variable temps
results_all_pred %>% filter(str_detect(temp,"V")) %>%
  ggplot(aes(x=exp_day,y=N)) +
  geom_line() +
  facet_grid(rows=vars(temp),cols=vars(rep),scales = "free") +
  geom_hline(aes(yintercept=K),color="darksalmon",
             linetype="dashed")+
  geom_point(data = filter(algae_all,str_detect(temp,"V")),
             mapping = aes(x=exp_day,y=count),
             size=2,alpha=0.2)+
  labs(y="Population size",x="Time (days)") +
  theme_bw()

# logistic model 1: summarized results ---- 
results_all %>% 
  filter(r <= 5) %>%
  ggplot(aes(x=temp,y=r)) +
  geom_boxplot() +
  theme_classic()

results_all %>% 
  filter(K <= 7500) %>%
  ggplot(aes(x=temp,y=K)) +
  geom_boxplot() +
  theme_classic()

results_all_sum <- results_all %>% 
  select(temp,rep,r,K) %>%
  unique() %>%
  group_by(temp) %>%
  summarize(mean_r = mean(r),
            mean_K = mean(K),
            se_r = sd(r)/sqrt(n()),
            se_K = sd(K)/sqrt(n()))

ggplot(results_all_sum,aes(x=temp,y=mean_r)) +
  geom_point() +
  geom_errorbar(aes(ymin=mean_r-se_r,ymax=mean_r+se_r), width = 0.2) +
  theme_classic()

# logistic model 2: summarized replicates ---- 
algae_sum <- algae %>% group_by(temp,exp_day) %>%
  summarize(count = mean(count,na.rm=T))

fit_logistic_model_sum <- function(algae_all) {
  time <- algae_all$exp_day
  population <- algae_all$count
  
  logistic_growth <- function(time, r, K) {
    N0 <- population[1]
    N <- N0 * K / (N0 + (K - N0) * exp(-r * time))
    return(N)
  }
  
  model <- optim(par = c(r = 1.5, K = 120),
                 fn = function(par) sum((population - logistic_growth(time, par[1], par[2]))^2),
                 method = "BFGS",
                 control = list(maxit = 200))
  
  r_hat <- model$par[1]
  K_hat <- model$par[2]
  
  return(data.frame(temp = unique(algae_all$temp),
                    exp_day = time,
                    observed_population = population,
                    r = r_hat, 
                    K = K_hat,
                    N0 = population[1]))
}

results_sum <- algae_sum %>%
  group_by(temp) %>%
  do(fit_logistic_model_sum(.))

results_sum_pred <- data.frame(
  exp_day = rep(sequence_values,length(temps)),
  temp = rep(temps,each=length(sequence_values)))

# get r and K for each temperature
randK_sum <- results_sum %>% 
  select(r,K,N0,temp) %>% 
  unique()

results_sum_pred <- merge(results_sum_pred,randK_sum,by="temp")

# use r, K, and N0 to get model values
results_sum_pred$N <- results_sum_pred$N0 * results_sum_pred$K / 
  (results_sum_pred$N0 + 
     (results_sum_pred$K - results_sum_pred$N0) * 
     exp(-results_sum_pred$r * results_sum_pred$exp_day))

# plot fitted model ----
ggplot(results_sum_pred,aes(x=exp_day,y=N)) +
  geom_line() +
  facet_wrap(~temp) +
  geom_hline(aes(yintercept=K),color="darksalmon",
             linetype="dashed")+
  geom_point(data = algae_all,
             mapping = aes(x=exp_day,y=count),
             size=2,alpha=0.2)+
  labs(y="Population size",x="Time (days)") +
  theme_classic()

# now compare r and K for each temperature treatment ----
# first look
ggplot(results_sum,aes(x=temp,y=r)) +
  geom_point() +
  theme_classic()

ggplot(results_sum,aes(x=temp,y=K)) +
  geom_point() +
  theme_classic()

# compare r's and K's for each method ----
randK_sum$rep <- NA
randK_sum$method <- "summarized"
randK_all$method <- "by replicate"

compare <- rbind(randK_all,randK_sum)

# plot
compare %>%
  filter(r <=5) %>%
  ggplot(aes(x=temp,y=r,color=method)) +
  geom_boxplot() +
  theme_classic()

compare %>%
  filter(K <= 7500) %>%
  ggplot(aes(x=temp,y=K,color=method)) +
  geom_boxplot() +
  theme_classic()

# real temp data ----
#compare with HOBO data ----
hobo<-read.csv("processed-data/hobo_summary_daily.csv") 
hobo<-hobo[,2:8]

colnames(compare)[4] <- "treatment_name"

merged_results_all <- merge(compare,hobo,by="treatment_name")

#graph slope (r) by mean temperature from HOBO data with temp range (daily) as the color    
merged_results_all %<>% mutate(fluct_size = case_when(
  str_detect(treatment_name,"S")~"small",
  str_detect(treatment_name,"M")~"medium",
  str_detect(treatment_name,"B")~"big",
  temp_type=="constant" ~ "constant"
),
  fluct_size_num = case_when(
  temp_desc == "temp_range" & temp_rounded <= 5 ~ 0,
  temp_desc == "temp_range" & temp_rounded == 7 ~ 7.5,
  temp_desc == "temp_range" & temp_rounded == 8 ~ 7.5,
  temp_desc == "temp_range" & temp_rounded == 15 ~ 15,
  temp_desc == "temp_range" & temp >= 20 ~ 22,
  TRUE ~ NA))

# try to look at size of fluctuation
merged_results_all %>%
  filter(temp_desc == "mean_temp") %>%
  ggplot(aes(x=temp,y=mean_r,
             color=fluct_size)) +
  scale_color_manual(values = c("gray40","black","grey60","grey80"))+
  stat_smooth(data = filter(merged_results_all,
                            temp_type=="constant" &
                              temp_desc=="mean_temp"),
              method = "lm", formula = y ~ x + I(x^2), size = 1,se=FALSE) +
  labs("x=Temperature (deg C)",y="Intrinsic growth rate, r") +
  geom_point() +
  theme_classic()

merged_results_all$temp_rounded <- round(merged_results_all$temp,digits=0)

merged_results_all %>%
  filter(temp_desc == "mean_temp") %>%
  ggplot(aes(x=temp,y=r,
             color=fluct_size)) +
  scale_color_manual(values = c("gray40","black","grey60","grey80"))+
  stat_smooth(data = filter(merged_results_all,
                            temp_type=="constant" &
                              temp_desc=="mean_temp"),
              method = "lm", formula = y ~ x + I(x^2), size = 1,se=FALSE) +
  labs("x=Temperature (deg C)",y="Intrinsic growth rate, r") +
  geom_point() +
  facet_wrap(~method) +
  theme_classic()

merged_results_all %>%
  filter(temp_desc == "mean_temp" &
           temp_type == "constant") %>%
  ggplot(aes(x=temp_rounded,y=r)) +
  geom_point() +
  stat_smooth(data = filter(merged_results_all,
                            temp_type=="constant" &
                              temp_desc=="mean_temp"),
              method = "lm", formula = y ~ x + I(x^2), 
              size = 1,se=FALSE, color="grey") +
  labs("x=Temperature (deg C)",y="Intrinsic growth rate, r") +
  facet_wrap(~method) +
  theme_classic()

# for K
merged_results_all %>%
  filter(temp_desc == "mean_temp" &
           temp_type == "constant") %>%
  filter(treatment_name != "34") %>%
  ggplot(aes(x=temp_rounded,y=K)) +
  geom_point() +
  stat_smooth(data = filter(merged_results_all,
                            temp_type=="constant" &
                              temp_desc=="mean_temp" & 
                              treatment_name != "34"),
              method = "lm", formula = y ~ x + I(x^2), 
              size = 1,se=FALSE, color="grey") +
  labs("x=Temperature (deg C)",y="Carrying capacity, K") +
  facet_wrap(~method) +
  theme_classic()

