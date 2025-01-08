#script purpose: take algal counts from haemocytometer taken over time and calculate algal population growth rates

setwd("/Users/katieschroeder/Documents/GitHub/Emilie-Dudgeon-algae-growth")
algae <- read.csv("raw-data/algae_counts.csv")

library(tidyverse)
library(magrittr)

#raw data ---- 
#convert diluted counts to real count
colnames(algae)[5] <- "raw_count"

algae$count <- ifelse(algae$dilution=="none",algae$raw_count,algae$raw_count*10)

#first look: algal growth over time in each treatment:
ggplot(algae,aes(x=exp_day,y=count,group=rep)) +
  geom_point() +
  geom_line() +
  facet_wrap(~temp) +
  theme_classic() +
  labs(y="Count of algal cells on haemocytometer",x="Day of experiment")

#look at just the maximum count (the last day)
#first, add some temp treatment information
algae %<>%
  mutate(mean_temp = ifelse(str_detect(temp,"V"),
                            as.numeric(gsub("\\D", "",temp)),
                            as.numeric(temp)),
         temp_var = case_when(
           str_detect(temp,"VM") ~ "medium",
           str_detect(temp,"VS") ~ "small",
           str_detect(temp,"VB") ~ "big",
           TRUE ~ "none"
         ))

algae$mean_temp <- ifelse(algae$mean_temp == 23.0, 23.5, algae$mean_temp)

#now plot max at last experiment day (19) as a boxplot for each with mean temp as the x-axis
algae %>% filter(exp_day == 19) %>%
  ggplot(aes(x=factor(mean_temp),y=count,fill=as.factor(temp_var))) +
  geom_boxplot() +
  #facet_wrap(~temp_var,nrow=1) +
  theme_classic()

ggsave("plots/exploratory/algae_last_day_count_boxplot.png",width=12,height=5,units="in",dpi=600)

#daily summary for each treatment
#summarize within each treatment and day
algae_sum <- algae %>% group_by(temp,exp_day) %>%
  dplyr::summarize(mean_count = mean(count,na.rm=T),
            sd_count = sd(count,na.rm=T),
            se_count = sd(count,na.rm=T)/sqrt(n()),
            mean_temp = unique(mean_temp),
            temp_var = unique(temp_var))

#calculate log N
algae_sum$log_mean_count <- log(algae_sum$mean_count)

#graph log N vs t for each treatment and add rough trendlines
algae_sum %>%
  ggplot(aes(x=exp_day,y=log_mean_count)) +
  geom_point() +
  facet_wrap(~temp) +
  geom_smooth(method='lm',se=FALSE)+
  theme_classic()

ggsave("plots/exploratory/all_trt_log_mean_count_time_lms.png",width=10,height=8,units="in",dpi=600)

#specific growth rate should be slope of the line of ln(density) vs time
#do a bunch of linear models ---- 
list_of_algae_dfs = split(algae_sum, algae_sum$temp)

results = lapply(list_of_algae_dfs, function(dat) lm(log_mean_count ~ exp_day, data = dat))
lapply(results, summary)

algae_pval <- sapply(results, function(x) summary(x)$coefficients)
algae_pval 

algae_pval <- as.data.frame(t(algae_pval))

algae_pval$treatment_name <- rownames(algae_pval)

colnames(algae_pval)[1:8] <- c("intercept","slope","intercept_se","slope_se",
                               "intercept_t_val","slope_t_val","intercept_p_value",
                               "slope_p_value")
algae_pval_ordered <- algae_pval[, c(9,1,2,3,4,5,6,7,8)]

write.csv(algae_pval_ordered,"processed-data/algae_growth_rates.csv")

#graph slopes (r) for each treatment ----
ggplot(algae_pval_ordered,aes(x=treatment_name,y=slope)) +
  geom_point() +
  theme_classic()

#merge with HOBO data to get actual temp response
hobo<-read.csv("processed-data/hobo_summary_daily.csv") 
hobo<-hobo[,2:8]

merged <- merge(algae_pval_ordered,hobo,by="treatment_name")

#graph slope (r) by mean temperature from HOBO data with temp range (daily) as the color    
merged %<>% mutate(fluct_size = case_when(
  str_detect(treatment_name,"S")~"small",
  str_detect(treatment_name,"M")~"medium",
  str_detect(treatment_name,"B")~"big",
  temp_type=="constant" ~ "constant"
))

merged$temp_rounded <- round(merged$temp,digits=0)

merged %>% filter(temp_var=="temp_range"|temp_var=="mean_temp") %>%
  pivot_wider(names_from= "temp_var" ,values_from = "temp_rounded") %>%
  mutate(temp_range2 = lead(temp_range)) %>%
  filter(!is.na(mean_temp)) %>%
  ggplot(aes(x=mean_temp,y=slope,color=fluct_size)) +
  geom_point() +
  #geom_errorbar(aes(ymin=slope-slope_se,ymax=slope+slope_se,alpha=0.3,width=0.1))+
  theme_classic() +
  geom_line(aes(group=fluct_size))
ggsave("plots/exploratory/all_trt_algal_growth_temp_corrected.png",width=10,height=8,units="in",dpi=600)

#just constant 
merged %>% filter(temp_var=="temp_range"|temp_var=="mean_temp") %>%
  pivot_wider(names_from= "temp_var" ,values_from = "temp_rounded") %>%
  filter(temp_type=="constant") %>%
  mutate(temp_range2 = lead(temp_range)) %>%
  filter(!is.na(mean_temp)) %>%
  ggplot(aes(x=mean_temp,y=slope,color=fluct_size)) +
  geom_point() +
  #geom_errorbar(aes(ymin=slope-slope_se,ymax=slope+slope_se,alpha=0.3,width=0.1))+
  theme_classic() +
  geom_line(aes(group=fluct_size))
ggsave("plots/exploratory/constant_algal_growth_temp_corrected.png",width=10,height=8,units="in",dpi=600)

#per capita growth rates - get r and K ---- 
#get t2 and N2 columns
algae <- algae %>% group_by(temp,rep) %>% 
  mutate(t2 = lead(exp_day),
         N2 = lead(count)) %>%
  ungroup()

#calculate population growth (dN/dt)
algae$pop_growth <- (algae$N2-algae$count)/(algae$t2-algae$exp_day)
algae$per_capita <- algae$pop_growth/algae$count
algae$ID <- paste(algae$temp,algae$rep,sep="_")

ggplot(algae,aes(x=count,y=per_capita,group=ID)) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  facet_wrap(~temp) +
  theme_classic() +
  labs(x="N",y="Per capita growth rate (1/N dN/dt)")

#do a bunch of lms
list_of_algae_dfs = split(algae, algae$ID)

results_algae = lapply(list_of_algae_dfs, function(dat) lm(per_capita ~ count, data = dat))
lapply(results_algae, summary)

algae_pval <- sapply(results_algae, function(x) summary(x)$coefficients)
algae_pval

algae_pval <- as.data.frame(t(algae_pval))

algae_pval$ID <- rownames(algae_pval)

colnames(algae_pval)[1:8] <- c("intercept","slope","intercept_se","slope_se",
                              "intercept_t_val","slope_t_val","intercept_p_value",
                              "slope_p_value")
algae_pval_ordered <- algae_pval[, c(9,1,2,3,4,5,6,7,8)]

algae_pval_ordered$r <- algae_pval_ordered$intercept
algae_pval_ordered$K <- algae_pval_ordered$intercept / algae_pval_ordered$slope * -1

#add in actual temperature info
hobo<-read.csv("processed-data/hobo_summary_daily.csv") 
hobo<-hobo[,2:8]

algae_pval_ordered$treatment_name <- gsub("_.","",algae_pval_ordered$ID)

merged <- merge(algae_pval_ordered,hobo,by="treatment_name")

#graph slope (r) by mean temperature from HOBO data with temp range (daily) as the color    
merged %<>% mutate(fluct_size = case_when(
  str_detect(treatment_name,"S")~"small",
  str_detect(treatment_name,"M")~"medium",
  str_detect(treatment_name,"B")~"big",
  temp_type=="constant" ~ "constant"
))

merged$temp_rounded <- round(merged$temp,digits=0)

#plot r and K, raw values
#constant treatments
merged %>%
  filter(slope_p_value<=0.5) %>%
  filter(r<=2) %>%
  filter(temp_type=="constant") %>%
  filter(temp_desc=="mean_temp") %>%
  ggplot(aes(x=temp_rounded,y=r,group=treatment_name)) +
  labs("x=Temperature (deg C)",y="Intrinsic growth rate, r") +
  geom_boxplot(varwidth=F) +
  #facet_wrap(~fluct_size) +
  theme_classic()

merged %>%
  filter(slope_p_value<=0.5) %>%
  filter(temp_type=="constant") %>%
  filter(temp_desc=="mean_temp") %>%
  ggplot(aes(x=temp_rounded,y=K,group=treatment_name)) +
  labs("x=Temperature (deg C)",y="Carrying capacity, K") +
  geom_boxplot(varwidth=F) +
  #facet_wrap(~fluct_size) +
  theme_classic()

#now for fluctuating
merged$fluct_size <- factor(merged$fluct_size,levels=c("constant","small","medium", "big"))
#small = 7C range, medium = 15C range, large = 22C range

merged %>%
  filter(slope_p_value<=0.5) %>%
  filter(r<2) %>%
  filter(temp_rounded >=15 & temp_rounded <=30) %>%
  filter(temp_desc=="mean_temp") %>%
  ggplot(aes(x=temp_rounded,y=r,group=treatment_name,fill=fluct_size)) +
  labs("x=Temperature (deg C)",y="Intrinsic growth rate, r") +
  scale_fill_manual(values=c("gray","lightblue","lightblue4","blue4"))+
  geom_boxplot(varwidth=F) +
  #facet_wrap(~fluct_size) +
  theme_classic()
  
merged %>%
  filter(temp_desc=="mean_temp") %>%
  filter(slope_p_value<=0.5) %>%
  ggplot(aes(x=temp_rounded,y=K,fill=fluct_size,group=treatment_name)) +
  geom_boxplot(varwidth=F) +
  facet_wrap(~fluct_size) +
  theme_classic()

#get mean values for each, use cut off of p=0.4 for the slope
algae_r_K_sum <- algae_pval_ordered %>%
  filter(slope_p_value <= 0.5) %>%
  group_by(treatment_name) %>%
  dplyr::summarize(mean_r = mean(r,na.rm=T),
                   se_r = sd(r)/sqrt(n()),
                   mean_K = mean(K),
                   se_K = sd(K)/sqrt(n()),
                   count = n())

#add temp data back in
algae_r_K_sum <- merge(algae_r_K_sum,hobo,by="treatment_name")
algae_r_K_sum %<>% mutate(fluct_size = case_when(
  str_detect(treatment_name,"S")~"small",
  str_detect(treatment_name,"M")~"medium",
  str_detect(treatment_name,"B")~"big",
  temp_type=="constant" ~ "constant"
))

algae_r_K_sum$temp_rounded <- round(algae_r_K_sum$temp,digits=0)

#r and K for constant
algae_r_K_sum %>%
  filter(temp_type=="constant") %>%
  filter(temp_desc=="mean_temp") %>%
  ggplot(aes(x=temp_rounded,y=mean_r,group=treatment_name)) +
  labs("x=Temperature (deg C)",y="Intrinsic growth rate, r") +
  geom_point() +
  geom_errorbar(aes(ymin=mean_r-se_r,ymax=mean_r+se_r),width=0.2)+
  theme_classic() 

algae_r_K_sum %>%
  filter(temp_type=="constant") %>%
  filter(temp_desc=="mean_temp") %>%
  ggplot(aes(x=temp_rounded,y=mean_K,group=treatment_name)) +
  labs("x=Temperature (deg C)",y="Carrying capacity, K") +
  geom_point() +
  geom_errorbar(aes(ymin=mean_K-se_K,ymax=mean_K+se_K),width=0.2)+
  theme_classic() 

#constant vs fluctuating comparisons
algae_r_K_sum$fluct_size <- factor(algae_r_K_sum$fluct_size,levels=c("constant","small","medium", "big"))

algae_r_K_sum %>%
  filter(temp_rounded >=15 & temp_rounded <=30) %>%
  filter(temp_desc=="mean_temp") %>%
  ggplot(aes(x=temp_rounded,y=mean_r,group=treatment_name,color=fluct_size)) +
  labs("x=Temperature (deg C)",y="Intrinsic growth rate, r") +
  geom_point() +
  geom_errorbar(aes(ymin=mean_r-se_r,ymax=mean_r+se_r),width=0.2)+
  scale_color_manual(values=c("lightblue","gray80","gray60","gray40")) +
  theme_classic()

algae_r_K_sum %<>%
  mutate(temp_comp = case_when(
    temp_rounded == 14 ~ 14,
    temp_rounded == 18 ~ 20,
    temp_rounded == 19 ~ 20,
    temp_rounded == 20 ~ 20,
    temp_rounded == 23 ~ 23.5,
    temp_rounded == 24 ~ 23.5,
    temp_rounded == 26 ~ 27,
    temp_rounded == 27 ~ 27,
    temp_rounded == 32 ~ 32,
    temp_rounded == 34 ~ 34
  ))

algae_r_K_sum %>%
  filter(temp_rounded >=15 & temp_rounded <=30) %>%
  filter(temp_desc=="mean_temp") %>%
  ggplot(aes(x=fluct_size,y=mean_r,group=treatment_name,color=fluct_size)) +
  labs("x=Temperature (deg C)",y="Intrinsic growth rate, r") +
  geom_point() +
  geom_errorbar(aes(ymin=mean_r-se_r,ymax=mean_r+se_r),width=0.2)+
  scale_color_manual(values=c("red3","gray80","gray60","gray40")) +
  facet_wrap(~temp_comp) +
  theme_classic()

algae_r_K_sum %>%
  filter(temp_rounded >=15 & temp_rounded <=30) %>%
  filter(temp_desc=="mean_temp") %>%
  ggplot(aes(x=fluct_size,y=mean_K,group=treatment_name,color=fluct_size)) +
  labs("x=Temperature (deg C)",y="Carrying capacity, K") +
  geom_point() +
  geom_errorbar(aes(ymin=mean_K-se_K,ymax=mean_K+se_K),width=0.2)+
  scale_color_manual(values=c("red3","gray80","gray60","gray40")) +
  facet_wrap(~temp_comp) +
  theme_classic()


# now logistic growth model ----
library(dplyr)
library(ggplot2)
library(nls2)

# Sample data (replace with your actual data)
test <- 

data <- data.frame(
  time = c(0, 3, 6, 9, 12, 15, 18, 21),
  population = c(10, 25, 50, 80, 120, 150, 170, 180)
)

# Fit the logistic growth model
model <- nls(population ~ K / (1 + A * exp(-r * time)), 
             data = data, 
             start = list(K = 200, r = 0.1, A = 0.5))

# Extract parameters
K_est <- coef(model)["K"]
r_est <- coef(model)["r"]

# Print the estimated parameters
cat("Estimated carrying capacity (K):", K_est, "\n")
cat("Estimated intrinsic growth rate (r):", r_est, "\n")

# Visualize the fit
ggplot(data, aes(x = time, y = population)) +
  geom_point() +
  geom_line(aes(y = fitted(model)), color = "red") +
  labs(x = "Time (days)", y = "Population")



library(nls2)
  # Function to fit the logistic growth model and extract parameters
  fit_logistic_model <- function(df) {
    model <- nls(population ~ K / (1 + A * exp(-r * time)), 
                 data = df, 
                 start = list(K = 200, r = 1, A = 0.5))
    
    K_est <- coef(model)["K"]
    r_est <- coef(model)["r"]
    
    return(data.frame(K = K_est, r = r_est))
  }
  
  # Group the data by replicate and group, and fit the model for each group
  test <- algae %>% rename(time = exp_day,
                           population = count)
  
  results <- test %>%
    group_by(rep, temp) %>%
    nest() %>%
    mutate(model_fit = map(data, fit_logistic_model)) %>%
    unnest(model_fit)
  
  # Print the results
  print(results)
  
  # Visualize the fits for a specific replicate and group (adjust as needed)
  replicate_num <- 1
  group_letter <- "A"
  
  data_subset <- data %>%
    filter(replicate == replicate_num, group == group_letter)
  
  ggplot(data_subset, aes(x = time, y = population)) +
    geom_point() +
    geom_line(aes(y = fitted(model)), color = "red") +
    labs(x = "Time (days)", y = "Population")
