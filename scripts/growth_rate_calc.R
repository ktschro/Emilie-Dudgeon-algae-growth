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
  theme_classic()

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
  summarize(mean_count = mean(count,na.rm=T),
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
  theme_classic()

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
ggplot(algae_pval_ordered,aes(x=temp_trt,y=slope)) +
  geom_point() +
  theme_classic()

#merge with HOBO data to get actual temp response
hobo<-read.csv("processed-data/hobo_summary_daily.csv") 
hobo<-hobo[,2:8]

merged <- merge(algae_pval_ordered,hobo,by="treatment_name")
                  