#this script takes all of the HOBO data, restricts it to the data from when the experiment actually occurred and binds it to the temperature data
#at the end, it also makes a dataframe with summary information for each logger to use with the counts data later

#set up ---- 
#load in the hobo data (binded version)
setwd("/Users/katieschroeder/Documents/GitHub/Emilie-Dudgeon-algae-growth")

hobo_all <- read.csv("edited-files/hobo_all.csv")

#load in the algae count data so that we know which dates we need from the loggers
algae <- read.csv("raw-data/algae_counts.csv")
unique(algae$date)

#experiment started 7/24/2024 and ended 8/12/2024
#use 7/25/2024 to 8/11/2024 just in case first couple days aren't representative

#date formatting ----
hobo_all$date_time <- strptime(as.character(hobo_all$date_time),"%m/%d/%Y %H:%M:%S")
hobo_all$date <- as.Date(hobo_all$date_time)

hobo_algae <- hobo_all %>% filter(date >= "2024-07-27" & date <= "2024-08-12") 

#graph raw data just to do some checking ----
hobo_algae$time <- format( 
  as.POSIXct(hobo_algae$date_time), format = "%H:%M:%S")
library(chron)
hobo_algae$time2 <- chron(times=hobo_algae$time)

hobo_algae %>% 
  filter(date >= "2024-07-27" & date <= "2024-08-10") %>% #filtered again here 
  #because there were some weird temp sequences during the early days of the experiment
  ggplot(aes(x=time2,y=temp_C,group=date)) +
  geom_point() +
  geom_line() +
  facet_wrap(~logger_ID) +
  theme_classic()


#summarize by each day for each logger ----
hobo_day_summary <- hobo_algae %>% group_by(logger_ID,date) %>%
  summarize(mean_temp = mean(temp_C),
            max_temp = max(temp_C),
            min_temp = min(temp_C),
            sd_temp = sd(temp_C),
            temp_range = max_temp - min_temp,
            mean_light = mean(light_lux))

#plot it out
library(ggrepel)

#mean temp
hobo_day_summary %>% 
  mutate(label = ifelse(date==max(date), as.character(logger_ID),NA_character_)) %>%
  ggplot(aes(x=date,y=mean_temp,color=logger_ID)) +
  geom_point() +
  geom_line(aes(group=logger_ID)) +
  theme_classic()+
  geom_label_repel(aes(label = label),
                                   nudge_x = 1,
                                   na.rm = TRUE) +
  theme(legend.position="none")

#look at light. It's probably kind of bad
hobo_day_summary %>% 
  mutate(label = ifelse(date==max(date), as.character(logger_ID),NA_character_)) %>%
  ggplot(aes(x=date,y=mean_light,color=logger_ID)) +
  geom_point() +
  geom_line(aes(group=logger_ID)) +
  theme_classic()+
  geom_label_repel(aes(label = label),
                   nudge_x = 1,
                   na.rm = TRUE) +
  theme(legend.position="none")

#now overall summary table ----
#summarize for each logger
hobo_summary <- hobo_algae %>% group_by(logger_ID) %>%
  summarize(mean_temp = mean(temp_C),
            max_temp = max(temp_C),
            min_temp = min(temp_C),
            sd_temp = sd(temp_C),
            temp_range = max_temp - min_temp,
            mean_light = mean(light_lux))

#planned vs actual ----
#bind with temp information to see how far off we were
planned_temp <- read.csv("raw-data/temp_trts.csv")
colnames(hobo_summary)[1] <- "treatment_name"
hobo_summary$temp_type <- ifelse(str_detect(hobo_summary$treatment_name,"V"),"fluctuating","constant")
colnames(planned_temp)[1:6] <- c("treatment_name","planned_mean_temp","temp_type","planned_temp_range","planned_min_temp","planned_max_temp")

hobo_compare <- merge(planned_temp,hobo_summary,by=c("treatment_name","temp_type"))

#reconfiguring to make graphing easier
hobo_compare_long <- hobo_compare %>% pivot_longer(cols=planned_mean_temp:temp_range,
                                                   names_to = "temp_desc",
                                                   values_to = "temp")
hobo_compare_long$planned <-
  ifelse(str_detect(hobo_compare_long$temp_desc,"plan"),"planned","actual") 

hobo_compare_long$temp_var <- gsub("planned_","",hobo_compare_long$temp_desc)

#visualize differences between planned and actual
hobo_compare_long %>% filter(temp_var!="min_temp"&temp_var!="max_temp") %>%
  ggplot(aes(x=planned,y=temp,color=planned)) +
  geom_point() +
  geom_line(aes(group=treatment_name)) +
  facet_wrap(~temp_var,scale="free_y") +
  geom_label_repel(aes(label = treatment_name),
                   nudge_x = 1,
                   na.rm = TRUE) +
  theme_classic()

#temp range is from min and max temp, probably not super representative of daily temp range. Take daily hobo summary and get average daily temp range
hobo_daily_temp_range <- hobo_day_summary %>% 
  select(logger_ID,temp_range,mean_light) %>%
  group_by(logger_ID) %>%
  summarize(temp_range_daily = mean(temp_range),
            mean_light = mean(mean_light))

colnames(hobo_daily_temp_range)[1:2] <- c("treatment_name","temp")
hobo_daily_temp_range$planned <-"actual_daily"
hobo_daily_temp_range$temp_desc <- "temp_range"
hobo_daily_temp_range$temp_var <- "temp_range"
hobo_daily_temp_range$temp_type <- ifelse(
  str_detect(hobo_daily_temp_range$treatment_name,"V"),
  "fluctuating","constant")

hobo_compare_long<-rbind(hobo_daily_temp_range,hobo_compare_long)

#graph it out
hobo_compare_long %>% filter(temp_var=="temp_range"|temp_var=="mean_temp") %>%
  ggplot(aes(x=planned,y=temp,color=planned)) +
  geom_point() +
  geom_line(aes(group=treatment_name)) +
  geom_label_repel(aes(label = treatment_name),
                   nudge_x = 1,
                   na.rm = TRUE) +
  facet_wrap(~temp_var) +
  theme_classic()

#doing it all again with daily means ---- 
hobo_summary_daily <- hobo_day_summary %>% group_by(logger_ID) %>%
  summarize(mean_temp_daily = mean(mean_temp),
            max_temp_daily = mean(max_temp),
            min_temp_daily = mean(min_temp),
            sd_temp_daily = sd(mean_temp),
            temp_range_daily = max_temp_daily - min_temp_daily,
            mean_light = mean(mean_light))

hobo_summary_daily %<>% 
  pivot_longer(cols=contains("temp"), names_to= "temp_var",values_to = "temp") %>%
  mutate(planned = "actual_daily",
         temp_var = gsub("_daily","",temp_var))

colnames(hobo_summary_daily)[1] <- "treatment_name"
hobo_summary_daily$temp_desc <- hobo_summary_daily$temp_var
hobo_summary_daily$temp_type <- ifelse(
  str_detect(hobo_summary_daily$treatment_name,"V"),
  "fluctuating","constant")

write.csv(hobo_summary_daily,"processed-data/hobo_summary_daily.csv")

compare_all <- rbind(hobo_compare_long,hobo_summary_daily)

#all graphs
compare_all %>% 
  mutate(label2 = ifelse(planned=="planned", as.character(treatment_name),NA_character_)) %>%
  filter(temp_var!="sd_temp") %>%
  ggplot(aes(x=planned,y=temp,color=temp_type)) +
  geom_point() +
  geom_line(aes(group=treatment_name)) +
  facet_wrap(~temp_var,scale="free_y",nrow=1) +
  geom_label_repel(aes(label = label2),
                   nudge_x = 1,
                   na.rm = TRUE) +
  theme_classic()
ggsave("plots/exploratory/HOBO_compare_min_max_mean_range.png",width=15,height=5,units="in",dpi=600)

#just mean and temp_range
compare_all %>% 
  filter(temp_var=="mean_temp"|temp_var=="temp_range") %>%
  mutate(label2 = ifelse(planned=="planned", as.character(treatment_name),NA_character_)) %>%
  ggplot(aes(x=planned,y=temp,color=temp_type)) +
  geom_point() +
  geom_line(aes(group=treatment_name)) +
  facet_wrap(~temp_var,scale="free_y") +
  geom_label_repel(aes(label = label2),
                   nudge_x = 1,
                   na.rm = TRUE) +
  theme_classic()
ggsave("plots/exploratory/HOBO_compare_mean_range.png",width=10,height=5,units="in",dpi=600)


