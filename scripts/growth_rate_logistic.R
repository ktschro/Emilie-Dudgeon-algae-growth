#script purpose: estimate r and K using logistic growth equation for populations of Ankistrodesmus falcatus algae grown under different temperatures

setwd("/Users/katieschroeder/Documents/GitHub/Emilie-Dudgeon-algae-growth")
algae <- read.csv("raw-data/algae_counts.csv")

library(tidyverse)
library(magrittr)

# Convert raw counts to actual counts
colnames(algae)[5] <- "raw_count"
algae$count <- ifelse(algae$dilution=="none",algae$raw_count,algae$raw_count*10)

# maybe do this at some point? I'm not sure how much it changes the
# initial parameters and I don't want to mess with it too much
# convert count on haemocytometer to count in flask
# algae$count <- algae$haemo_count / 8 * 10000 * 250

# summarize within each temperature on each day before fitting logistic function
algae_sum <- algae %>% group_by(temp,exp_day) %>%
  summarize(count = mean(count,na.rm=T))

algae_all <- algae %>% filter(!is.na(count))

# Summarized model: define model ----
fit_logistic_model <- function(algae_sum) {
  time <- algae_sum$exp_day
  population <- algae_sum$count
  
  logistic_growth <- function(time, r, K) {
    N0 <- population[1]
    N <- N0 * K / (N0 + (K - N0) * exp(-r * time))
    return(N)
  }
  
  model <- optim(par = c(r = 1, K = 120),
                 fn = function(par) sum((population - logistic_growth(time, par[1], par[2]))^2),
                 method = "BFGS",
                 control = list(maxit = 200))
  
  r_hat <- model$par[1]
  K_hat <- model$par[2]
  
  return(data.frame(temp = unique(algae_sum$temp),
                    exp_day = time,
                    observed_population = population,
                    r = r_hat, 
                    K = K_hat,
                    N0 = population[1]))
}

# Group the data by temperature and run
results <- algae_sum %>%
  group_by(temp) %>%
  do(fit_logistic_model(.))

# Summarized model: predictions ---- 
# Define the 17 groups
temps <- unique(algae$temp)
sequence_values <- seq(1,50,by=0.25)

results_pred <- data.frame(
  exp_day = rep(sequence_values,length(temps)),
  temp = rep(temps,each=length(sequence_values)))

# get r and K for each temperature
randK <- results %>% 
  select(r,K,N0,temp) %>% 
  unique()

results_pred <- merge(results_pred,randK,by="temp")

# use r, K, and N0 to get model values
results_pred$N <- results_pred$N0 * results_pred$K / 
  (results_pred$N0 + 
     (results_pred$K - results_pred$N0) * 
     exp(-results_pred$r * results_pred$exp_day))

# Summarized model: plotting ----
ggplot(results_pred,aes(x=exp_day,y=N)) +
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
ggplot(results,aes(x=temp,y=r)) +
  geom_point() +
  theme_classic()

ggplot(results,aes(x=temp,y=K)) +
  geom_point() +
  theme_classic()

# # free scales so see more of plot shape
# ggplot(results,aes(x=exp_day,y=observed_population)) +
#   geom_point() +
#   geom_line(aes(x=exp_day,y=predicted_population),color="skyblue4") +
#   geom_hline(aes(yintercept=K),color="darksalmon",linetype="dashed")+
#   facet_wrap(~temp,scales="free") +
#   labs(y="Population size",x="Time (days)") +
#   theme_classic()
# 
# #set scales to compare between temperatures
# ggplot(results,aes(x=exp_day,y=observed_population)) +
#   geom_point() +
#   geom_line(aes(x=exp_day,y=predicted_population),color="skyblue4") +
#   geom_hline(aes(yintercept=K),color="darksalmon",linetype="dashed")+
#   facet_wrap(~temp) +
#   labs(y="Population size",x="Time (days)") +
#   theme_classic()

# adding uncertainty to r and K estimates ---- 
library(boot)

fit_logistic_model_boot <- function(algae_all, R = 1000) {
  time <- algae_all$exp_day
  population <- algae_all$count
  
  logistic_growth <- function(time, r, K) {
    N0 <- population[1]
    N <- N0 * K / (N0 + (K - N0) * exp(-r * time))
    return(N)
  }
  
  # Bootstrap function
  boot_func <- function(data, indices) {
    data_boot <- data[indices,]
    model_boot <- optim(par = c(r = 1, K = 200),
                        fn = function(par) sum((data_boot$count - logistic_growth(
                          data_boot$exp_day, par[1], par[2]))^2),
                        method = "Nelder-Mead",
                        control = list(maxit = 200))
    return(model_boot$par)
  }
  
  # Perform bootstrapping
  boot_results <- boot(algae_all, boot_func, R = R)
  
  # Calculate confidence intervals
  ci <- boot.ci(boot_results, type = "norm")
  
  # Extract parameter estimates and confidence intervals
  r_hat <- boot_results$t0[1]
  K_hat <- boot_results$t0[2]
  ci_r <- ci$normal[2:3] 
  ci_K <- ci$normal[4:5] 
  
  return(data.frame(temp = unique(algae_all$temp),
                    r = r_hat,
                    K = K_hat,
                    r_lower_ci = ci_r[1],
                    r_upper_ci = ci_r[2],
                    K_lower_ci = ci_K[1],
                    K_upper_ci = ci_K[2]))
}


# Group the data by temperature and apply the function
results_boot <- algae_all %>%
  group_by(temp) %>%
  do(fit_logistic_model_boot(.))

# plot
ggplot(results_boot,aes(x=temp,y=r)) +
  geom_point() +
  geom_errorbar(aes(ymin=r_lower_ci,ymax=r_upper_ci),width=0.2) +
  theme_classic()

ggplot(results_boot,aes(x=temp,y=K)) +
  geom_point() +
 #geom_errorbar(aes(ymin=K_lower_ci,ymax=K_upper_ci),width=0.2) +
  theme_classic()

# fit with unsummarized data ----
fit_logistic_model_all <- function(algae_all) {
  time <- algae_all$exp_day
  population <- algae_all$count
  
  logistic_growth <- function(time, r, K) {
    N0 <- population[1]
    N <- N0 * K / (N0 + (K - N0) * exp(-r * time))
    return(N)
  }
  
  model <- optim(par = c(r = 2, K = 150),
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

results_all_sum <- results_all %>% 
  select(temp,rep,r,K) %>%
  unique() %>%
  group_by(temp) %>%
  summarize(mean_r = mean(r),
            mean_K = mean(K),
            se_r = sd(r)/sqrt(n()),
            se_K = sd(K)/sqrt(n()))

# plot
ggplot(results_all,aes(x=exp_day,y=observed_population)) +
  geom_point() +
  geom_line(aes(x=exp_day,y=predicted_population),color="skyblue4") +
  geom_hline(aes(yintercept=K),color="darksalmon",linetype="dashed")+
  facet_wrap(~temp,scales="free") +
  labs(y="Population size",x="Time (days)") +
  theme_classic()


#compare with HOBO data ----
hobo<-read.csv("processed-data/hobo_summary_daily.csv") 
hobo<-hobo[,2:8]

colnames(results_all_sum)[1] <- "treatment_name"

merged_results_all <- merge(results_all_sum,hobo,by="treatment_name")

#graph slope (r) by mean temperature from HOBO data with temp range (daily) as the color    
merged_results_all %<>% mutate(fluct_size = case_when(
  str_detect(treatment_name,"S")~"small",
  str_detect(treatment_name,"M")~"medium",
  str_detect(treatment_name,"B")~"big",
  temp_type=="constant" ~ "constant"
))

#merged_results_all$temp_rounded <- round(merged_results_all$temp,digits=0)

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

merged_results_all %<>% mutate(fluct_size_num = case_when(
  str_detect(treatment_name,"S")~"small",
  str_detect(treatment_name,"M")~"medium",
  str_detect(treatment_name,"B")~"big",
  temp_type=="constant" ~ "constant"
))

merged_results_all %>%
  filter(temp_desc == "mean_temp") %>%
  ggplot(aes(x=fluct_size,y=r)) +
  labs("x=Temperature (deg C)",y="Intrinsic growth rate, r") +
  geom_point() +
  facet_wrap(~temp_rounded) +
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
  theme_classic()

# Arrhenius curve ---- 
constant_data <- merged_results_all %>% 
  filter(temp_type == "constant") %>%
  select(temp,temp_rounded,r,K,temp_desc) %>%
  filter(temp_desc == "mean_temp") %>% 
  mutate(temp_K = temp + 273.15) %>%
  unique()

# fit TPC with rTPC
library(rTPC)
library(nls.multstart)

mod = 'briere2_1999'

# get start vals
start_vals <- get_start_vals(constant_data$temp, constant_data$r,
                             model_name = 'briere2_1999')
# fit model
briere_mod <- nls_multstart(r~briere2_1999(temp = temp, tmin, tmax, a, b),
                                    data = constant_data,
                                    iter = c(4,4,4,4),
                                    start_lower = start_vals - 10,
                                    start_upper = start_vals + 10,
                                    lower = get_lower_lims(
                                      constant_data$temp, d$rate,
                                      model_name = 'briere2_1999'),
                                    upper = get_upper_lims(
                                      constant_data$temp, d$rate,
                                      model_name = 'briere2_1999'),
                                    supp_errors = 'Y',
                                    convergence_count = FALSE)

# look at model fit
summary(briere_mod)

# predictions from model fit
briere_pred <- data.frame(temp = seq(min(constant_data$temp), 
                                  max(constant_data$temp), 0.5))
briere_pred <- augment(fit, newdata = briere_pred)

# plot data and model fit
ggplot(constant_data, aes(temp, r)) +
  geom_point() +
  geom_line(aes(temp, .fitted), briere_pred, col = 'lightblue') +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Intrinsic growth rate, r')

# repeat with K
start_vals_K <- get_start_vals(constant_data$temp, constant_data$K, 
                             model_name = 'boatman_2017')

# get limits
low_lims_K <- get_lower_lims(constant_data$temp, constant_data$K, 
                           model_name = 'boatman_2017')
upper_lims_K <- get_upper_lims(constant_data$temp, constant_data$K, 
                             model_name = 'boatman_2017')

fit_K <- nls_multstart(r~boatman_2017(temp = temp, rmax, tmin, tmax, a, b),
                                            data = constant_data,
                                            iter = c(4,4,4,4,4),
                                            start_lower = start_vals - 10,
                                            start_upper = start_vals + 10,
                                            lower = lower_lims_K,
                                            upper = upper_lims_K,
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)
fit_K

~nls_multstart(rate~boatman_2017(temp = temp, rmax, tmin, tmax, a,b),
               data = .x,
               iter = c(4,4,4,4,4),
               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') - 10,
               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') + 10,
               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
               supp_errors = 'Y',
               convergence_count = FALSE)),

# predictions from model fit
preds_K <- augment(fit_K, newdata = new_data)

# plot data and model fit
ggplot(constant_data, aes(temp, K)) +
  geom_point() +
  geom_line(aes(temp, .fitted), preds_K, col = 'blue') +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Carrying capacity, K')

# trying to fit the modified Briere equation ----

# Define the modified Briere equation function
modified_briere <- function(temp, a, m, xmax) {
  y <- a * ((temp^m) / ((temp^m) + xmax^m))^(1/m)
  y[temp < 0 | temp > xmax] <- 0  # Set growth rate to 0 outside valid range
}

# Fit the model to your data
mod_briere_model <- nls(r ~ modified_briere(temp, a, m, xmax), 
             data = constant_data, 
             start = list(a = 1, m = 2, xmax = 10))

# Extract fitted parameters
fitted_a <- coef(model)[["a"]]
fitted_m <- coef(model)[["m"]]
fitted_xmax <- coef(model)[["xmax"]]

# approach with optim
your_data = constant_data %>% select(temp,r) %>%
  rename(temperature=temp,growth_rate=r)

# Define the objective function to minimize
objective_function <- function(params, temperature, growth_rate) {
  a <- params[1]
  m <- params[2]
  xmax <- params[3]

  
  predicted_growth_rate <- modified_briere(temperature, a, m, xmax)
  sum((predicted_growth_rate - growth_rate)^2)
}

# Initial parameter values
initial_params <- c(a = 0.1, m = 20, xmax = 100)

# Optimization using `optim()` with different methods and bounds
methods <- c("Nelder-Mead", "BFGS", "L-BFGS-B")
best_fit <- NULL
best_obj_value <- Inf

for (method in methods) {
  result <- optim(par = initial_params, fn = objective_function, 
                  method = method, 
                  lower = c(0, 0, 0, -Inf), upper = c(Inf, Inf, Inf, Inf),
                  temperature = your_data$temperature, 
                  growth_rate = your_data$growth_rate)
  
  if (result$value < best_obj_value) {
    best_fit <- result
    best_obj_value <- result$value
  }
}

# Extract fitted parameters from the best fit
fitted_params <- best_fit$par

# old code ----
# boltzmann_arrhenius_func <- function(T, A, Ea) {
#   Rk_B <- 1.380649e-23  # J/K Boltzmann constant
#   return(A * exp(-Ea / (k_B * T)))
# }

# # Fit the model
# boltz_arrhenius_fit <- nls(r ~ boltzmann_arrhenius_func(temp_K, A, Ea), 
#                            data = constant_data,
#                            start = list(A = 1, Ea = 1000000))
# A <- coef(arrhenius_fit)[1]
# Ea <- coef(arrhenius_fit)[2]
# 
# # add model predictions into data frame
# R <- 8.314
# arr_preds <- data.frame(
#   temp_pred_K = 
#     seq(min(constant_data$temp_K)-20,
#         max(constant_data$temp_K)+20,by=2))
# 
# arr_preds$rate_pred <- A*exp(-Ea / (R*arr_preds$temp_pred_K))
# 
# plot(arr_preds$temp_pred_K,arr_preds$rate_pred)
# 
# # plot
# merged_results_all %>%
#   filter(temp_desc == "mean_temp" &
#            temp_type == "constant") %>%
#   ggplot(aes(x=temp_rounded,y=r)) +
#   geom_point() +
#   geom_line(data = arr_preds, aes(x=temp_pred_C,y=rate_pred)) +
#   labs("x=Temperature (deg C)",y="Intrinsic growth rate, r") +
#   theme_classic()

