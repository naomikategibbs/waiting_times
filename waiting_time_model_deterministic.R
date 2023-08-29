
# This script is for a decision model estimating the health impact of waiting 
# for elective procedures in the NHS in England

# The model accounts for age/sex/IMD quintile
# The model can be used for various health conditions
# The week that proc is carried out is set exogenously

#### load packages ####
library(tidyverse)
library(flexsurv)
library(stats)
library(readxl)

#### CHOOSE PROCEDURE (e.g. cabg, hip, cataract etc) ####

# set the working directory to whichever procedure you want to use by changing
# the final part of the path (e.g. cabg, hip, cataract etc)
setwd("G:/Shared drives/ESHCRU EEPRU Waiting Times project/Waiting times model/hip")


#### load models and data and functions ####
# all of the inputs needed for the model (parameters and regression models) are 
# kept in the input folder for the specific condition of interest

# three regression models related to survival, created using HES data linked with
# ONS mortality data, repeated for females and males
periop_model_female <- readRDS("inputs/survival/periop_model_female.Rds")
days_less_30_model_female <- readRDS("inputs/survival/days_less_30_model_female.Rds")
survival_model_female <- readRDS("inputs/survival/survival_model_female.Rds")

periop_model_male <- readRDS("inputs/survival/periop_model_male.Rds")
days_less_30_model_male <- readRDS("inputs/survival/days_less_30_model_male.Rds")
survival_model_male <- readRDS("inputs/survival/survival_model_male.Rds")

# the variance covariance matrix corresponding to the models, these had to be 
# saved as separate objects to allow for extraction from the safehaven
periop_model_vcov_female <- readRDS("inputs/survival/periop_model_vcov_female.Rds")
days_less_30_model_vcov_female <- readRDS("inputs/survival/days_less_30_model_vcov_female.Rds")
survival_model_vcov_female <- readRDS("inputs/survival/survival_model_vcov_female.Rds")

periop_model_vcov_male <- readRDS("inputs/survival/periop_model_vcov_male.Rds")
days_less_30_model_vcov_male <- readRDS("inputs/survival/days_less_30_model_vcov_male.Rds")
survival_model_vcov_male <- readRDS("inputs/survival/survival_model_vcov_male.Rds")

# subgroup data drawn from the patients receiving the procedure
# between 2010 and 2019 in HES, it gives mean age at proc and charlson index by sex/imd
subgroup_data <- readRDS("inputs/demographics/subgroup_data.Rds")

# pre procedure utility measured using the EQ-5D, this relates to the condition 
# which the procedure relates to
utilities <- read_excel("inputs/utilities/pre_proc_utility.xlsx", 
                        sheet = "data")

# IMD adjustments for utility values
EQ5D_IMD <- read_excel("inputs/utilities/EQ-5D IMD gradients.xlsx", 
                                                    sheet = "Combined")

# post procedure utility, this is an EQ-5D improvement 
utility_improvement <- read_excel("inputs/utilities/Utility_improvement.xlsx", 
                                  sheet = "Scores")

# post procedure utility decrements based on waiting
postproc_HRQoL_decrements <- read_excel("inputs/utilities/utility_decrement.xlsx", sheet = "results")
postproc_HRQoL_decrements_vcov <- read_excel("inputs/utilities/utility_decrement.xlsx", sheet = "vcov")

# the probability of exiting the waiting list for emergency procedure
prob_exit_acute_df <-  read_excel("inputs/survival/prob_emergency_admission.xlsx", 
                               sheet = "data")

# the increase in perioperative mortality if the admission is an emergency
uplift_mortality_emergency_admission <-  read_excel("inputs/survival/uplift_mortality_emergency_admission.xlsx", 
                                                    sheet = "data")

# the R function which creates the mortality transition probabilities from the
# survival functions
source("inputs/survival/generating_transition_probabilities_function_deterministic.R")

# probability of switching to private whilst waiting
prob_exit_private_df <- read_excel("inputs/survival/prob_switch_private.xlsx", sheet = "data")

# general population mortality
life_tables_uk_imd_both <- readRDS("inputs/survival/life_tables_uk_imd.Rds")

#### setting deterministic model parameters ####

# The patient can either be waiting for proc, exit for private care, exit for 
# emergency admission, post proc or dead, either from a perioperative death 
# (within 30 days) or after 30 days. There are many reasons for exit whilst 
# on the waiting list but we model just three: dead, private, acute admission
state_names <- c("on_waiting_list","exit_acute","exit_private", "post_proc","peri_dead", "dead")
n_states <- length(state_names)

#  allocate where people start in the model
starting_population <- c(1,0,0,0,0,0) ## all people/person starts in the first state, waiting for proc
## order defined inline with state_names order

dr_o <- 0.035 ## set the discount rate for outcomes (3.5%)

sim_runs <- 1 ## sets the number of simulations for the PSA (we don't need 
# this line when we run it on viking)

#### the model function #####
waiting_times_model <- function(day_of_proc, age, female, imd_quintile,
                charlsonindex, mean_wait){

  # this is a function to run the waiting times model
  # INPUTS: day_of_proc = the day on which the procedure takes place, set by the modeller
  #         age = a continuous numeric value, age when added to wait list
  #         female = 1 if female, 0 if male
  #         imd_quintile = 1, 2, 3, 4, 5 (1 = poorest)
  #         charlson_index = a numeric value
  #         mean_wait = the mean wait time for the subgroup stays fixed
  #         even when the day of proc changes, this is because we are not able
  #         to get at the causal effect of waiting in the survival analysis and 
  #         so waiting appears to increase LYs for certain conditions which we know
  #         is incorrect
  # OUTPUTS:Vector of 12 values: QALYs, discounted QALYs, LYs, discounted LYs
  #         split into "normal", acute and private

# to run the code inside the function comment out the function
# and comment in the below parameters
# 
# day_of_proc <- 10
# age <-  64
# female <- 0
# imd_quintile <- 3
# charlsonindex <- 4.3
# mean_wait <- 63.045


#### 0: Setting the number of cycles for the model ####

# the model will run to age 101, in order to make it a lifetime model
# we need the age at the point they enter the waiting list
cycles <- (101 - age)*365   ## cycles in days
cycles_v <- 1:cycles ## a vector of days 
cycles_v_years <- c(1:(101 - (age))) ## a vector of years which will be used for discounting


#### 1: setting up utility parameters ####
  
## mean utility for pre proc utility, filtered for female yes or no
pre_proc_utility_base <- utilities$EQ5D_UK[utilities$female == female] 

# adding a coefficient for age allows for the fact that HRQoL deteriorates
# as the patient ages within the model, this is hard coded as it is the same 
# across all conditions. Sullivan et al. 2011

pre_proc_utility_age_decrement <- 0.0002747

# adding a coefficient allowing for the fact that HRQoL deteriorates
# as the patient waits within the model, this is only possible where there is
# evidence and is zero in most cases, it also only occurs up to the time
# period there is evidence for
pre_proc_utility_wait_decrement <- utilities$daily_decrement[utilities$female == female] ## mean coefficient to adjust utility by wait, pre proc

# creating a vector of pre procedure utility which deteriorates with age
# and daily decrements for waiting
pre_proc_utility <- pre_proc_utility_base - 
  pre_proc_utility_age_decrement * (age + cycles_v/365) -
           if_else(utilities$end_daily_decrement[utilities$female == female] > cycles_v, 
                   pre_proc_utility_wait_decrement * cycles_v, 
                   pre_proc_utility_wait_decrement * utilities$end_daily_decrement[utilities$female == female])

# we need to create a pre_proc_utility with no daily decrements (nodd) to add the uplift
# to otherwise we assume the decrements whilst waiting have post procedural impact,
# which we do not have evidence for at this point. 
pre_proc_utility_nodd <- pre_proc_utility_base - 
  pre_proc_utility_age_decrement * (age + cycles_v/365)

# post procedure utility, adding the pre procedure utility to a utility increase,
# based on age at day of procedure and sex
utility_imp <- utility_improvement$Difference[utility_improvement$Age == 
           round(age + day_of_proc/365) & utility_improvement$Female == female]

# adding the pre procedure utility and the utility improvement
post_proc_utility <- pre_proc_utility_nodd + utility_imp

# adding the post procedure decrements due to waiting, these are only available
# for hip and knee, for all other conditions they are zero.
m <- as.matrix(postproc_HRQoL_decrements[,2])

post_proc_utility <- post_proc_utility + day_of_proc*m[1,1] + (1-female)*day_of_proc*m[3,1]

# applying a gradient for IMD
# select the row based on age and female and then
# select the column that corresponds to the correct IMD_quintile. We use the IMD
# multiplication factor from the age the person enters the model
imd_multiplication_factor <- as.numeric(EQ5D_IMD[EQ5D_IMD$Age == 
         round(age) & EQ5D_IMD$Female == female, (imd_quintile + 2)])
  
# multiply the estimated health utility by the imd_multiplication_factor
pre_proc_utility <- pre_proc_utility * imd_multiplication_factor
post_proc_utility <- post_proc_utility * imd_multiplication_factor

# if dead then utility is zero. If patient exits for acute admission or private
# care then their payoffs are estimated separately, hence keep utility at zero here
exit_acute_utility <- 0
exit_private_utility <- 0
peridead_utility <- 0
dead_utility <- 0

# health state utilities
state_utilities <- data.frame(pre_proc_utility, exit_acute_utility, 
        exit_private_utility, post_proc_utility, peridead_utility, dead_utility) 

#### 2: setting up peri-op mortality parameters ####

# To use the regression models first create a dataframe of model inputs which 
# is in the same format as the imported regression models
df1 <- data.frame(day_of_proc)
df1$mean_wait <- mean_wait
df1 <- rename(df1, days_wait = mean_wait) # we are fixing mean wait for the subgroup
df1$agecont <- age + day_of_proc/365 # age at the point they have their procedure
df1$charlsonindex <- charlsonindex
df1$female <- female
df1$female <- factor(df1$female,c(0,1)) # to match the original regression model 
# female and IMD need to be factors
df1$imd_quintile <- imd_quintile
df1$imd_quintile <- factor(df1$imd_quintile,c(1,2,3,4,5)) # to match the 
# structure of the original regression model

# create a list with both the male and female models
periop_model_list <- list(periop_model_male, periop_model_female)
# then select the one you want
periop_model <- periop_model_list[[female + 1]]

# predicting the probability of dying within 30 days
periop_mortality <- predict(periop_model, newdata = df1, type = "response")

#### 3: setting up peri-days survived parameters ####

# create a list with both the male and female models
days_less_30_model_list <- list(days_less_30_model_male, days_less_30_model_female)
# then select the one you want
days_less_30_model <- days_less_30_model_list[[female + 1]]

# predict the number of days survived
peri_days_survived <- predict(days_less_30_model, newdata = df1, type = "response")  
  
#### 4: setting up post proc survival parameters ####

# Using the R function we loaded above, and inputting the survival model and 
# patient characteristics this outputs two vectors of transition 
# probabilities to death, one for pre and one for post procedure

# create a list with both the male and female models
survival_model_list <- list(survival_model_male, survival_model_female)
# then select the one you want
survival_model <- survival_model_list[[female + 1]]

# selecting for the female/male
life_tables_uk_imd <- life_tables_uk_imd_both[life_tables_uk_imd_both$female == female,]

# running the function to generate the pre and post transition probabilities to 
# death
tps_list <-  mortality_tps(df1, life_tables_uk_imd, survival_model, 
                           survival_model_vcov, cycles_v)

pre_proc_mortality <- tps_list[[1]]
post_proc_mortality <- tps_list[[2]]

#### 5: creating the transition probabilities ####
# this section generates all the transition probabilities and gathers them into
# a data frame

## probability of having procedure 
# depends on the week we set the proc at above
prob_proc <- ifelse(day_of_proc <= cycles_v, 1, 0)

## probability of pre procedure mortality 
# putting the pre proc vector we generated earlier into a data frame with the 
# probability of having the procedure
time_dependent_probabilities <- data.frame(pre_proc_mortality, prob_proc)

## probability of exit to acute admission

# currently this is constant with time. This equals zero for some conditions.
prob_exit_acute <- prob_exit_acute_df$mean
time_dependent_probabilities$prob_exit_acute <- prob_exit_acute
# once the procedure has happened for all those on the waiting list then 
# probability of acute admission becomes zero
time_dependent_probabilities$prob_exit_acute <- 
  ifelse(time_dependent_probabilities$prob_proc == 1, 0, time_dependent_probabilities$prob_exit_acute)

## probability of exit to private
# currently this is constant with time.
prob_exit_private <- prob_exit_private_df$mean[prob_exit_private_df$imd == imd_quintile]
                           
time_dependent_probabilities$prob_exit_private <- prob_exit_private
# once the procedure has happened for all those on the waiting list then the 
# probability of switching to private care becomes zero
time_dependent_probabilities$prob_exit_private <- 
  ifelse(time_dependent_probabilities$prob_proc == 1, 0, time_dependent_probabilities$prob_exit_private)

## probability of perioperative mortality
# creating a column for periop mortality
time_dependent_probabilities$periop_mortality <- 0
# adding in the periop mortality at the time of procedure
time_dependent_probabilities$periop_mortality[day_of_proc] <- periop_mortality

## probability of post procedure mortality 
# adding the post_proc_mortality into the data frame
time_dependent_probabilities$post_proc_mortality <- post_proc_mortality
# remove any post procedure mortality for 30 days post procedure as the periop
# risk has been dealt with separately
time_dependent_probabilities$post_proc_mortality[1:(day_of_proc+30)] <- 0

#### 6a: creating a transition matrix ####
  
  # start with a three dimensional empty array
  # (number of states by number of states by number of cycles)
  # and we name the dimensions using state_names, set earlier
  transition_matrix <- array(data=0,dim=c(n_states, n_states, cycles),
                             dimnames= list(state_names, state_names, 1:cycles)) 
  
  ### create a loop that creates a time dependent transition matrix for each cycle
  for (i in 1:cycles) {
    
    ## transitions from being on the waiting list
    # staying on the waiting list...
    transition_matrix["on_waiting_list","on_waiting_list",i] <- 
      ifelse(time_dependent_probabilities$prob_proc[i] == 1, 0,
             1 - time_dependent_probabilities$pre_proc_mortality[i] - 
               time_dependent_probabilities$prob_exit_acute[i] - 
               time_dependent_probabilities$prob_exit_private[i])
    
    ## die whilst on the waiting list...
    transition_matrix["on_waiting_list","dead",i] <-  
      ifelse(time_dependent_probabilities$prob_proc[i] == 1, 0,
             time_dependent_probabilities$pre_proc_mortality[i])             
    
    ## exit the waiting list for acute op  ...
    transition_matrix["on_waiting_list","exit_acute",i] <- 
      ifelse(time_dependent_probabilities$prob_proc[i] == 1, 0, time_dependent_probabilities$prob_exit_acute[i])
    
    ## exit the waiting list for private  ...
    transition_matrix["on_waiting_list","exit_private",i] <- 
      ifelse(time_dependent_probabilities$prob_proc[i] == 1, 0, 
             time_dependent_probabilities$prob_exit_private[i])
    
    ## have proc and move to the peri dead state 
    transition_matrix["on_waiting_list","peri_dead",i] <-  
      ifelse(time_dependent_probabilities$prob_proc[i] == 1, 
             time_dependent_probabilities$periop_mortality[i], 0)    
    
    ## have proc and move to the post proc state 
    transition_matrix["on_waiting_list","post_proc",i] <-  
      ifelse(time_dependent_probabilities$prob_proc[i] == 1, 
             1 - time_dependent_probabilities$periop_mortality[i], 0)    
    
    
    ## transitions from exiting for acute care
    transition_matrix["exit_acute","exit_acute",i] <- 1 ## no transitions out of exit
    
    ## transitions from exiting to get private care
    transition_matrix["exit_private","exit_private",i] <- 1 ## no transitions out of exit
    
    ## transitions out of post_proc
    transition_matrix["post_proc","post_proc",i] <- 1 - time_dependent_probabilities$post_proc_mortality[i]
    transition_matrix["post_proc","dead",i] <- time_dependent_probabilities$post_proc_mortality[i]
    
    
    ## no transitions out of peri dead
    transition_matrix["peri_dead","peri_dead",i] <- 1 ## no transitions out of death
    
    
    ## no transitions out of dead
    transition_matrix["dead","dead",i] <- 1 ## no transitions out of death
  }
  
  # check
  transition_matrix[,,41]
  
  
#### 6b: creating a second transition matrix to cap private procedures ####

# We do not want a situation where an unrealistic proportion of the patients 
# go private. No evidence was found to inform the correct proportion so we use a 
# "rule of thumb" of no more than twice the proportion at baseline
  
# first of all the maximum
max_private_patients <- starting_population[1] * (2 * prob_exit_private_df$prob_private[prob_exit_private_df$imd == imd_quintile])

# secondly we need two transition matrices, the one we have created above and one
# where nobody flows off into private. We will switch to this when creating the 
# trace if the max proportion is reached

transition_matrix_2 <- array(data=0,dim=c(n_states, n_states, cycles),
                           dimnames= list(state_names, state_names, 1:cycles)) 

### create a loop that creates a time dependent transition matrix for each cycle
for (i in 1:cycles) {
  
  ## transitions from being on the waiting list
  # staying on the waiting list...
  transition_matrix_2["on_waiting_list","on_waiting_list",i] <- 
    ifelse(time_dependent_probabilities$prob_proc[i] == 1, 0,
           1 - time_dependent_probabilities$pre_proc_mortality[i] - 
             time_dependent_probabilities$prob_exit_acute[i])
  
  ## die whilst on the waiting list...
  transition_matrix_2["on_waiting_list","dead",i] <-  
    ifelse(time_dependent_probabilities$prob_proc[i] == 1, 0,
           time_dependent_probabilities$pre_proc_mortality[i])             
  
  ## exit the waiting list for acute op  ...
  transition_matrix_2["on_waiting_list","exit_acute",i] <- 
    ifelse(time_dependent_probabilities$prob_proc[i] == 1, 0, time_dependent_probabilities$prob_exit_acute[i])
  
  ## exit the waiting list for private  ...
  transition_matrix_2["on_waiting_list","exit_private",i] <- 0
  
  ## have proc and move to the peri dead state 
  transition_matrix_2["on_waiting_list","peri_dead",i] <-  
    ifelse(time_dependent_probabilities$prob_proc[i] == 1, 
           time_dependent_probabilities$periop_mortality[i], 0)    
  
  ## have proc and move to the post proc state 
  transition_matrix_2["on_waiting_list","post_proc",i] <-  
    ifelse(time_dependent_probabilities$prob_proc[i] == 1, 
           1 - time_dependent_probabilities$periop_mortality[i], 0)    
  
  
  ## transitions from exiting for acute care
  transition_matrix_2["exit_acute","exit_acute",i] <- 1 ## no transitions out of exit
  
  ## transitions from exiting to get private care
  transition_matrix_2["exit_private","exit_private",i] <- 1 ## no transitions out of exit
  
  ## transitions out of post_proc
  transition_matrix_2["post_proc","post_proc",i] <- 1 - time_dependent_probabilities$post_proc_mortality[i]
  transition_matrix_2["post_proc","dead",i] <- time_dependent_probabilities$post_proc_mortality[i]
  
  
  ## no transitions out of peri dead
  transition_matrix_2["peri_dead","peri_dead",i] <- 1 ## no transitions out of death
  
  
  ## no transitions out of dead
  transition_matrix_2["dead","dead",i] <- 1 ## no transitions out of death
}

# check
transition_matrix_2[,,40]


# now we have everything we need to set up the trace
  
#### 7: creating a trace ####

# create a trace which is an empty matrix with one column for each state
# and one row for each cycle
trace_waiting_times <- matrix(data = 0, nrow = cycles, ncol = n_states)
colnames(trace_waiting_times) <- state_names

# putting the population in the first row of the trace, all patients
# start at time 0 on the waiting list
trace_waiting_times[1,] <- starting_population%*%transition_matrix[,,1]

# this selects the matrix to multiply by depending on whether the max
# proportion of private patients has been reached
for (i in 2:cycles) {
  
  trace_waiting_times[i,] <- 
    if (trace_waiting_times[i-1,3] < max_private_patients){
      trace_waiting_times[i-1,]%*%transition_matrix[,,i]
    } 
  else {trace_waiting_times[i-1,]%*%transition_matrix_2[,,i]
}
}

rowSums(trace_waiting_times) # check each row sums to the number who started 

# in the model
  
#### 8: summing days of LDs and QALDs ####

# multiplying the daily trace by the daily utilities
# Calling them Quality Adjusted Life Days at this point as we are still in days.

QALDs <- trace_waiting_times * state_utilities

# Adding in a few QALDs for those who died perioperatively
# Using the estimate of days lived less than 30 apply the pre proc QALY,
# this assumes they are a less well group (otherwise we would use post proc QALY) 
periop_QALDS <- peri_days_survived * pre_proc_utility[day_of_proc] * 
  trace_waiting_times[day_of_proc, "peri_dead"]

# Add these QALDs into the post_procedure utility before discounting
QALDs[ceiling(day_of_proc), "post_proc_utility"] <- 
  QALDs[ceiling(day_of_proc), "post_proc_utility"] + periop_QALDS

# calculating the discount rate for each day
discount_factor <- 1/(1 + dr_o)^(cycles_v/365)

# summing all QALDs from across the different states
QALDs$total <- rowSums(QALDs) 

# creating a new data frame with just the QALDs total and renaming
outputs <- QALDs %>% select(total) %>% rename(QALDs = total) 

# adding a column for discounted QALDs
outputs$discounted_QALDs <- outputs$QALDs * discount_factor 

# summing days for those on_waiting_list and post_proc 
# those who exited for acute or private care are counted elsewhere
outputs$LDs <- trace_waiting_times[,"on_waiting_list"] + trace_waiting_times[,"post_proc"] 
outputs$discounted_LDs <- outputs$LDs * discount_factor

# summing the columns into one total
output_markov <- apply(outputs, 2, sum)

# changing from days to years by dividing by 365 and changing the names
output_markov <- output_markov/365
names(output_markov) <- c("QALYs", "Discounted_QALYs","LYs", "Discounted_LYs")


#### 9: predicting periop mortality if exited the wait list ####
 
# If people exit the waiting list (for either private or acute admission) at 
# day 1, 2, 3.... they can expect a payoff based on the three regression 
# models computed using the HES data. Namely the probability of perioperative 
# mortality, number of days survived less than 30, and the survival model for 
# post procedure survival. Their payoff will differ depending on when they 
# exit the list, as their wait days are different which is an independent 
# variable in all three regression models. The payoff also depends on their 
# reason for exit as acute admission may be subject to higher periop mortality.

# prepare the dataframe which has all the input variables needed for the 
# models
payoff_dataframe <- do.call(rbind, replicate(day_of_proc, df1, simplify = FALSE))

# create a vector for days but only up to the day of procedure, as there will be 
# no one waiting after that point
payoff_dataframe$day_of_proc <- 1:(day_of_proc) 
# create a year vector which we will need later on to select the appropriate QALY
payoff_dataframe$year <- round(ceiling(payoff_dataframe$day_of_proc/365))

# change age so you are start at the point of referral
payoff_dataframe$agecont <- payoff_dataframe$agecont - day_of_proc/365
# add on days wait, because the transition prob script will take them off
payoff_dataframe$agecont <- payoff_dataframe$agecont + payoff_dataframe$day_of_proc/365


## estimate the probability of perioperative mortality for each day

# perioperative mortality for exit to acute care is higher than for electives.
# this will be 0 for some procedures
or_periop <- uplift_mortality_emergency_admission$mean

# increasing the perioperative mortality for acute admissions using the OR
payoff_dataframe$periop_mortality_acute <- 
  predict(periop_model, newdata = payoff_dataframe, type = "response") * or_periop

# and for exit to private care it is the same as those who receive the treatment
# from the waiting list
payoff_dataframe$periop_mortality_private <- 
  predict(periop_model, newdata = payoff_dataframe, type = "response")

#### 10: predicting LDs and QALDs if exited wait list and dies perioperatively ####

## estimate the number of days survived if you die perioperatively, this does not
# differ for private or acute
payoff_dataframe$LDs_peri <- 
    predict(days_less_30_model, 
            newdata = payoff_dataframe, type = "response")

# estimate discounted LDs, we are discounting using the rate they left the wait
# list, rather than cycling through the 10 days or whatever they lived
payoff_dataframe$Discounted_LDs_peri <- payoff_dataframe$LDs_peri * 
  discount_factor[payoff_dataframe$days_wait]

# estimate QALDs (pre proc utility is applied as patients who die peri 
# operatively are assumed to be less well)
payoff_dataframe$QALDs_peri <- payoff_dataframe$LDs_peri * 
  pre_proc_utility[payoff_dataframe$days_wait]

# estimate discounted QALDs
payoff_dataframe$Discounted_QALDs_peri <- payoff_dataframe$QALDs_peri * 
  discount_factor[payoff_dataframe$days_wait]

#### 11: predicting LDs and QALDs if exited the wait list and peri survived ####

## estimate LDs and QALDs if patient does not die perioperatively
# for this we need a mini markov model for each day of exit with just post proc 
# and dead in it. If we used the regression model to predict mean survival 
# we would over estimate it, as we know survival is right skewed with a long tail

# creating empty columns in the payoff dataframe
payoff_dataframe$QALDs_post <- 0
payoff_dataframe$Discounted_QALDs_post <- 0
payoff_dataframe$LDs_post <- 0
payoff_dataframe$Discounted_LDs_post <- 0

for(i in 1:nrow(payoff_dataframe)){
  
# first create an empty mini trace with just post_proc and dead in it
mini_trace <- matrix(data = 0, nrow = cycles-i, ncol = 2)
colnames(mini_trace) <- c("post_proc", "dead")

# generate the transition probabilities, which currently differ for each day 
# people exit because wait days is linked to mortality.
tps_list_mini <- mortality_tps(payoff_dataframe[i,], life_tables_uk_imd, survival_model, survival_model_vcov, cycles_v)
post_proc_mortality_mini <- tps_list_mini[[2]] # select the post procedural transition probabilities

# this is going into a mini matrix that starts after exit to proc, so we need 
# to remove the first i observations of the vector. Remember our hes survival 
# model starts from entering the waiting list. But they also need 30 days not
# being subject to the risk as this would count as a peri operative death which 
# has been counted separately
post_proc_mortality_mini <- post_proc_mortality_mini[(i+1):cycles]
post_proc_mortality_mini[1:30] <- 0


# then create an empty mini transition matrix with just post_proc and dead in it
mini_transition_matrix <- array(data = 0, dim = c(2, 2, cycles-i), 
                                dimnames = list(colnames(mini_trace), colnames(mini_trace), 1:(cycles-i)))

# populating the transition matrix
for(j in 1:(cycles-i)){

## transitions out of post_proc
mini_transition_matrix["post_proc","post_proc",j] <- 1 - post_proc_mortality_mini[j]
mini_transition_matrix["post_proc","dead",j] <- post_proc_mortality_mini[j]

## no transitions out of dead
mini_transition_matrix["dead","dead",j] <- 1 ## no transitions out of death

}
    
# check
mini_transition_matrix[,,10]
  
# now we need to put the population in the first row of the trace
starting_population_mini <- c(1, 0)

# enter the first row of the matrix
mini_trace[1,] <- starting_population_mini %*% mini_transition_matrix[,,1]

# run the trace through all the cycles
for (k in 2:(cycles-i)) {
  mini_trace[k,] <- mini_trace[k-1,] %*% mini_transition_matrix[,,k]
}
    
rowSums(mini_trace) # check each row sums to the starting population 
# in the model
    
# now I need to calculate LDs and QALDs and discounted versions of both
outputs_mini <- as.data.frame(mini_trace)

# calculating QALDS, 
# post procedure utility has to be recalculated as the difference depends on the
# age at which you have the procedure which could change if you have to wait for
# a long time

# utility increase, selected by age and sex from the appropriate data frame
utility_imp <- utility_improvement$Difference[utility_improvement$Age == 
  round(payoff_dataframe$agecont[i] + payoff_dataframe$days_wait[i]/365) & utility_improvement$Female == female]

# adding the pre procedure utility and the utility improvement
post_proc_utility <- pre_proc_utility_nodd + utility_imp

# adding the post procedure utility decrement from waiting (we only have data
# on this for hips and knees)
post_proc_utility <- post_proc_utility + day_of_proc*m[1,1] +
                      (1-female)*day_of_proc*m[3,1]

# then don't forget the multiplication factor which was selected based on the age
# the person entered the model and their sex
post_proc_utility <- post_proc_utility * imd_multiplication_factor

# then change the vector length so it only lasts for the duration of the model run
post_proc_utility <- post_proc_utility[(i+1):cycles] 

# multiply the post proc_utility by the mini trace which we have renamed outputs mini 
# correct point
outputs_mini$QALDs <- outputs_mini$post_proc * post_proc_utility
  
# and discounted QALDs, starting the discount factor at the correct point 
outputs_mini$Discounted_QALDs <- outputs_mini$QALDs * discount_factor[(i+1):cycles]

# and LDs
outputs_mini$LDs <- outputs_mini$post_proc
# and discounted LDs
outputs_mini$Discounted_LDs <- outputs_mini$LDs * discount_factor[(i+1):cycles]

# tidy it up, keep everything in days at this point
outputs_mini <- outputs_mini[ ,3:6]
output_markov_mini <- apply(outputs_mini, 2, sum)
  
# creating empty columns in the payoff dataframe
payoff_dataframe$QALDs_post[i] <- output_markov_mini[1]
payoff_dataframe$Discounted_QALDs_post[i] <- output_markov_mini[2]
payoff_dataframe$LDs_post[i] <- output_markov_mini[3]
payoff_dataframe$Discounted_LDs_post[i] <- output_markov_mini[4]
  
}

#### 12: combining and calculating expected payoffs for those who exit the model ####

# create separate dataframes for acute and private

## acute
expected_payoff_acute <- payoff_dataframe %>% select(days_wait, periop_mortality_acute,
                                                     LDs_peri, Discounted_LDs_peri,
                                                     QALDs_peri, Discounted_QALDs_peri,
                                                     LDs_post, Discounted_LDs_post, 
                                                     QALDs_post, 
                                                     Discounted_QALDs_post)
  
# calculate the expected payoffs using the probability of dying peri operatively
expected_payoff_acute$QALDs <- 
  expected_payoff_acute$periop_mortality_acute * expected_payoff_acute$QALDs_peri +
  (1 - expected_payoff_acute$periop_mortality_acute)*expected_payoff_acute$QALDs_post
  
expected_payoff_acute$Discounted_QALDs <- 
  expected_payoff_acute$periop_mortality_acute * expected_payoff_acute$Discounted_QALDs_peri +
  (1 - expected_payoff_acute$periop_mortality_acute)*expected_payoff_acute$Discounted_QALDs_post
  
expected_payoff_acute$LDs <- 
  expected_payoff_acute$periop_mortality_acute * expected_payoff_acute$LDs_peri +
  (1 - expected_payoff_acute$periop_mortality_acute)*expected_payoff_acute$LDs_post

expected_payoff_acute$Discounted_LDs <- 
  expected_payoff_acute$periop_mortality_acute * expected_payoff_acute$Discounted_LDs_peri +
  (1 - expected_payoff_acute$periop_mortality_acute)*expected_payoff_acute$Discounted_LDs_post

# just keep the variables needed
expected_payoff_acute <- expected_payoff_acute %>% select(days_wait, QALDs, Discounted_QALDs, LDs, Discounted_LDs)

## private
expected_payoff_private <- payoff_dataframe %>% select(days_wait, periop_mortality_private,
                                                         LDs_peri, Discounted_LDs_peri,
                                                         QALDs_peri, Discounted_QALDs_peri,
                                                         LDs_post, Discounted_LDs_post, 
                                                         QALDs_post, 
                                                         Discounted_QALDs_post)
  
# calculate the expected payoffs using the probability of dying peri operatively
expected_payoff_private$QALDs <- 
  expected_payoff_private$periop_mortality_private * expected_payoff_private$QALDs_peri +
  (1 - expected_payoff_private$periop_mortality_private)*expected_payoff_private$QALDs_post
  
expected_payoff_private$Discounted_QALDs <- 
  expected_payoff_private$periop_mortality_private * expected_payoff_private$Discounted_QALDs_peri +
  (1 - expected_payoff_private$periop_mortality_private)*expected_payoff_private$Discounted_QALDs_post

expected_payoff_private$LDs <- 
  expected_payoff_private$periop_mortality_private * expected_payoff_private$LDs_peri +
  (1 - expected_payoff_private$periop_mortality_private)*expected_payoff_private$LDs_post

expected_payoff_private$Discounted_LDs <- 
  expected_payoff_private$periop_mortality_private * expected_payoff_private$Discounted_LDs_peri +
  (1 - expected_payoff_private$periop_mortality_private)*expected_payoff_private$Discounted_LDs_post

# just keep the variables needed
expected_payoff_private <- expected_payoff_private %>% select(days_wait, QALDs, Discounted_QALDs, LDs, Discounted_LDs)
  
#### 13: multiplying exit payoffs by numbers who exit ####

# finally working out how many people exit on each day
trace_waiting_times <- as.data.frame(trace_waiting_times)
trace_waiting_times$daily_exits_private <- trace_waiting_times$exit_private -
  dplyr::lag(trace_waiting_times$exit_private,1)
trace_waiting_times$daily_exits_acute <- trace_waiting_times$exit_acute -
  dplyr::lag(trace_waiting_times$exit_acute,1)
  
# replace the first row from NA to the number who left
trace_waiting_times$daily_exits_private[1] <- trace_waiting_times$exit_private[1]
trace_waiting_times$daily_exits_acute[1] <- trace_waiting_times$exit_acute[1]

sum(trace_waiting_times$daily_exits_private)

# then multiply the number who exit by the payoffs
# create the empty matrix to put them in
payoff_acute <- array(data = 0, dim = c(day_of_proc, 4),
                      dimnames = list(1:day_of_proc, c("QALYs_A", "Discounted_QALYs_A",
                                                       "LYs_A", "Discounted_LYs_A")))

payoff_private <- array(data = 0, dim = c(day_of_proc, 4),
                        dimnames = list(1:day_of_proc, c("QALYs_P", "Discounted_QALYs_P",
                                                         "LYs_P", "Discounted_LYs_P")))

# acute
payoff_acute[,1] <- trace_waiting_times$daily_exits_acute[1:day_of_proc] * expected_payoff_acute$QALDs
payoff_acute[,2] <- trace_waiting_times$daily_exits_acute[1:day_of_proc] * expected_payoff_acute$Discounted_QALDs
payoff_acute[,3] <- trace_waiting_times$daily_exits_acute[1:day_of_proc] * expected_payoff_acute$LDs
payoff_acute[,4] <- trace_waiting_times$daily_exits_acute[1:day_of_proc] * expected_payoff_acute$Discounted_LDs
  
# private
payoff_private[,1] <- trace_waiting_times$daily_exits_private[1:day_of_proc] * expected_payoff_private$QALDs
payoff_private[,2] <- trace_waiting_times$daily_exits_private[1:day_of_proc] * expected_payoff_private$Discounted_QALDs
payoff_private[,3] <- trace_waiting_times$daily_exits_private[1:day_of_proc] * expected_payoff_private$LDs
payoff_private[,4] <- trace_waiting_times$daily_exits_private[1:day_of_proc] * expected_payoff_private$Discounted_LDs

# remember at this last stage to divide by 365
payoff_acute <- apply(payoff_acute, 2, sum)
payoff_acute <- payoff_acute/365

# remember at this last stage to divide by 365
payoff_private <- apply(payoff_private, 2, sum)
payoff_private <- payoff_private/365

#### 14: adding utility from the NHS waiting period to private and acute ####

# the pre proc utility whilst waiting for an NHS procedure, before switching
# to private/acute needs to be calculated and added on to the private/acute
# payoff and taken away from the NHS payoffs.

# first create a vector of cumulative utility, and a discounted one

cumulative_utility  <- matrix(data = 0, nrow = day_of_proc, ncol = 4)
cumulative_utility[,1] <- pre_proc_utility[1:day_of_proc]
cumulative_utility[,2] <- pre_proc_utility[1:day_of_proc] * discount_factor[1:day_of_proc]
cumulative_utility[,3] <- 1
cumulative_utility[,4] <- discount_factor[1:day_of_proc]
colnames(cumulative_utility) <- c("utility", "discounted_utility", "LDs", "dLDs")

cumulative_utility <- apply(cumulative_utility, 2, cumsum)

# then multiply this by the numbers who exit on each day

# first for private
pre_proc_private_utility <- cumulative_utility * trace_waiting_times$daily_exits_private[1:day_of_proc]
# then I need to sum it and divide by 365
pre_proc_private_utility <- ifelse( day_of_proc > 1, apply(pre_proc_private_utility, 2, sum), 0)
pre_proc_private_utility <- pre_proc_private_utility/365

# first for private
pre_proc_acute_utility <- cumulative_utility * trace_waiting_times$daily_exits_acute[1:day_of_proc]
# then I need to sum it and divide by 365
pre_proc_acute_utility <- ifelse(day_of_proc > 1, apply(pre_proc_acute_utility, 2, sum),0)
pre_proc_acute_utility <- pre_proc_acute_utility/365

# finally add this on to the private/acute payoffs and take it away from the NHS payoff
# check
#sum(output_markov + payoff_private + payoff_acute)

payoff_private <- payoff_private + pre_proc_private_utility

payoff_acute <- payoff_acute + pre_proc_acute_utility

# and remember to take it off the outcomes for those who receive their proc on
# the NHS

output_markov <- output_markov - pre_proc_private_utility - pre_proc_acute_utility

# check
#sum(output_markov + payoff_private + payoff_acute)

#### 15: adding all the model outcomes together ####

# payoffs for all those who exited to acute admission
payoff_acute

# payoffs for all those who exited to go to private care
payoff_private

# payoff for those who have the operation on the NHS
output_markov

# total payoff
output_1 <- c(output_markov, payoff_acute, payoff_private)

# we also want some clinical outcomes for the report
output_clinical <- trace_waiting_times[day_of_proc,1:6] # select the appropriate row
output_clinical_names <- names(trace_waiting_times[1:6]) # select the names

output_2 <- setNames(as.numeric(output_clinical), output_clinical_names) # we want it as a named number

# these are the prob of dying peri op if you have acute or private admission
prob_peri_acute_private <- as.numeric( c(
  sum(payoff_dataframe$periop_mortality_acute * trace_waiting_times$daily_exits_acute[1:day_of_proc]),
sum(payoff_dataframe$periop_mortality_private * trace_waiting_times$daily_exits_private[1:day_of_proc]) 
))

output_3 <- setNames(as.numeric(prob_peri_acute_private), c("peri_acute", "peri_private")) # we want it as a named number

output <- c(output_1, output_2, output_3) # one group of named numbers

output 

  
}  

#### running the function ####


waiting_times_model(day_of_proc = 7*18, age = 60, female = 0, imd_quintile = 3,
                    charlsonindex = 4.5, mean_wait = 60)

