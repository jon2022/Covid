#######
#covid simulation model JAM April 2020

#open required libraries
library(data.table)
library(tidyverse)
library(gridExtra)

#to prevent timezone error messages
options(tz = "Europe/London")
Sys.setenv(TZ="Europe/London")

'%ni%' <- Negate('%in%')  # define 'not in' func
options(scipen=999)  # prevents printing scientific notations.


#function to take median (test_med) and nth quantile (default = 0.95) (test_quant) and return sd parameter for lognormal distribution
sdlog <- function(test_med, test_quant, nth = 0.95)
{result<-((sqrt(2)*(qnorm(nth))/sqrt(2))^-1)*log(test_quant/test_med)
return(result)}

#function to calculate mean of lognormal distribution from meanlog and sdlog
test_mean <- function(x,y) exp(x+(y^2/2))

#function to make new cases
make_cases <- function(incm, inc95, infm, inf95, ncases = 1, asy = 0.179, mort = 0.0066){
  incubation <- round(rlnorm(ncases, log(incm), sdlog(incm, inc95)),0)
  infective <- round(rlnorm(ncases, log(infm), sdlog(infm, inf95)),0)
  asymptomatic <- runif(ncases) < asy
  died <- runif(ncases) < mort
  result <- data.table(incubation, infective, asymptomatic, died)
  return(result)
}
#function to simulate population and calculate SEIR results
simulate_population <- function(pop_size = 100000,   #population size
                                inc_med = 4, inc_95 = 7, #incubation period, median and 95th centile
                                inf_med = 5, inf_95 = 10,  #infectious period, median and 95th centile
                                asymp = 0.179,  #proportion asymptomatic
                                ifr = 0.0066,  #IFR
                                r0 = 2.5,  #R0
                                r0_psup = 1.5,  #R0 with partial suppression (e.g. key workers)
                                r0_sup = 0.5,  #R0_suppressed
                                r0_group = 3,  #R0 within groups
                                gr_size = 4,  #set group size
                                prop_sup = 0.9,  #proportion of population fully suppressing
                                prop_psup = 0.1,  #proportion of partially suppressing
                                last_day = 90,  # how many days to simulate
                                seed_numb = 1000  #starting population (seeding)
                                )
{
  #generate populations - all, unsupressed, partially supressed
  if((prop_psup+prop_sup)>1){
    stop("sum of proportion supressed and partially suppressed must be <= 1")
  }
  all_pop <- 1:pop_size
  psup_pop <- sample(all_pop,pop_size*prop_psup)
  unsup_pop <- sample(all_pop[all_pop %ni% psup_pop],(1-(prop_sup+prop_psup))*pop_size)
  sup_pop <- all_pop[all_pop %ni% c(psup_pop,unsup_pop)]
  
  #mean infectious period
  inf_mean <- test_mean(log(inf_med),sdlog(inf_med, inf_95))
  
  #daily infection transmission
  r0_daily <- r0/inf_mean
  
  #daily infection transmission - suppressed group
  r0_daily_psup <- r0_psup/inf_mean
  
  #daily infection transmission - suppressed group
  r0_daily_sup <- r0_sup/inf_mean
  
  #daily infection transmission - within group
  r0_daily_group <- r0_group/inf_mean
  
  #vector of case id's for new cases
  seed_cases <- sample(1:pop_size,seed_numb)
  
  #make data.table of case numbers and initial cases with start day as exponentially distributed prior to day 0 and result of make cases
  cases <- make_cases(inc_med,inc_95,inf_med,inf_95,seed_numb)
  days_into = rexp(seed_numb, 0.1) 
  incubating = days_into < cases$incubation 
  infectious = (days_into >= cases$incubation) & (days_into < (cases$incubation + cases$infective)) 
  rec_or_dead = days_into >= (cases$incubation + cases$infective) 
  id = seed_cases
  cases <- data.table(days_into, incubating, infectious, rec_or_dead, id, cases)
  
  #set up data.table for results of daily simulation
  out_table <- data.table(day_number = 0L,
                          new = sum(cases$days_into==0),
                          cases %>%
                            summarise(n_inc_sup = sum(incubating[id %in% sup_pop]),
                                      n_inc_psup = sum(incubating[id %in% psup_pop]),
                                      n_inc_unsup = sum(incubating[id %in% unsup_pop]),
                                      n_inf_sup = sum(infectious[id %in% sup_pop]),
                                      n_inf_psup = sum(infectious[id %in% psup_pop]),
                                      n_inf_unsup = sum(infectious[id %in% unsup_pop]),
                                      dead_sup = sum(died[id %in% sup_pop] & rec_or_dead[id %in% sup_pop]),
                                      dead_psup = sum(died[id %in% psup_pop] & rec_or_dead[id %in% psup_pop]),
                                      dead_unsup = sum(died[id %in% unsup_pop] & rec_or_dead[id %in% unsup_pop]),
                                      recovered_sup = sum(!died[id %in% sup_pop] & rec_or_dead[id %in% sup_pop]),
                                      recovered_psup = sum(!died[id %in% psup_pop] & rec_or_dead[id %in% psup_pop]),
                                      recovered_unsup = sum(!died[id %in% unsup_pop] & rec_or_dead[id %in% unsup_pop])))
  out_table$currentR0 <- r0
  



  # loop to populate table of outcomes
  for(i in 1:last_day){
    
    # number of new cases infected - in suppressed population
    new_cases_sup <- round(r0_daily_sup * (out_table$n_inf_sup[nrow(out_table)]+out_table$n_inf_unsup[nrow(out_table)]))
    # number of new cases infected - unsupressed
    new_cases_unsup <- round((r0_daily - r0_daily_sup) * out_table$n_inf_unsup[nrow(out_table)])
    
    #identify groups in which there are infectious cases and find population in these groups
    inf_groups <- round(cases$id/gr_size)
    group_pop <- all_pop[round(all_pop/gr_size) %in% inf_groups]
    
    #new cases within groups
    new_cases_grp <- round(r0_daily_group * (out_table$n_inf_sup[nrow(out_table)]+out_table$n_inf_unsup[nrow(out_table)]))
    
    #sample whole population to find additional newly infected cases from suppressed - replace = TRUE to allow that the same person may potentially contact more than one infected person
    added_cases <- sample(all_pop,new_cases_sup, replace = T)
    
    #add sample from unsupressed population
    added_cases <- append(added_cases,sample(unsup_pop,new_cases_unsup, replace = T))
    
    #add sample within groups
    added_cases <- append(added_cases,sample(group_pop,new_cases_grp, replace = T))
    
    #identify groups in which there are infectious cases
    inf_groups <- round(cases$id/gr_size)
    
    #exclude duplicates
    added_cases <- unique(added_cases)
    
    #find cases that are not already infected, dead or recovered
    added_cases <- added_cases[added_cases %ni% cases$id]
    
    #add the new cases with current days since infected (days-into) of 0 and random distribution of incubation period, infective period, mortality and asymptomatic defined by make-cases function
    cases <- rbind(cases,data.table(days_into = rep(0,length(added_cases)), incubating = rep(T,length(added_cases)), infectious = rep(F,length(added_cases)), rec_or_dead = rep(F,length(added_cases)), id = added_cases, make_cases(inc_med,inc_95,inf_med,inf_95,length(added_cases))))
    
    #update days into and current state and summarise result for both suppressed and unsuppressed cases.
    cases <- cases %>% 
      mutate(days_into = days_into+1,
             incubating = days_into < incubation,
             infectious = (days_into >= incubation) & (days_into < (incubation + infective)),
             rec_or_dead = days_into >= (incubation + infective))
    
    
    new_data <- data.table(day_number = i, 
                           new = length(added_cases),
                           cases %>%
                             summarise(n_inc_sup = sum(incubating[id %in% sup_pop]),
                                       n_inc_psup = sum(incubating[id %in% psup_pop]),
                                       n_inc_unsup = sum(incubating[id %in% unsup_pop]),
                                       n_inf_sup = sum(infectious[id %in% sup_pop]),
                                       n_inf_psup = sum(infectious[id %in% psup_pop]),
                                       n_inf_unsup = sum(infectious[id %in% unsup_pop]),
                                       dead_sup = sum(died[id %in% sup_pop] & rec_or_dead[id %in% sup_pop]),
                                       dead_psup = sum(died[id %in% psup_pop] & rec_or_dead[id %in% psup_pop]),
                                       dead_unsup = sum(died[id %in% unsup_pop] & rec_or_dead[id %in% unsup_pop]),
                                       recovered_sup = sum(!died[id %in% sup_pop] & rec_or_dead[id %in% sup_pop]),
                                       recovered_psup = sum(!died[id %in% psup_pop] & rec_or_dead[id %in% psup_pop]),
                                       recovered_unsup = sum(!died[id %in% unsup_pop] & rec_or_dead[id %in% unsup_pop])))
    new_data$currentR0 <- new_data$new*inf_mean/(new_data$n_inf_unsup+new_data$n_inf_psup+new_data$n_inf_sup)
    out_table <- rbind(out_table, new_data)
    }
  return(out_table)
}

