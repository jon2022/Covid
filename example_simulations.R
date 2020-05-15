###############################
#a few comparisons of parameters

# various proportions suppressed

ptm <- proc.time() #time process
sim_result1 <- simulate_population(prop_psup = 0.05, prop_sup = 0.95, seed_numb = 3000, r0_psup = 3, inc_med = 3, inc_95 = 6, inf_med = 3, inf_95 = 8)
sim_result1$dead <- sim_result1$dead_sup + sim_result1$dead_unsup + sim_result1$dead_unsup
sim_result2 <- simulate_population(prop_psup = 0.1, prop_sup = 0.9, seed_numb = 3000, r0_psup = 3, inc_med = 3, inc_95 = 6, inf_med = 3, inf_95 = 8)
sim_result2$dead <- sim_result2$dead_sup + sim_result2$dead_unsup + sim_result2$dead_unsup
sim_result3 <- simulate_population(prop_psup = 0.15, prop_sup = 0.85, seed_numb = 3000, r0_psup = 3, inc_med = 3, inc_95 = 6, inf_med = 3, inf_95 = 8)
sim_result3$dead <- sim_result3$dead_sup + sim_result3$dead_unsup + sim_result3$dead_unsup
sim_result4 <- simulate_population(prop_psup = 0.2, prop_sup = 0.8, seed_numb = 3000, r0_psup = 3, inc_med = 3, inc_95 = 6, inf_med = 3, inf_95 = 8)
sim_result4$dead <- sim_result4$dead_sup + sim_result4$dead_unsup + sim_result4$dead_unsup
proc.time() - ptm

p1<-ggplot() + 
  geom_point(mapping = aes(x=sim_result1$day_number, y=sim_result1$currentR0, colour = "95% supressed"))+
  geom_point(mapping = aes(x=sim_result2$day_number, y=sim_result2$currentR0, colour = "90% supressed"))+
  geom_point(mapping = aes(x=sim_result3$day_number, y=sim_result3$currentR0, colour = "85% supressed"))+
  geom_point(mapping = aes(x=sim_result4$day_number, y=sim_result4$currentR0, colour = "80% supressed"))+
  xlab("days from seeding")+
  ylab("current R0")+
  ggtitle(paste0("Apparent R0"))+
  theme(legend.position = c(0.7,0.7))+
  labs(colour = "Parameters")


p2<-ggplot() + 
  geom_point(mapping = aes(x=sim_result1$day_number, y=sim_result1$dead, colour = "1"))+
  geom_point(mapping = aes(x=sim_result2$day_number, y=sim_result2$dead, colour = "2"))+
  geom_point(mapping = aes(x=sim_result3$day_number, y=sim_result3$dead, colour = "3"))+
  geom_point(mapping = aes(x=sim_result4$day_number, y=sim_result4$dead, colour = "4"))+
  xlab("days from seeding")+
  ylab("Deaths")+
  ggtitle(paste0("Cumulative deaths"))+
  theme(legend.position = "none")

p3<-ggplot() + 
  geom_point(mapping = aes(x=sim_result1$day_number, y=sim_result1$new, colour = "1"))+
  geom_point(mapping = aes(x=sim_result2$day_number, y=sim_result2$new, colour = "2"))+
  geom_point(mapping = aes(x=sim_result3$day_number, y=sim_result3$new, colour = "3"))+
  geom_point(mapping = aes(x=sim_result4$day_number, y=sim_result4$new, colour = "4"))+
  xlab("days from seeding")+
  ylab("Cases")+
  ggtitle(paste0("Rate of new cases"))+
  theme(legend.position = "none")

p4<-ggplot() + 
  geom_point(mapping = aes(x=sim_result1$day_number, y=sim_result1$n_inf_sup+sim_result1$n_inf_psup+sim_result1$n_inf_unsup, colour = "1"))+
  geom_point(mapping = aes(x=sim_result2$day_number, y=sim_result2$n_inf_sup+sim_result2$n_inf_psup+sim_result2$n_inf_unsup, colour = "2"))+
  geom_point(mapping = aes(x=sim_result3$day_number, y=sim_result3$n_inf_sup+sim_result3$n_inf_psup+sim_result3$n_inf_unsup, colour = "3"))+
  geom_point(mapping = aes(x=sim_result4$day_number, y=sim_result4$n_inf_sup+sim_result4$n_inf_psup+sim_result4$n_inf_unsup, colour = "4"))+
  xlab("days from seeding")+
  ylab("Infectious cases")+
  ggtitle(paste0("Number of infectious cases"))+
  theme(legend.position = "none")

grid.arrange(p1,p2,p3,p4,ncol = 2, nrow = 2)
