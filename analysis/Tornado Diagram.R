#### RHTA LMICs - Workshop September ####

### Tornado Diagram ###

devtools::source_url("https://github.com/mbounthavong/Decision_Analysis/blob/master/tornado_diagram_code.R?raw=TRUE")


### 1. Create a new function for the base case Markov model ----

markov <- function(params){
  with(
    as.list(params), 
    {
      
      cycle_length <- 1 # cycle length equal to one year (use 1/12 for monthly)
      n_age_init <- 25  # age at baseline
      n_age_max  <- 100 # maximum age of follow up
      n_cycles <- (n_age_max - n_age_init)/cycle_length # time horizon, number of cycles
      #* Age labels 
      v_age_names <- paste(rep(n_age_init:(n_age_max-1), each = 1/cycle_length), 
                           1:(1/cycle_length), 
                           sep = ".")
      #* the 4 health states of the model:
      v_names_states <- c("H",  # Healthy (H)
                          "S1", # Sick (S1)
                          "S2", # Sicker (S2)
                          "D")  # Dead (D)
      
      n_states <- length(v_names_states)   # number of health states 
      
      ### Discounting factors ----
      d_c <- 0.03 # annual discount rate for costs 
      d_e <- 0.03 # annual discount rate for QALYs
      
      ### Strategies ----
      v_names_str <- c("Standard of care",      # store the strategy names
                       "Strategy A", 
                       "Strategy B",
                       "Strategy AB") 
      n_str       <- length(v_names_str)        # number of strategies
      
      ## Within-cycle correction (WCC) using Simpson's 1/3 rule ----
      v_wcc <- gen_wcc(n_cycles = n_cycles,  # Function included in "R/Functions.R". The latest version can be found in `darthtools` package
                       method = "Simpson1/3") # vector of wcc
      
      ### Transition rates (annual), and hazard ratios (HRs) ----
      r_HS1  <- 0.15  # constant annual rate of becoming Sick when Healthy
      r_S1H  <- 0.5   # constant annual rate of becoming Healthy when Sick
      r_S1S2 <- 0.105 # constant annual rate of becoming Sicker when Sick
      hr_S1  # hazard ratio of death in Sick vs Healthy 
      hr_S2  <- 10    # hazard ratio of death in Sicker vs Healthy 
      
      ### Effectiveness of treatment B ----
      hr_S1S2_trtB  # hazard ratio of becoming Sicker when Sick under treatment B
      
      ## Age-dependent mortality rates ----
      lt_usa_2015 <- read.csv("data/LifeTable_USA_Mx_2015.csv")
      #* Extract age-specific all-cause mortality for ages in model time horizon
      v_r_mort_by_age <- lt_usa_2015 %>% 
        dplyr::filter(Age >= n_age_init & Age < n_age_max) %>%
        dplyr::select(Total) %>%
        as.matrix()
      
      ### State rewards ----
      #### Costs ----
      c_H    <- 2000  # annual cost of being Healthy
      c_S1   <- 4000  # annual cost of being Sick
      c_S2   <- 15000 # annual cost of being Sicker
      c_D    <- 0     # annual cost of being dead
      c_trtA # annual cost of receiving treatment A
      c_trtB # annual cost of receiving treatment B
      #### Utilities ----
      u_H    <- 1     # annual utility of being Healthy
      u_S1   <- 0.75  # annual utility of being Sick
      u_S2   <- 0.5   # annual utility of being Sicker
      u_D    <- 0     # annual utility of being dead
      u_trtA # annual utility when receiving treatment A
      
      ### Transition rewards ----
      du_HS1 <- 0.01  # disutility when transitioning from Healthy to Sick
      ic_HS1 <- 1000  # increase in cost when transitioning from Healthy to Sick
      ic_D   <- 2000  # increase in cost when dying
      
      ### Discount weight for costs and effects ----
      v_dwc  <- 1 / ((1 + (d_e * cycle_length)) ^ (0:n_cycles))
      v_dwe  <- 1 / ((1 + (d_c * cycle_length)) ^ (0:n_cycles))
      
      # Process model inputs ----
      ## Age-specific transition rates to the Dead state for all cycles ----
      v_r_HDage  <- rep(v_r_mort_by_age, each = 1/cycle_length)
      #* Name age-specific mortality vector 
      names(v_r_HDage) <- v_age_names
      
      #* compute mortality rates
      v_r_S1Dage <- v_r_HDage * hr_S1 # Age-specific mortality rate in the Sick state 
      v_r_S2Dage <- v_r_HDage * hr_S2 # Age-specific mortality rate in the Sicker state 
      #* transform rates to probabilities adjusting by cycle length
      #* Function included in "R/Functions.R". The latest version can be found in `darthtools` package
      p_HS1  <- rate_to_prob(r = r_HS1, t = cycle_length) # constant annual probability of becoming Sick when Healthy conditional on surviving 
      p_S1H  <- rate_to_prob(r = r_S1H, t = cycle_length) # constant annual probability of becoming Healthy when Sick conditional on surviving
      p_S1S2 <- rate_to_prob(r = r_S1S2, t = cycle_length)# constant annual probability of becoming Sicker when Sick conditional on surviving
      v_p_HDage  <- rate_to_prob(v_r_HDage, t = cycle_length)  # Age-specific mortality risk in the Healthy state 
      v_p_S1Dage <- rate_to_prob(v_r_S1Dage, t = cycle_length) # Age-specific mortality risk in the Sick state
      v_p_S2Dage <- rate_to_prob(v_r_S2Dage, t = cycle_length) # Age-specific mortality risk in the Sicker state
      
      ## Annual transition probability of becoming Sicker when Sick for treatment B ----
      #* Apply hazard ratio to rate to obtain transition rate of becoming Sicker when 
      #* Sick for treatment B
      r_S1S2_trtB <- r_S1S2 * hr_S1S2_trtB
      #* Transform rate to probability to become Sicker when Sick under treatment B 
      #* adjusting by cycle length conditional on surviving
      #* (Function included in "R/Functions.R". The latest version can be found in 
      #* `darthtools` package)
      p_S1S2_trtB <- rate_to_prob(r = r_S1S2_trtB, t = cycle_length)
      
      # Construct state-transition models ----
      ## Initial state vector ----
      #* All starting healthy
      v_m_init <- c(H = 1, S1 = 0, S2 = 0, D = 0) # initial state vector
      
      ## Initialize cohort traces ----
      ### Initialize cohort trace under SoC ----
      m_M_SoC <- matrix(NA, 
                        nrow = (n_cycles + 1), ncol = n_states, 
                        dimnames = list(0:n_cycles, v_names_states))
      #* Store the initial state vector in the first row of the cohort trace
      m_M_SoC[1, ] <- v_m_init
      
      ### Initialize cohort trace for strategies A, B, and AB ----
      #* Structure and initial states are the same as for SoC
      m_M_strA  <- m_M_SoC # Strategy A
      m_M_strB  <- m_M_SoC # Strategy B
      m_M_strAB <- m_M_SoC # Strategy AB
      
      ## Create transition probability arrays for strategy SoC ----
      ### Initialize transition probability array for strategy SoC ----
      #* All transitions to a non-death state are assumed to be conditional on survival
      a_P_SoC <- array(0,
                       dim  = c(n_states, n_states, n_cycles),
                       dimnames = list(v_names_states, 
                                       v_names_states, 
                                       0:(n_cycles - 1)))
      ### Fill in array
      ## From H
      a_P_SoC["H", "H", ]   <- (1 - v_p_HDage) * (1 - p_HS1)
      a_P_SoC["H", "S1", ]  <- (1 - v_p_HDage) * p_HS1
      a_P_SoC["H", "D", ]   <- v_p_HDage
      ## From S1
      a_P_SoC["S1", "H", ]  <- (1 - v_p_S1Dage) * p_S1H
      a_P_SoC["S1", "S1", ] <- (1 - v_p_S1Dage) * (1 - (p_S1H + p_S1S2))
      a_P_SoC["S1", "S2", ] <- (1 - v_p_S1Dage) * p_S1S2
      a_P_SoC["S1", "D", ]  <- v_p_S1Dage
      ## From S2
      a_P_SoC["S2", "S2", ] <- 1 - v_p_S2Dage
      a_P_SoC["S2", "D", ]  <- v_p_S2Dage
      ## From D
      a_P_SoC["D", "D", ]   <- 1
      
      ### Initialize transition probability array for strategy A as a copy of SoC's ----
      a_P_strA <- a_P_SoC
      
      ### Initialize transition probability array for strategy B ----
      a_P_strB <- a_P_SoC
      #* Update only transition probabilities from S1 involving p_S1S2
      a_P_strB["S1", "S1", ] <- (1 - v_p_S1Dage) * (1 - (p_S1H + p_S1S2_trtB))
      a_P_strB["S1", "S2", ] <- (1 - v_p_S1Dage) * p_S1S2_trtB
      
      ### Initialize transition probability array for strategy AB as a copy of B's ----
      a_P_strAB <- a_P_strB
      
      ## Check if transition probability arrays are valid ----
      #* Functions included in "R/Functions.R". The latest version can be found in `darthtools` package
      ### Check that transition probabilities are [0, 1] ----
      check_transition_probability(a_P_SoC, verbose = TRUE)
      check_transition_probability(a_P_strA, verbose = TRUE)
      check_transition_probability(a_P_strB, verbose = TRUE)
      check_transition_probability(a_P_strAB, verbose = TRUE)
      ### Check that all rows for each slice of the array sum to 1 ----
      check_sum_of_transition_array(a_P_SoC, n_states = n_states, n_cycles = n_cycles, verbose = TRUE)
      check_sum_of_transition_array(a_P_strA, n_states = n_states, n_cycles = n_cycles, verbose = TRUE)
      check_sum_of_transition_array(a_P_strB, n_states = n_states, n_cycles = n_cycles, verbose = TRUE)
      check_sum_of_transition_array(a_P_strAB, n_states = n_states, n_cycles = n_cycles, verbose = TRUE)
      
      ## Create transition dynamics arrays ----
      #* These arrays will capture transitions from each state to another over time 
      ### Initialize transition dynamics array for strategy SoC ----
      a_A_SoC <- array(0,
                       dim      = c(n_states, n_states, n_cycles + 1),
                       dimnames = list(v_names_states, v_names_states, 0:n_cycles))
      #* Set first slice of a_A_SoC with the initial state vector in its diagonal
      diag(a_A_SoC[, , 1]) <- v_m_init
      ### Initialize transition-dynamics array for strategies A, B, and AB ----
      #* Structure and initial states are the same as for SoC
      a_A_strA  <- a_A_SoC
      a_A_strB  <- a_A_SoC
      a_A_strAB <- a_A_SoC
      
      #  Run Markov model ----
      #* Iterative solution of age-dependent cSTM
      for(t in 1:n_cycles){
        ## Fill in cohort trace
        # For SoC
        m_M_SoC[t + 1, ]  <- m_M_SoC[t, ]  %*% a_P_SoC[, , t]
        # For strategy A
        m_M_strA[t + 1, ] <- m_M_strA[t, ] %*% a_P_strA[, , t]
        # For strategy B 
        m_M_strB[t + 1, ] <- m_M_strB[t, ] %*% a_P_strB[, , t]
        # For strategy ZB 
        m_M_strAB[t + 1, ] <- m_M_strAB[t, ] %*% a_P_strAB[, , t]
        
        ## Fill in transition-dynamics array
        # For SoC
        a_A_SoC[, , t + 1]  <- diag(m_M_SoC[t, ]) %*% a_P_SoC[, , t]
        # For strategy A
        a_A_strA[, , t + 1] <- diag(m_M_strA[t, ]) %*% a_P_strA[, , t]
        # For strategy B
        a_A_strB[, , t + 1] <- diag(m_M_strB[t, ]) %*% a_P_strB[, , t]
        # For strategy AB
        a_A_strAB[, , t + 1] <- diag(m_M_strAB[t, ]) %*% a_P_strAB[, , t]
      }
      
      ## Store the cohort traces in a list ----
      l_m_M <- list(SoC =  m_M_SoC,
                    A   =  m_M_strA,
                    B   =  m_M_strB,
                    AB  =  m_M_strAB)
      names(l_m_M) <- v_names_str
      
      ## Store the transition dynamics array for each strategy in a list ----
      l_a_A <- list(SoC =  a_A_SoC,
                    A   =  a_A_strA,
                    B   =  a_A_strB,
                    AB  =  a_A_strAB)
      names(l_a_A) <- v_names_str
      
      # Plot Outputs ----
      #* (Functions included in "R/Functions.R"; depends on the `ggplot2` package)
      
      ## Plot the cohort trace for strategy SoC ----
      plot_trace(m_M_SoC)
      ## Plot the cohort trace for all strategies ----
      plot_trace_strategy(l_m_M)
      
      ## Plot the epidemiology outcomes ----
      ### Survival ----
      survival_plot <- plot_surv(l_m_M, v_names_death_states = "D") +
        theme(legend.position = "bottom")
      survival_plot
      ### Prevalence ----
      prevalence_S1_plot   <- plot_prevalence(l_m_M, 
                                              v_names_sick_states = c("S1"), 
                                              v_names_dead_states = "D")  +
        theme(legend.position = "")
      prevalence_S1_plot
      prevalence_S2_plot   <- plot_prevalence(l_m_M, 
                                              v_names_sick_states = c("S2"), 
                                              v_names_dead_states = "D")  +
        theme(legend.position = "")
      prevalence_S2_plot
      prevalence_S1S2_plot <- plot_prevalence(l_m_M, 
                                              v_names_sick_states = c("S1", "S2"), 
                                              v_names_dead_states = "D") +
        theme(legend.position = "")
      prevalence_S1S2_plot
      prop_sicker_plot     <- plot_proportion_sicker(l_m_M, 
                                                     v_names_sick_states = c("S1", "S2"), 
                                                     v_names_sicker_states = c("S2")) +
        theme(legend.position = "bottom")
      prop_sicker_plot
      
      ## Combine plots ----
      gridExtra::grid.arrange(survival_plot, 
                              prevalence_S1_plot, 
                              prevalence_S2_plot, 
                              prevalence_S1S2_plot, 
                              prop_sicker_plot, 
                              ncol = 1, heights = c(0.75, 0.75, 0.75, 0.75, 1))
      
      # State Rewards ----
      ## Scale by the cycle length ----
      #* Vector of state utilities under strategy SoC
      v_u_SoC    <- c(H  = u_H, 
                      S1 = u_S1, 
                      S2 = u_S2, 
                      D  = u_D) * cycle_length
      #* Vector of state costs under strategy SoC
      v_c_SoC    <- c(H  = c_H, 
                      S1 = c_S1,
                      S2 = c_S2, 
                      D  = c_D) * cycle_length
      #* Vector of state utilities under strategy A
      v_u_strA   <- c(H  = u_H, 
                      S1 = u_trtA, 
                      S2 = u_S2, 
                      D  = u_D) * cycle_length
      #* Vector of state costs under strategy A
      v_c_strA   <- c(H  = c_H, 
                      S1 = c_S1 + c_trtA,
                      S2 = c_S2 + c_trtA, 
                      D  = c_D) * cycle_length
      #* Vector of state utilities under strategy B
      v_u_strB   <- c(H  = u_H, 
                      S1 = u_S1, 
                      S2 = u_S2, 
                      D  = u_D) * cycle_length
      #* Vector of state costs under strategy B
      v_c_strB   <- c(H  = c_H, 
                      S1 = c_S1 + c_trtB, 
                      S2 = c_S2 + c_trtB, 
                      D  = c_D) * cycle_length
      #* Vector of state utilities under strategy AB
      v_u_strAB  <- c(H  = u_H, 
                      S1 = u_trtA, 
                      S2 = u_S2, 
                      D  = u_D) * cycle_length
      #* Vector of state costs under strategy AB
      v_c_strAB  <- c(H  = c_H, 
                      S1 = c_S1 + (c_trtA + c_trtB), 
                      S2 = c_S2 + (c_trtA + c_trtB), 
                      D  = c_D) * cycle_length
      
      ## Store state rewards ----
      #* Store the vectors of state utilities for each strategy in a list 
      l_u <- list(SoC = v_u_SoC,
                  A   = v_u_strA,
                  B   = v_u_strB,
                  AB  = v_u_strAB)
      #* Store the vectors of state cost for each strategy in a list 
      l_c <- list(SoC =  v_c_SoC,
                  A   =  v_c_strA,
                  B   =  v_c_strB,
                  AB  =  v_c_strAB)
      
      #* assign strategy names to matching items in the lists
      names(l_u) <- names(l_c) <- v_names_str
      
      # Compute expected outcomes ----
      #* Create empty vectors to store total utilities and costs 
      v_tot_qaly <- v_tot_cost <- vector(mode = "numeric", length = n_str)
      names(v_tot_qaly) <- names(v_tot_cost) <- v_names_str
      
      ## Loop through each strategy and calculate total utilities and costs ----
      for (i in 1:n_str) { # i <- 1
        v_u_str <- l_u[[i]]   # select the vector of state utilities for the i-th strategy
        v_c_str <- l_c[[i]]   # select the vector of state costs for the i-th strategy
        a_A_str <- l_a_A[[i]] # select the transition array for the i-th strategy
        
        ##* Array of state rewards 
        #* Create transition matrices of state utilities and state costs for the i-th strategy 
        m_u_str   <- matrix(v_u_str, nrow = n_states, ncol = n_states, byrow = T)
        m_c_str   <- matrix(v_c_str, nrow = n_states, ncol = n_states, byrow = T)
        #* Expand the transition matrix of state utilities across cycles to form a transition array of state utilities
        a_R_u_str <- array(m_u_str, 
                           dim      = c(n_states, n_states, n_cycles + 1),
                           dimnames = list(v_names_states, v_names_states, 0:n_cycles))
        # Expand the transition matrix of state costs across cycles to form a transition array of state costs
        a_R_c_str <- array(m_c_str, 
                           dim      = c(n_states, n_states, n_cycles + 1),
                           dimnames = list(v_names_states, v_names_states, 0:n_cycles))
        
        ##* Apply transition rewards
        #* Apply disutility due to transition from H to S1
        a_R_u_str["H", "S1", ]      <- a_R_u_str["H", "S1", ]       - du_HS1
        #* Add transition cost per cycle due to transition from H to S1
        a_R_c_str["H", "S1", ]      <- a_R_c_str["H", "S1", ]       + ic_HS1
        #* Add transition cost  per cycle of dying from all non-dead states
        a_R_c_str[-n_states, "D", ] <- a_R_c_str[-n_states, "D", ] + ic_D
        
        ###* Expected QALYs and costs for all transitions per cycle
        #* QALYs = life years x QoL
        #* Note: all parameters are annual in our example. In case your own case example is different make sure you correctly apply.
        a_Y_c_str <- a_A_str * a_R_c_str
        a_Y_u_str <- a_A_str * a_R_u_str 
        
        ###* Expected QALYs and costs per cycle
        ##* Vector of QALYs and costs
        v_qaly_str <- apply(a_Y_u_str, 3, sum) # sum the proportion of the cohort across transitions 
        v_cost_str <- apply(a_Y_c_str, 3, sum) # sum the proportion of the cohort across transitions
        
        ##* Discounted total expected QALYs and Costs per strategy and apply within-cycle correction if applicable
        #* QALYs
        v_tot_qaly[i] <- t(v_qaly_str) %*% (v_dwe * v_wcc)
        #* Costs
        v_tot_cost[i] <- t(v_cost_str) %*% (v_dwc * v_wcc)
      }
      
      # Cost-effectiveness analysis (CEA) ----
      ## Incremental cost-effectiveness ratios (ICERs) ----
      #* Function included in "R/Functions.R"; depends on the `dplyr` package
      #* The latest version can be found in `dampack` package
      df_cea <- calculate_icers(cost       = v_tot_cost, 
                                effect     = v_tot_qaly,
                                strategies = v_names_str)
      df_cea
      a<-df_cea[3,2]
      b<-df_cea[1,2]
      c<-df_cea[3,3]
      d<-df_cea[1,3]
  return((a-b)/(c-d))
    }
)
      }



### 2. Define the parameters to be plotted on the Tornado diagram ----

inputs<- data.frame( 
  hr_S1  <- 3,     
  hr_S1S2_trtB <- 0.6,
  c_trtA <- 12000,
  c_trtB<- 13000,
  u_trtA <- 0.95)



### 3. Test the function ----

markov(inputs)


### 4. Define possible variation ranges for each parameter considered. ----

hr_S1_range  <- c(BaseCase = hr_S1,    low = 2,  high = 4)      
hr_S1S2_trtB_range <- c(BaseCase = hr_S1S2_trtB,    low = 0.3,  high = 0.9)  
c_trtA_range <- c(BaseCase = c_trtA,    low = 8000,  high = 18000) 
c_trtB_range<- c(BaseCase = c_trtB, low = 8000, high = 18000)
u_trtA_range <- c(BaseCase = u_trtA,    low = 0.85, high = 1) 


paramNames <-   c( "hr_S1",     
                   "hr_S1S2_trtB",
                   "c_trtA",
                   "c_trtB",
                    "u_trtA")

l.tor.in <- vector("list", 5)
names(l.tor.in) <- paramNames



### 5. Run the model n times according to the number of parameters included. ----
l.tor.in$hr_S1  <- cbind(hr_S1  = hr_S1_range, inputs[-1])
l.tor.in$hr_S1S2_trtB  <- cbind(hr_S1S2_trtB  = hr_S1S2_trtB_range, inputs[-2])
l.tor.in$c_trtA  <- cbind(c_trtA  = c_trtA_range, inputs[-3])
l.tor.in$c_trtB  <- cbind(c_trtB  = c_trtB_range, inputs[-4])
l.tor.in$u_trtA  <- cbind(u_trtA  = u_trtA_range, inputs[-5])

l.tor.out <- vector("list", 5)
names(l.tor.out) <- paramNames


n_time_init_psa_parallel <- Sys.time()
for(i in 1:5){
  l.tor.out[[i]] <- t(apply(l.tor.in[[i]], 1, markov)) 
}
n_time_end_psa_parallel <- Sys.time()
time<-n_time_end_psa_parallel-n_time_init_psa_parallel
time

### 6. Save the ICERs obtained. ----
m.tor <- matrix(unlist(l.tor.out), nrow = 5, ncol = 3, byrow = TRUE, 
                dimnames = list(paramNames, c("basecase", "low", "high")))

### 7. Run the Tornado Diagram. ----
TornadoPlot(main_title = "Tornado Plot", Parms = paramNames, Outcomes = m.tor, 
            outcomeName = "Incremental Cost-Effectiveness Ratio (ICER)", 
            xlab = "ICER", 
            ylab = "Parameters", 
            col1="#3182bd", col2="#6baed6")
