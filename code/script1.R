library(tidyverse)
library(mvtnorm)
library(deSolve)

# 100 species i.e. parameter sets for thermal performance curves of gross 
# photosynthetic rate and respiration rate 
# these were generate from a model fitted on data in 10.1111/ele.13469
TPCs = read_csv("data/TPCpairs100.csv") |> 
  rename_with(~ .x |> 
                str_replace("r_", "x_") |> 
                str_replace("gp_", "r_"))

# the thermal performance curve 
schoolfield_high = function(c, Ea, Eh, Th, temp) {
  Tc = 293.15
  k = 8.62e-5
  
  boltzmann_term = 24 * c * exp(Ea / k * (1 / Tc - 1 / temp))
  inactivation_term = 1 + exp(Eh / k * (1 / Th - 1 / temp))
  
  boltzmann_term / inactivation_term
}

# to generate matrices with a desired niche difference
competition_matrix = function(n, niche_diff) {
  A = matrix(1 - niche_diff, ncol = n, nrow = n)
  diag(A) = 1
  A
}
competition_matrix(3, .5)

# calculate equilibrium biomasses
equilibrium = function(community, temp, niche_diff) {
  
  n = nrow(community)
  r = numeric(n)
  x = numeric(n)
  
  for(i in seq_len(n)) {
    
    sp = community[i,]
    r[i] = schoolfield_high(sp$r_c, sp$r_Ea, sp$r_Eh, sp$r_Th, temp)
    x[i] = schoolfield_high(sp$x_c, sp$x_Ea, sp$x_Eh, sp$x_Th, temp)
    
  }
  
  K = 1 * (1 - x / r)
  A = competition_matrix(n, niche_diff)
  Bstar = solve(A, K)
  
  list(Bstar = Bstar, A = A)
  
}

# is the equilibrium feasible?
feasible = function(Bstar) {
  all(Bstar > 0)
}

# based on https://doi.org/10.5061/dryad.v9f5s but returning Omega not log(Omega)
feasibility_domain = function(A) {
  
  n = nrow(A)
  Sigma = solve(t(A) %*% A)
  p = pmvnorm(mean = rep(0,n), sigma = Sigma, lower = rep(0,n), upper = rep(Inf,n))
  as.numeric(p) * 2^n
  
}


Tmu = 12 + 273.15                                
target = 100
richness_levels = 2:4
niche_differences = seq(.1, .9, length.out = 5)

set.seed(2026)

feasible_communities = list()

for(n in richness_levels) {
  
  for(nd in niche_differences) {
    
    # avoid duplicates
    checked = list()
    count = 0
    
    while(count < target) {
      
      sp = sort(sample(nrow(TPCs), n))
      
      if(any(map_lgl(checked, ~ identical(.x, sp)))) {
        message("duplicate")
        next
      }
      
      checked = append(checked, list(sp))
      
      community = TPCs[sp,]
      
      eq = equilibrium(community, Tmu, nd)
      
      if(feasible(eq$Bstar)) {
        
        Omega = feasibility_domain(eq$A)
        
        feasible_communities[[length(feasible_communities) + 1]] =
          list(species = sp,
               richness = n,
               niche_diff = nd,
               Bstar = eq$Bstar,
               A = eq$A,
               Omega = Omega)
        
        count = count + 1
      }
    }
    
    message("richness = ", n, ", 
            niche difference = ", nd, ": found ", count, " communities.")
  }
}
feasible_communities[[111]]


pd = do.call(rbind, lapply(feasible_communities, function(x) {
  data.frame(richness = x$richness,
             niche_difference = x$niche_diff,
             Omega = x$Omega)
  })) |> 
  mutate(community = row_number(),
         .before = "richness")

#check
ggplot(pd, aes(x = niche_difference, y = Omega, colour = factor(richness))) +
  geom_point()


dynamics = function(t, B, pars) {
  
  with(pars, {
    
    T = Tmean + Tamp * sin(2 * pi * t / period)
    
    r = schoolfield_high(TPCs$r_c, TPCs$r_Ea, TPCs$r_Eh, TPCs$r_Th, T)
    x = schoolfield_high(TPCs$x_c, TPCs$x_Ea, TPCs$x_Eh, TPCs$x_Th, T)
    
    B[B < ext] = 0
    dB = B * (r * (1 - (A %*% B) / K0) - x)
    
    list(dB)
  })
}
# extinction handling
root = function(t, B, pars) {
  with(pars, { B - ext })
}
event = function(t, B, pars) {
  with(pars, { extinct = B <= ext
               B[extinct] = 0
               return(B) })
}
simulate_community = function(community, 
                              Tmean = 12 + 273.15, 
                              Tamp = 7, 
                              period = 365,
                              ext = 1e-2) {
  pars = list(A = community$A,
              K0 = 1,
              TPCs = TPCs[community$species, ],
              Tmean = Tmean,
              Tamp = Tamp,
              period = period,
              ext = ext)
  
  B0 = rep(0.1, length(community$species))
  names(B0) = paste0("sp_", community$species)
  times = seq(0, 365 * 3, by = 1)
  
  out = ode(y = B0, times = times, func = dynamics, parms = pars,
            rootfun = root, 
            events = list(func = event, root = TRUE))
  
  out = as.data.frame(out) |>
    pivot_longer(cols = -time,
                 names_to = "species",
                 values_to = "biomass") |> 
    filter(biomass > 0)
  
  out
}


all_dynamics = lapply(feasible_communities, 
                      simulate_community) |>
  bind_rows(.id = "community") |>
  mutate(community = as.integer(community))

# a couple of examples
all_dynamics |> filter(community == 645) |> 
  ggplot(aes(x = time, y = biomass, 
             colour = species)) +
  geom_line()
all_dynamics |> filter(community == 741) |> 
  ggplot(aes(x = time, y = biomass, 
             colour = species)) +
  geom_line()

# which communities did not lose species
intact = all_dynamics |> 
  filter(time %in% c(0, 3*365)) |> 
  group_by(community, time) |>
  summarise(N = n(), .groups = "drop") |> 
  pivot_wider(names_from = time, values_from = N) |> 
  filter(`0` == `1095`) |> 
  pull(community) 
  

all_intact = all_dynamics |> filter(community %in% intact)

temporal = all_intact |> 
  filter(time > 365) |> 
  group_by(community, time) |> 
  summarise(N = n(),
            Btot = sum(biomass), .groups = "drop") |> 
  group_by(community) |>
  summarise(N = first(N),
            meanBtot = mean(Btot),
            cvBtot = sd(Btot)/meanBtot, 
            .groups = "drop")

# include niche and Omega
stuctural_temporal = temporal |>
  left_join(pd |> select(community, 
                         niche_difference, 
                         Omega), 
            by = "community")

stuctural_temporal |> 
  ggplot(aes(x = niche_difference, 
             y = cvBtot, colour = as.factor(N))) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(~N)

stuctural_temporal |> 
  ggplot(aes(x = Omega, 
             y = cvBtot, colour = as.factor(N))) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(~N)
