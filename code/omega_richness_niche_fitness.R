library(mvtnorm)
library(tidyverse)
library(patchwork)

# Showing the effect of fitness differences (competitive asymmetry) and 
# niche differences on the size of the feasible domain Omega 
# for three different levels of species richness


# based on https://doi.org/10.5061/dryad.v9f5s but returning Omega not log(Omega)
feasibility_domain = function(A) {
  
  n = nrow(A)
  Sigma = solve(t(A) %*% A)
  p = pmvnorm(mean = rep(0,n), sigma = Sigma, lower = rep(0,n), upper = rep(Inf,n))
  as.numeric(p) * 2^n
  
}

make_matrix_fitness = function(n, base, asymmetry) {
  
  g = outer(1:n, 1:n, FUN = function(i,j) (j-i)/(n-1))
  A = base + asymmetry * g
  diag(A) = 1
  A
  
}
make_matrix_fitness(3, .5, .4)
#      [,1] [,2] [,3]
# [1,]  1.0  0.7  0.9
# [2,]  0.3  1.0  0.7
# [3,]  0.1  0.3  1.0

comb = expand.grid(richness = 2:4,
                   asym = seq(0, .5, length.out = 51))

d1 = vector("list", nrow(comb))

for (i in 1:nrow(comb)) {
  
  n = comb$richness[i]
  j = comb$asym[i]
  
  A = make_matrix_fitness(.5, j, n)
  
  d1[[i]] = data.frame(competitive_asymmetry = j,
                       Omega = feasibility_domain(A),
                       richness = n)
}

p1 = bind_rows(d1) |> 
  ggplot(aes(competitive_asymmetry, Omega, colour = factor(richness))) +
  geom_line()


make_matrix_niche = function(n, niche_diff) {
  
  A = matrix(1 - niche_diff, nrow = n, ncol = n)
  diag(A) = 1
  A
}
make_matrix_niche(3, .5)
#      [,1] [,2] [,3]
# [1,]  1.0  0.5  0.5
# [2,]  0.5  1.0  0.5
# [3,]  0.5  0.5  1.0

comb = expand.grid(richness = 2:4,
                   niche_difference = seq(.01, .99, length.out = 51))
d2 = vector("list", nrow(comb))

for (i in seq_len(nrow(comb))) {
  
  n = comb$richness[i]
  nd = comb$niche_difference[i]
  
  A = make_matrix_niche(n, nd)
  
  d2[[i]] = data.frame(niche_difference = nd,
                       Omega = feasibility_domain(A),
                       richness = n)
}

p2 = bind_rows(d2) |> 
  ggplot(aes(niche_difference, Omega, colour = factor(richness))) +
  geom_line()

p1 / p2
