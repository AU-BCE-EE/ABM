# Concept for stoich()
# S. Hafner

# Substrate state variables
y <- c(NDF = 10, CP = 5, CF = 2, VSd = 0)

# Coef matrix, set in parameters, or calculated before lsoda() call e.g., if composition of VSd varies
cc <- matrix(rep(1, 8), nrow = 4, dimnames = list(names(y), c('CO2', 'VFA')))

# Hydrolysis rates, all 0.1 here
alpha <- y 
alpha[] <- 0.1

# Derivatives 
# Substrate state variables
dydt <- alpha * y # + influent

# Products state variables
dcdt <- colSums(cc * dydt)

dydt
dcdt
