estimate_R0 <- function(beta, kappa, mu, rho, gamma_1, gamma_2, alpha, delta, q){
  
  var1 <- (beta * kappa) / (kappa + mu)
  
  var6 <- gamma_1 + alpha + mu
  var7 <- gamma_2 + delta + mu
  var5 <- alpha / (var6 * var7)
  var3 <- ((gamma_1 + alpha + mu) ** -1) + var5
  var4 <- q / (gamma_1 + mu)
  
  var2 <- rho * var3 + (1 - rho) * var4 
  
  var1 * var2
}