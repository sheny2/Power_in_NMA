

# Helper function

calculate_pi <- function(odds_ratio, p_A) {
  odds_A <- p_A / (1 - p_A)
  odds_B <- odds_A * odds_ratio
  p_B <- odds_B / (1 + odds_B)
  
  return(p_B)
}

# from prob to log odds
logit <- function(p) {
  return(log(p/(1-p)))
}

# from log odds to prob
inverse_logit <- function(x) {
  odds <- exp(x)  
  probability <- odds / (1 + odds) 
  return(probability)  
}
