get_ANOVA_f_stat <- function(v1, v2){
  ## input two vectors v1, v2
  ## output the f-statistics of ANOVA test
  v = c(v1, v2)
  ss1 = sum((v - mean(v))^2)
  ss2 = (sum((v1 - mean(v1))^2) + sum((v2 - mean(v2))^2))
  tmp = mean(c(v1,v2))
  ss3 = (mean(v1) - tmp)^2*length(v1) + (mean(v2) -tmp)^2*length(v2)
  
  ss1 == ss2 + ss3
  a = (ss2/(length(v)-2))
  b = (ss3/(2-1))
  f = b/a
  return(f)
}

