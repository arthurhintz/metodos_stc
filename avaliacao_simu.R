library(dplyr)

result <- read.table("simulation_results.txt", header = TRUE)

resu <- result |> 
  group_by(n, mu_true, method) |>
  summarise(
    media = mean(mu_est),                     
    bies = mean(mu_est) - mu_true,            
    rb = ((mean(mu_est) - mu_true) / mu_true) * 100,  
    vari = var(mu_est),                       
    eqm = (mean(mu_est) - mu_true)^2 + var(mu_est), 
    mae = mean(abs(mu_est - mu_true)),    
    cv = sqrt(var(mu_est)) / mean(mu_est) * 100,
    mape = mean(abs((mu_est - mu_true) / mu_true)) * 100,
    .groups = "drop"  
  ) |> 
  distinct() |> 
  mutate(across(c(media, bies, rb, vari, eqm, mae, cv, mape), round, 4))


write.table(
  resu,
  file = "summary_results.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

