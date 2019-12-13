##############################################################################
# Code for estimating treatment effect in a real dataset

# Source the code
source('code_sigma_Z.R')

#####################################################################################
estimate_real_data <- function(delta = 0){

  dt = read.csv('generated_data.csv')

  # Scale continuous covariates to mean 0 and SD = 0.5
  dt$bthwght = (dt$bthwght - mean(dt$bthwght))/(2*sd(dt$bthwght))
  dt$gestage_weeks = (dt$gestage_weeks - mean(dt$gestage_weeks))/(2*sd(dt$gestage_weeks))
  dt$precare = (dt$precare - mean(dt$precare))/(2*sd(dt$precare))
  dt$parity = (dt$parity - mean(dt$parity))/(2*sd(dt$parity))
  dt$agem = (dt$agem - mean(dt$agem))/(2*sd(dt$agem))
  dt$medu = (dt$medu - mean(dt$medu))/(2*sd(dt$medu))
  dt$income = (dt$income - mean(dt$income))/(2*sd(dt$income))
  dt$below_poverty = (dt$below_poverty - mean(dt$below_poverty))/(2*sd(dt$below_poverty))
  dt$home_value = (dt$home_value - mean(dt$home_value))/(2*sd(dt$home_value))
  dt$highschool = (dt$highschool - mean(dt$highschool))/(2*sd(dt$highschool))
  dt$college = (dt$college - mean(dt$college))/(2*sd(dt$college))
  dt$rent = (dt$rent - mean(dt$rent))/(2*sd(dt$rent))
  dt$num_anomaly = (dt$num_anomaly - mean(dt$num_anomaly))/(2*sd(dt$num_anomaly))
  dt$diff_travel_time = (dt$diff_travel_time - mean(dt$diff_travel_time))/(2*sd(dt$diff_travel_time))

  dt$losI = dt$losI - delta * dt$diff_travel_time

  dt = dt[order(dt$diff_travel_time),]
  # Estimate E[L|Z]
  dt$mean_bthwght = ksmooth(dt$diff_travel_time, dt$bthwght, n.points = dim(dt)[1], x.points = dt$diff_travel_time)$y
  dt$mean_gestage_weeks = ksmooth(dt$diff_travel_time, dt$gestage_weeks, n.points = dim(dt)[1], x.points = dt$diff_travel_time)$y
  dt$mean_gest_diabetes = ksmooth(dt$diff_travel_time, dt$Gestational_DiabetesM, n.points = dim(dt)[1], x.points = dt$diff_travel_time)$y
  dt$mean_precare = ksmooth(dt$diff_travel_time, dt$precare, n.points = dim(dt)[1], x.points = dt$diff_travel_time)$y
  dt$mean_parity = ksmooth(dt$diff_travel_time, dt$parity, n.points = dim(dt)[1], x.points = dt$diff_travel_time)$y
  dt$mean_agem = ksmooth(dt$diff_travel_time, dt$agem, n.points = dim(dt)[1], x.points = dt$diff_travel_time)$y
  dt$mean_medu = ksmooth(dt$diff_travel_time, dt$medu, n.points = dim(dt)[1], x.points = dt$diff_travel_time)$y
  dt$mean_income = ksmooth(dt$diff_travel_time, dt$income, n.points = dim(dt)[1], x.points = dt$diff_travel_time)$y
  dt$mean_below_poverty = ksmooth(dt$diff_travel_time, dt$below_poverty, n.points = dim(dt)[1], x.points = dt$diff_travel_time)$y
  dt$mean_home_value = ksmooth(dt$diff_travel_time, dt$home_value, n.points = dim(dt)[1], x.points = dt$diff_travel_time)$y
  dt$mean_highschool = ksmooth(dt$diff_travel_time, dt$highschool, n.points = dim(dt)[1], x.points = dt$diff_travel_time)$y
  dt$mean_college = ksmooth(dt$diff_travel_time, dt$college, n.points = dim(dt)[1], x.points = dt$diff_travel_time)$y
  dt$mean_rent = ksmooth(dt$diff_travel_time, dt$rent, n.points = dim(dt)[1], x.points = dt$diff_travel_time)$y
  dt$mean_num_anomaly = ksmooth(dt$diff_travel_time, dt$num_anomaly, n.points = dim(dt)[1], x.points = dt$diff_travel_time)$y
  dt$mean_raceM = ksmooth(dt$diff_travel_time, dt$raceM, n.points = dim(dt)[1], x.points = dt$diff_travel_time)$y
  dt$mean_insuranceM_HMO = ksmooth(dt$diff_travel_time, dt$insuranceM_HMO, n.points = dim(dt)[1], x.points = dt$diff_travel_time)$y
  dt$mean_insuranceM_FS = ksmooth(dt$diff_travel_time, dt$insuranceM_FS, n.points = dim(dt)[1], x.points = dt$diff_travel_time)$y

  # Estimate A and B
  md_A = lm(dt$losI ~ (dt$bthwght - dt$mean_bthwght) + (dt$gestage_weeks - dt$mean_gestage_weeks) + (dt$Gestational_DiabetesM - dt$mean_gest_diabetes) +
              (dt$precare - dt$mean_precare) + (dt$parity - dt$mean_parity) + (dt$agem - dt$mean_agem) + (dt$medu - dt$mean_medu) +
              (dt$income - dt$mean_income) + (dt$below_poverty - dt$mean_below_poverty) + (dt$home_value - dt$mean_home_value) +
              (dt$highschool - dt$mean_highschool) + (dt$college - dt$mean_college) + (dt$rent - dt$mean_rent) + (dt$num_anomaly - dt$mean_num_anomaly) +
              (dt$raceM - dt$mean_raceM) + (dt$insuranceM_HMO - dt$mean_insuranceM_HMO) + (dt$insuranceM_FS - dt$mean_insuranceM_FS))

  md_B = lm(dt$nicu_level_ind ~ (dt$bthwght - dt$mean_bthwght) + (dt$gestage_weeks - dt$mean_gestage_weeks) + (dt$Gestational_DiabetesM - dt$mean_gest_diabetes) +
              (dt$precare - dt$mean_precare) + (dt$parity - dt$mean_parity) + (dt$agem - dt$mean_agem) + (dt$medu - dt$mean_medu) +
              (dt$income - dt$mean_income) + (dt$below_poverty - dt$mean_below_poverty) + (dt$home_value - dt$mean_home_value) +
              (dt$highschool - dt$mean_highschool) + (dt$college - dt$mean_college) + (dt$rent - dt$mean_rent) + (dt$num_anomaly - dt$mean_num_anomaly) +
              (dt$raceM - dt$mean_raceM) + (dt$insuranceM_HMO - dt$mean_insuranceM_HMO) + (dt$insuranceM_FS - dt$mean_insuranceM_FS))


  A_hat_L = predict(md_A, newdata = dt) - coef(md_A)[1]
  B_hat_L = predict(md_B, newdata = dt) - coef(md_B)[1]

  T = dt$nicu_level_ind - B_hat_L
  Y = dt$losI - A_hat_L

  data1 = data.frame(Z = dt$diff_travel_time, T = T, Y = Y)

  start_0 = rev(unname(coef(ivreg(Y~T | .-T+Z, data = data1))))
  temp = estimate_beta_estimate_var(data1, start_0)
  point_est = temp[1]
  se_est = temp[2]

  return(list(point_est = point_est,
         estimate_se = se_est,
         CI = c(temp[1] - 1.96*temp[2], temp[1] + 1.96*temp[2])))
}

# Run
result = estimate_real_data(delta = 0)
