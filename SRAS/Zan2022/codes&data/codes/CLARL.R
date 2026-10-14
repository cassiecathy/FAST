library(Rcpp)
sourceCpp("ARL.cpp")
library(openxlsx)
set.seed(0)

# here we list an example for Table 3

#parameter setting
IC_ARL <- 370
lambda <- 0.1
n <- 1000 # running time periods
q <- 5000 # experiment times
m0 <- 200 # number of historical observations
shift <- 0 # mean shift (0 for IC ARL)

p <- 10 # total number of data streams
r_list <- c(8, 6, 4) # list of different numbers of observable data streams 

#determine IC mean and covariance matrix of X
dis <- "norm" # determine distribution type
h <- 1.1
mu <- rep(0,p) # mean vector
cov <- diag(1, nrow = p, ncol = p) # covariance vector
# cov <- matrix(0, nr = p, nc = p) # covariance vector
for (i in 1:p){
  for (j in 1:p){
    cov[i,j] = 0.5 # Table 3: exchangeable
    # cov[i,j] = 0.5 ^ abs(i - j) # Table 4: unexchangeable
  }
}


# bootstrap Function of obtaining control limit according to different delta
CL <- function(Desired, lambda, r, delta){
  L_min <- 8
  L_max <- 16
  arl <- ARL(q, m0, n, h, p, r, lambda, L_max, dis, mu, cov, shift, delta)[1]
  while(arl < Desired){
    L_min <- L_max
    L_max <- L_min * 2
    arl <- ARL(q, m0, n, h, p, r, lambda, L_max, dis, mu, cov, shift, delta)[1]
  }
  
  L_mid <- (L_min + L_max) / 2
  arl <- ARL(q, m0, n, h, p, r, lambda, L_mid, dis, mu, cov, shift, delta)[1]
  
  while(abs(arl - Desired) > 0.5 & (L_max - L_min) > 0.05){
    if (arl < Desired){
      L_min <- L_mid
    }else{
      L_max <- L_mid
    }
    L_mid <- (L_min + L_max) / 2
    arl <- ARL(q, m0, n, h, p, r, lambda, L_mid, dis, mu, cov, shift, delta)[1]
  }
  return(L_mid)
}


for (k in 1:length(r_list)){

  r <- r_list[k]
  
  ## below is for Table 1&2 to select best delta:
  # Set of delta
  # Delta_list <- c(0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3)
  # Delta_list <- seq(0.02, 0.3, 0.02)
  # 
  # CL_list <- rep(0, length(Delta_list))
  # 
  # 
  # for (l in 1:length(Delta_list)){
  #   delta <- Delta_list[l]
  #   CL_list[l] <- CL(IC_ARL, lambda, r, delta)
  #   print(delta)
  #   print(CL_list[l])
  # }
  # 
  # df_in <- data.frame(Delta_list, CL_list)
  # colnames(df_in) <- c("Delta", "CL")
  # write.xlsx(df_in, sprintf("normcor_CL_p=%s_r=%s.xlsx",  p, r), sheetName = "lambda=0.1", col.names = TRUE, row.names = FALSE, append = FALSE)
  
  
  # shift_mean_List <- c(-3, -2, -1, 1, 2, 3)
  shift_mean_List <- c(1, 2, 3)
  
  ARL1 <- matrix(nrow = 1, ncol = length(shift_mean_List))
  
  for (i in 1:nrow(df_in)){
    # delta <- df_in[i, "Delta"]
    # print(delta)
    # L <- df_in[i, "CL"]
    delta <- 0.1
    L <- CL(IC_ARL, lambda, r, delta)
    
    for (l in 1:length(shift_mean_List)){
      shift <- shift_mean_List[l]
      RL_Vec <- ARL(q, m0, n, h, p, r, lambda, L, dis, mu, cov, shift, delta)
      ARL1[i, l] <- sprintf("%s(%s)", RL_Vec[1], RL_Vec[2])
      # print(shift_mean_List[l])
    }
  }
  
  # df_out <- data.frame(df_in["Delta"], ARL1)
  # colnames(df_out) <- c("Delta", "-3", "-2", "-1", "1", "2", "3")
  df_out <- data.frame(ARL1)
  colnames(df_out) <- c("1", "2", "3")
  
  write.xlsx(df_out, sprintf("normcor_ARL1_p=%s_r=%s_delta=%s.xlsx",  p, r, delta), sheetName = "singleshift", col.names = TRUE, row.names = FALSE, append = FALSE)
  
  
}

