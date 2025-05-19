search_l <- function(p, adj, start = 0.01, end = 1000, tol = 0.01, max_run = 100){
  run <- 0
  p_low <- calculate_p(adj, start)
  p_high <- calculate_p(adj, end)
  if (p_low > p + tol){
    print("l not found, try smaller start point.")
    return(NULL)
  } else if (p_high < p - tol){
    print("l not found, try bigger end point.")
    return(NULL)
  } else if (abs(p_low - p) <= tol){
    #print(paste("recommended l = ", start))
    return(start)
  } else if (abs(p_high - p) <= tol){
    #print(paste("recommended l = ", end))
    return(end)
  }
  
  while (findInterval(p, c(p_low + tol,p_high - tol)) == 1L){
  #while ((p_low + tol)< p < (p_high - tol)){
    run <- run + 1
   # print(paste("Run ", run, ": l [", start, ", ", end, "], p [", p_low, ", ", p_high, "]"))
    if (run > max_run){
      print(paste("Exact l not found, closest values are:\nl=", start, ": p=", p_low, "\nl=", end, ": p=", p_high))
      return(NULL)
    }
    mid <- (start + end)/2
    p_mid <- calculate_p(adj, mid)
    if (abs(p_mid - p) <= tol){
      #print(paste("recommended l = ", mid))
      return(mid)
    }
    if (p_mid <= p){
      start <- mid
      p_low <- p_mid
    } else {
      end <- mid
      p_high <- p_mid
    }
  }
}

calculate_p <- function(adj, l){
  adj_exp <- exp(-1*(adj^2)/(2*(l^2)))
  return(mean(rowSums(adj_exp))-1)
}

test_l <- function(adj, list_l){
  for (l in list_l){
    print(paste("l is ",str(l),"Percentage of total expression contributed by neighborhoods:",calculate_p(adj, l)))
  }
}

find_l <- function(p, adj, start=0.5, end=2,sep=0.01, tol=0.01){
  for (l in seq(start, end, sep)){
    q <- calculate_p(adj, l)
    print(paste("L=", str(l), "P=", str(round(q,5))))
    if (abs(p-q)<=tol){
      return(l)
    }
  }
  print("l not found, try bigger range or smaller sep!")
}


