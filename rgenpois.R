rgenpois <-
  function(n, mu = stop("no mu arg"), phi = stop("no phi arg"), omega = 0)
  {
    # check if parameters are valid
    if(omega < 0) {return("omega has to be in [0,1]!")}
    if(omega > 1) {return("omega has to be in [0,1]!")}
    
    # inversion method
    x <- double(n)
    u <- runif(n, 0, 1)
    upper <- max(u)
    s <- double(1000)
    #P(X=0)
    p <- omega + (1-omega) * exp(-mu/phi)
    s[1] <- p
    if (upper > 0) {
      rekursive <- FALSE
      i <- 1
      while (s[i] < upper) {
        #P(X=x)
        if (rekursive==FALSE) {
          p <- (1-omega)*mu*(mu+(phi-1)*i)^(i-1)/exp(lgamma(i+1))*
            phi^(-i)*exp(-1/phi*(mu+(phi-1)*i))}
        if (p==Inf) { 
          rekursive <- TRUE 
          log.p.alt <- log( (1-omega)*mu*(mu+(phi-1)*(i-1))^(i-2)/exp(lgamma(i-1+1))*
                              phi^(-(i-1))*exp(-1/phi*(mu+(phi-1)*(i-1))))
        }
        if (rekursive==TRUE) {
          log.p <- log( (mu+(i-1)*(phi-1))/(phi*i)*
                          (1+(phi-1)/(mu+(i-1)*(phi-1)))^(i-1)*
                          exp(1/phi-1) ) + log.p.alt
          log.p.alt <- log.p
          p <- exp(log.p)
        }
        if (ceiling(i/1000)==floor(i/1000)) {
          temp <- double(1000)
          s <- c(s,temp)
        }
        s[i+1] <- s[i] + p
        i <- i+1
      }
    }
    for (j in 1:length(u)) {
      i <- 1
      while (u[j] > s[i]) {   
        i <- i+1
      }
      x[j] <- i-1
    }
    return(x)
  }