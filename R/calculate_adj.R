euclid_dist <- function(t1, t2) {
  sum <- 0
  for (i in 1:length(t1)) {
    sum <- sum + (t1[i] - t2[i])^2
  }
  return(sqrt(sum))
}

pairwise_distance <- function(X) {
  n <- dim(X)[1]
  adj <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      adj[i, j] <- euclid_dist(X[i, ], X[j, ])
    }
  }
  return(adj)
}

extract_color <- function(x_pixel, y_pixel, image, beta = 49) {
  # beta to control the range of neighbourhood when calculate grey vale for one spot
  beta_half <- round(beta / 2)
  g <- vector()
  for (i in seq_along(x_pixel)) {
    max_x <- dim(image)[1]
    max_y <- dim(image)[2]
    nbs <- image[max(1, x_pixel[i] - beta_half):min(max_x, x_pixel[i] + beta_half + 1), 
                 max(1, y_pixel[i] - beta_half):min(max_y, y_pixel[i] + beta_half + 1), ]
    g[i] <- mean(mean(nbs, na.rm = TRUE), na.rm = TRUE)
  }
  c0 <- c1 <- c2 <- vector()
  for (i in g) {
    c0 <- c(c0, i[1])
    c1 <- c(c1, i[2])
    c2 <- c(c2, i[3])
  }
  c3 <- (c0 * var(c0) + c1 * var(c1) + c2 * var(c2)) / (var(c0) + var(c1) + var(c2))
  return(c3)
}


calculate_adj_matrix <- function(x, y, x_pixel=NULL, y_pixel=NULL, image=NULL, beta=49, alpha=1, histology=TRUE) {
  # x, y, x_pixel, y_pixel are vectors
  if (histology) {
    stopifnot(!is.null(x_pixel) & !is.null(y_pixel) & !is.null(image))
    stopifnot(length(x) == length(x_pixel) & length(y) == length(y_pixel))
   # print("Calculating adj matrix using histology image...")
    # beta to control the range of neighbourhood when calculating gray value for one spot
    # alpha to control the color scale
    beta_half <- round(beta/2)
    g <- list()
    for (i in seq_along(x_pixel)) {
      max_x <- dim(image)[1]
      max_y <- dim(image)[2]
      nbs <- image[max(1, x_pixel[i]-beta_half):min(max_x, x_pixel[i]+beta_half), max(1, y_pixel[i]-beta_half):min(max_y, y_pixel[i]+beta_half), ]
      nbs_rgb = col2rgb(nbs)
      
      g[[i]] <- rowMeans(nbs_rgb)
    }
    c0 <- numeric()
    c1 <- numeric()
    c2 <- numeric()
    for (i in g) {
      c0 <- c(c0, i[1])
      c1 <- c(c1, i[2])
      c2 <- c(c2, i[3])
    }
    c3 <- (c0*var(c0)+c1*var(c1)+c2*var(c2))/(var(c0)+var(c1)+var(c2))
    c4 <- (c3-mean(c3))/sd(c3)
    z_scale <- max(sd(x), sd(y))*alpha
    z <- c4*z_scale
    z <- as.numeric(z)
    X <- cbind(x, y, z)
  } else {
   # print("Calculating adj matrix using xy only...")
    X <- cbind(x, y)
  }
  return(pairwise_distance(X))
}


