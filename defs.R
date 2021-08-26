emission_prob <- function(genotype, state_tuple, source_a, source_d, a_prime=0.1, d_prime=0.1) {
  # Homozygous state
  if(state_tuple[1] == state_tuple[2]) {
    k <- state_tuple[1]
    a <- source_a[k]
    d <- source_d[k]
    numerator <- beta(genotype + d + d_prime, 2 - genotype + a + a_prime)
    denominator <- beta(d + d_prime, a + a_prime)
    binomial_coefficient <- if(genotype == 1) 2
                            else if(genotype == 0 || genotype == 2) 1
                            else stop(paste0("unknown genotype value ", genotype))
    return(binomial_coefficient * numerator / denominator)
  }
  # Heterozygous state
  else {
    k1 <- state_tuple[1]
    k2 <- state_tuple[2]
    a1 <- source_a[k1]
    a2 <- source_a[k2]
    d1 <- source_d[k1]
    d2 <- source_d[k2]
    denominator <- (d1 + d_prime + a1 + a_prime) * (d2 + d_prime + a2 + a_prime)
    if(genotype == 0) {
      numerator <- (a1 + a_prime) * (a2 + a_prime)
    }
    else if(genotype == 1) {
      numerator <- (a1 + a_prime) * (d2 + d_prime) + (a2 + a_prime) * (d1 + d_prime)
    }
    else if(genotype == 2) {
      numerator <- (d1 + d_prime) * (d2 + d_prime)
    } else {
      stop(paste0("unkown genotype value ", genotype))
    }
    return(numerator / denominator)
  }
}

equilibrium_probs <- function(prop2) {
  prop1 = 1.0 - prop2
  c(prop1 ^ 2, 2 * prop1 * prop2, prop2 ^ 2)
}

transition_probs <- function(prop2, length2) {
  prop1 = 1.0 - prop2
  length1 = prop1 * length2 / prop2
  t12 = 2.0 / length1
  t13 = 1.0 / length1 ^ 2
  t11 = 1.0 - t12 - t13
  t21 = 1.0 / length2
  t23 = 1.0 / length1
  t22 = 1.0 - t21 - t23
  t31 = 1.0 / length2 ^ 2
  t32 = 2.0 / length2
  t33 = 1.0 - t31 - t32
  matrix(c(t11, t12, t13, t21, t22, t23, t31, t32, t33), nrow = 3, byrow = TRUE)
}


# Expect data to be a data-frame with the following column-names:
# chrom
# pos
# genotype
# source_a (list column)
# source_d (list column)

# Currently this is hard-coded for 2 sources.
run_forward <- function(data, transition_matrix, a_prime, d_prime) {
  state_tuples <- list(c(1, 1), c(1, 2), c(2, 2))
  L <- nrow(data)
  forward_vec <- matrix(0.0, ncol=L, nrow=3)
  scaling_factors <- rep(1.0, L)
  for(i in 1:L) {
    g <- data$genotype[i]
    a <- data$source_a[[i]]
    d <- data$source_d[[i]]
    emission_vec <- purrr::map_dbl(state_tuples, function(s) {emission_prob(g, s, a, d)})
    previous_forward_vec <- if(i == 1) rep(1.0, 3) else forward_vec[,i-1]
    new_forward_vec <- (t(transition_matrix) %*% previous_forward_vec) * emission_vec
    scaling_factors[i] <- sum(new_forward_vec)
    forward_vec[,i] <- new_forward_vec / scaling_factors[i]
  }
  return(list(forward_vec, scaling_factors))
}

run_backward <- function(data, scaling_factors, equilibrium_probs, transition_matrix, a_prime, d_prime) {
  state_tuples <- list(c(1, 1), c(1, 2), c(2, 2))
  L <- nrow(data)
  backward_vec <- matrix(0.0, ncol=L, nrow=3)
  for(i in L:1) {
    g <- data$genotype[i]
    a <- data$source_a[[i]]
    d <- data$source_d[[i]]
    emission_vec <- purrr::map_dbl(state_tuples, function(s) {emission_prob(g, s, a, d)})
    previous_backward_vec <- if(i == L) equilibrium_probs else backward_vec[,i+1]
    new_backward_vec <- transition_matrix %*% (previous_backward_vec * emission_vec)
    backward_vec[,i] <- new_backward_vec / scaling_factors[i]
  }
  return(backward_vec)
}

