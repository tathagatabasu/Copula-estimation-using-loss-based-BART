# numbers of trees with x terminal nodes
catalan_n <- function(x) {
  
  exp(lfactorial(2 * (x - 1)) - lfactorial(x) - lfactorial(x - 1))
}


# function to calculate the number of possible trees with given depth and number of terminal nodes
enumerate_trees <- function(depth, nterm) {
  # Base cases
  if (depth == 0 || depth == 1) {
    return(1)
  }
  
  max_terminal <- 2^depth
  if (nterm == max_terminal) {
    return(1)
  }
  
  if (nterm < (depth + 1) || nterm > max_terminal) {
    stop('Number of terminal nodes outside valid range: (depth + 1) to 2^depth.')
  }
  
  # Candidate splits for left subtree terminal nodes
  min_left <- max(1, nterm - 2^(depth - 1))
  max_left <- min(2^(depth - 1), nterm - 1)
  poss_left_nodes <- seq(min_left, max_left)
  
  # Refine candidate splits based on depth vs terminal nodes
  if (depth > (nterm - depth)) {
    poss_left_nodes <- poss_left_nodes[poss_left_nodes <= (nterm - depth) | poss_left_nodes >= depth]
  }
  
  count <- 0L
  for (nl in poss_left_nodes) {
    nr <- nterm - nl
    
    left_large <- (nl >= depth)
    right_large <- (nr >= depth)
    
    left_count <- if (left_large) enumerate_trees(depth - 1, nl) else catalan_n(nl)
    right_count <- if (right_large) enumerate_trees(depth - 1, nr) else catalan_n(nr)
    
    count <- count + left_count * right_count
  }
  
  return(count)
}
