#' Evaluate a tree 
#' 
#' The function evaluate a tree for a single observation
#'
#' @param tree_top a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch, 
#' `cond` (only for internal nodes) the splitting rule, 
#' `value` (only for terminal nodes) the value at the terminal node
#' @param x 
#' A row of a `data.frame` representing an observation on which the tree is evaluated
#' @return a scalar representing the value at the terminal node associated with `x`
#' @export
#'
#' @examples 
#' tree_top <- generate_random_binary_tree(4,5)
#' tree_top <- assign_node_idx(tree_top)
#' XX <- data.frame(x1 = rnorm(100), x2 = sample(c('A', 'B'), 500, replace = TRUE))
#' tree_top <- assign_split_rules(tree_top, XX)
#' tree_top <- assign_term_node_values(tree_top, 0, 10)
#' get_tree_plot(tree_top)
#' cbind(XX[1,], tree_value = g.T(tree_top, XX[1,]))
#' cbind(XX[10,], tree_value = g.T(tree_top, XX[10,]))

g.T <- function(tree_top, x) {
  cond <- tree_top$cond
  
  if (is.null(cond)) {
    return(tree_top$value)
  }
  
  if (eval_cond(cond, x)) {
    return(g.T(tree_top$left, x))
  }
  
  g.T(tree_top$right, x)
}

#' Evaluate a tree on multiple observations
#' 
#' The function evaluates the tree for each row of `X` 
#'
#' @param tree_top a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch, 
#' `cond` (only for internal nodes) the splitting rule, 
#' `value` (only for terminal nodes) the value at the terminal node
#' @param X
#' A `data.frame` representing the observations on which the tree is evaluated
#'
#' @return a numeric vector with length equal to the number of rows of `X`. 
#' Each element represents the value of the tree for the corresponding row of `X`
#' @export
#'
#' @examples
#' tree_top <- generate_random_binary_tree(4,5)
#' tree_top <- assign_node_idx(tree_top)
#' XX <- data.frame(x1 = rnorm(100), x2 = sample(c('A', 'B'), 500, replace = TRUE))
#' tree_top <- assign_split_rules(tree_top, XX)
#' tree_top <- assign_term_node_values(tree_top, 0, 10)
#' get_tree_plot(tree_top)
#' cbind(XX[1:5,], tree_value = get_value_tree(tree_top, XX[1:5,]))
get_value_tree <- function(tree_top, X){
  vapply(1:nrow(X), \(idx) g.T(tree_top, X[idx,]), 0)
}
