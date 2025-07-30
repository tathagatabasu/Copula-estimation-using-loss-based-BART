#' @title Chipman Prior Tree Generation and Analysis Functions
#' @description Suite of functions for generating and analyzing trees with the Chipman prior,
#'   as used in Bayesian Additive Regression Trees (BART).
#' @author Your Name

#' @param n Number of samples or nodes.
#' @param alpha, beta Parameters for the Chipman prior.
#' @param tree_top The root node of the tree.
#' @param node.idx Index of the node for depth calculation.
#' @return Various outputs depending on function.
#' @examples
#' generate.tree.chip(alpha. = 0.95, beta. = 2)
#' chipman_prior_tree(tree, alpha. = 0.95, beta. = 2)

#' Generate a tree using the Chipman prior
#' @param alpha. Numeric, root-splitting probability (0 < alpha. <= 1).
#' @param beta. Numeric, depth penalty (beta. >= 0).
#' @return A list representing the tree.
generate.tree.chip <- function(alpha., beta.) {
  stopifnot(is.numeric(alpha.), alpha. > 0, alpha. <= 1,
            is.numeric(beta.), beta. >= 0)
  tree <- list(node.idx = 1, depth = 0, left = NULL, right = NULL)
  nodes <- list(tree)
  n.split <- 1
  max_depth <- 100
  split_count <- 0
  
  while (n.split > 0 && split_count < max_depth) {
    split_count <- split_count + 1
    to_split <- nodes[sapply(nodes, function(node) is.null(node$left) && is.null(node$right))]
    for (node in to_split) {
      p_split <- alpha. / (1 + node$depth)^beta.
      if (runif(1) < p_split) {
        # Split node
        node$left <- list(node.idx = length(nodes) + 1, depth = node$depth + 1, left = NULL, right = NULL)
        node$right <- list(node.idx = length(nodes) + 2, depth = node$depth + 1, left = NULL, right = NULL)
        nodes <- c(nodes, list(node$left, node$right))
        n.split <- n.split + 2
      }
    }
    n.split <- length(nodes[sapply(nodes, function(node) is.null(node$left) && is.null(node$right))])
  }
  if (split_count == max_depth) warning("Max tree depth reached; stopping splits.")
  return(tree)
}

#' Calculate the depth of a node in a tree
#' @param tree_top The root node.
#' @param node.idx The index of the node to find.
#' @return Integer depth or NA if not found.
get_depth_at_node <- function(tree_top, node.idx) {
  if (tree_top$node.idx == node.idx) return(tree_top$depth)
  if (!is.null(tree_top$left)) {
    result <- get_depth_at_node(tree_top$left, node.idx)
    if (!is.na(result)) return(result)
  }
  if (!is.null(tree_top$right)) {
    result <- get_depth_at_node(tree_top$right, node.idx)
    if (!is.na(result)) return(result)
  }
  return(NA)
}

#' Calculate the prior probability of a tree under the Chipman prior
#' Note: Requires helper functions get_terminal_nodes_idx and get_internal_nodes_idx
#' @param tree The tree object.
#' @param alpha. Numeric, root-splitting probability.
#' @param beta. Numeric, depth penalty.
#' @return Numeric prior probability.
chipman_prior_tree <- function(tree, alpha., beta.) {
  stopifnot(is.numeric(alpha.), alpha. > 0, alpha. <= 1,
            is.numeric(beta.), beta. >= 0)
  terminal_nodes <- get_terminal_nodes_idx(tree)
  internal_nodes <- get_internal_nodes_idx(tree)
  all_nodes <- c(terminal_nodes, internal_nodes)
  depth.node <- vapply(all_nodes, function(idx) get_depth_at_node(tree, idx), numeric(1))
  prior <- prod(
    (alpha. / ((1 + depth.node[internal_nodes])^beta.)),
    (1 - alpha. / ((1 + depth.node[terminal_nodes])^beta.))
  )
  return(prior)
}

#' Simulate depths of root nodes under the Chipman prior
#' @param n Number of samples.
#' @param alpha. Numeric, root-splitting probability.
#' @param beta. Numeric, depth penalty.
#' @return Integer vector of depths.
rdepth.chip <- function(n, alpha., beta.) {
  depths <- integer(n)
  for (i in seq_len(n)) {
    tree <- generate.tree.chip(alpha., beta.)
    depths[i] <- get_depth_at_node(tree, 1)
  }
  return(depths)
}

#' Simulate number of terminal nodes under the Chipman prior
#' @param n Number of samples.
#' @param alpha. Numeric, root-splitting probability.
#' @param beta. Numeric, depth penalty.
#' @return Integer vector of number of terminal nodes.
rnterm.chip <- function(n, alpha., beta.) {
  nterm <- integer(n)
  for (i in seq_len(n)) {
    tree <- generate.tree.chip(alpha., beta.)
    nterm[i] <- length(get_terminal_nodes_idx(tree))
  }
  return(nterm)
}

#' Simulate joint (depth, nterm) under Chipman prior
#' @param n Number of samples.
#' @param alpha. Numeric, root-splitting probability.
#' @param beta. Numeric, depth penalty.
#' @return Data frame with columns depth and nterm.
rdepth.nterm.chip <- function(n, alpha., beta.) {
  depths <- integer(n)
  nterms <- integer(n)
  for (i in seq_len(n)) {
    tree <- generate.tree.chip(alpha., beta.)
    depths[i] <- get_depth_at_node(tree, 1)
    nterms[i] <- length(get_terminal_nodes_idx(tree))
  }
  return(data.frame(depth = depths, nterm = nterms))
}

# Note: get_terminal_nodes_idx and get_internal_nodes_idx must be defined elsewhere for this file to work!
