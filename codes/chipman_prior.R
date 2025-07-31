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
  # Initialize tree root node
  tree <- list(node.idx = 1, depth = 0, left = NULL, right = NULL)
  
  # Maintain queue of nodes to possibly split (store references as paths)
  nodes_to_split <- list(tree)
  
  split_count <- 0
  max_depth <- 100
  
  while (length(nodes_to_split) > 0 && split_count < max_depth) {
    # Pop first node to split
    current_node <- nodes_to_split[[1]]
    nodes_to_split <- nodes_to_split[-1]
    
    split_count <- split_count + 1
    p_split <- alpha. / (1 + current_node$depth)^beta.
    
    if (runif(1) < p_split) {
      # Create children
      left_node <- list(node.idx = current_node$node.idx * 2, depth = current_node$depth + 1, left = NULL, right = NULL)
      right_node <- list(node.idx = current_node$node.idx * 2 + 1, depth = current_node$depth + 1, left = NULL, right = NULL)
      
      current_node$left <- left_node
      current_node$right <- right_node
      
      # Add children to queue to split later
      nodes_to_split <- c(nodes_to_split, list(left_node, right_node))
    }
    
    # Rebuild tree from updated nodes: Because R uses copy-on-modify,
    # we need a helper function to update nodes correctly
    update_tree_node <- function(tree, node) {
      if (tree$node.idx == node$node.idx) {
        return(node)
      } else if (!is.null(tree$left) && node$node.idx %in% c(tree$left$node.idx, tree$left$node.idx * 2, tree$left$node.idx * 2 + 1)) {
        tree$left <- update_tree_node(tree$left, node)
      } else if (!is.null(tree$right)) {
        tree$right <- update_tree_node(tree$right, node)
      }
      tree
    }
    
    tree <- update_tree_node(tree, current_node)
  }
  
  if (split_count == max_depth) warning("Max tree depth reached; stopping splits.")
  
  tree
}

#' Calculate the depth of a node in a tree
#' @param tree_top The root node.
#' @param node.idx The index of the node to find.
#' @return Integer depth or NA if not found.
get_depth_at_node <- function(tree_top, node.idx) {
  
  if (tree_top$node.idx == node.idx) return(tree_top$depth)
  
  if (!is.null(tree_top$left)) {
    depth_left <- get_depth_at_node(tree_top$left, node.idx)
    if (!is.na(depth_left)) return(depth_left)
  }
  
  if (!is.null(tree_top$right)) {
    depth_right <- get_depth_at_node(tree_top$right, node.idx)
    if (!is.na(depth_right)) return(depth_right)
  }
  
  NA
}

#' Calculate the prior probability of a tree under the Chipman prior
#' Note: Requires helper functions get_terminal_nodes_idx and get_internal_nodes_idx
#' @param tree The tree object.
#' @param alpha. Numeric, root-splitting probability.
#' @param beta. Numeric, depth penalty.
#' @return Numeric prior probability.
chipman_prior_tree <- function(tree, alpha., beta.) {
  
  terminal_nodes <- get_terminal_nodes_idx(tree)
  internal_nodes <- get_internal_nodes_idx(tree)
  all_nodes <- c(terminal_nodes, internal_nodes)
  
  depth.node <- vapply(all_nodes, function(idx) get_depth_at_node(tree, idx), numeric(1))
  names(depth.node) <- all_nodes  # Safe mapping
  
  # Compute prior components
  internal_part <- alpha. / ((1 + depth.node[as.character(internal_nodes)])^beta.)
  terminal_part <- 1 - alpha. / ((1 + depth.node[as.character(terminal_nodes)])^beta.)
  
  # Final prior
  prior <- prod(c(internal_part, terminal_part))
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
    depths[i] <- get_depth_at_node(tree, 1L)
  }
  
  depths
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
  
  nterm
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
    depths[i] <- get_depth_at_node(tree, 1L)
    nterms[i] <- length(get_terminal_nodes_idx(tree))
  }
  
  data.frame(depth = depths, nterm = nterms)
}

# Note: get_terminal_nodes_idx and get_internal_nodes_idx must be defined elsewhere for this file to work!
