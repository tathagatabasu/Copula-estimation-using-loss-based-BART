# check if an idx is good for swap
check_swap <- function(node.idx, term.node.idx) {
  # Check if both children (node.idx + 1 and node.idx + 2) are terminal nodes
  children_are_terminal <- ((node.idx + 1) %in% term.node.idx) && ((node.idx + 2) %in% term.node.idx)
  # Return TRUE if swap is allowed (i.e., not both children terminal), FALSE otherwise
  return(!children_are_terminal)
}

# function to get valid set of swap nodes
get_valid_swap_idx <- function(tree_top) {
  # Get terminal and internal node indices
  term.node.idx <- get_terminal_nodes_idx(tree_top)
  all.node.idx <- seq_len(max(term.node.idx))
  intern.node.idx <- setdiff(all.node.idx, term.node.idx)
  
  # Filter internal nodes where swap is allowed by check_swap
  intern.node.idx.check <- vapply(intern.node.idx, function(x) check_swap(x, term.node.idx), logical(1))
  intern.node.idx.swap <- intern.node.idx[intern.node.idx.check]
  
  node.idx.df <- data.frame(idx.1 = integer(), idx.2 = integer())
  
  for (node.idx in intern.node.idx.swap) {
    subtree <- get_offsprings(node.idx, tree_top)
    parent_cond_idx <- subtree$cond$x.idx
    
    # Check left child
    if (!is.null(subtree$left) && subtree$left$node.idx %in% intern.node.idx) {
      left_cond_idx <- subtree$left$cond$x.idx
      if (left_cond_idx != parent_cond_idx) {
        node.idx.df <- rbind(node.idx.df, data.frame(idx.1 = node.idx, idx.2 = subtree$left$node.idx))
      }
    }
    
    # Check right child
    if (!is.null(subtree$right) && subtree$right$node.idx %in% intern.node.idx) {
      right_cond_idx <- subtree$right$cond$x.idx
      if (right_cond_idx != parent_cond_idx) {
        node.idx.df <- rbind(node.idx.df, data.frame(idx.1 = node.idx, idx.2 = subtree$right$node.idx))
      }
    }
  }
  
  return(node.idx.df)
}

# swap two nodes conditions
swap_node_condition <- function(first.node.idx, second.node.idx, tree_top){
  cond1 <- get_node_condition(first.node.idx, tree_top)
  cond2 <- get_node_condition(second.node.idx, tree_top)
  tree_top1 <- set_node_condition(first.node.idx, tree_top, cond2)
  tree_top2 <- set_node_condition(second.node.idx, tree_top1, cond1)
  return(tree_top2)
}

# given an idx of terminal, create a split
grow_terminal <- function(node.idx, tree_top, new.cond, value.l = 1, value.r = 1) {
  # If current node is terminal
  if (is.null(tree_top$left) && is.null(tree_top$right)) {
    if (tree_top$node.idx == node.idx) {
      # Replace terminal node with new internal node and two terminal children
      return(list(
        left = list(left = NULL, right = NULL, value = value.l),
        right = list(left = NULL, right = NULL, value = value.r),
        cond = new.cond,
        node.idx = tree_top$node.idx  # Preserve node index
      ))
    } else {
      return(tree_top)
    }
  } else {
    # Recursively descend left and right subtrees
    tree.left <- grow_terminal(node.idx, tree_top$left, new.cond, value.l, value.r)
    tree.right <- grow_terminal(node.idx, tree_top$right, new.cond, value.l, value.r)
    return(list(
      left = tree.left,
      right = tree.right,
      node.idx = tree_top$node.idx,
      cond = tree_top$cond
    ))
  }
}

# this is assuming obs.per.term = 1
check_grow_valid <- function(node.idx, tree_top, X){
  obs.at.node <- get_obs_at_node(node.idx = node.idx,
                                 X = X,
                                 tree_top = tree_top, X.orig = X)
  valid.pred <- get_set_valid_predictor(X.sel = obs.at.node,
                                        n.term.left = 1,
                                        n.term.right = 1,
                                        obs.per.term = 1)
  return(length(valid.pred) > 0)
}

# check if a node can be pruned 
check_prune_valid <- function(node.idx, tree_top){
  node.offspring <- get_offsprings(node.idx, tree_top)
  return(is.null(node.offspring$left$left) & 
           is.null(node.offspring$left$right) & 
           is.null(node.offspring$right$left) & 
           is.null(node.offspring$right$right))
}

# get idx of nodes that can be pruned
get_prune_idx <- function(tree_top){
  intern.node.idx <- get_internal_nodes_idx(tree_top)
  check.prune <- vapply(intern.node.idx, \(x) check_prune_valid(x, tree_top), TRUE)
  return(intern.node.idx[check.prune])
}


# prune a node given the idx
prune_terminal <- function(node.idx, tree_top, value.node) {
  # If current node is terminal, just return it (nothing to prune here)
  if (is.null(tree_top$left) && is.null(tree_top$right)) {
    return(tree_top)
  }
  
  # If current node is the one to prune, replace it with a terminal node
  if (tree_top$node.idx == node.idx) {
    return(list(left = NULL, right = NULL, value = value.node, node.idx = tree_top$node.idx))
  } else {
    # Otherwise, recursively prune left and right subtrees
    tree.left <- prune_terminal(node.idx, tree_top$left, value.node)
    tree.right <- prune_terminal(node.idx, tree_top$right, value.node)
    
    # Return updated internal node with its (possibly pruned) subtrees
    return(list(
      left = tree.left,
      right = tree.right,
      node.idx = tree_top$node.idx,
      cond = tree_top$cond
    ))
  }
}
