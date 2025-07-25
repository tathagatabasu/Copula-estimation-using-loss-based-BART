# check if an idx is good for swap
check_swap <- function(node.idx, term.node.idx){
  ! (( (node.idx + 1) %in% term.node.idx ) & ((node.idx + 2) %in% term.node.idx) )
}

# function to get valid set of swap nodes
get_valid_swap_idx <- function(tree_top){
  # take terminal node idx
  term.node.idx <- get_terminal_nodes_idx(tree_top)
  all.node.idx <- 1:max(term.node.idx)
  # take internal node idx
  intern.node.idx <- all.node.idx[-term.node.idx]
  intern.node.idx.check <- vapply(intern.node.idx, \(x) check_swap(x, term.node.idx), TRUE)
  intern.node.idx.swap <- intern.node.idx[intern.node.idx.check]
  node.idx.df <- data.frame(idx.1 = NA, idx.2 = NA)[-1,]
  for(node.idx in intern.node.idx.swap){
    subtree <- get_offsprings(node.idx, tree_top)
    cond.idx <- subtree$cond$x.idx
    if(subtree$left$node.idx %in% intern.node.idx){
      if(subtree$left$cond$x.idx != cond.idx){
        node.idx.df <- rbind(node.idx.df, 
                             data.frame(idx.1 = node.idx, 
                                        idx.2 = subtree$left$node.idx) )  
      }
    } 
    if(subtree$right$node.idx %in% intern.node.idx){
      if(subtree$right$cond$x.idx != cond.idx){
        node.idx.df <- rbind(node.idx.df, 
                             data.frame(idx.1 = node.idx, 
                                        idx.2 = subtree$right$node.idx) )  
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
grow_terminal <- function(node.idx, tree_top, new.cond, value.l = 1, value.r = 1){
  if(is.null(tree_top$left) & is.null(tree_top$right)){
    if(tree_top$node.idx == node.idx){
      return(list(left = list(left = NULL, right = NULL, value = value.l),
                  right = list(left = NULL, right = NULL, value = value.r),
                  cond = new.cond))
    } else {
      return(tree_top)
    }
  } else {
    tree.left <- grow_terminal(node.idx, tree_top$left, new.cond, value.l, value.r)
    tree.right <- grow_terminal(node.idx, tree_top$right, new.cond, value.l, value.r)
    return(list(left = tree.left, right = tree.right, node.idx = tree_top$node.idx,
                cond = tree_top$cond))
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
prune_terminal <- function(node.idx, tree_top, value.node){
  if(is.null(tree_top$left) & is.null(tree_top$right)){
    return(tree_top) 
  }
  if(tree_top$node.idx == node.idx){
    return(list(left = NULL, right = NULL, value = value.node))
  }  else {
    tree.left <- prune_terminal(node.idx, tree_top$left, value.node)
    tree.right <- prune_terminal(node.idx, tree_top$right, value.node)
    return(list(left = tree.left, right = tree.right, node.idx = tree_top$node.idx,
                cond = tree_top$cond))
  }
}
