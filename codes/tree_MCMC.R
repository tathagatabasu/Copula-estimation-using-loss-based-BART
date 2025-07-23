# the transition probabilities are independent from the nodes that are actually swapped
transition.p.swap <- function(tree_top){
  # valid.swap.idx <- get_valid_swap_idx(tree_top)
  # 1/nrow(valid.swap.idx)
  1
}

# here we assume it is always continuous uniform
transition.p.change <- function(change_move_list){
  # n.int <- length(get_internal_nodes_idx(change_move_list$tree))
  # if(is.character(change_move_list$valid.split)){
  #   unique.valid.split <- unique(change_move_list$valid.split)
  #   prob.split <- 1/length(unique.valid.split)
  # } else {
  #   if(length(change_move_list$valid.split) == 1){
  #     prob.split <- 1
  #   } else {
  #     valid.set.extremes <- range(change_move_list$valid.split)
  #     prob.split <- (1/(valid.set.extremes[2] - valid.set.extremes[1]))
  #   }
  # }
  # (1/n.int)*(1/length(change_move_list$valid.pred))*prob.split
  1
}

transition.p.change.reverse <- function(change_move_list, X, old.tree){
  1
  # n.int <- length(get_internal_nodes_idx(change_move_list$tree))
  # old.pred <- get_node_condition(change_move_list$node.idx, old.tree)$x.idx
  # x.sel <- get_obs_at_node(change_move_list$node.idx, X, change_move_list$tree)
  # n.left <- change_move_list$n.left
  # n.right <- change_move_list$n.right
  # valid.pred <- get_set_valid_predictor(X.sel = x.sel,
  #                                       n.term.left = n.left,
  #                                       n.term.right = n.right,
  #                                       obs.per.term = 1)
  # valid.split <- get_set_valid_split(X.pred = x.sel[,old.pred],
  #                                    n.term.left = n.left,
  #                                    n.term.right = n.right,
  #                                    obs.per.term = 1)
  # if(is.character(valid.split)){
  #   unique.valid.split <- unique(valid.split)
  #   prob.split <- 1/length(unique.valid.split)
  # } else {
  #   if(length(valid.split) == 1){
  #     prob.split <- 1
  #   } else {
  #     valid.set.extremes <- range(valid.split)
  #     prob.split <- (1/(valid.set.extremes[2] - valid.set.extremes[1]))
  #   }
  # }
  # 
  # (1/n.int)*(1/length(valid.pred))*prob.split
}

transition.p.grow <- function(grow_move_list, include.split, cont.unif = TRUE){
  n.term <- get_num_terminal_nodes(grow_move_list$tree) - 1
  if(include.split){
    n.valid.pred <- length(grow_move_list$valid.pred)
    if(is.character(grow_move_list$valid.split)){
      unique.valid.split <- unique(grow_move_list$valid.split)
      prob.split_ <- 1/length(unique.valid.split)
    } else {
      if(length(grow_move_list$valid.split) == 1){
        prob.split_ <- 1
      } else {
        if(cont.unif){
          valid.set.extremes <- range(grow_move_list$valid.split)
          prob.split_ <- (1/(valid.set.extremes[2] - valid.set.extremes[1]))  
        } else {
          prob.split_ <- 1/length(grow_move_list$valid.split)
        }
        
      }
    }
    return((1/n.term)*(1/n.valid.pred)*prob.split_)
  } else{
    return(1/n.term)
  }
  
}

transition.p.grow.reverse <- function(grow_move_list){
  valid.prune.idx <- get_prune_idx(grow_move_list$tree)
  1/length(valid.prune.idx)
}

transition.p.prune <- function(prune_move_list){
  1/length(prune_move_list$prune.node.idx)
}

transition.p.prune.reverse <- function(prune_move_list, old.tree, X, obs.per.term,
                                       cont.unif = TRUE, include.split){
  n.term <- get_num_terminal_nodes(prune_move_list$tree)
  if(include.split){
    node.idx <- prune_move_list$node.idx
    obs.at.node <- get_obs_at_node(node.idx = node.idx, 
                                   X = X, 
                                   tree_top = prune_move_list$tree, X.orig = X)
    valid.pred <- get_set_valid_predictor(X.sel = obs.at.node, 
                                          n.term.left = 1, n.term.right = 1, 
                                          obs.per.term = 1)
    old.cond <- get_node_condition(node.idx, old.tree)
    valid.split <- get_set_valid_split(X.pred = obs.at.node[,old.cond$x.idx], 
                                       n.term.left = 1, n.term.right = 1, 
                                       obs.per.term = 1)
    if(is.character(valid.split)){
      unique.valid.split <- unique(valid.split)
      prob.split_ <- 1/length(unique.valid.split)
    } else {
      if(length(valid.split) == 1){
        prob.split_ <- 1
      } else {
        if(cont.unif){
          valid.set.extremes <- range(valid.split)
          prob.split_ <- (1/(valid.set.extremes[2] - valid.set.extremes[1]))  
        } else{
          prob.split_ <- 1/length(valid.split)
        }
      }
    }
    return((1/n.term)*(1/length(valid.pred))*prob.split_) 
    } else{
      return(1/n.term)
  }
 
}