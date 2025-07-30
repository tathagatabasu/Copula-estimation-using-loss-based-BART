#' Evaluate a splitting rule
#' 
#' @param cond 
#' `list` representing a splitting rule in the form `list(x.idx = .., x.val = ..)` 
#' where `x.idx` represents the index of the value of `x` on which the condition is evaluated, 
#' and `x.val` is the splitting value of the condition. 
#' If `x.val` is numeric the condition is `x[x.idx] <= x.val`; if `x.val` is a character or a  
#' vectore of characters the condition is `x[x.idx] %in% x.val`.
#' @param x 
#' A row of a `data.frame` representing an observation on which the condition is evaluated
#'
#' @return `logical` TRUE if the condition is verified and FALSE otherwise
#'
#' @examples 
#' x = data.frame(x1 = 1, x2 = 'A')
#' cond = list(x.idx = 1, x.val = 0)
#' eval_cond(cond, x)
#' x = data.frame(x1 = 1, x2 = 'A')
#' cond = list(x.idx = 1, x.val = 2)
#' eval_cond(cond, x)
#' x = data.frame(x1 = 1, x2 = 'A')
#' cond = list(x.idx = 2, x.val = c('A', 'B')
#' eval_cond(cond, x)
#' x = data.frame(x1 = 1, x2 = 'A')
#' cond = list(x.idx = 2, x.val = 'C')
#' eval_cond(cond, x)

eval_cond <- function(cond, x) {
  val <- cond$x.val
  idx <- cond$x.idx
  xi <- x[idx]
  
  if (is.character(val)) {
    xi %in% val
  } else if (is.numeric(val)) {
    xi <= val
  } else {
    stop("Unsupported type in cond$x.val")
  }
}

#' Number of leaves for a given binary tree
#'
#' The function returns the number of leaves of the binary tree given as input 
#' @param tree a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch, 
#' `cond` (only for internal nodes) the splitting rule, 
#' `value` (only for terminal nodes) the value at the terminal node
#'
#' @return a scalar representing the number of leaves of `tree`
#' @export
#'
#' @examples 
#' tree_top <- generate_random_binary_tree(4,5)
#' get_num_terminal_nodes(tree_top)
get_num_terminal_nodes <- function(tree) {
  if (is.null(tree)) return(0)
  if (is.null(tree$left) && is.null(tree$right)) return(1)
  get_num_terminal_nodes(tree$left) + get_num_terminal_nodes(tree$right)
}

#' Depth plus 1 of a binary tree
#'
#' The function returns the depth + 1 of the binary tree given as input
#' @param tree a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch, 
#' `cond` (only for internal nodes) the splitting rule, 
#' `value` (only for terminal nodes) the value at the terminal node
#'
#' @return the value of the depth of the tree plus 1
#' @export
#'
#' @examples
#' tree_top <- generate_random_binary_tree(4,5)
#' get_depth_plus1(tree_top)
get_depth_plus1 <- function(tree) {
  if (is.null(tree)) return(0)
  max(get_depth_plus1(tree$left), get_depth_plus1(tree$right)) + 1
}

#' Depth of a binary tree
#'
#' The function returns the depth of the binary tree given as input
#' @param tree a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch, 
#' `cond` (only for internal nodes) the splitting rule, 
#' `value` (only for terminal nodes) the value at the terminal node
#'
#' @return the value of the depth of the tree
#' @export
#'
#' @examples
#' tree_top <- generate_random_binary_tree(4,5)
#' get_depth(tree_top)
get_depth <- function(tree) {
  get_depth_plus1(tree) - 1
}

#' Observations on a node of a binary tree
#'
#' The function retrieve the observations at a given node of a binary tree provided as input
#'
#' @param node.idx a scalar representing the index of the node. The nodes are indexed going from left to right 
#' @param X a `data.frame` representing the full set of observations 
#' @param tree_top a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch, 
#' `cond` (only for internal nodes) the splitting rule, 
#' `value` (only for terminal nodes) the value at the terminal node
#' @param node.count a scalar used for internal calculations, the default value is 1. It does not need to be changed.
#'
#' @return a `data.frame` representing the observations corresponding to a node
#' @export
#'
#' @examples
#' XX <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
#' tree_top <- generate_random_binary_tree(4,5)
#' tree_top <- assign_node_idx(tree_top)
#' get_tree_plot.idx(tree_top)
#' tree_top <- assign_split_rules(tree_top, XX)
#' tree_top <- assign_term_node_values(tree_top, 0, 2)
#' get_tree_plot(tree_top)
#' XX.at.node <- get_obs_at_node(node.idx = 2, X = XX, tree_top = tree_top)
#' nrow(XX.at.node)
get_obs_at_node <- function(node.idx = 2, X, tree_top, node.count = 1, X.orig) {
  # Handle edge case where X is a numeric vector (1-row matrix)
  if (is.numeric(X) && is.null(dim(X))) {
    X <- matrix(X, nrow = 1)
    colnames(X) <- names(X.orig)
    idx.obs <- which(apply(X.orig, 1, function(row) all(row == X)))
    rownames(X) <- idx.obs
  }
  
  # Early exit if no observations
  if (nrow(X) == 0) return(X)
  
  # Terminal node
  if (is.null(tree_top$left) && is.null(tree_top$right)) {
    return(if (node.idx == node.count) X else NULL)
  }
  
  # Internal node
  if (node.idx == node.count) return(X)
  
  # Evaluate condition
  cond.sel.idx <- vapply(seq_len(nrow(X)), function(i) eval_cond(tree_top$cond, X[i, ]), logical(1))
  
  X.left <- X[cond.sel.idx, , drop = FALSE]
  X.right <- X[!cond.sel.idx, , drop = FALSE]
  
  node.count.left <- node.count + 1
  n.term.left <- get_num_terminal_nodes(tree_top$left)
  node.count.right <- node.count + 2 * n.term.left
  
  if (node.idx < node.count.right) {
    return(get_obs_at_node(node.idx, X.left, tree_top$left, node.count.left, X.orig))
  } else {
    return(get_obs_at_node(node.idx, X.right, tree_top$right, node.count.right, X.orig))
  }
}

#' Terminal node index values
#'
#' The function returns the index corresponding to the terminal nodes of 
#' the binary tree provided as input. The indexes are placed going from left
#' to right
#' @param tree_top a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch, 
#' `cond` (only for internal nodes) the splitting rule, 
#' `value` (only for terminal nodes) the value at the terminal node
#' @param counter a scalar used for internal calculations, the default value is 1. It does not need to be changed.
#'
#' @return the node indexes corresponding to the leaves of the binary tree
#' @export
#'
#' @examples
#' tree_top <- generate_random_binary_tree(4,5)
#' tree_top <- assign_node_idx(tree_top)
#' get_tree_plot.idx(tree_top)
#' get_terminal_nodes_idx(tree_top)
get_terminal_nodes_idx <- function(tree_top, counter = 1) {
  if (is.null(tree_top$left) && is.null(tree_top$right)) {
    return(counter)
  }
  
  count.left <- counter + 1
  n.term.left <- get_num_terminal_nodes(tree_top$left)
  count.right <- counter + 2 * n.term.left
  
  idx.left <- get_terminal_nodes_idx(tree_top$left, counter = count.left)
  idx.right <- get_terminal_nodes_idx(tree_top$right, counter = count.right)
  
  c(idx.left, idx.right)
}

#' Internal node index values
#'
#' The function returns the index corresponding to the internal nodes of 
#' the binary tree provided as input. The indexes are placed going from left
#' to right
#' @param tree_top a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch, 
#' `cond` (only for internal nodes) the splitting rule, 
#' `value` (only for terminal nodes) the value at the terminal node
#' @param counter a scalar used for internal calculations, the default value is 1. It does not need to be changed.
#'
#' @return the node indexes corresponding to the internal nodes of the binary tree
#' @export
#'
#' @examples
#' tree_top <- generate_random_binary_tree(4,5)
#' tree_top <- assign_node_idx(tree_top)
#' get_tree_plot.idx(tree_top)
#' get_internal_nodes_idx(tree_top)
get_internal_nodes_idx <- function(tree_top) {
  n.nodes <- 2 * get_num_terminal_nodes(tree_top) - 1
  setdiff(seq_len(n.nodes), get_terminal_nodes_idx(tree_top))
}

#' Create name from splitting rule
#'
#' The function returns a string representing the condition provided as input. 
#' It is used by the `get_tree_plot()` function to fill the internal nodes.
#' @param cond `list` representing a splitting rule in the form `list(x.idx = .., x.val = ..)` 
#' where `x.idx` represents the index of the value of `x` on which the condition is evaluated, 
#' and `x.val` is the splitting value of the condition. 
#' If `x.val` is numeric the condition is `x[x.idx] <= x.val`; if `x.val` is a character or a  
#' vectore of characters the condition is `x[x.idx] %in% x.val`.
#'
#' @return a string representing the condition
#' @export
#'
#' @examples
#' cond <- list(x.idx = 1, x.val = 2)
#' node.name.from.cond(cond)
#' cond <- list(x.idx = 1, x.val = c('A', 'B'))
#' node.name.from.cond(cond)
node.name.from.cond <- function(cond) {
  if (is.character(cond$x.val)) {
    paste0("X", cond$x.idx, " in {", paste(cond$x.val, collapse = ", "), "}")
  } else {
    paste0("X", cond$x.idx, " <= ", format(round(cond$x.val, 3), trim = TRUE))
  }
}

#' Reshape tree for plotting
#' 
#' The function reshapes the tree topology in a format suitable to be used by 
#' the `data.tree` R-pacakge. It is used by the `get_tree_plot` function.
#'
#' @param tree_top a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch, 
#' `cond` (only for internal nodes) the splitting rule, 
#' `value` (only for terminal nodes) the value at the terminal node
#'
#' @return a nested `list` representing a binary tree 
#' @export
#'
#' @examples 
#' XX <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
#' tree_top <- generate_random_binary_tree(4,5)
#' tree_top <- assign_node_idx(tree_top)
#' tree_top <- assign_split_rules(tree_top, XX)
#' tree_top <- assign_term_node_values(tree_top, 0, 2)
#' tree.to.plot(tree_top)
tree.to.plot <- function(tree_top) {
  if (is.null(tree_top$left) && is.null(tree_top$right)) {
    list(
      left = NULL,
      right = NULL,
      node.name = round(tree_top$value, 3)
    )
  } else {
    list(
      node.name = node.name.from.cond(tree_top$cond),
      children = list(
        left = tree.to.plot(tree_top$left),
        right = tree.to.plot(tree_top$right)
      )
    )
  }
}

#' Reshape tree for plotting
#' 
#' The function reshapes the tree topology in a format suitable to be used by 
#' the `data.tree` R-pacakge. It is used by the `get_tree_plot.idx` function.
#'
#' @param tree_top a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch,
#' `node.idx` the index of the node.
#'
#' @return a nested `list` representing a binary tree 
#' @export
#'
#' @examples 
#' XX <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
#' tree_top <- generate_random_binary_tree(4,5)
#' tree_top <- assign_node_idx(tree_top)
#' tree.to.plot.idx(tree_top)
tree.to.plot.idx <- function(tree_top) {
  if (is.null(tree_top$left) && is.null(tree_top$right)) {
    list(node.name = tree_top$node.idx)
  } else {
    list(
      node.name = tree_top$node.idx,
      children = list(
        left = tree.to.plot.idx(tree_top$left),
        right = tree.to.plot.idx(tree_top$right)
      )
    )
  }
}

#' Plot a binary tree with splitting rules and values at terminal nodes
#' 
#' The function plots a binary tree given as input. For the internal nodes, 
#' the splitting rules are shown, assuming that if they are true the left branch is activated.
#' For the terminal nodes, the value of the leaves is shown.
#' @param tree_top a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch,
#' `cond` (only for internal nodes) the splitting rule, 
#' `value` (only for terminal nodes) the value at the terminal node
#'
#' @return a plot of the tree showing at each internal node the splitting rule and at each terminal node a value.
#' @export
#'
#' @examples
#' XX <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
#' tree_top <- generate_random_binary_tree(4,5)
#' tree_top <- assign_node_idx(tree_top)
#' tree_top <- assign_split_rules(tree_top, XX)
#' tree_top <- assign_term_node_values(tree_top, 0, 2)
#' get_tree_plot(tree_top)
get_tree_plot <- function(tree_top){
  mod.tree.top <- tree.to.plot(tree_top)
  tree.d <- as.Node(mod.tree.top, 
                    mode = 'explicit', 
                    childrenName = 'children', 
                    nameName = 'node.name')
  plot(tree.d)
  
}

#' Plot a binary tree with nodes index value
#' 
#' The function plots a binary tree given as input. The index of each node is shown.
#' 
#' @param tree_top a tree topology provided as nested `list` with arguments 
#' `left` the left branch, `right` the right branch,
#' `node.idx` the index of the node.
#' 
#' @return a plot of the binary tree showing the index value of each node
#' @export
#'
#' @examples
#' tree_top <- generate_random_binary_tree(4,5)
#' tree_top <- assign_node_idx(tree_top)
#' get_tree_plot.idx(tree_top)
get_tree_plot.idx <- function(tree_top){
  mod.tree.top <- tree.to.plot.idx(tree_top)
  tree.d <- as.Node(mod.tree.top, 
                    mode = 'explicit', 
                    childrenName = 'children', 
                    nameName = 'node.name')
  plot(tree.d)
  
}

# extract offsprings of a given node
get_offsprings <- function(node.idx, tree_top) {
  if (is.null(tree_top$left) && is.null(tree_top$right)) return(NULL)
  
  if (tree_top$node.idx == node.idx) {
    return(tree_top)
  }
  
  left_result <- get_offsprings(node.idx, tree_top$left)
  if (!is.null(left_result)) return(left_result)
  
  get_offsprings(node.idx, tree_top$right)
}

# extract condition at node for switch
get_node_condition <- function(node.idx, tree_top) {
  if (is.null(tree_top$left) && is.null(tree_top$right)) return(NULL)
  
  if (tree_top$node.idx == node.idx) {
    return(list(x.idx = tree_top$cond$x.idx, x.val = tree_top$cond$x.val))
  }
  
  cond.left <- get_node_condition(node.idx, tree_top$left)
  if (!is.null(cond.left)) return(cond.left)
  
  get_node_condition(node.idx, tree_top$right)
}
