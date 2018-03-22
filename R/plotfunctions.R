##  Plot functions

#' Heatmap Transition Probabilities
#'
#' Plots a heatmap of an HMM's initial and transition probabilities.
#'
#' @param transition A matrix with transition probabilities (see also
#'   \code{\link[seqHMM]{build_hmm}}).
#' @param initial An (optional) vector of initial probabilities.
#' @param base_size Numerical, a size parameter for the plots made using \code{ggplot2}
#'   (see \link[ggplot2]{theme}), default = 10.
#' @seealso See \code{\link{hmm.clust}} for an example.
#' @export
transition.heatmap <- function(transition, initial = NULL, base_size = 10) {
  if(!is.null(initial)) {
    if(is.numeric(initial)) {
      if(sum(initial) > 1) {
        # Vector initial is provided in percentages, not in decimals
        initial <- initial/100
      }
    } else {
      stop("initial must be a numerical vector.")
    }
  }
  if(class(transition) != "matrix") {
    stop("transition must be a matrix.")
  }
  if(!is.numeric(transition)) {
    stop("transition must be numerical.")
  }
  transition.m <- melt(transition)
  colnames(transition.m) <- c("from", "to", "value")
  if(!is.null(initial)) {
    initial.part <- data.frame(from = matrix("Initial", nrow = nrow(transition)),
                               to = names(initial),
                               value = initial)
    transition.m <- rbind(initial.part, transition.m)
    transition.m$from <- factor(transition.m$from, levels = c(rev(rownames(transition)), "Initial"))
    transition.m$to <- factor(transition.m$to, levels = colnames(transition))
  }

  ggplot(transition.m, aes(transition.m$to, transition.m$from)) +
    geom_tile(aes(fill = transition.m$value), colour = "#272B30") +
    scale_fill_gradient2(name = "Scale", low = "white", high = "black") +
    theme_grey(base_size = base_size) + labs(x = "", y = "") +
    scale_x_discrete(expand = c(0, 0), position = "top") +
    scale_y_discrete(expand = c(0, 0)) +
    coord_fixed(ratio = 1) +
    theme(legend.position = "none", axis.ticks = element_blank(),
          axis.text.x = element_text(size = base_size * 1, angle = 90,
                                     hjust = 0, colour = "black"),
          axis.text.y = element_text(size = base_size, colour = "black"))

}

#' Heatmap Emission Probabilities
#'
#' Plots a heatmap of an HMM's emission probabilities.
#'
#' @param emission A matrix with emission probabilities (see also
#'   \code{\link[seqHMM]{build_hmm}}).
#' @param base_size Numerical, a size parameter for the plots made using \code{ggplot2}
#'   (see \link[ggplot2]{theme}), default = 10.
#' @seealso See \code{\link{hmm.clust}} for an example.
#' @export
emission.heatmap <- function(emission, base_size = 10) {
  if(class(emission) != "matrix") {
    stop("emission must be a matrix.")
  }
  if(!is.numeric(emission)) {
    stop("emission must be numerical.")
  }
  emission.m <- melt(emission)
  colnames(emission.m) <- c("state_names", "symbol_names", "value")
  emission.m$state_names <- factor(emission.m$state_names,
                                   levels = rev(rownames(emission)))
  emission.m$symbol_names <- factor(emission.m$symbol_names,
                                    levels = colnames(emission))
  ggplot(emission.m, aes(emission.m$symbol_names, emission.m$state_names)) +
    geom_tile(aes(fill = emission.m$value), colour = "#272B30") +
    scale_fill_gradient2(name = "Scale", low = "white", high = "black") +
    theme_grey(base_size = base_size) + labs(x = "", y = "") +
    scale_x_discrete(expand = c(0, 0), position = "top") +
    scale_y_discrete(expand = c(0, 0)) +
    coord_fixed(ratio = 1) +
    theme(legend.position = "none", axis.ticks = element_blank(),
          axis.text.x = element_text(size = base_size * 1, angle = 90,
                                     hjust = 0, colour = "black"),
          axis.text.y = element_text(size = base_size, colour = "black"))

}
