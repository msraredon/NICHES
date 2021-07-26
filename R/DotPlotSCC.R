#' @importFrom rlang %||%
DotPlotSCC <- function (object, assay = NULL, features, cols = c("lightgrey", 
                                                   "blue"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
          idents = NULL, group.by = NULL, split.by = NULL, cluster.idents = FALSE, 
          scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA) 
{
  assay <- assay %||% Seurat::DefaultAssay(object = object)
  Seurat::DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in% 
                                                   rownames(x = RColorBrewer::brewer.pal.info))
  scale.func <- switch(EXPR = scale.by, size = ggplot2::scale_size, 
                       radius = ggplot2::scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(X = 1:length(features), 
                                        FUN = function(x) {
                                          return(rep(x = names(x = features)[x], each = length(features[[x]])))
                                        }))
    if (any(is.na(x = feature.groups))) {
      warning("Some feature groups are unnamed.", call. = FALSE, 
              immediate. = TRUE)
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = Seurat::CellsByIdentities(object = object, idents = idents))
  data.features <- Seurat::FetchData(object = object, vars = features, 
                             cells = cells)
  data.features$id <- if (is.null(x = group.by)) {
    Seurat::Idents(object = object)[cells, drop = TRUE]
  }
  else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop("Not enough colors for the number of groups")
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits, sep = "_")
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), 
                        "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident, 
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      #return(mean(x = expm1(x = x))) #Problematic line for SCC values, which can be very large
      return(mean(x=x))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = Seurat:::PercentAbove,  #Edited to reflect Seurat:::
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(what = rbind, args = lapply(X = data.plot, 
                                               FUN = unlist))
    mat <- scale(x = mat)
    id.levels <- id.levels[stats::hclust(d = stats::dist(x = mat))$order]
  }
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  if (length(x = levels(x = data.plot$id)) == 1) {
    scale <- FALSE
    warning("Only one identity present, the expression values will be not scaled", 
            call. = FALSE, immediate. = TRUE)
  }
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot == 
                                                     x, "avg.exp"]
                             if (scale) {
                               data.use <- scale(x = data.use)
                               data.use <- Seurat::MinMax(data = data.use, min = col.min, 
                                                  max = col.max)
                             }
                             else {
                               data.use <- log(x = data.use)
                             }
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, 
                                         breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                    levels = features)
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- vapply(X = as.character(x = data.plot$id), 
                         FUN = gsub, FUN.VALUE = character(length = 1L), pattern = paste0("^((", 
                                                                                          paste(sort(x = levels(x = object), decreasing = TRUE), 
                                                                                                collapse = "|"), ")_)"), replacement = "", 
                         USE.NAMES = FALSE)
    data.plot$colors <- mapply(FUN = function(color, value) {
      return(colorRampPalette(colors = c("grey", color))(20)[value])
    }, color = cols[splits.use], value = avg.exp.scaled)
  }
  color.by <- ifelse(test = split.colors, yes = "colors", no = "avg.exp.scaled")
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(x = feature.groups[data.plot$features.plot], 
                                       levels = unique(x = feature.groups))
  }
  plot <- ggplot2::ggplot(data = data.plot, mapping = ggplot2::aes_string(x = "features.plot", 
                                                        y = "id")) + 
    ggplot2::geom_point(mapping = aes_string(size = "pct.exp", color = color.by)) + 
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) + 
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),axis.title.y = element_blank()) + 
    ggplot2::guides(size = guide_legend(title = "Percent Expressed")) + 
    ggplot2::labs(x = "Features", y = ifelse(test = is.null(x = split.by), 
                                    yes = "Identity", no = "Split Identity")) + 
    cowplot::theme_cowplot()
  if (!is.null(x = feature.groups)) {
    plot <- plot + ggplot2::facet_grid(facets = ~feature.groups, scales = "free_x", 
                              space = "free_x", switch = "y") + 
      ggplot2::theme(panel.spacing = ggplot2::unit(x = 1,units = "lines"), strip.background = element_blank())
  }
  if (split.colors) {
    plot <- plot + ggplot2::scale_color_identity()
  }
  else if (length(x = cols) == 1) {
    plot <- plot + ggplot2::scale_color_distiller(palette = cols)
  }
  else {
    plot <- plot + ggplot2::scale_color_gradient(low = cols[1], high = cols[2])
  }
  if (!split.colors) {
    plot <- plot + ggplot2::guides(color = ggplot2::guide_colorbar(title = "Average Expression"))
  }
  return(plot)
}
