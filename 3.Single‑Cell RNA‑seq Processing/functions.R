
down_stream= function(obj){
    
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj, nfeatures = 3000)
    obj <- ScaleData(obj, features = rownames(obj))
    
    obj <- RunPCA(obj, features = VariableFeatures(object = obj))
    pct <- obj [["pca"]]@stdev / sum(obj [["pca"]]@stdev) * 100
    cumu <- cumsum(pct)
    co1 <- which(cumu > 90 & pct < 5)[1]
    co1
    co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
    co2
    pcs <- min(co1, co2)
    print(pcs)
    
    obj <- FindNeighbors(obj, dims = 1:pcs, reduction = "pca")
    obj <- FindClusters(obj, resolution =0.1, cluster.name = "unintegrated_clusters")
    
    obj <- RunUMAP(obj, dims = 1:pcs, reduction = "pca", reduction.name = "umap.unintegrated")
    return(obj)
}

auto_cluster_resolution <- function(
  obj,
  reduction = "pca",
  dims = 1:24,
  resolutions = seq(0.1, 1, 0.1),
  graph.name = "autores_snn",
  distance = "euclidean",
  max_cells_for_sil = 10000,  
  min_cluster_size = 10,     
  prefer_simpler = TRUE,      
  eps = 0.02,               
  palette = NULL,             
  verbose = TRUE
) {
  stopifnot(inherits(obj, "Seurat"))
  suppressPackageStartupMessages({
    require(Seurat)
    require(dplyr)
    require(ggplot2)
    require(clustree)
    require(cluster)
    require(tidyr)
    require(purrr)
  })


  if (verbose) message("Finding neighbors...")
  obj <- Seurat::FindNeighbors(obj, reduction = reduction, dims = dims,
                               graph.name = graph.name, verbose = FALSE)
  if (verbose) message("Finding clusters at multiple resolutions: ",
                       paste(resolutions, collapse = ", "))
  obj <- Seurat::FindClusters(obj, graph.name = graph.name,
                              resolution = resolutions, verbose = FALSE)

  prefix <- paste0(graph.name, "_res.")
  p_tree <- clustree::clustree(obj, prefix = prefix)


  emb <- Seurat::Embeddings(obj[[reduction]])[, dims, drop = FALSE]
  n_cells <- nrow(emb)
  if (n_cells > max_cells_for_sil) {
    if (verbose) message("Sampling cells for silhouette: ",
                         max_cells_for_sil, " / ", n_cells)
    set.seed(1)
    sel_idx <- sample(seq_len(n_cells), max_cells_for_sil)
  } else {
    sel_idx <- seq_len(n_cells)
  }
  emb_sel <- emb[sel_idx, , drop = FALSE]
  if (verbose) message("Computing distance matrix on ", nrow(emb_sel), " cells...")
  distance_matrix <- dist(emb_sel, method = distance)


  meta <- obj@meta.data
  cluster_cols <- grep(paste0("^", prefix), colnames(meta), value = TRUE)
  if (length(cluster_cols) == 0) stop("找不到聚类列（前缀：", prefix, "）。")

  res_vals <- suppressWarnings(as.numeric(sub(paste0("^", prefix), "", cluster_cols)))
  ord <- order(res_vals)
  cluster_cols <- cluster_cols[ord]
  res_vals <- res_vals[ord]


  if (verbose) message("Evaluating silhouette by resolution...")
  sil_summ <- purrr::map2_dfr(cluster_cols, res_vals, function(col, resv) {
    cl_all <- meta[[col]]
    cl_all <- as.factor(cl_all)

    if (nlevels(cl_all) < 2) {
      return(tibble::tibble(
        resolution = resv, k = nlevels(cl_all),
        mean_sil = NA_real_, median_sil = NA_real_, min_sil = NA_real_,
        min_size = min(table(cl_all)),
        pass_min_size = FALSE
      ))
    }
    cl_sel <- droplevels(cl_all[sel_idx])

    if (any(table(cl_sel) < 2) || nlevels(cl_sel) < 2) {
      return(tibble::tibble(
        resolution = resv, k = nlevels(cl_all),
        mean_sil = NA_real_, median_sil = NA_real_, min_sil = NA_real_,
        min_size = min(table(cl_all)),
        pass_min_size = all(table(cl_all) >= min_cluster_size)
      ))
    }
    sil <- cluster::silhouette(as.numeric(cl_sel), dist = distance_matrix)
    tibble::tibble(
      resolution   = resv,
      k            = nlevels(cl_all),
      mean_sil     = mean(sil[, 3]),
      median_sil   = median(sil[, 3]),
      min_sil      = min(sil[, 3]),
      min_size     = min(table(cl_all)),
      pass_min_size = all(table(cl_all) >= min_cluster_size)
    )
  })

  sil_ok <- sil_summ %>% filter(pass_min_size, !is.na(mean_sil))
  base_tbl <- if (nrow(sil_ok) > 0) sil_ok else sil_summ %>% filter(!is.na(mean_sil))
  if (nrow(base_tbl) == 0) stop("The contour coefficients for all resolutions are unavailable (possibly due to too small sampling or too few clusters)")
  max_mean <- max(base_tbl$mean_sil, na.rm = TRUE)
  cand <- base_tbl %>% filter(mean_sil >= max_mean - eps)
  rec_row <- if (prefer_simpler) cand %>% arrange(resolution) %>% slice(1)
             else cand %>% arrange(desc(resolution)) %>% slice(1)
  rec_res <- rec_row$resolution

  if (verbose) {
    message(sprintf("Recommended resolution = %.3g  (k=%d, mean silhouette=%.3f)",
                    rec_res, rec_row$k, rec_row$mean_sil))
  }


  p_sil_curve <- ggplot(sil_summ, aes(x = resolution, y = mean_sil)) +
    geom_line() + geom_point() +
    geom_vline(xintercept = rec_res, linetype = 2, color = "red") +
    labs(x = "Resolution", y = "Mean silhouette (sampled cells)",
         title = "Silhouette vs. resolution") +
    theme_bw()


  col_rec <- paste0(prefix, rec_res)
  col_rec <- cluster_cols[which.min(abs(res_vals - rec_res))]
  cl_rec  <- as.factor(meta[[col_rec]])
  k_rec <- nlevels(cl_rec)
  if (is.null(palette)) {
    pal <- scales::hue_pal()(k_rec)
  } else {
    if (length(palette) < k_rec) {
      warning("the length of the palette is less than the number of clusters, it will be reused cyclically.")
    }
    pal <- rep(palette, length.out = k_rec)
  }

  sil_cell <- rep(NA_real_, n_cells)
  cl_sel_rec <- droplevels(cl_rec[sel_idx])
  p_bar <- NULL
  if (nlevels(cl_sel_rec) >= 2 && all(table(cl_sel_rec) >= 2)) {
    sil_rec <- cluster::silhouette(as.numeric(cl_sel_rec), dist = distance_matrix)
    sil_cell[sel_idx] <- sil_rec[, 3]
    mean_sil_overall <- mean(sil_rec[, 3])

    df_bar <- obj@meta.data %>%
      mutate(
        .barcode = rownames(.),
        .cluster = cl_rec,
        .sil     = sil_cell
      ) %>%
      arrange(.cluster, desc(.sil)) %>%
      mutate(.barcode = factor(.barcode, levels = .barcode))

    p_bar <- ggplot(df_bar, aes(x = .barcode, y = .sil, fill = .cluster)) +
      geom_col(show.legend = FALSE) +
      geom_hline(yintercept = mean_sil_overall, color = "red", linetype = "dashed") +
      scale_x_discrete(name = "Cells") +
      scale_y_continuous(name = "Silhouette score") +
      scale_fill_manual(values = pal) +
      ggtitle(paste0(col_rec)) +
      theme_bw() +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
  } else {
    if (verbose) message("The recommended resolution is to skip the histogram for certain clusters where the number of cells is less than 2 after sampling.")
  }

  out <- list(
    object = obj,
    recommended_resolution = rec_res,
    summary = sil_summ,
    plots = list(
      clustree = p_tree,
      sil_curve = p_sil_curve,
      bar = p_bar
    ),
    prefix = prefix,
    graph.name = graph.name
  )
  return(out)
}


add_percent_spe <- function(obj, spe_n){
    spe_percent_vector = c()

    for (i in 1:length(spe_n)) {
      spe_percent <- paste0("percent.", spe_n[i], "_genes")
  
      spe_percent_vector <- append(spe_percent_vector, spe_percent)
  
  obj[[spe_percent]] <- PercentageFeatureSet(obj, pattern = paste0("^", spe_n[i]))
}
    print(spe_percent_vector)

    return(obj)
    
}

plot_spe_density <- function(obj, spe_percent_vector){

    for (i in 1:length(spe_percent_vector)){

        select = spe_percent_vector[[i]]
        # print(select)
        
        spe_data = obj@meta.data[ , select]

        density_data <- density(spe_data)

        turnpoints_result <- turnpoints(density_data$y)

        plot(density_data, main = select)
        points(density_data$x[turnpoints_result$peaks], 
       density_data$y[turnpoints_result$peaks], 
       col = "red", pch = 19)
        points(density_data$x[turnpoints_result$pits], 
       density_data$y[turnpoints_result$pits], 
       col = "blue", pch = 19)

        print(density_data$x[turnpoints_result$pits])
    }
}


library(dplyr)
library(ggplot2)
.cap_sentence <- function(x) {
  if (is.null(x)) return(x)
  sub("^\\s*([[:alpha:]])", "\\U\\1", x, perl = TRUE)
}

.pick_top_simple <- function(x, n = 10, direction = c("up","down")) {
  direction <- match.arg(direction)
  df <- if (inherits(x, "enrichResult")) as.data.frame(x) else as.data.frame(x)
  if (nrow(df) == 0) return(dplyr::tibble())
  pcol <- if ("pvalue" %in% names(df)) "pvalue" else "p.adjust"

  df %>%
    mutate(
      logp = -log10(.data[[pcol]]),
      direction = direction
    ) %>%
    arrange(.data[[pcol]], desc(Count)) %>%
    slice_head(n = n) %>%
    rename(term = Description)
}

.wrap_labels <- function(x, width) {
  if (!is.finite(width)) return(x)
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
}

plot_go_lollipop_divergent <- function(
  ego_up, ego_down,
  top_n = 10,
  left_group = c("down","up"),
  term_text_size = 11,
  colors = c("#6AA84F", "white", "#D66075"),
  title = NULL,
  wrap_width = 28,        
  term_lineheight = 1.05  
) {
  left_group <- match.arg(left_group)

  up_df   <- .pick_top_simple(ego_up,   n = top_n, direction = "up")
  down_df <- .pick_top_simple(ego_down, n = top_n, direction = "down")
  dat <- bind_rows(up_df, down_df)
  if (nrow(dat) == 0) stop("There are no items to draw.")

  dat <- dat %>% mutate(signed_logp = ifelse(direction == left_group, -logp, logp))

  lev_left  <- dat %>% filter(direction == left_group) %>% arrange(desc(logp)) %>% pull(term)
  lev_right <- dat %>% filter(direction != left_group) %>% arrange(desc(logp)) %>% pull(term)
  dat <- dat %>% mutate(term = factor(term, levels = rev(c(lev_left, lev_right))))

  max_abs <- ceiling(max(abs(dat$signed_logp), na.rm = TRUE))
  legend_breaks <- c(-10, -5, 0, 5, 10)

  legend_breaks <- legend_breaks[legend_breaks >= -max_abs & legend_breaks <=  max_abs]
  legend_labels <- abs(legend_breaks)
  
  brks <- pretty(c(-max_abs, max_abs))
  brks <- unique(c(-rev(brks[brks > 0]), 0, brks[brks > 0]))

  ggplot(dat, aes(y = term)) +
    geom_segment(aes(x = 0, xend = signed_logp, yend = term),
                 linewidth = 0.8, colour = "grey75", lineend = "round") +
    geom_point(aes(x = signed_logp, size = Count, fill = signed_logp),
               shape = 21, colour = "black", stroke = 0.2) +
    scale_size_continuous(
      name = "Count",
      range = c(2.5, 7),
      guide = guide_legend(
        override.aes = list(fill = "white", colour = "black", stroke = 0.2),
        order = 1
      )
    ) +
    scale_fill_gradient2(
      name = "-log10pvalue",
      low = colors[1], mid = colors[2], high = colors[3], midpoint = 0,
      breaks = legend_breaks,
      labels = legend_labels,
      limits = c(-max(max_abs, max(abs(legend_breaks))),
                max(max_abs, max(abs(legend_breaks)))),
      guide = guide_colorbar(order = 2)
    ) +
    scale_x_continuous(
      breaks = brks, labels = abs(brks),
      expand = expansion(mult = c(0.05, 0.05))
    ) +

    scale_y_discrete(labels = function(x) .wrap_labels(.cap_sentence(x), wrap_width)) +
    labs(x = "-log10pvalue", y = NULL, title = title) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      axis.text.y  = element_text(size = term_text_size, lineheight = term_lineheight),
      axis.line.x  = element_line(linewidth = 0.6, colour = "black"),
      axis.line.y  = element_line(linewidth = 0.8, colour = "black"),
      plot.margin  = margin(10, 12, 10, 12)  
    )
}