# functions/Run_CellChat.R
Run_CellChat <- function(seurat_obj,
                         Run_CellChat_output_path,
                         Load_QC_species,
                         Run_CellChat_group_by,
                         Run_CellChat_source_celltype,
                         Run_CellChat_target_celltype,
                         Run_CellChat_plot_heatmap,
                         Run_CellChat_ntop_signaling,
                         Run_CellChat_MaxGroup) {
  suppressPackageStartupMessages({
    library(CellChat)
    library(patchwork)
    library(scales)
    library(ComplexHeatmap)
  })
  
  # Check species
  if (Load_QC_species == "human") {
    CellChatDB.use <- CellChatDB.human
  } else if (Load_QC_species == "mouse") {
    CellChatDB.use <- CellChatDB.mouse
  } else {
    stop("[Run_CellChat] Unsupported species: ", Load_QC_species)
  }

  dir.create(Run_CellChat_output_path, showWarnings = FALSE, recursive = TRUE)
  assay_to_use <- DefaultAssay(seurat_obj)
  
  # Check integrated or SCT
  if (assay_to_use %in% c("integrated", "SCT")) {
  message("[Run_CellChat] Detected integrated/SCT assay. Switching to RNA for CellChat input")
  assay_to_use <- "RNA"
  DefaultAssay(seurat_obj) <- "RNA"
  }
  # Check harmony
  harmony_detected <- any(grepl("^harmony", names(seurat_obj@reductions))) |
                      any(grepl("^harmony", colnames(seurat_obj@meta.data)))

  if (harmony_detected) {
    # Using RNA assay as CellChat input
    if ("RNA" %in% names(seurat_obj@assays)) {
      assay_to_use <- "RNA"
      DefaultAssay(seurat_obj) <- "RNA"
      message("[Run_CellChat] Harmony detected, using RNA assay for CellChat input")
    } else {
      stop("[Run_CellChat] Harmony detected but RNA assay not found in Seurat object")
    }
  }
  
  # Flatten if Seurat Assay5 structure detected
  if (inherits(seurat_obj[[assay_to_use]], "Assay5")) {
    message("[Run_CellChat] Detected multi-layer ", assay_to_use, " assay")

    layer_names <- SeuratObject::Layers(seurat_obj[[assay_to_use]])
    message("[Run_CellChat] Found layers: ", paste(layer_names, collapse = ", "))

    # Using data layer as default
    if ("data" %in% layer_names) {
      mat <- GetAssayData(seurat_obj, assay = assay_to_use, layer = "data")
      chosen_layer <- "data"
    } else if ("counts" %in% layer_names) {
      mat <- GetAssayData(seurat_obj, assay = assay_to_use, layer = "counts")
      chosen_layer <- "counts"
    } else {
      stop("[Run_CellChat] No usable data layer found in assay ", assay_to_use)
    }

    flat_assay_name <- paste0(assay_to_use, "_flat")
    message("[Run_CellChat] Using layer '", chosen_layer, "' to create assay: ", flat_assay_name)

    seurat_obj[[flat_assay_name]] <- Seurat::CreateAssayObject(data = mat)
    DefaultAssay(seurat_obj) <- flat_assay_name
    assay_to_use <- flat_assay_name
  }

  # Check grouping
  is_multi_group <- FALSE
  groups <- NULL
  cellchat <- NULL
  cellchat_list <- NULL

  if (Run_CellChat_group_by %in% colnames(seurat_obj@meta.data) &&
      length(unique(seurat_obj@meta.data[[Run_CellChat_group_by]])) > 1 &&
      Run_CellChat_group_by == "sample_group") {
    
    # Multiple groups
    is_multi_group <- TRUE
    groups <- unique(seurat_obj@meta.data[[Run_CellChat_group_by]])
    message("[Run_CellChat] Multi-group CellChat analysis based on: ", Run_CellChat_group_by, "; groups: ", paste(groups, collapse = ", "))

    if (!is.null(Run_CellChat_MaxGroup)) {
      if (!all(Run_CellChat_MaxGroup %in% groups)) {
        stop("[Run_CellChat] Run_CellChat_MaxGroup contains invalid group names: ", paste(setdiff(Run_CellChat_MaxGroup, groups), collapse = ", "))
      }
    } else {
      Run_CellChat_MaxGroup <- groups
    }

    # Run CellChat for each group
    cellchat_list <- list()
    for (grp in groups) {
      message("[Run_CellChat] Running CellChat for group: ", grp)
      sub_obj <- subset(seurat_obj,
                        cells = rownames(seurat_obj@meta.data[seurat_obj@meta.data[[Run_CellChat_group_by]] == grp, ]))
      mat <- SeuratObject::LayerData(sub_obj[[assay_to_use]], layer = "data")

      sub_cellchat <- createCellChat(object = mat,
                                     meta = sub_obj@meta.data,
                                     group.by = "Celltype")
      sub_cellchat@DB <- CellChatDB.use
      sub_cellchat <- subsetData(sub_cellchat)
      sub_cellchat <- identifyOverExpressedGenes(sub_cellchat)
      sub_cellchat <- identifyOverExpressedInteractions(sub_cellchat)
      sub_cellchat <- computeCommunProb(sub_cellchat)
      sub_cellchat <- filterCommunication(sub_cellchat, min.cells = 10)
      sub_cellchat <- computeCommunProbPathway(sub_cellchat)
      sub_cellchat <- aggregateNet(sub_cellchat)
      sub_cellchat <- netAnalysis_computeCentrality(sub_cellchat, slot.name = "netP")
      cellchat_list[[grp]] <- sub_cellchat
    }

    # Merge
    cellchat <- mergeCellChat(cellchat_list, add.names = groups)
    cellchat@DB <- CellChatDB.use
    cellchat <- updateCellChat(cellchat)
    saveRDS(cellchat, file = file.path(Run_CellChat_output_path, "cellchat_merged.rds"))
    message("[Run_CellChat] Merged multi-group cellchat object saved to: ", Run_CellChat_output_path)

  } else {
    # Single group
    message("[Run_CellChat] Single-group CellChat analysis, grouping by: ", Run_CellChat_group_by)

    mat <- SeuratObject::LayerData(seurat_obj[[assay_to_use]], layer = "data")
    if (!is(mat, "AnyMatrix") || nrow(mat) == 0) {
      warning("[Run_CellChat] Data slot invalid, switching to counts.")
      mat <- SeuratObject::LayerData(seurat_obj[[assay_to_use]], layer = "counts")
    }
    if (!is(mat, "AnyMatrix") || nrow(mat) == 0) {
      stop("[Run_CellChat] Expression matrix invalid — cannot run CellChat.")
    }

    cellchat <- createCellChat(object = mat,
                               meta = seurat_obj@meta.data,
                               group.by = Run_CellChat_group_by)
    cellchat@DB <- CellChatDB.use
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- computeCommunProb(cellchat)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    saveRDS(cellchat, file = file.path(Run_CellChat_output_path, "cellchat.rds"))
    message("[Run_CellChat] Single-group cellchat object saved to: ", Run_CellChat_output_path)
  }

  # Plot signaling role heatmap
  if (Run_CellChat_plot_heatmap) {
    message("[Run_CellChat] Drawing netAnalysis_signalingRole_heatmap ...")

    # Multiple groups
    if (is_multi_group) {
      message("[Run_CellChat] Multi-group mode: drawing heatmaps for each group with unified pathway set")

      # Get the union of all pathways
      pathway.union <- unique(unlist(lapply(cellchat_list, function(x) x@netP$pathways)))
      n_rows <- length(pathway.union)
      n_cols <- length(unique(cellchat@idents))
      
      # Plot unit: increase 0.5/0.25 cm per pathway
      heatmap_pdf_height <- max(20, n_rows * 0.25)
      bubble_pdf_height <- max(20, n_rows * 0.5)
      # Plot unit: increase 2 / 2 * group number cm per cell type
      heatmap_pdf_width <- max(20, n_cols * 2)
      bubble_pdf_width <- max(20, n_cols * 2) * length(cellchat_list)

      # Modify font size based on pathway and Celltype number
      row_font_size <- ifelse(n_rows > 60, 8.5, ifelse(n_rows > 30, 9, 10))
      col_font_size <- ifelse(n_cols > 20, 9, 8)

      # Outgoing signaling heatmap
      ht_list_out <- list()
      for (i in seq_along(cellchat_list)) {
        grp_name <- names(cellchat_list)[i]
        ht_list_out[[i]] <- netAnalysis_signalingRole_heatmap(
          cellchat_list[[i]],
          pattern = "outgoing",
          signaling = pathway.union,
          title = grp_name,
          width = heatmap_pdf_width, height = heatmap_pdf_height
        )
        ht_list_out[[i]]@row_names_param$gp <- grid::gpar(fontsize = row_font_size)
        ht_list_out[[i]]@column_names_param$gp <- grid::gpar(fontsize = col_font_size)
      }

      pdf(file = file.path(Run_CellChat_output_path, "Outgoing_signaling.pdf"),
          width = heatmap_pdf_width, 
          height = heatmap_pdf_height)
      draw(Reduce(`+`, ht_list_out), ht_gap = unit(0.5, "cm"))
      dev.off()

      # Incoming signaling heatmap
      ht_list_in <- list()
      for (i in seq_along(cellchat_list)) {
        grp_name <- names(cellchat_list)[i]
        ht_list_in[[i]] <- netAnalysis_signalingRole_heatmap(
          cellchat_list[[i]],
          pattern = "incoming",
          signaling = pathway.union,
          title = grp_name,
          width = heatmap_pdf_width, height = heatmap_pdf_height,
          color.heatmap = "GnBu"
        )
        ht_list_in[[i]]@row_names_param$gp <- grid::gpar(fontsize = row_font_size)
        ht_list_in[[i]]@column_names_param$gp <- grid::gpar(fontsize = col_font_size)
      }

      pdf(file = file.path(Run_CellChat_output_path, "Incoming_signaling.pdf"),
          width = heatmap_pdf_width, 
          height = heatmap_pdf_height)
      draw(Reduce(`+`, ht_list_in), ht_gap = unit(0.5, "cm"))
      dev.off()

      message("[Run_CellChat] signalingRole_heatmap (multi-group, union pathways) saved: incoming & outgoing")

    } else {
      # Single group
      # Get the union of all pathways
      n_rows <- length(cellchat@netP$pathways)
      n_cols <- length(unique(cellchat@idents))
      
      # Plot unit: increase 0.25/0.5 cm per pathway
      heatmap_pdf_height <- max(20, n_rows * 0.25)
      bubble_pdf_height <- max(20, n_rows * 0.5)
      # Plot unit: increase 0.5 cm per cell type
      heatmap_pdf_width <- max(20, n_cols * 0.5)
      bubble_pdf_width <- max(20, n_cols * 0.5)

      # Modify font size based on pathway and Celltype number
      row_font_size <- ifelse(n_rows > 60, 8.5, ifelse(n_rows > 30, 9, 10))
      col_font_size <- ifelse(n_cols > 20, 9, 8)

      # Outgoing signaling heatmap
      ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", width = heatmap_pdf_width, height = heatmap_pdf_height)
      ht1@row_names_param$gp <- grid::gpar(fontsize = row_font_size)
      ht1@column_names_param$gp <- grid::gpar(fontsize = col_font_size)
      # Incoming signaling heatmap
      ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", width = heatmap_pdf_width, height = heatmap_pdf_height)
      ht2@row_names_param$gp <- grid::gpar(fontsize = row_font_size)
      ht2@column_names_param$gp <- grid::gpar(fontsize = col_font_size)
      pdf(file = file.path(Run_CellChat_output_path, "signalingRole_heatmap.pdf"),
          width = heatmap_pdf_width, height = heatmap_pdf_height)
      print(ht1 + ht2)
      dev.off()
      message("[Run_CellChat] signalingRole_heatmap (single-group) saved")
    }
  }
  # Bubble plot
  if (!is_multi_group) {
    # Single group
    comm <- subsetCommunication(cellchat, slot.name = "net")
    pcol <- if ("pval" %in% colnames(comm)) "pval" else if ("p.value" %in% colnames(comm)) "p.value" else NULL
    if (!is.null(pcol)) comm <- comm[comm[[pcol]] < 0.05, ]
    if ("prob" %in% colnames(comm)) comm <- comm[order(-comm$prob), ]
    if (nrow(comm) > 0) {
      top_comm <- head(comm, Run_CellChat_ntop_signaling)
      csv_file <- file.path(Run_CellChat_output_path, paste0("top", Run_CellChat_ntop_signaling, "_LR_singleGroup.csv"))
      write.csv(top_comm, file = csv_file, row.names = FALSE)
      message("[Run_CellChat] Saved top L-R pairs to: ", csv_file)

      pdf(file = file.path(Run_CellChat_output_path, "Bubble_plot_singleGroup.pdf"),
          width = bubble_pdf_width,
          height = bubble_pdf_height)
      print(netVisual_bubble(cellchat,
                             signaling = unique(as.character(top_comm$pathway_name)),
                             sources.use = Run_CellChat_source_celltype,
                             targets.use = Run_CellChat_target_celltype,
                             remove.isolate = TRUE))
      dev.off()
      message("[Run_CellChat] Single-group bubble plot saved")
    }
  } else {
    # Multiple groups
    message("[Run_CellChat] Performing multi-group bubble plots based on Run_CellChat_MaxGroup ...")
    for (max_grp in Run_CellChat_MaxGroup) {
      max_idx <- which(groups == max_grp)
      comparison_vec <- c(1:length(groups))
      n_groups <- length(comparison_vec)
      color_text <- RColorBrewer::brewer.pal(min(n_groups, 8), "Set2")
      if (n_groups > 8) {
        color_text <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_groups)
      }
      
      pdf(file = file.path(Run_CellChat_output_path, paste0("Bubble_plot_maxGroup_", max_grp, ".pdf")),
          width = bubble_pdf_width, height = bubble_pdf_height)
      print(netVisual_bubble(cellchat,
                             comparison = comparison_vec,
                             max.dataset = max_idx,
                             sources.use = Run_CellChat_source_celltype,
                             targets.use = Run_CellChat_target_celltype,
                             remove.isolate = TRUE,
                             title.name = paste("Max signaling in", max_grp),
                             color.text = color_text))
      dev.off()
      message("[Run_CellChat] Bubble plot saved for max group: ", max_grp)
    }

    # Calculate the difference of ligand-receptor pairs between multiple groups
    comm_list <- lapply(cellchat_list, function(cc) {
      df <- subsetCommunication(cc, slot.name = "net")
      if ("p.value" %in% colnames(df) && !"pval" %in% colnames(df)) df$pval <- df[["p.value"]]
      if (!("interaction_name" %in% colnames(df)) && "interaction" %in% colnames(df)) df$interaction_name <- df$interaction
      if (!("source" %in% colnames(df)) && "from" %in% colnames(df)) df$source <- df$from
      if (!("target" %in% colnames(df)) && "to" %in% colnames(df)) df$target <- df$to
      if ("pval" %in% colnames(df)) df <- df[df$pval < 0.05, , drop = FALSE]
      if (nrow(df) > 0) df$key <- paste(df$interaction_name, df$source, df$target, sep = "|")
      df
    })
    names(comm_list) <- groups

    sig_details <- list()
    for (i in seq_along(groups)) {
      for (j in seq_along(groups)) {
        if (i >= j) next
        g1 <- groups[i]; g2 <- groups[j]
        comm1 <- comm_list[[g1]]; comm2 <- comm_list[[g2]]
        keys_to_check <- unique(c(comm1$key, comm2$key))
        if (length(keys_to_check) == 0) next
        for (k in keys_to_check) {
          row1 <- comm1[comm1$key == k, , drop = FALSE]
          row2 <- comm2[comm2$key == k, , drop = FALSE]
          prob1 <- if (nrow(row1) > 0) as.numeric(row1$prob[1]) else 0
          prob2 <- if (nrow(row2) > 0) as.numeric(row2$prob[1]) else 0
          pval1 <- if (nrow(row1) > 0) as.numeric(row1$pval[1]) else NA
          pval2 <- if (nrow(row2) > 0) as.numeric(row2$pval[1]) else NA
          parts <- strsplit(k, "\\|")[[1]]
          sig_details[[length(sig_details) + 1]] <- data.frame(
            interaction = parts[1],
            source = parts[2],
            target = parts[3],
            key = k,
            group1 = g1,
            group2 = g2,
            prob1 = prob1,
            prob2 = prob2,
            pval1 = pval1,
            pval2 = pval2,
            prob_diff = prob1 - prob2,
            direction = paste(g1, "vs", g2),
            higher_group = ifelse(prob1 > prob2, g1, g2),
            stringsAsFactors = FALSE
          )
        }
      }
    }

    if (length(sig_details) > 0) {
      diff_df <- do.call(rbind, sig_details)
  
      # Only remain the results of higher_group (Run_CellChat_MaxGroup)
      diff_df <- diff_df[diff_df$higher_group %in% Run_CellChat_MaxGroup, , drop = FALSE]
  
      # If arguments "Run_CellChat_source_celltype" and "Run_CellChat_target_celltype" are set, filter the results.
      if (!is.null(Run_CellChat_source_celltype)) {
        diff_df <- diff_df[diff_df$source %in% Run_CellChat_source_celltype, , drop = FALSE]
      }
      if (!is.null(Run_CellChat_target_celltype)) {
        diff_df <- diff_df[diff_df$target %in% Run_CellChat_target_celltype, , drop = FALSE]
      }
  
      # Rank the absolute probability difference 
      diff_df <- diff_df[order(-abs(diff_df$prob_diff)), ]
  
      # Get the top signalings
      if (nrow(diff_df) > Run_CellChat_ntop_signaling) {
        diff_df <- diff_df[1:Run_CellChat_ntop_signaling, ]
      }

      # Get the ligand and receptor information from interaction column
      if ("interaction" %in% colnames(diff_df)) {
        diff_df$ligand <- sub("_.*$", "", diff_df$interaction)
        diff_df$receptor <- sub("^[^_]+_", "", diff_df$interaction)
      }
  
      # Output CSV file
      out_csv <- file.path(Run_CellChat_output_path, "multiGroup_significant_LR.csv")
      write.csv(diff_df, file = out_csv, row.names = FALSE)
      message("[Run_CellChat] Multi-group significant L-R pairs saved to: ", out_csv)
    } else {
      message("[Run_CellChat] No significant L-R pairs found across group comparisons.")
    }
}
  return(cellchat)
}
