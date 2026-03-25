# export_CellChatDB.R
#
# Export CellChatDB ligand-receptor interactions to tab-separated text files
# compatible with scGEAToolbox_dev and sc_dock_ccc.m.
#
# Outputs (written to the same directory as this script):
#   CellChatDB_human_LR.txt  - human LR pairs with pathway annotation
#   CellChatDB_mouse_LR.txt  - mouse LR pairs with pathway annotation
#
# Requirements:
#   install.packages("BiocManager")
#   BiocManager::install("CellChat")   # or: devtools::install_github("jinworks/CellChat")
#
# Usage:
#   Rscript export_CellChatDB.R
#   -- or run interactively in RStudio --

library(CellChat)

out_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
# Fallback if not in RStudio
if (!nzchar(out_dir) || out_dir == ".") out_dir <- getwd()

export_lr <- function(db, species_tag) {
  int_df <- db$interaction

  # Core columns always present in CellChatDB
  cols_required <- c("ligand", "receptor", "pathway_name", "annotation")
  cols_keep <- intersect(cols_required, colnames(int_df))

  # Some versions store multi-subunit receptors as receptor.1 / receptor.2 etc.
  # Flatten to a single representative receptor (first subunit)
  if (!"receptor" %in% colnames(int_df) && "receptor.1" %in% colnames(int_df)) {
    int_df$receptor <- int_df$receptor.1
    cols_keep <- c("ligand", "receptor",
                   intersect(c("pathway_name","annotation"), colnames(int_df)))
  }

  out <- int_df[, cols_keep, drop = FALSE]

  # Standardise column names to match sc_dock_ccc expectation
  names(out)[names(out) == "pathway_name"] <- "pathway"
  names(out)[names(out) == "annotation"]   <- "interaction_type"

  # Remove rows with missing ligand or receptor
  out <- out[nzchar(out$ligand) & nzchar(out$receptor), ]

  fname <- file.path(out_dir, paste0("CellChatDB_", species_tag, "_LR.txt"))
  write.table(out, fname, sep = "\t", quote = FALSE, row.names = FALSE)
  message(sprintf("Wrote %d LR pairs to %s", nrow(out), fname))
  invisible(out)
}

# Human
export_lr(CellChatDB.human, "human")

# Mouse
export_lr(CellChatDB.mouse, "mouse")

message("Done. Load in MATLAB with:")
message("  T = readtable('CellChatDB_human_LR.txt', 'FileType','text', 'Delimiter','\\t', 'TextType','string');")
