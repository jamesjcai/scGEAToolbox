# functions/Vina_Docking.R

#' Run AutoDock Vina Molecular Docking (Module 3 of scDock)
#'
#' Bridges cell-cell communication (CCC) results from CellChat (Module 2) to
#' structure-based virtual screening. It extracts the top ligand-receptor gene
#' pairs, maps them to 3D protein structures (RCSB PDB or AlphaFold), downloads
#' and preprocesses those structures, downloads a small-molecule drug library
#' (FDA-approved or user-supplied CAS numbers), and runs AutoDock Vina for every
#' (protein target, drug) pair. Results are ranked by binding affinity and
#' written to CSV files alongside the docked pose PDBQT files.
#'
#' @section Pipeline stages:
#' \enumerate{
#'   \item Parse LR pairs from CellChat output CSVs.
#'   \item Convert mouse gene symbols to human via \pkg{homologene}.
#'   \item Match gene names to PDB/AlphaFold IDs using reference CSVs.
#'   \item Download protein structures (RCSB or AlphaFold) and preprocess to
#'         PDBQT using AutoDockTools (\code{prepare_receptor.py}).
#'   \item Download small-molecule structures from PubChem via CAS numbers
#'         (\code{download_cas_pubchem.py}).
#'   \item Run Vina all-vs-all (each drug against each protein target) and
#'         collect per-receptor ranked score tables.
#' }
#'
#' @section Grid box requirement:
#' Each receptor PDBQT file must have a companion \code{<name>_grid.txt} in the
#' same directory, containing key=value lines for \code{center_x}, \code{center_y},
#' \code{center_z}, \code{size_x}, \code{size_y}, \code{size_z}. This file
#' defines the docking search box and must be prepared by the user in advance.
#'
#' @section Terminology note:
#' "Ligand" is used in two senses throughout: CellChat \emph{signaling ligands}
#' (proteins, e.g. VEGFA) are treated as docking \emph{protein targets}, alongside
#' the signaling receptors. Small-molecule \emph{drug ligands} are the compounds
#' docked against those protein targets.
#'
#' @param Run_CellChat_output_path Character. Path to the directory containing
#'   CellChat output CSVs. Expects files matching \code{top\d+_LR_.*\.csv}
#'   (single-group) or \code{multiGroup_significant_LR.csv} (multi-group).
#' @param Vina_Docking_output_path Character. Path to the root output directory.
#'   Created if it does not exist.
#' @param Vina_Docking_ligand_ref_file Character. Path to a CSV with columns
#'   \code{protein_name} and \code{PDB_model} mapping signaling-ligand gene
#'   names to PDB/AlphaFold IDs (semicolon-separated for multiple structures).
#' @param Vina_Docking_receptor_ref_file Character. Same format as
#'   \code{Vina_Docking_ligand_ref_file} but for signaling-receptor proteins.
#' @param Vina_Docking_cas_txt_file Character or \code{NULL}. Path to a plain-text
#'   file of CAS numbers (one per line) for the user-supplied compound library.
#'   Ignored when \code{Vina_Docking_use_fda = TRUE}.
#' @param Vina_Docking_use_fda Logical. If \code{TRUE}, use the FDA-approved drug
#'   library defined in \code{Vina_Docking_fda_txt} instead of CAS numbers.
#' @param Vina_Docking_fda_txt Character or \code{NULL}. Path to a plain-text file
#'   of FDA drug CAS numbers. Required when \code{Vina_Docking_use_fda = TRUE}.
#' @param Vina_Docking_docking_ligand_dir Character or \code{NULL}. Path to a
#'   directory of pre-built drug PDBQT files. If provided, the compound download
#'   step is skipped and these files are used directly.
#' @param Vina_Docking_docking_receptor_dir Character or \code{NULL}. Path to a
#'   directory of pre-built protein target PDBQT files (with companion
#'   \code{_grid.txt} files). If provided, protein download and preprocessing are
#'   skipped entirely.
#' @param Vina_Docking_vina_exhaustiveness Integer. Vina \code{--exhaustiveness}
#'   parameter (search thoroughness; default in Vina is 8).
#' @param Vina_Docking_vina_num_modes Integer. Vina \code{--num_modes} parameter
#'   (maximum number of binding modes to output).
#' @param Vina_Docking_vina_seed Integer. Vina \code{--seed} for reproducibility.
#' @param Vina_Docking_vina_cpu Integer. Number of CPU cores for Vina
#'   (\code{--cpu}).
#'
#' @return Called for side effects. Writes to \code{Vina_Docking_output_path}:
#' \describe{
#'   \item{\code{ligands_with_PDB.csv}}{Signaling-ligand proteins matched to PDB IDs.}
#'   \item{\code{receptors_with_PDB.csv}}{Signaling-receptor proteins matched to PDB IDs.}
#'   \item{\code{CellChat_ligand/}}{Downloaded and preprocessed protein-ligand PDBQT files.}
#'   \item{\code{CellChat_receptor/}}{Downloaded and preprocessed protein-receptor PDBQT files.}
#'   \item{\code{<receptor_subdir>/<drug>_result_structure.pdbqt}}{Docked poses.}
#'   \item{\code{<receptor_subdir>/<drug>_result_score.txt}}{Full Vina scoring log.}
#'   \item{\code{<receptor_subdir>/AutoDockVina_score.csv}}{Ranked drug candidates by
#'         best binding affinity (kcal/mol, ascending) for each protein target.}
#' }
#'
#' @section Dependencies:
#' R packages: \pkg{dplyr}, \pkg{homologene}.
#' System tools (must be on PATH): \code{vina} (AutoDock Vina), \code{python3}.
#' Python scripts (in \code{functions/}): \code{download_pdb_chains_from_csv.py},
#' \code{download_alphafold.py}, \code{prepare_receptor.py},
#' \code{download_cas_pubchem.py}.
#'
#' @examples
#' \dontrun{
#' Vina_Docking(
#'   Run_CellChat_output_path       = "output/cellchat/",
#'   Vina_Docking_output_path       = "output/docking/",
#'   Vina_Docking_ligand_ref_file   = "ref/ligand_pdb_ref.csv",
#'   Vina_Docking_receptor_ref_file = "ref/receptor_pdb_ref.csv",
#'   Vina_Docking_cas_txt_file      = NULL,
#'   Vina_Docking_use_fda           = TRUE,
#'   Vina_Docking_fda_txt           = "ref/fda_cas.txt",
#'   Vina_Docking_docking_ligand_dir  = NULL,
#'   Vina_Docking_docking_receptor_dir = NULL,
#'   Vina_Docking_vina_exhaustiveness = 8,
#'   Vina_Docking_vina_num_modes      = 9,
#'   Vina_Docking_vina_seed           = 42,
#'   Vina_Docking_vina_cpu            = 4
#' )
#' }
Vina_Docking <- function(Run_CellChat_output_path,
                         Vina_Docking_output_path,
                         Vina_Docking_ligand_ref_file,
                         Vina_Docking_receptor_ref_file,
                         Vina_Docking_cas_txt_file,
                         Vina_Docking_use_fda,
                         Vina_Docking_fda_txt,
                         Vina_Docking_docking_ligand_dir,
                         Vina_Docking_docking_receptor_dir,
                         Vina_Docking_vina_exhaustiveness,
                         Vina_Docking_vina_num_modes,
                         Vina_Docking_vina_seed,
                         Vina_Docking_vina_cpu) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("[Vina_Docking] Package 'dplyr' is required.")
  suppressPackageStartupMessages({
  library(dplyr)
})

  # grid.txt read function
  read_grid_file <- function(grid_file) {
    lines <- readLines(grid_file)
    vals <- sapply(lines, function(x) as.numeric(strsplit(x, "=")[[1]][2]))
    names(vals) <- sapply(lines, function(x) trimws(strsplit(x, "=")[[1]][1]))
    return(vals)
  }

  # Vina docking running function
  run_vina_docking <- function(ligand_file, receptor_file, out_dir, params) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    ligand_name <- tools::file_path_sans_ext(basename(ligand_file))
    output_pdbqt <- file.path(out_dir, paste0(ligand_name, "_result_structure.pdbqt"))
    log_file <- file.path(out_dir, paste0(ligand_name, "_result_score.txt"))

    grid_file <- sub("\\.pdbqt$", "_grid.txt", receptor_file)
    if (!file.exists(grid_file)) {
      stop("[Vina] Missing grid file: ", grid_file, " for receptor: ", receptor_file)
    }
    grid <- read_grid_file(grid_file)

    cmd <- c("--ligand", ligand_file,
             "--receptor", receptor_file,
             "--out", output_pdbqt,
             "--log", log_file,
             "--exhaustiveness", as.character(params$exhaustiveness),
             "--num_modes", as.character(params$num_modes),
             "--seed", as.character(params$seed),
             "--cpu", as.character(params$cpu),
             "--center_x", as.character(grid["center_x"]),
             "--center_y", as.character(grid["center_y"]),
             "--center_z", as.character(grid["center_z"]),
             "--size_x", as.character(grid["size_x"]),
             "--size_y", as.character(grid["size_y"]),
             "--size_z", as.character(grid["size_z"])
    )
    message("")
    message("[Vina] Docking ligand ", ligand_name, " with receptor ", basename(receptor_file))
    status <- system2("vina", args = cmd)
    if (status != 0) {
      warning("[Vina] Docking failed for ligand ", ligand_file, " with receptor ", receptor_file)
    }
  }

  # Set Vina arguments
  vina_params <- list(
    exhaustiveness = Vina_Docking_vina_exhaustiveness,
    num_modes = Vina_Docking_vina_num_modes,
    seed = Vina_Docking_vina_seed,
    cpu = Vina_Docking_vina_cpu
  )

  # Read LR CSV
  csv_files <- list.files(Run_CellChat_output_path, pattern = "^top\\d+_LR_.*\\.csv$", full.names = TRUE)

  if (length(csv_files) == 0) {
    # fallback: try multi-group file
    multi_csv <- file.path(Run_CellChat_output_path, "multiGroup_significant_LR.csv")
    if (file.exists(multi_csv)) {
      message("[INFO] No top*_LR_*.csv found. Using multiGroup_significant_LR.csv instead.")
      csv_files <- c(multi_csv)
    } else {
      stop("No LR CSV files found (neither top*_LR_*.csv nor multiGroup_significant_LR.csv) in: ",
           Run_CellChat_output_path)
    }
  }

  ligands <- receptors <- c()
  for (csv in csv_files) {
    df <- tryCatch(read.csv(csv, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(df)) next

    # Ensure ligand/receptor columns exist
    if (!all(c("ligand","receptor") %in% colnames(df)) && "interaction" %in% colnames(df)) {
      parts <- strsplit(df$interaction, "_")
      df$ligand <- sapply(parts, `[`, 1)
      df$receptor <- sapply(parts, `[`, 2)
    }

    if (all(c("ligand","receptor") %in% colnames(df))) {
      ligands <- c(ligands, df$ligand)
      receptors <- c(receptors, df$receptor)
    } else {
      warning("[Vina] CSV ", basename(csv), " does not contain ligand/receptor information.")
    }
  }
  ligands <- unique(ligands)
  receptors <- unique(receptors)
  
  # Gene name conversion: Mouse -> Human (offline via homologene)
  if (!requireNamespace("homologene", quietly = TRUE)) {
    stop("Package 'homologene' is required for gene name conversion.")
  }
  suppressPackageStartupMessages({
  library(homologene)
})

  convert_mouse_to_human <- function(gene_vec) {
    gene_vec <- unique(gene_vec)
    gene_map <- homologene::homologene(gene_vec, inTax = 10090, outTax = 9606)
    colnames(gene_map) <- c("MouseGene", "HumanGene")
    map_dict <- setNames(gene_map$HumanGene, gene_map$MouseGene)

    converted <- ifelse(gene_vec %in% names(map_dict),
                        map_dict[gene_vec],
                        gene_vec)
    converted <- converted[!is.na(converted) & converted != ""]
    return(unique(converted))
  }

  ligands <- convert_mouse_to_human(ligands)
  receptors <- convert_mouse_to_human(receptors)
  
  # Load reference
  ligand_ref <- read.csv(Vina_Docking_ligand_ref_file, stringsAsFactors = FALSE)
  receptor_ref <- read.csv(Vina_Docking_receptor_ref_file, stringsAsFactors = FALSE)
  ligand_match <- ligand_ref %>% filter(protein_name %in% ligands)
  receptor_match <- receptor_ref %>% filter(protein_name %in% receptors)

  if (!dir.exists(Vina_Docking_output_path)) dir.create(Vina_Docking_output_path, recursive = TRUE)
  write.csv(ligand_match, file.path(Vina_Docking_output_path, "ligands_with_PDB.csv"), row.names = FALSE)
  write.csv(receptor_match, file.path(Vina_Docking_output_path, "receptors_with_PDB.csv"), row.names = FALSE)

  # Create PDB directory
  ligand_dir <- file.path(Vina_Docking_output_path, "CellChat_ligand")
  receptor_dir <- file.path(Vina_Docking_output_path, "CellChat_receptor")
  if (!dir.exists(ligand_dir)) dir.create(ligand_dir, recursive = TRUE)
  if (!dir.exists(receptor_dir)) dir.create(receptor_dir, recursive = TRUE)

  message("[INFO] Downloading ligand structures from RCSB...")
  system2("python3", args = c(
    normalizePath("functions/download_pdb_chains_from_csv.py"),
    file.path(Vina_Docking_output_path, "ligands_with_PDB.csv"),
    ligand_dir
    ))

  if (is.null(Vina_Docking_docking_receptor_dir)) {
  message("[INFO] Downloading receptor structures from RCSB...")
  system2("python3", args = c(
    normalizePath("functions/download_pdb_chains_from_csv.py"),
    file.path(Vina_Docking_output_path, "receptors_with_PDB.csv"),
    receptor_dir
    ))
  }

  # Download protein structure
  download_structures <- function(df, outdir, label) {
  for (raw_ids in df$PDB_model) {
    if (is.na(raw_ids) || raw_ids == "") next
    ids <- trimws(unlist(strsplit(raw_ids, ";")))
    for (id in ids) {
      if (id == "") next
      subdir <- file.path(outdir, id)
      if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)

      outfile <- NULL

      if (grepl("^AF-", id)) {
        # AlphaFold structure
        outfile <- file.path(subdir, paste0(id, ".pdb"))
        if (!file.exists(outfile)) {
          message("[", label, "] Downloading AlphaFold model: ", id)
          system2("python3", args = c(
            normalizePath("functions/download_alphafold.py"),
            id, subdir
          ))
        }
      } else if (grepl("^[0-9][A-Za-z0-9]{3}\\.[A-Z]$", id)) {
        parts <- strsplit(id, "\\.")[[1]]
        pdb_id <- toupper(parts[1])
        chain <- toupper(parts[2])
        outfile <- file.path(subdir, paste0(pdb_id, "_", chain, ".pdb"))

        if (!file.exists(outfile)) {
          message("[", label, "] Downloading PDB chain from RCSB: ", pdb_id, " chain ", chain)
          status <- system2("python3", args = c(
            normalizePath("functions/download_pdb_chains_from_csv.py"),
            tempfile_csv <- tempfile(fileext = ".csv"),
            subdir
          ))

          write.csv(data.frame(PDB_model = id), tempfile_csv, row.names = FALSE)

          system2("python3", args = c(normalizePath("functions/download_pdb_chains_from_csv.py"),
                                      tempfile_csv, subdir))

          if (!file.exists(outfile)) {
            af_id <- paste0("AF-", pdb_id, "-", chain)
            message("[", label, "] RCSB not found, fallback to AlphaFold: ", af_id)
            outfile <- file.path(subdir, paste0(af_id, ".pdb"))
            system2("python3", args = c(
              normalizePath("functions/download_alphafold.py"),
              af_id, subdir
            ))
          }
        }
      } else warning("[", label, "] Unrecognized ID format: ", id)

      # preprocessing with AutoDockTools
      if (!is.null(outfile) && file.exists(outfile)) {
        message("[", label, "] Preprocessing with AutoDockTools: ", outfile)
        system2("python3", args = c(
          normalizePath("functions/prepare_receptor.py"),
          outfile, subdir
        ))
      }
    }
  }
}


  download_structures(ligand_match, ligand_dir, "ligand")

  # Only download receptors if user did NOT provide Vina_Docking_docking_receptor_dir
  if (is.null(Vina_Docking_docking_receptor_dir)) {
    download_structures(receptor_match, receptor_dir, "receptor")
  }

  # Ligand structure source
  cas_dir <- NULL

  if (Vina_Docking_use_fda) {
    # FDA compounds support
    if (is.null(Vina_Docking_fda_txt) || !file.exists(Vina_Docking_fda_txt)) {
      stop("[FDA] Vina_Docking_use_fda = TRUE but Vina_Docking_fda_txt is missing or does not exist!")
    }

    fda_dir <- file.path(Vina_Docking_output_path, "fda_results")
    if (!dir.exists(fda_dir)) dir.create(fda_dir, recursive = TRUE)

    fda_cas <- trimws(readLines(Vina_Docking_fda_txt))
    if (length(fda_cas) > 0) {
      message("[FDA] Downloading structures for ", length(fda_cas), " FDA compounds...")
      system2("python3", args = c(
        normalizePath("functions/download_cas_pubchem.py"),
        Vina_Docking_fda_txt, fda_dir
      ))
      Vina_Docking_docking_ligand_dir <- fda_dir
    } else {
      stop("[FDA] fda.txt is empty! Please provide CAS numbers.")
    }

  } else {
    # CAS compounds support
    if (!is.null(Vina_Docking_cas_txt_file) && file.exists(Vina_Docking_cas_txt_file)) {
      cas_dir <- file.path(Vina_Docking_output_path, "ligand_structures_from_CAStxt_for_AutoDockVina")
      if (!dir.exists(cas_dir)) dir.create(cas_dir)
      cas_numbers <- trimws(readLines(Vina_Docking_cas_txt_file))
      message("[CAS] Downloading structures for ", length(cas_numbers), " compounds...")
      system2("python3", args = c(
        normalizePath("functions/download_cas_pubchem.py"),
        Vina_Docking_cas_txt_file, cas_dir
      ))
      Vina_Docking_docking_ligand_dir <- cas_dir
    }
  }

  # AutoDock Vina docking
  Vina_Docking_docking_ligand_dir <- if (is.null(Vina_Docking_docking_ligand_dir)) cas_dir else Vina_Docking_docking_ligand_dir

  if (is.null(Vina_Docking_docking_receptor_dir)) {
    Vina_Docking_docking_receptor_dir <- list(ligand_dir, receptor_dir)
  } else {
    message("[INFO] Using user-provided receptor directory: ", Vina_Docking_docking_receptor_dir)
  }

  receptor_files <- unlist(lapply(Vina_Docking_docking_receptor_dir, function(d) {
    list.files(d, pattern = "_prepared\\.pdbqt$", full.names = TRUE, recursive = TRUE)
  }))

  ligand_files <- list.files(Vina_Docking_docking_ligand_dir, pattern = "\\.pdbqt$", full.names = TRUE, recursive = TRUE)
  ligand_files <- ligand_files[file.size(ligand_files) > 0]

  for (receptor in receptor_files) {
    rec_subdir <- dirname(receptor)
    for (ligand in ligand_files) {
      run_vina_docking(
        ligand_file = ligand,
        receptor_file = receptor,
        out_dir = rec_subdir,
        params = vina_params
      )
    }

    # Collect docking results per receptor
    message("[Vina] Collecting docking results for receptor folder: ", rec_subdir)
    log_files <- list.files(rec_subdir,
                            pattern = "_result_score.txt$",
                            full.names = TRUE,
                            recursive = FALSE)

    results <- lapply(log_files, function(logf) {
      lines <- readLines(logf, warn = FALSE)
      affinity <- NA
      # Get the first mode 
      line1 <- grep("^\\s*1\\s", lines, value = TRUE)
      if (length(line1) > 0) {
        affinity <- as.numeric(strsplit(trimws(line1), "\\s+")[[1]][2])
      }
      ligand <- sub("_result_score.txt$", "", basename(logf))
      data.frame(ligand = ligand,
                 affinity = affinity,
                 stringsAsFactors = FALSE)
    })

    if (length(results) > 0) {
      results_df <- do.call(rbind, results) %>%
        arrange(affinity)
      out_csv <- file.path(rec_subdir, "AutoDockVina_score.csv")
      write.csv(results_df, out_csv, row.names = FALSE)
      message("[Vina] Receptor score summary saved to: ", out_csv)
    }
  }
}