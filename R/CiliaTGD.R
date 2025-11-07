#' CiliaTGD Pipeline
#'
#' Automates the segmentation of nuclei and cilia, interactive editing, and CiliaQ processing.
#'
#' @param MainDirectory Path to the main directory containing image subdirectories.
#' @param pixel_size Numeric. Microns per pixel.
#' @param nuc_length Numeric. Expected diameter of a nucleus in microns.
#' @param watershed Numeric. Watershed splitting parameter.
#' @param NucleusBorderRemoval Should nuclei touching the right/bottom border be removed from counting?
#' @param cilium_threshold_value Numeric. Threshold for initial cilium segmentation.
#' @param min_cilium_area Numeric. Minimum area for cilium filter.
#' @param background_subtraction Numeric. Brush size for background subtraction.
#' @param cilium_aspectratio Numeric. Aspect ratio threshold for cilia.
#' @param brightnessfactor Numeric. Value to increase brightness by for cilium mask image
#' @param CiliaQEditor Should CiliaQ editor step be skipped entirely?
#' @param CiliaQPreparator Should CiliaQPreparator be used for segmentation?
#' @param subtract_radius_val Radius value in pixels for subtraction brush if CiliaQPreparator is used
#' @param smooth_radius_val Radius value in pixels for smoothing brush if CiliaQPreparator is used
#' @param seg_method_val Segmentation method if CiliaQPreparator is used
#' @param SegmentationTester Should CiliaTGD iterate over various segmentation methods to help you identify the best approach for your data?
#' @param fiji_app_path Path to the Fiji.app bundle.
#' @export

  CiliaTGD <- function(
    MainDirectory = NULL,          # required
    pixel_size = NULL,             # required
    nuc_length = 12,               # defaulted
    watershed = 3,                 # defaulted
    NucleusBorderRemoval = TRUE,   # defaulted
    cilium_threshold_value = NULL, # required only if Preparator = FALSE 
    min_cilium_area = NULL,        # required
    background_subtraction = NULL, # required only if Preparator = FALSE
    cilium_aspectratio = NULL,      # required only if Preparator = FALSE
    brightnessfactor = 2.5,        # defaulted
    CiliaQPreparator = NULL,       # must be set TRUE/FALSE in normal runs
    subtract_radius_val = NULL,    # required always
    smooth_radius_val = NULL,      # required only if Preparator = TRUE
    seg_method_val = NULL,         # required only if Preparator = TRUE
    CiliaQEditor = TRUE,           # defaulted
    SegmentationTester = FALSE,    # controls behavior below
    fiji_app_path = NULL           # required
  ) {

    # ---------- validation & normalization ----------
    fail <- function(msg) stop(msg, call. = FALSE)
    need <- function(x, name) if (is.null(x)) fail(sprintf("`%s` must be set by the user.", name))
    
    # 1️⃣ Always required (for all modes)
    need(MainDirectory, "MainDirectory")
    need(pixel_size, "pixel_size")
    need(min_cilium_area, "min_cilium_area")
    need(fiji_app_path, "fiji_app_path")
    need(NucleusBorderRemoval, "NucleusBorderRemoval")
    
    # 2️⃣ Normalize paths
    MainDirectory <- normalizePath(path.expand(MainDirectory), winslash = "/", mustWork = TRUE)
    fiji_app_path <- normalizePath(path.expand(fiji_app_path), winslash = "/", mustWork = TRUE)
    
    # 3️⃣ Mode-specific validation
    if (isTRUE(SegmentationTester)) {
      message("ℹ️  SegmentationTester mode detected — skipping parameter completeness checks, as these will be iterated across.")
    } else {
      # Normal run: user must define whether Preparator is used
      if (is.null(CiliaQPreparator))
        fail("`CiliaQPreparator` must be explicitly set to TRUE or FALSE for normal runs. How do you want to segement cilia?")
      
      if (isTRUE(CiliaQPreparator)) {
        # Preparator mode → require these three
        need(subtract_radius_val, "subtract_radius_val")
        need(smooth_radius_val,   "smooth_radius_val")
        need(seg_method_val,      "seg_method_val")
      } else {
        # Non-preparator mode → require these three
        need(cilium_threshold_value, "cilium_threshold_value")
        need(background_subtraction, "background_subtraction")
        need(cilium_aspectratio,      "cilium_aspectratio")
      }
    }
    
    # 4️⃣ Numeric sanity checks
    num1 <- function(x, nm) if (!(is.numeric(x) && is.finite(x) && length(x)==1)) fail(sprintf("`%s` must be a single finite number.", nm))
    num1(pixel_size, "pixel_size")
    num1(nuc_length, "nuc_length"); if (nuc_length <= 0) fail("`nuc_length` must be > 0.")
    num1(watershed, "watershed"); if (watershed < 1) fail("`watershed` must be >= 1.")
    num1(min_cilium_area, "min_cilium_area"); if (min_cilium_area < 1) fail("`min_cilium_area` must be >= 1.")
    
    if (!isTRUE(SegmentationTester)) {
      if (isTRUE(CiliaQPreparator)) {
        num1(subtract_radius_val, "subtract_radius_val"); if (subtract_radius_val < 0) fail("`subtract_radius_val` must be >= 0.")
        num1(smooth_radius_val, "smooth_radius_val"); if (smooth_radius_val < 0) fail("`smooth_radius_val` must be >= 0.")
        if (!is.null(seg_method_val)) {
          ok <- c("RenyiEntropy","Otsu","Triangle","Default","Huang","Li","Mean","Minimum","Moments","Percentile","Yen")
          if (!(seg_method_val %in% ok))
            fail(sprintf("`seg_method_val` must be one of: %s", paste(ok, collapse=", ")))
        }
      } else {
        num1(cilium_threshold_value, "cilium_threshold_value")
        if (cilium_threshold_value < 0 || cilium_threshold_value > 255)
          fail("`cilium_threshold_value` must be in [0,255].")
        num1(background_subtraction, "background_subtraction")
        if (background_subtraction < 1)
          fail("`background_subtraction` must be >= 1.")
        num1(cilium_aspectratio, "cilium_aspectratio")
        if (cilium_aspectratio < 1)
          fail("`cilium_aspectratio` must be >= 1.")
      }
    }
    
    # 5️⃣ Integer coercions
    watershed <- as.integer(round(watershed))
    background_subtraction <- if (!is.null(background_subtraction)) as.integer(round(background_subtraction)) else background_subtraction
    subtract_radius_val <- if (!is.null(subtract_radius_val)) as.integer(round(subtract_radius_val)) else subtract_radius_val
    
    # ---------- end validation ----------
    
    
    
    
    # Ensure stringr functions are available
  library(stringr)
  # Redirect output to a log file
  logfile <- file.path(MainDirectory, "CiliaTGD_RunLog.txt")
  zz <- file(logfile, open = "wt")
  sink(zz, split = TRUE)
  on.exit({
    sink()
    close(zz)
  })
  
  
  
  # Start pipeline
  cat("\n############################\n")
  cat("# CiliaTGD Pipeline starting...\n")
  cat("############################\n")
  # Print parameters
  cat("\nParameters used are as follows:\n")
  cat("  - MainDirectory =", MainDirectory, "\n")
  cat("  - pixel_size =", pixel_size, "\n")
  cat("  - nuc_length =", nuc_length, "\n")
  cat("  - watershed =", watershed, "\n")
  cat("  - NucleusBorderRemoval =", NucleusBorderRemoval, "\n")
  cat("  - cilium_threshold_value =", cilium_threshold_value, "\n")
  cat("  - min_cilium_area =", min_cilium_area, "\n")
  cat("  - background_subtraction =", background_subtraction, "\n")
  cat("  - cilium_aspectratio =", cilium_aspectratio, "\n")
  cat("  - fiji_app_path =", fiji_app_path, "\n")
  cat("  - brightnessfactor =", brightnessfactor, "\n")
  cat("  - CiliaQPreparator =", CiliaQPreparator, "\n")
  cat("  - subtract_radius_val =", subtract_radius_val, "\n")
  cat("  - smooth_radius_val =", smooth_radius_val, "\n")
  cat("  - CiliaQEditor =", CiliaQEditor, "\n")
  cat("  - SegmentationTester =", SegmentationTester, "\n")
  cat("  - seg_method_val =", seg_method_val, "\n")
  
  
  # Script directory
  script_dir <- system.file("scripts", package = "CiliaTGD")
 
  RUN_NUCLEI <- TRUE

  # SegmentationTester toggles the fixer OFF (to keep the tester fast)
  SEGTEST_ACTIVE <- isTRUE(get0("SegmentationTester", ifnotfound = FALSE))
  
  # Default: ON for normal pipeline runs
  if (!exists("FIX_TIFF_METADATA", inherits = TRUE)) FIX_TIFF_METADATA <- TRUE
  if (SEGTEST_ACTIVE) FIX_TIFF_METADATA <- FALSE
  
  
  ################################################################################
  ################################################################################
  
  ################################################################################
  ### SegmentationTester (runs only if SegmentationTester = TRUE)
  ################################################################################
  
  
  if (isTRUE(SegmentationTester)) {
    cat("\n================ SegmentationTester ENABLED ================\n")
    
    segtests_dir_main <- file.path(MainDirectory, "SegmentationTests")
    dir.create(segtests_dir_main, recursive = TRUE, showWarnings = FALSE)
    
    .make_label <- function(x) {
      x <- lapply(x, function(v) if (is.factor(v)) as.character(v) else v)
      x <- lapply(x, as.character)
      paste(names(x), x, sep = "=", collapse = "__")
    }
    
    .sanitize_label <- function(x) gsub("[^A-Za-z0-9_=.-]+", "_", x)
    
    .purge_outputs <- function(MainDirectory) {
      subs <- list.dirs(MainDirectory, full.names=TRUE, recursive=FALSE)
      for (d in subs) {
        outd <- file.path(d, "output")
        if (dir.exists(outd)) unlink(outd, recursive=TRUE, force=TRUE)
      }
      unlink(file.path(MainDirectory, "CiliaQ-Editor_Files"), recursive=TRUE, force=TRUE)
      unlink(file.path(MainDirectory, "CiliaQ-Preparator_Files"), recursive=TRUE, force=TRUE)
      unlink(file.path(MainDirectory, "ProcessedImages.txt"), force=TRUE)
    }
    
    .collect_masks <- function(test_root, label, main_dir) {
      # Always copy curated outputs back to the real project SegmentationTests dir
      dest_root <- file.path(main_dir, "SegmentationTests")
      dir.create(dest_root, showWarnings = FALSE, recursive = TRUE)
      
      files <- list.files(
        test_root,
        pattern = "_ShiftedCiliumMask\\.tif$",
        recursive = TRUE,
        full.names = TRUE
      )
      if (!length(files)) {
        message("⚠️  No ShiftedCiliumMask files found for label: ", label)
        return(invisible(0L))
      }
      
      info <- file.info(files)
      latest <- files[which.max(info$mtime)]
      
      clean_label <- .sanitize_label(label)
      out_path <- file.path(dest_root, paste0(clean_label, ".tif"))
      
      ok <- file.copy(latest, out_path, overwrite = TRUE)
      if (!ok) warning("Copy failed for: ", latest, " -> ", out_path)
      
      cat(" 📦 Saved mask → ", out_path, "\n", sep = "")
      invisible(as.integer(ok))
    }
    
    .run_seg_test <- function(MainDirectory, overrides, note = NULL) {
      # choose first subdirectory from the real MainDirectory
      subdirs <- list.dirs(MainDirectory, full.names = TRUE, recursive = FALSE)
      if (!length(subdirs)) stop("No subdirectories found in MainDirectory")
      first_subdir <- subdirs[1]
      
      # temp root with exactly one subdir (what the pipeline expects)
      test_root <- file.path(tempdir(), "CiliaTGD_SegmentationTester_Temp")
      if (dir.exists(test_root)) unlink(test_root, recursive = TRUE, force = TRUE)
      dir.create(test_root, recursive = TRUE)
      test_subdir <- file.path(test_root, basename(first_subdir))
      dir.create(test_subdir, recursive = TRUE)
      file.copy(list.files(first_subdir, full.names = TRUE), test_subdir, recursive = TRUE)
      
      cat("\n🧪 Running test:", if (is.null(note)) .make_label(overrides) else note, "\n")
      
      .purge_outputs(test_root)
      
      # --- build a clean environment for the sourced script
      env <- new.env(parent = .GlobalEnv)
      env$plan           <- future::plan
      env$multisession   <- future::multisession
      env$future_lapply  <- future.apply::future_lapply
      env$with_progress  <- progressr::with_progress
      env$handlers       <- progressr::handlers
      env$RUN_NUCLEI <- FALSE
      env$FIX_TIFF_METADATA <- FALSE
      env$SEGTEST_ACTIVE <- TRUE
      
      # script dir
      env$script_dir <- system.file("scripts", package = "CiliaTGD")
      
      # inject parameters the pipeline reads (defaults from current scope)
      inherit <- parent.frame()
      param_names <- c(
        "MainDirectory","pixel_size","nuc_length","watershed", "NucleusBorderRemoval",
        "cilium_threshold_value","min_cilium_area","background_subtraction",
        "cilium_aspectratio","brightnessfactor","CiliaQPreparator",
        "subtract_radius_val","smooth_radius_val","seg_method_val",
        "CiliaQEditor","fiji_app_path"
      )
      for (nm in param_names) {
        if (exists(nm, envir = inherit, inherits = TRUE)) {
          assign(nm, get(nm, envir = inherit, inherits = TRUE), envir = env)
        }
      }
      
      # force test-specific values
      env$MainDirectory <- test_root
      env$CiliaQEditor  <- FALSE  # suppress Editor during tests
      
      # apply grid overrides
      for (nm in names(overrides)) assign(nm, overrides[[nm]], envir = env)
      
      # <<< NEW: load required packages inside the test env so plan(), handlers(), etc. exist
      evalq({
        suppressPackageStartupMessages({
          library(future)
          library(future.apply)
          library(progressr)
          library(dplyr)
          library(EBImage)
          library(abind)
        })
      }, envir = env)
      
      # run segmentation ONLY, inside this env
      processed_manifest <- file.path(test_root, "ProcessedImages.txt")
      source(file.path(env$script_dir, "TGD_detect_cilia_nuclei_pipeline.R"), local = env)
      
      # collect outputs from the test root
      label <- .make_label(overrides)
      .collect_masks(test_root, label, MainDirectory)
      
      # cleanup this test root
      unlink(test_root, recursive = TRUE, force = TRUE)
      invisible(label)
    }
    
    # GROUP A: normal (non-Preparator) parameters
    groupA <- expand.grid(
      cilium_threshold_value=c(5,10,15,20,25), 
      background_subtraction=c(5,10,15,20,25),
      cilium_aspectratio=c(1,1.5,2),
      KEEP.OUT.ATTRS=FALSE
    )
    
    baseA <- list(CiliaQPreparator=FALSE)
    for (i in seq_len(nrow(groupA))) {
      o <- c(baseA, as.list(groupA[i,]))
      .run_seg_test(MainDirectory, o, note=paste("GROUP A", .make_label(o)))
    }
    
    # GROUP B: Preparator parameters
    groupB <- expand.grid(
      seg_method_val=c("RenyiEntropy","Otsu","Triangle"), 
      smooth_radius_val=c(1,2,3), 
      subtract_radius_val=c(5,10,20,30),
      KEEP.OUT.ATTRS=FALSE
    )

    
    baseB <- list(CiliaQPreparator=TRUE)
    for (i in seq_len(nrow(groupB))) {
      o <- c(baseB, as.list(groupB[i,]))
      .run_seg_test(MainDirectory, o, note=paste("GROUP B", .make_label(o)))
    }
    
    cat("\n================ SegmentationTester Iterations COMPLETE ================\n")
    cat("All outputs saved under:\n  ", file.path(MainDirectory, "SegmentationTests"), "\n")
    cat("\nPreparing Segmentation Iteration Graph\n")
    
    
    # libs needed for plotting
    suppressPackageStartupMessages({
      if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
      library(ggplot2)
    })
    
    segtests_dir <- segtests_dir_main
    
    mask_files <- list.files(segtests_dir, pattern = "\\.tif$", full.names = TRUE)
    
    if (!length(mask_files)) {
      warning("No ShiftedCiliumMask files found in: ", segtests_dir)
    } else {
      # helper: choose the slice that most looks like a binary mask and count objects
      .count_cilia_in_file <- function(fp) {
        img <- EBImage::readImage(fp)
        d <- dim(img)
        
        # Prefer channel 2 explicitly; fall back if not present
        if (length(d) >= 3 && d[3] >= 2) {
          plane <- img[,,2]          # <-- mask channel
        } else if (length(d) >= 3 && d[3] == 1) {
          plane <- img[,,1]
        } else {
          # single-plane image fallback
          plane <- img
        }
        
        # binarize (your masks are 0/1 in EBImage scale 0..1)
        bin <- plane > 0.5
        
        # count connected components (exclude background 0)
        lab <- EBImage::bwlabel(bin)
        tbl <- table(lab)
        if ("0" %in% names(tbl)) tbl <- tbl[names(tbl) != "0"]
        count <- length(tbl)
        
        data.frame(image = basename(fp), cilia_count = count, stringsAsFactors = FALSE)
      }
      
      counts_df_list <- lapply(mask_files, .count_cilia_in_file)
      counts_df <- do.call(rbind, counts_df_list)
      
      # write text/CSV
      out_txt  <- file.path(segtests_dir, "cilia_counts.txt")
      out_csv  <- file.path(segtests_dir, "cilia_counts.csv")
      write.table(counts_df, file = out_txt, sep = "\t", row.names = FALSE, quote = FALSE)
      utils::write.csv(counts_df, file = out_csv, row.names = FALSE)
      
      # plot: bar chart, ordered by count desc, rotated x labels
      counts_df$image <- factor(counts_df$image, levels = counts_df$image[order(counts_df$cilia_count, decreasing = TRUE)])
      
      p <- ggplot(counts_df, aes(x = image, y = cilia_count)) +
        geom_col() +
        labs(
          x = "Image (parameter-labeled)",
          y = "Number of cilia",
          title = "Cilia counts per mask"
        ) +
        theme_minimal() +
        coord_cartesian(clip = "off") + 
        theme(
          axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 1),
          axis.title.x = element_text(margin = margin(t = 14))
        )
      
      out_png <- file.path(segtests_dir, "cilia_counts_barplot.png")
      ggsave(out_png, p, width = 20, height = 20, dpi = 150, bg = "white", limitsize = FALSE)
      cat(" 📊 Wrote counts to:\n  -", out_txt, "\n  -", out_csv, "\n")
      cat(" 🖼️  Saved bar plot to:\n  -", out_png, "\n")
    }
    
    cat("All outputs saved under:\n  ", segtests_dir_main, "\n")
    cat("✅ SegmentationTester finished. Inspect results and rerun normally.\n")
    return(invisible(NULL))
    
    cat("\n================ SegmentationTester COMPLETE ================\n")
    
  }
  

  
  ################################################################################
  ################################################################################

  
  
  
  
  
  # Path to manifest file
  processed_manifest <- file.path(MainDirectory, "ProcessedImages.txt")
  
  # Source the internal scripts
  
  # Start nucleus and cilium segmentation
  if (!file.exists(processed_manifest)) {
    cat("\n🟢🟢🟢 Running nucleus and cilium segmentation script\n")
    source(file.path(script_dir, "TGD_detect_cilia_nuclei_pipeline.R"), local = TRUE)
    cat(" ✔️ Completed nucleus and cilium segmentation\n")
  } else {
    cat("\n🔁 ProcessedImages.txt found — skipping segmentation\n")
  }
  source(file.path(script_dir, "CiliaQ_Editor_Script_pipeline.R"), local = TRUE)
  cat(" ✔️ Completed CiliaQ Editor\n")
  source(file.path(script_dir, "CiliaQ_Processor_Script_pipeline.R"), local = TRUE)
  cat(" ✔️ Completed CiliaQ Analysis\n")
  
  cat("############################\n")
  cat("# CiliaTGD Pipeline complete!\n")
  cat("############################\n")
}
