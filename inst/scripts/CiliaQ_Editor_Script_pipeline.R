################################################################################
###   Launch ImageJ to process images interactively (or bypass if requested)
################################################################################

# Expect these to be defined upstream:
#   MainDirectory  <- "/path/to/project"
#   fiji_app_path  <- "/Applications/Fiji.app"
#   CiliaQEditor   <- TRUE/FALSE

fiji_app_path <- path.expand(fiji_app_path)

################################################################################
# Quit Fiji (cross-platform)
quitFiji <- function() {
  sys <- Sys.info()[["sysname"]]
  if (sys == "Darwin") {
    system2("osascript", args = c("-e", 'tell application "Fiji" to quit'), stdout = FALSE, stderr = FALSE)
    system2("pkill", args = c("-f", "Fiji.app"), wait = FALSE)
    system2("pkill", args = c("-f", "Fiji"), wait = FALSE)
  } else if (sys == "Windows") {
    system2("taskkill", args = c("/IM", "Fiji.exe", "/F"), stdout = FALSE, stderr = FALSE)
  } else {
    system("pkill -f Fiji", wait = FALSE)
  }
}

################################################################################
# Set up directories and paths
input_dir <- file.path(MainDirectory, "CiliaQ-Editor_Files")
if (!dir.exists(input_dir)) {
  dir.create(input_dir, recursive = TRUE)
  cat("Cannot find", input_dir, "containing segmented images\n")
  stop("Exiting")
} else {
  cat(" - Directory containing segmented images is", input_dir, "\n")
}

# Full list of images to process
files_to_process <- list.files(
  path      = input_dir,
  pattern   = "_ShiftedCiliumMask\\.tif$",
  recursive = TRUE,
  full.names = TRUE
)

# Define checkpoint log files
processed_log  <- file.path(MainDirectory, "ProcessedImages.txt")
file_list_path <- file.path(input_dir, "ImagesToEdit.txt")

# Remove already-processed images
if (file.exists(processed_log)) {
  processed_images <- readLines(processed_log, warn = FALSE)
  files_to_process <- setdiff(files_to_process, processed_images)
  writeLines(files_to_process, file_list_path)
  writeLines(processed_images, processed_log)
  cat(" 🔁 Resuming: Skipping", length(processed_images), "already-processed images\n")
} else {
  processed_images <- character(0)
  writeLines(files_to_process, file_list_path)
  writeLines(processed_images, processed_log)
}

# Predefine flag path (used in both branches)
done_flag <- file.path(input_dir, "ciliaq_editor_complete.flag")

# Check if anything remains
if (length(files_to_process) == 0) {
  cat(" ✅ All images already processed. Nothing to do.\n")
  
} else {
  
  if (CiliaQEditor == FALSE) {
    ############################################################################
    # BYPASS FIJI: copy masks as-is but with `_ed` suffix into output dir
    ############################################################################
    dest_dir <- file.path(MainDirectory, "CiliaQ_Files")
    if (!dir.exists(dest_dir)) {
      dir.create(dest_dir, recursive = TRUE)
      cat(" - Created output directory:", dest_dir, "\n")
    } else {
      cat(" - Output directory already exists:", dest_dir, "\n")
    }
    
    # Build destination filenames with _ed suffix
    dest_basenames <- sub("_ShiftedCiliumMask\\.tif$",
                          "_ShiftedCiliumMask_ed.tif",
                          basename(files_to_process))
    dest_paths <- file.path(dest_dir, dest_basenames)
    
    # Copy (not move) into destination with new names
    copied_ok <- file.copy(from = files_to_process,
                           to   = dest_paths,
                           overwrite = FALSE)
    
    # Report
    cat(" - Attempted to copy", length(files_to_process), "files\n")
    if (any(copied_ok)) {
      cat("   • Copied", sum(copied_ok), "files to", dest_dir, "\n")
    }
    if (any(!copied_ok)) {
      warning("Some files were not copied (already exist or error):\n",
              paste(basename(dest_paths[!copied_ok]), collapse = "\n"))
    }
    
    # Mark as processed in the log (append)
    try({
      con <- file(processed_log, open = "a")
      writeLines(files_to_process, con)
      close(con)
    }, silent = TRUE)
    
    # Clear the edit list since we're done
    suppressWarnings(file.remove(file_list_path))
    
    # Create a flag for downstream cleanup (mirrors the Fiji path)
    file.create(done_flag)
    
    cat(" ✔️ Skipped Fiji; copied masks with _ed suffix.\n")
    
  } else {
    ############################################################################
    # RUN FIJI INTERACTIVELY
    ############################################################################
    cat(" - Found", length(files_to_process), "image files to process\n")
    
    # Prepare the macro to batch process files in Fiji
    input_dir_norm <- path.expand(input_dir)
    ijm_dir  <- gsub("\\\\", "/", normalizePath(input_dir_norm, winslash = "/"))
    macro_path <- file.path(input_dir_norm, "run_ciliaq_editor.ijm")
    
    # Macro script (quotes/escapes balanced for R + IJM)
    macro_code <- sprintf(
      'dir = "%s";
listPath = dir + "/ImagesToEdit.txt";
logRoot  = "%s";

lines = split(File.openAsString(listPath), "\\n");

for (i = 0; i < lines.length; i++) {
    file = lines[i];
    if (endsWith(file, "ShiftedCiliumMask.tif") && File.exists(file)) {
        print("Opening: " + file);
        open(file);
        run("CiliaQ Editor v0.0.3", "channel=2 channel_0=3 output=[save as filename + suffix \'_ed\']");
        waitForUser("Edit complete. Click OK to continue.");
        close();
        // Append processed original path to ProcessedImages.txt at project root
        File.append(file + "\\n", logRoot + "/ProcessedImages.txt");
    }
}
File.saveString("done", dir + "/ciliaq_editor_complete.flag");',
  gsub("\\\\", "/", ijm_dir),
  gsub("\\\\", "/", normalizePath(MainDirectory, winslash = "/"))
    )
    
    writeLines(macro_code, con = macro_path)
    
    # Launch Fiji with macro
    cat(
      " 🚀 Launching Fiji CiliaQ Processor macro\n",
      "    (If the script fails at any point during editing,\n",
      "     abort the script in R and rerun it using the\n",
      "     same parameters to pick up where you left off.)\n",
      sep = ""
    )
    
    cmd <- sprintf('open -n -a "%s" --args -macro "%s"', fiji_app_path, macro_path)
    system(cmd)
    
    # Wait for completion flag
    cat(" ⏳ Waiting for you to finish editing cilia in Fiji...\n")
    while (!file.exists(done_flag)) {
      Sys.sleep(2)
    }
    quitFiji()
    cat(" ✔️ Fiji CiliaQ Editor macro complete.\n")
  }
}

################################################################################
# Move edited files to output directory (only when Fiji was used)
if (isTRUE(CiliaQEditor)) {
  dest_dir <- file.path(MainDirectory, "CiliaQ_Files")
  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE)
    cat(" - Created output directory:", dest_dir, "\n")
  } else {
    cat(" - Output directory already exists:", dest_dir, "\n")
  }
  
  edited_files <- list.files(input_dir, pattern = "_ed\\.tif$", full.names = TRUE)
  if (length(edited_files) > 0) {
    moved <- file.rename(edited_files, file.path(dest_dir, basename(edited_files)))
    if (any(!moved)) {
      warning("Some edited files failed to move:\n", paste(edited_files[!moved], collapse = "\n"))
    } else {
      cat(" - Moved", length(edited_files), "files to", dest_dir, "\n")
    }
  } else {
    cat(" - No _ed files found to move (did you save edits?).\n")
  }
} else {
  cat(" - Move step skipped (no Fiji run; files were already copied to output).\n")
}

################################################################################
# Optional cleanup
suppressWarnings(file.remove(file_list_path))
suppressWarnings(file.remove(processed_log))
if (file.exists(done_flag)) file.remove(done_flag)

cat(" ✅ Done.\n")

  