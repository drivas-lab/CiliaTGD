
################################################################################
###    SECTION 2: Define Functions:
################################################################################
  
  # Function to shift ciliary mask a bit
  shift_image <- function(img, x_shift = 0, y_shift = 0) {
    shifted <- EBImage::Image(0, dim = dim(img))
    
    # Calculate source and destination ranges
    x_src <- (1 + max(0, x_shift)):(ncol(img) + min(0, x_shift))
    y_src <- (1 + max(0, y_shift)):(nrow(img) + min(0, y_shift))
    
    x_dest <- (1 + max(0, -x_shift)):(ncol(img) + min(0, -x_shift))
    y_dest <- (1 + max(0, -y_shift)):(nrow(img) + min(0, -y_shift))
    
    # Copy pixels
    shifted[y_dest, x_dest] <- img[y_src, x_src]
    
    return(shifted)
  }
  
  
  .safe_dir <- function(path) {
    d <- dirname(path)
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
    invisible(d)
  }
  
  .safe_write_tiff <- function(img, path, bits = 8L, type = "tiff") {
    .safe_dir(path)
    
    # Try #1: direct write
    ok <- tryCatch({
      if (file.exists(path)) unlink(path, force = TRUE)
      EBImage::writeImage(img, files = path, bits.per.sample = bits, type = type)
      TRUE
    }, error = function(e) FALSE)
    
    if (ok) return(invisible(path))
    
    # Try #2: write to a temp file next to target, then rename
    tmp_local <- paste0(path, ".tmpwrite.tif")
    ok <- tryCatch({
      if (file.exists(tmp_local)) unlink(tmp_local, force = TRUE)
      EBImage::writeImage(img, files = tmp_local, bits.per.sample = bits, type = type)
      if (!file.rename(tmp_local, path)) file.copy(tmp_local, path, overwrite = TRUE)
      if (file.exists(tmp_local)) unlink(tmp_local, force = TRUE)
      TRUE
    }, error = function(e) FALSE)
    
    if (ok) return(invisible(path))
    
    # Try #3: write to OS tempdir, then copy over
    tmp_os <- file.path(tempdir(), paste0("cilia_", basename(path)))
    ok <- tryCatch({
      if (file.exists(tmp_os)) unlink(tmp_os, force = TRUE)
      EBImage::writeImage(img, files = tmp_os, bits.per.sample = bits, type = type)
      file.copy(tmp_os, path, overwrite = TRUE)
      if (file.exists(tmp_os)) unlink(tmp_os, force = TRUE)
      TRUE
    }, error = function(e) FALSE)
    
    if (!ok) stop("Failed to write TIFF after 3 attempts: ", path)
    invisible(path)
  }

  
  # Ensure Fiji is truly not running
  waitForFijiToExit <- function(timeout_sec = 8) {
    sys <- Sys.info()[["sysname"]]
    t0 <- Sys.time()
    repeat {
      alive <- switch(sys,
                      "Darwin"  = system2("pgrep", c("-f", "ImageJ-macosx"), stdout = TRUE, stderr = FALSE),
                      "Linux"   = system("pgrep -f ImageJ", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0,
                      "Windows" = FALSE  # crude; ImageJ.exe check could be added for Windows
      )
      if (identical(alive, character(0)) || isFALSE(alive)) break
      if (as.numeric(difftime(Sys.time(), t0, units = "secs")) > timeout_sec) break
      Sys.sleep(0.25)
    }
  }
  
  # Write a file and wait until its size stops changing (reduces race conditions)
  writeStableFile <- function(path, lines, settle_checks = 3, wait_ms = 100L) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(lines, con = path, useBytes = TRUE)
    Sys.chmod(path, mode = "0644")
    last <- -1L; stable <- 0L
    repeat {
      if (!file.exists(path)) { Sys.sleep(wait_ms/1000); next }
      sz <- file.info(path)$size
      if (is.na(sz)) { Sys.sleep(wait_ms/1000); next }
      if (sz == last) {
        stable <- stable + 1L
        if (stable >= settle_checks) break
      } else {
        stable <- 0L
        last <- sz
      }
      Sys.sleep(wait_ms/1000)
    }
    normalizePath(path, winslash = "/", mustWork = TRUE)
  }
  
  
  
  # Safer launcher that returns TRUE if the completion flag appears
  runFijiMacroOnce <- function(fiji_bin, macro_path, done_flag, timeout_sec = 3600) {
    # Clear any stale flag
    if (file.exists(done_flag)) file.remove(done_flag)
    
    # Launch
    status <- tryCatch(
      {
        system2(fiji_bin, args = c("-macro", macro_path),
                wait = FALSE, stdout = "", stderr = "")
        0L
      },
      error = function(e) 1L
    )
    
    if (!is.null(attr(status, "status")) && attr(status, "status") != 0L) {
      return(FALSE)
    }
    
    # Wait until the macro writes the completion flag (or timeout)
    t0 <- Sys.time()
    repeat {
      if (file.exists(done_flag)) break
      if (as.numeric(difftime(Sys.time(), t0, units = "secs")) > timeout_sec) {
        return(FALSE)
      }
      Sys.sleep(1)
    }
    TRUE
  }
  
  
  FIX_TIFF_METADATA <- FALSE 
  
################################################################################
###    SECTION 3: Script:
################################################################################

## 1. Load Libraries
require(parallel)
require(future.apply)
require(progressr)
require(dplyr)
  
if(!("EBImage" %in% utils::installed.packages())){
  print("Installing EBImage.")
  BiocManager::install("EBImage")
}

if (!exists("RUN_NUCLEI", inherits = TRUE)) RUN_NUCLEI <- TRUE  
if (!exists("FIX_TIFF_METADATA", inherits = TRUE)) FIX_TIFF_METADATA <- TRUE
if (!exists("SEGTEST_ACTIVE", inherits = TRUE)) SEGTEST_ACTIVE <- FALSE  
  
  
  
  
  
  
################################################################################
# If CiliaPreparator is true, preparate files
################################################################################
  
  if (isTRUE(CiliaQPreparator)) {
    
    # Move input files to new directory
    
    cat(" Copying raw cilium files to CiliaQ-Preparator directory\n")
    
    # Create staging folder
    dest_dir_prep <- file.path(MainDirectory, "CiliaQ-Preparator_Files")
    dir.create(dest_dir_prep, showWarnings = FALSE, recursive = TRUE)
    
    # Find all C01 inputs under MainDirectory
    files_to_copy <- list.files(
      path = MainDirectory,
      pattern = "--C01\\.tif$",
      recursive = TRUE,
      full.names = TRUE
    )
    
    # Build a flat, unique staged path for each file (subdirs collapsed with "__")
    root <- normalizePath(MainDirectory, winslash = "/")
    stage_map <- lapply(files_to_copy, function(f) {
      rel <- sub(paste0("^", root, "/?"), "", normalizePath(f, winslash = "/"))
      dst <- file.path(dest_dir_prep, gsub("/", "__", rel))   # <-- FLATTENED NAME
      list(src = f, dst_staged = dst)
    })
    
    # Copy using the flattened destination (matches the map)
    invisible(vapply(stage_map, function(m) {
      dir.create(dirname(m$dst_staged), showWarnings = FALSE, recursive = TRUE)
      file.copy(m$src, m$dst_staged, overwrite = TRUE)
    }, logical(1)))
    
    # Write map + list for the macro
    map_df <- do.call(rbind, lapply(stage_map, as.data.frame))
    write.csv(map_df, file.path(dest_dir_prep, "staging_map.csv"), row.names = FALSE)
    
    # Macro will read exactly these staged, flattened paths
    files_to_process <- vapply(stage_map, `[[`, "", "dst_staged")
    writeLines(files_to_process, file.path(dest_dir_prep, "ImagesToEdit.txt"))
    
    
    ################################################################################
    # Quit Fiji (cross-platform)
    fiji_app_path <- path.expand(fiji_app_path)
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
    
    
    # Set up directories and paths
    input_dir <- dest_dir_prep
    if (!dir.exists(input_dir)) {
      dir.create(input_dir)
      cat("Cannot find", input_dir, "containing raw cilium images\n")
      stop("Exiting")
    } else {
      cat(" - Directory containing raw cilium images", input_dir, "\n")
    }
    
    # Full list of images to process
    files_to_process <- list.files(
      path = input_dir,
      pattern = "--C01.tif$",
      recursive = TRUE,
      full.names = TRUE
    )
    
    # Define checkpoint log files
    processed_log <- file.path(dest_dir_prep, "CiliaQPreparatorProcessedImages.txt")
    file_list_path <- file.path(input_dir, "ImagesToEdit.txt")
    
    # Remove already-processed images
    if (file.exists(processed_log)) {
      processed_images <- readLines(processed_log)
      files_to_process <- setdiff(files_to_process, processed_images)
      writeLines(files_to_process, file_list_path)
      writeLines(processed_images, processed_log)
      cat(" 🔁 Resuming: Skipping", length(processed_images), "already-processed images\n")
    } else {
      processed_images <- ""
      writeLines(files_to_process, file_list_path)
      writeLines(processed_images, processed_log)
    }
    
    
    # Check if anything remains
    if (length(files_to_process) == 0) {
      cat(" ✅ All images already processed. Nothing to do.\n")
    } else {
      cat(" - Found", length(files_to_process), "image files to process\n")
      
      
      # Prepare the macro to batch process files in Fiji
      input_dir <- path.expand(input_dir)
      ijm_dir <- gsub("\\\\", "/", normalizePath(input_dir, winslash = "/"))
      macro_path <- file.path(input_dir, "run_ciliaq_preparator.ijm")
     
      
      # Macro script
      macro_code <- sprintf(
        'setBatchMode(true);

// params
subtract_radius = %g;
smooth_radius   = %g;
seg_method      = "%s";

dir = "%s";
listPath = dir + "/ImagesToEdit.txt";
logPath  = "%s";

lines = split(File.openAsString(listPath), "\\n");
for (i = 0; i < lines.length; i++) {
    file = trim(lines[i]);
    if (file == "") continue;
    if (endsWith(file, "--C01.tif") && File.exists(file)) {
        print("\\n[+] Opening: " + file);
        open(file);

        // preprocess
        if (subtract_radius > 0) run("Subtract Background...", "rolling=" + subtract_radius);
        if (smooth_radius > 0)   run("Gaussian Blur...",       "sigma="   + smooth_radius);

        // save temp, reopen
        tmpPath = replace(file, "\\\\.[Tt][Ii][Ff]{1,2}$", "_preproc.tif");
        if (File.exists(tmpPath)) File.delete(tmpPath);
        saveAs("Tiff", tmpPath);
        close();
        open(tmpPath);

        // run plugin (note: output based on tmp basename)
        opts  = "process=[active image in FIJI] channel=1 segmentation=" + seg_method + " ";
        opts += "stack=[apply threshold determined in a maximum-intensity-projection] segment=0 ";
        opts += "output=[save as filename + suffix \'_CQP\']";
        run("CiliaQ Preparator v0.1.2", opts);

        // expected output paths
        tmpCQP = replace(tmpPath, "\\\\.[Tt][Ii][Ff]{1,2}$", "_CQP.tif");
        outCQP = replace(file,    "\\\\.[Tt][Ii][Ff]{1,2}$", "_CQP.tif");

        // normalize name: move temp-based output to original-based output name
        if (File.exists(tmpCQP)) {
            if (outCQP != tmpCQP) {
                File.copy(tmpCQP, outCQP);
                File.delete(tmpCQP);
            }
        } else {
            print("[!] Expected temp CQP not found: " + tmpCQP);
        }

        // cleanup temp image
        close();
        if (File.exists(tmpPath)) File.delete(tmpPath);

        File.append(file + "\\n", logPath + "/CiliaQPreparatorProcessedImages.txt");
    }
}
File.saveString("done", dir + "/ciliaq_preparator_complete.flag");
setBatchMode(false);',
    subtract_radius_val,
    smooth_radius_val,
    seg_method_val,
    gsub("\\\\", "/", ijm_dir),
    gsub("\\\\", "/", normalizePath(dest_dir_prep, winslash = "/"))
      )
      
      
      ###################
      # --- Write macro to a unique, safe temp file ---
      macro_name <- sprintf("cqp_%s.ijm", format(Sys.time(), "%Y%m%d_%H%M%S_%03d"))
      macro_path <- file.path(tempdir(), macro_name)
      macro_path <- writeStableFile(macro_path, macro_code)
      cat(" 🧩 Macro path: ", macro_path, "\n")
      
      # --- Find the Fiji binary inside the app bundle (macOS) ---
      fiji_app_path <- path.expand(fiji_app_path)
      fiji_bin <- file.path(fiji_app_path, "Contents", "MacOS", "ImageJ-macosx")
      if (!file.exists(fiji_bin)) {
        stop("Fiji binary not found at: ", fiji_bin,
             "\n(If Fiji moved, pass the new fiji_app_path to CiliaTGD.)")
      }
      
      # --- Make sure previous Fiji is fully closed before we start ---
      quitFiji()
      waitForFijiToExit(8)
      
      # --- Clear any old completion flag and launch ---
      done_flag <- file.path(input_dir, "ciliaq_preparator_complete.flag")
      if (file.exists(done_flag)) file.remove(done_flag)
      
      cat(" 🚀 Launching Fiji via binary with -macro\n")
      
      ok <- runFijiMacroOnce(fiji_bin, macro_path, done_flag, timeout_sec = 3600)
      if (!ok) {
        cat(" 🔁 First attempt didn’t complete. Retrying once…\n")
        quitFiji(); waitForFijiToExit(8)
        ok <- runFijiMacroOnce(fiji_bin, macro_path, done_flag, timeout_sec = 3600)
      }
      quitFiji(); waitForFijiToExit(8)
      
      if (!ok) stop("CiliaQ Preparator did not complete or macro not found. Last macro path: ", macro_path)
      cat(" ✔️ Fiji CiliaQ Preparator macro complete.\n")
      
    }
  
  
  
  ################################################
  # === Move CQP outputs back to original folders (map if present, else derive) ===
  cqp_files <- list.files(dest_dir_prep,
                          pattern = "--C01(_preproc)?_CQP\\.tif$",
                          full.names = TRUE, recursive = FALSE)
  
  map_csv <- file.path(dest_dir_prep, "staging_map.csv")
  
  if (file.exists(map_csv)) {
    # Use saved map (fast path)
    map <- read.csv(map_csv, stringsAsFactors = FALSE)
    key_base <- function(p) sub("\\.tif$", "", basename(p))
    map$key <- key_base(map$dst_staged)
    
    for (cqp in cqp_files) {
      base <- sub("(_preproc)?_CQP\\.tif$", "", basename(cqp))  # "Pos1__Pos001--C01"
      hit  <- map[map$key == base, , drop = FALSE]
      if (nrow(hit) != 1L) {
        warning("Could not uniquely map CQP file back: ", basename(cqp))
        next
      }
      orig_src <- hit$src                                   # ".../Pos1/Pos001--C01.tif"
      target   <- sub("\\.tif$", "_CQP.tif", orig_src)      # ".../Pos1/Pos001--C01_CQP.tif"
      
      # Ensure destination folder exists BEFORE copying
      dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
      ok <- file.copy(cqp, target, overwrite = TRUE)
      if (!ok) warning("Copy failed for: ", cqp, " -> ", target)
    }
    
  } else {
    # Map file missing (e.g., intermittent FS glitch). Derive path by reversing "__" -> "/".
    for (cqp in cqp_files) {
      base_flat <- sub("(_preproc)?_CQP\\.tif$", "", basename(cqp))  # "Pos1__Pos001--C01"
      rel       <- gsub("__", "/", base_flat, fixed = TRUE)           # "Pos1/Pos001--C01"
      orig_src  <- file.path(MainDirectory, rel)                      # ".../SegmentationTester_Temp/Pos1/Pos001--C01.tif"
      target    <- sub("\\.tif$", "_CQP.tif", orig_src)
      
      dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
      ok <- file.copy(cqp, target, overwrite = TRUE)
      if (!ok) warning("Copy failed (fallback) for: ", cqp, " -> ", target)
    }
  }
  ################################################
  
  # Optional cleanup
  file.remove(file_list_path)
  file.remove(processed_log)
  file.remove(macro_path)
  file.remove(file.path(dest_dir_prep, "staging_map.csv"))
  if (file.exists(done_flag)) file.remove(done_flag)
  
  ################################################
  
  }




  
  
  
  
## 2. Define loop to run across directories ------------------------------------
  
  # Get list of directories within MainDirectory
  directories <- list.dirs(MainDirectory, full.names = TRUE, recursive = FALSE)
  
  # Only keep subdirs that actually contain our expected inputs
  directories <- Filter(function(d) {
    any(file.exists(list.files(d, pattern = "--C0[01]\\.tif$", full.names = TRUE)))
  }, directories)
  
  # Add a trailing slash "/" to each directory if it doesn't already have one
  directories <- sapply(directories, function(dir) {
    if (substr(dir, nchar(dir), nchar(dir)) != "/") {
      return(paste0(dir, "/"))
    } else {
      return(dir)
    }
  })
  
  
## 3. Define the  detectNuclei function ----------------------------------------
    process_directory <- function(input_dir_tif) {
      cat("Processing directory:", input_dir_tif, "\n")
    
      
      # 0. Basics --------------------------------------------------------------
      .old.options <- options()
      on.exit(options(.old.options))
      
      options(stringsAsFactors = FALSE, warn=-1)
      
      
      
      # 1. Image preparation --------------------------------------------------
      
      # Get list of tiff files in input_dir_tif
      tif_files <- list.files(input_dir_tif, pattern = "\\.tif$", full.names = TRUE)
      
      # Find the file that ends with "--C00.tif," the nucleus image
      nuc_tif <- grep("--C00\\.tif$", tif_files, value = TRUE)
      
      # Find the file that ends with "--C01.tif," the cilium image
      cil_tif <- grep("--C01\\.tif$", tif_files, value = TRUE)
      
      # Check if nucleus image exists
      if (length(nuc_tif) == 0 || !file.exists(nuc_tif)) {
        cat("No nucleus file found in:", input_dir_tif, "\n")
        return(NULL)
      }
      
      # Check if cilium image exists
      if (length(cil_tif) == 0 || !file.exists(cil_tif)) {
        cat("No cilium file found in:", input_dir_tif, "\n")
        return(NULL)
      }
      
      # Make output directory
      output_dir <- file.path(input_dir_tif, "output")
      dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
      
      # Read the nucleus image
      image_data <- EBImage::readImage(nuc_tif, type = "tiff")
      
      # Make sure dimensions of image are multiples of 8 pixels, otherwise crop
      
          # Get image dimensions
          h <- dim(image_data)[1]  # Height (y-dimension)
          w <- dim(image_data)[2]  # Width (x-dimension)
          
          # Compute nearest lower multiple of 8
          new_h <- floor(h / 8) * 8
          new_w <- floor(w / 8) * 8
          
          # Crop image to new dimensions
          image_data <- image_data[1:new_h, 1:new_w]
          
      
          
          
              
      
      # 2. Nuclei detection ----------------------------------------------------
      
      if (isTRUE(RUN_NUCLEI)) {
            
      # Define nucleus image
      Image_nuclei <- image_data
      
      # Enhance image contrast
      Image_nuclei <- EBImage::clahe(x = Image_nuclei, nx = 8)
    
      # Compute the 20th percentile (bottom 20%) intensity value for all pixels
      threshold <- quantile(EBImage::imageData(Image_nuclei), probs = 0.20, na.rm = TRUE)
      
      # Set pixels below or equal to that threshold to 0 (removes dim background)
      EBImage::imageData(Image_nuclei)[EBImage::imageData(Image_nuclei) <= threshold] <- 0
      
      # Smooth the image
      filterSize <- floor(5*(0.22/pixel_size))
      Image_nuclei <-  EBImage::medianFilter(x = Image_nuclei, size = filterSize)
      
      # Make the image brighter if needed
      if(mean(Image_nuclei) < 0.05){
        Image_nuclei <- 10*Image_nuclei
      }else if(mean(Image_nuclei) < 0.1){
        Image_nuclei <- 5*Image_nuclei
      }else{
        Image_nuclei <- 2*Image_nuclei
      }
      
      # Define the window size within which to binarize the image (3 times the supplied expected nuclear diameter)
      nuc_mask_width_height_in_pixels <- 3*ceiling(nuc_length/pixel_size)
      
      # Generate binary nucleus mask
      nmask <- EBImage::thresh(x = Image_nuclei,
                               w = nuc_mask_width_height_in_pixels,
                               h = nuc_mask_width_height_in_pixels,
                               offset = 0.05)

      # Remove objects smaller than the expected nucleus length
      
          # Define the brush diameter (set to 1/2 the supplied expected nuclear diameter, in pixels)
          discSize <- floor(nuc_length/pixel_size)/2
          
          # Use this brush size to minimize objects smaller than half the expetced nuclear size
          nmask <- EBImage::opening(nmask, EBImage::makeBrush(discSize, shape='disc'))
      
      # Fill holes in nucleus mask
      nmask <- EBImage::fillHull(nmask)
      
      # Label each connected set of pixels with a distinct ID
      nmask <- EBImage::bwlabel(nmask)
      
      # Delete all  "nuclei" that are smaller than 35% of the median nucleus 
      # size or less than half the area expected by the nuc_length argument
      
            # First watershed the nuclei roughly, not the final watershed
            nmask2 <- EBImage::watershed(EBImage::distmap(nmask), 4)
            table_nmask <- table(nmask2)
          
            # Determine the median and mean "nucleus" size in the image
            if(length(table_nmask) >=2) {
            nuc_min_area_median <- median(table_nmask[-1])/3
            nuc_min_area_mean <- mean(table_nmask[-1])/3
            nuc_min_area <- if(nuc_min_area_mean > (2*nuc_min_area_median)){
              0.5 * (nuc_min_area_mean + nuc_min_area_median)
            } else {
              nuc_min_area_median
            }
        
            # Determine the expected nucleus area
            expected_nuc_area <- pi*(((nuc_length/pixel_size)/2)^2)
            
            # Remove nuclei less than half the expetced area or less than 1/3 the mean/median area
            to_be_removed <- c(as.integer(names(which(table_nmask < nuc_min_area))),
                               as.integer(names(which(table_nmask < 0.5*(expected_nuc_area) ))))

            to_be_removed <- sort(to_be_removed[to_be_removed != 0])
            
            if(length(to_be_removed) > 0){
              for(i in to_be_removed){
                EBImage::imageData(nmask2)[EBImage::imageData(nmask2) == i] <- 0
              }
            }
            }
      
      # Recount nuclei
      expected_nuc_area <- pi*(((nuc_length/pixel_size)/2)^2)
      nmask <- EBImage::bwlabel(nmask2)
      total_nuclei_area <- sum(nmask != 0)
      
      # Watershed the large nuclei only

            # Measure size of each object
            object_sizes <- table(nmask)
            object_sizes <- object_sizes[names(object_sizes) != "0"]
            object_labels <- as.integer(names(object_sizes))
            
            # Define size threshold to watershed as any object larger than 150% of the expected nucleus area (in pixels)
            watershed_threshold <- expected_nuc_area * 1.5
            
            # Identify labels to keep in small and large masks
            small_labels <- object_labels[object_sizes <= watershed_threshold]
            large_labels <- object_labels[object_sizes > watershed_threshold]
            
            # Create masks
            small_mask <- nmask
            large_mask <- nmask
            
            # Filter masks
            small_mask[!(nmask %in% small_labels)] <- 0
            large_mask[!(nmask %in% large_labels)] <- 0
            
            # Watershed only large objects
            # Convert large mask to binary
            large_mask_bin <- large_mask > 0
            
            # Apply distance transform and watershed to large object mask
            watershed_large <- EBImage::watershed(EBImage::distmap(large_mask_bin), watershed)
            
            # Mask the watershed result to the large object area only
            watershed_large[large_mask == 0] <- 0
            
            # Combine with small objects
            combined_mask <- small_mask
            
            # To avoid label overlap between large and small masks, shift watershed labels in large mask
            watershed_large_shifted <- watershed_large
            watershed_large_shifted[watershed_large > 0] <- watershed_large[watershed_large > 0] + max(small_mask)
            
            # Add them
            nmask_watershed <- combined_mask + watershed_large_shifted


      # Remove nuclei that touch the right or bottom borders
      right <- table(nmask_watershed[dim(nmask_watershed)[1], 1:dim(nmask_watershed)[2]])
      bottom <- table(nmask_watershed[1:dim(nmask_watershed)[1],dim(nmask_watershed)[2]])
      right <- as.integer(names(right))
      bottom <- as.integer(names(bottom))
      nuclei_at_borders <- unique(c(right, bottom))
      nuclei_at_borders <- nuclei_at_borders[nuclei_at_borders != 0]
            
      if(length(nuclei_at_borders) > 0){
        for(i in nuclei_at_borders){
          EBImage::imageData(nmask_watershed)[EBImage::imageData(nmask_watershed) == i] <- 0
        }
      }
        
      # Count number of nuclei and calculate mean area
      nucNo <- length(setdiff(unique(nmask_watershed), 0))
      mean_nucleus_area_in_pixels <- round(total_nuclei_area / nucNo)
      
      # Normalize and scale image for visualization
      Image_with_nuclei_borders <- EBImage::normalize(image_data)
      
      Image_with_nuclei_borders <- Image_with_nuclei_borders * 2
      Image_with_nuclei_borders[Image_with_nuclei_borders > 1] <- 1 
      
      # Create an image of the nuclei in blue
      Image_with_nuclei_borders <- EBImage::rgbImage(
        red = 0 * Image_with_nuclei_borders,
        green = 0 * Image_with_nuclei_borders,
        blue = Image_with_nuclei_borders 
      )
      
      # Paint the nuclei outlines in pink (#ff00ff)
      Image_with_nuclei_borders <- EBImage::paintObjects(
        x = nmask_watershed,
        tgt = Image_with_nuclei_borders,
        col = "magenta"
      )
      
      # Set image to color image
      EBImage::colorMode(Image_with_nuclei_borders) <- "color"
  
      # Save image
      if (isTRUE(RUN_NUCLEI)) {
        nuc_overlay_path <- file.path(
          output_dir,
          paste(basename(nuc_tif), "_nucleimask.tif", sep = "")
        )
        .safe_write_tiff(Image_with_nuclei_borders, nuc_overlay_path, bits = 8, type = "tiff")
      }
      
      
      
      
      
      # 5. Save results --------------------------------------------------------
      mean_nuclei_projection_diameter_in_um = sqrt(4 * mean_nucleus_area_in_pixels * pixel_size * pixel_size / pi)
      mean_nuclei_projection_diameter_in_um = round(x = mean_nuclei_projection_diameter_in_um, digits = 1)
      mean_nucleus_area_in_uM2 = mean_nucleus_area_in_pixels * (pixel_size^2)
      
      df_number_nuclei <- data.frame(
        "number_of_nuclei" = nucNo,
        "mean_nuclei_projection_area_in_pixels" = mean_nucleus_area_in_pixels,
        "mean_nuclei_projection_diameter_in_um" = mean_nuclei_projection_diameter_in_um,
        "mean_nucleus_area_in_um2" = mean_nucleus_area_in_uM2)
      
      df_number_nuclei$dirName <- nuc_tif
      df_number_nuclei$fileName <- basename(nuc_tif)
      df_number_nuclei <- df_number_nuclei %>%  
        dplyr::relocate(.data$fileName) %>%
        dplyr::relocate(.data$dirName)
      
      readr::write_csv(df_number_nuclei,
                       file = file.path(output_dir, "nuclei_number.csv"))
      
      
      } else {
        cat(" [SegmentationTester] Skipping nuclei to speed up.\n")
      }
      
      
      
      
      # ---- 6. Make cilium mask ---------------------------------------------------
      
      img <- EBImage::readImage(cil_tif, type = "tiff")  # expected 0..1 range
      
      if (isTRUE(CiliaQPreparator)) {
        
      # Retrieve CiliaQ Preparator Masks
      tif_files_CQP <- list.files(input_dir_tif,
                              pattern = "_CQP\\.tif$",
                              full.names = TRUE,
                              recursive = FALSE)
      
      # prefer the canonical (no _preproc) if present
      cil_tif_mask_CQ <- grep("--C01_CQP\\.tif$", tif_files_CQP, value = TRUE)
      if (!length(cil_tif_mask_CQ)) {
        # fall back to temp-based output name if rename didn’t happen
        cil_tif_mask_CQ <- grep("--C01_preproc_CQP\\.tif$", tif_files_CQP, value = TRUE)
      }
      stopifnot(length(cil_tif_mask_CQ) == 1)
      cil_tif_mask <- EBImage::readImage(cil_tif_mask_CQ, type = "tiff")
      
      
        # --- Binarize-only path ---------------------------------------------------
        th_otsu   <- EBImage::otsu(cil_tif_mask)
        floor_thr <- cilium_threshold_value / 255  # keep your manual floor (0..1)
        binary_img <- cil_tif_mask > max(th_otsu, floor_thr)
        
        # Label & remove small objects
        cmask <- EBImage::bwlabel(binary_img)
        
        tbl <- table(cmask)
        if ("0" %in% names(tbl)) tbl <- tbl[names(tbl) != "0"]
        
        # If min_cilium_area is µm^2, convert to pixels; else treat as pixels
        # min_pixels <- ceiling(min_cilium_area / (pixel_size^2))
        min_pixels <- ceiling(min_cilium_area)
        
        to_rm <- as.integer(names(tbl)[tbl < min_pixels])
        cmask2 <- if (length(to_rm)) EBImage::rmObjects(cmask, to_rm) else cmask
        
        final_mask <- cmask2 > 0
        
      } else {
        # --- Full processing path -------------------------------------------------
        # 1) Background subtraction (white top-hat)
        se <- EBImage::makeBrush(background_subtraction, shape = "box")
        bg <- EBImage::opening(img, se)
        img1 <- pmax(img - bg, 0)
        
        # 2) Gentle denoise
        img_blurred <- EBImage::medianFilter(img1, size = 3)
        
        # (optional) second background estimate
        bg2 <- EBImage::opening(img_blurred, se)
        img2 <- pmax(img_blurred - bg2, 0)
        
        # 3) Threshold
        th_otsu   <- EBImage::otsu(img2)
        floor_thr <- cilium_threshold_value / 255
        binary_img <- img2 > max(th_otsu, floor_thr)
        
        # 4) Small closing to bridge tiny gaps
        brush <- EBImage::makeBrush(3, shape = "disc")
        binary_img <- EBImage::closing(binary_img, brush)
        
        # 5) Remove small objects
        cmask <- EBImage::bwlabel(binary_img)
        tbl <- table(cmask)
        if ("0" %in% names(tbl)) tbl <- tbl[names(tbl) != "0"]
        
        # min_pixels <- ceiling(min_cilium_area / (pixel_size^2))  # if µm^2
        min_pixels <- ceiling(min_cilium_area)                      # if pixels
        
        to_rm <- as.integer(names(tbl)[tbl < min_pixels])
        cmask2 <- if (length(to_rm)) EBImage::rmObjects(cmask, to_rm) else cmask
        
        # 6) Aspect ratio / blobiness filters
        labeled <- cmask2
        if (max(labeled) == 0) {
          final_mask <- labeled > 0
        } else {
          fs <- EBImage::computeFeatures.shape(labeled)
          fm <- EBImage::computeFeatures.moment(labeled)
          
          df_shape  <- transform(as.data.frame(fs),  label = as.integer(rownames(fs)))
          df_moment <- transform(as.data.frame(fm),  label = as.integer(rownames(fm)))
          
          features <- merge(df_shape, df_moment, by = "label", all = TRUE)
          
          if (nrow(features) >= 1) {
            minor <- features$m.majoraxis * sqrt(pmax(1 - features$m.eccentricity^2, 0))
            aspect_ratio <- features$m.majoraxis / pmax(minor, .Machine$double.eps)
            perim_area_ratio <- features$s.perimeter / pmax(features$s.area, 1)
            blobiness <- features$s.perimeter / pmax(features$m.majoraxis, .Machine$double.eps)
            
            keep <- features$label[ aspect_ratio > cilium_aspectratio &
                                      perim_area_ratio < 0.45 &
                                      blobiness < 3 ]
            
            cleaned <- if (length(keep)) {
              EBImage::rmObjects(labeled, setdiff(seq_len(max(labeled)), keep))
            } else labeled
            
            final_mask <- cleaned > 0
          } else {
            final_mask <- labeled > 0
          }
        }
      }
          
      # Shift cilium mask
      shifted_mask <- shift_image(final_mask, x_shift = -20, y_shift = 20)
      
      # Create an image of the cilia in blue
      Image_with_cilium_mask <- EBImage::rgbImage(
        red = shifted_mask*0.5,
        green = 0 * final_mask,
         blue = img *5 
      )
    
      # Set image to color image
      EBImage::colorMode(Image_with_cilium_mask) <- "color"
      
      # Save cilium mask image
      overlay_path <- file.path(output_dir, paste0(basename(nuc_tif), "_CiliumOverlay.tif"))
      
      .safe_write_tiff(
        Image_with_cilium_mask,
        overlay_path,
        bits = 8
      )
      
      
      # Create an image of cilium mask for CiliaQ
          
          # Create an empty mask
          empty_mask <- final_mask *0
          
          # Stack the three channels along the 3rd dimension
          multi_channel_array <- abind::abind(empty_mask, final_mask, img * 5, final_mask, along = 3)
          
          # Create an EBImage object with 3 channels
          Image_with_cilium_mask <- EBImage::Image(multi_channel_array, colormode = 0)
          
          # Save the image as a multi-channel TIFF
          ciliq_path <- file.path(
            output_dir,
            paste(basename(nuc_tif), "_ciliummask_forCiliaQ.tif", sep = "")
          )
          .safe_write_tiff(Image_with_cilium_mask, ciliq_path, bits = 8, type = "tiff")
          
  ##################
  ### May need to install tiffile at system level:
  #
  # python3 -m venv ~/tiffenv
  # source ~/tiffenv/bin/activate
  # pip install tifffile
  #
  ##################
      
          # Define function to update image metadata
          fix_tiff_metadata_with_tifffile <- function(input_tif, output_tif, venv_path = "~/tiffenv") {
          
          py_script <- tempfile(fileext = ".py")
          
          script_contents <- sprintf(
            '\
import tifffile
import numpy as np

# Hardcoded ImageJ-style metadata
metadata = {
    "ImageJ": "1.54p",
    "images": 2,
    "channels": 2,
    "mode": "composite",
    "unit": "um",
    "loop": False
}

resolution = (13352921, 1000000)
resolutionunit = 1

# Load the image from R
original = tifffile.imread(r"%s")

if original.ndim != 3 or original.shape[0] < 2:
    raise ValueError("Expected at least 2 slices (channels) as first dimension, got shape: {}".format(original.shape))

# Save with copied metadata
tifffile.imwrite(
    r"%s",
    original,
    photometric="minisblack",
    metadata=metadata,
    imagej=True,
    resolution=resolution,
    resolutionunit=resolutionunit
)
', normalizePath(input_tif), normalizePath(output_tif))
          
          writeLines(script_contents, con = py_script)
          
          py_bin <- file.path(path.expand(venv_path), "bin", "python")
          cmd <- paste(shQuote(py_bin), shQuote(py_script))
          message("Running command:\n", cmd)
          system(cmd)
        }

          
        # Create an image of cilium mask for CiliaQ-Editor
  
            # Shift starting image
            shifted_image <- shift_image(img, x_shift = -20, y_shift = -20)
          
            # Create an empty mask
            empty_mask <- final_mask *0
            
            # Stack the three channels along the 3rd dimension
            multi_channel_array <- abind::abind(shifted_image * brightnessfactor, final_mask, along = 3)
            
            # Create an EBImage object with 2 channels
            Image_with_cilium_mask <- EBImage::Image(multi_channel_array, colormode = 0)
            
            # Save the image as a multi-channel TIFF
            shifted_path <- file.path(
              output_dir,
              paste(basename(nuc_tif), "_ShiftedCiliumMask.tif", sep = "")
            )
            .safe_write_tiff(Image_with_cilium_mask, shifted_path, bits = 8, type = "tiff")
            
      # Run function to update Shifted Cilium Mask
      
            # Define the shifted output path
            shifted_path <- file.path(output_dir,
                                      paste(basename(nuc_tif),
                                            "_ShiftedCiliumMask.tif", sep = ""))
            
            # Only try to fix metadata if enabled and file exists
            if (isTRUE(FIX_TIFF_METADATA) && file.exists(shifted_path)) {
              tryCatch(
                fix_tiff_metadata_with_tifffile(shifted_path, shifted_path),
                error = function(e) message("⚠️  Skipping metadata fix for ", shifted_path, ": ", conditionMessage(e))
              )
            }
      
      
          
          
}

  
  
  
################################################################################
###    SECTION 4: Run script:
################################################################################
  
    cat(" Segmenting nuclei and cilia in each image\n")

    num_cores <- max(1, parallel::detectCores() - 4)

    ## pick a supported backend (multisession) and a console handler
    plan(multisession, workers = num_cores)
    handlers("progress")    # simple percentage bar
    
    ## wrap creation of the progressor AND your loop in with_progress()
    with_progress({
      p <- progressor(steps = length(directories))
      
      results <- future_lapply(directories, function(dir) {
        p()    # tick the bar
        
        # === muffling EVERYTHING from process_directory() ===
        withCallingHandlers(
          {
            capture.output(
              process_directory(dir),
              file = NULL   # drops all cat()/print() output
            )
          },
          message = function(m) invokeRestart("muffleMessage"),
          warning = function(w) invokeRestart("muffleWarning")
        )
      })
    })
    
    cat(" Done segmenting nuclei and cilia in each image\n")    
    

################################################################################
###    SECTION 5: Combine nucleus output files into a single summary file:
################################################################################
    
    if (isTRUE(RUN_NUCLEI)) {
      
    cat(" Combining nuclei counts\n")
    
  # Initialize empty lists to hold data
  all_Summary_Nuclei <- list()
    
  # Loop over directories and read in the files  
  for (input_dir in directories) {
    results_file <- file.path(input_dir, "output/nuclei_number.csv")
 
    # Check if the files exist before reading
    if (file.exists(results_file)) {
      nucleus_data <- read.csv(results_file)
      all_Summary_Nuclei[[length(all_Summary_Nuclei) + 1]] <- nucleus_data
    } else {
      cat("Nucleus file missing in:", input_dir, "\n")
    }
  }
  
  # Combine all the data frames using rbind
  combined_Summary_Nuclei <- do.call(rbind, all_Summary_Nuclei)

  # Write the combined data to CSV files
  write.csv(combined_Summary_Nuclei, 
            file = paste(MainDirectory, "/Combined_Summary_Nuclei.csv", sep = ""), 
            row.names = FALSE)

    } else {
      cat(" Skipping nuclei summary (RUN_NUCLEI == FALSE)\n")
    }
  
  
  
################################################################################
###    SECTION 7: Move CiliaQ-Editor Masks to a new directory:
################################################################################
  
  cat(" Copying cilium masks to CiliaQ-Editor directory\n")
  
  # Create a new subdirectory "CiliaQ_Files" inside MainDirectory
  dest_dir <- file.path(MainDirectory, "CiliaQ-Editor_Files")
  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir)
    cat(" - Created directory:", dest_dir, "\n")
  } else {
    cat("Directory already exists:", dest_dir, "\n")
  }
  
  # Recursively list all files under MainDirectory that end with '_ShiftedCiliumMask.tif'
  # Using full.names = TRUE gives absolute or relative paths that include the subdirectory names
  files_to_copy <- list.files(path = MainDirectory, 
                              pattern = "_ShiftedCiliumMask\\.tif$", 
                              recursive = TRUE, 
                              full.names = TRUE)
  
  # Copy each file to the destination directory
  # Using file.copy() with basename() ensures the file name is preserved
  sapply(files_to_copy, function(f) {
    dest_file <- file.path(dest_dir, basename(f))
    file.copy(from = f, to = dest_file, overwrite = TRUE)
  })
  
 
  
  