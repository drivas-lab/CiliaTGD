################################################################################
###   Launch ImageJ to process CiliaQ Edited Files (from CiliaQ_Files)
###   and combine all output .CQ.txt files
################################################################################

suppressPackageStartupMessages({
  library(stringr)
  library(dplyr)
  library(readr)
})

# Ensure Fiji path expanded
fiji_app_path <- path.expand(fiji_app_path)

# ---------- Helper: quit Fiji safely ----------
quitFiji <- function() {
  sys <- Sys.info()[["sysname"]]
  if (sys == "Darwin") {
    system2("osascript", args = c("-e", 'tell application "Fiji" to quit'),
            stdout = FALSE, stderr = FALSE)
    system2("pkill", args = c("-f", "Fiji.app"), wait = FALSE)
    system2("pkill", args = c("-f", "Fiji"), wait = FALSE)
  } else if (sys == "Windows") {
    system2("taskkill", args = c("/IM", "Fiji.exe", "/F"),
            stdout = FALSE, stderr = FALSE)
  } else {
    system("pkill -f Fiji", wait = FALSE)
  }
}

# ---------- Define working folder ----------
input_dir <- file.path(MainDirectory, "CiliaQ_Files")  # always process here
if (!dir.exists(input_dir)) {
  stop("❌ Cannot find ", input_dir, " — expected folder containing *_ShiftedCiliumMask_ed.tif files")
} else {
  cat(" - Processing all edited masks from:", input_dir, "\n")
}

# ---------- Find all edited mask files ----------
files_to_process <- list.files(
  path = input_dir,
  pattern = "_ShiftedCiliumMask_ed\\.tif$",
  recursive = TRUE,
  full.names = TRUE
)

if (!length(files_to_process)) {
  stop("❌ No *_ShiftedCiliumMask_ed.tif files found in: ", input_dir)
} else {
  cat(" - Found", length(files_to_process), "image files to process\n")
}

# ---------- Build and run macro ----------
ijm_dir <- gsub("\\\\", "/", normalizePath(input_dir, winslash = "/"))
macro_path <- file.path(input_dir, "run_ciliaq_processor.ijm")

macro_code <- sprintf(
  'dir = "%s";
list = getFileList(dir);
for (i = 0; i < list.length; i++) {
    file = dir + "/" + list[i];
    if (endsWith(file, "ShiftedCiliumMask_ed.tif")) {
        print("Opening: " + file);
        open(file);
        run("CiliaQ v0.1.7", "process=[active image in FIJI] manually ->=%.2f ->_0=1.000000 ->_1=um time=0.50 time_0=min preferences=[manually enter preferences] output=[US (0.00...)] channels=2 =1 =3 =4 minimum=10 increase additionally=[cilia touching x or y borders] minimum_0=1 increase_0 determine before=2 =0 reference=1");
        close();
    }
}
File.saveString("done", dir + "/ciliaq_complete.flag");
', ijm_dir, pixel_size)

writeLines(macro_code, con = macro_path)
cmd <- sprintf('open -n -a "%s" --args -macro "%s"', fiji_app_path, macro_path)
cat(" 🚀 Launching Fiji CiliaQ Processor macro\n")
system(cmd)

# ---------- Wait for Fiji completion ----------
done_flag <- file.path(input_dir, "ciliaq_complete.flag")
cat(" ⏳ Waiting for CiliaQ to finish running...\n")
while (!file.exists(done_flag)) Sys.sleep(2)
quitFiji()
cat(" ✔️️ Fiji CiliaQ Processor macro complete.\n")

# ---------- Collect .CQ.txt results ----------
txt_files <- list.files(path = input_dir, pattern = "CQ\\.txt$", full.names = TRUE, recursive = TRUE)
if (!length(txt_files)) stop("❌ No .CQ.txt files produced in ", input_dir)

cat(" - Found", length(txt_files), ".CQ.txt result files\n")

all_data <- list()
for (file_path in txt_files) {
  lines <- readLines(file_path, warn = FALSE)
  
  # Extract image name
  image_name_line <- lines[grepl("image name:", lines)]
  image_name <- if (length(image_name_line) > 0) {
    str_trim(strsplit(image_name_line, "\t", fixed = TRUE)[[1]][2])
  } else {
    basename(file_path)
  }
  
  # Locate results section
  results_idx <- which(grepl("^\\s*Results:", lines))
  if (!length(results_idx)) next
  header_line <- lines[results_idx + 1]
  col_names <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
  
  # Column indices
  length_idx <- which(str_detect(col_names, "cilia length.*shortest path"))
  bend_idx   <- which(str_detect(col_names, "bending index"))
  if (length(length_idx) != 1 || length(bend_idx) != 1) {
    message("Skipping file (columns missing): ", basename(file_path))
    next
  }
  
  # Parse results
  data_lines <- lines[(results_idx + 2):length(lines)]
  data_lines <- data_lines[grepl("^\\s*\\d+\\s+", data_lines)]
  parsed_data <- lapply(data_lines, function(line) strsplit(line, "\t", fixed = TRUE)[[1]])
  if (!length(parsed_data)) next
  
  max_len <- length(col_names)
  parsed_data <- lapply(parsed_data, function(x) { length(x) <- max_len; x })
  data_matrix <- do.call(rbind, parsed_data)
  
  cilium_length_um <- suppressWarnings(as.numeric(data_matrix[, length_idx]))
  bending_index    <- suppressWarnings(as.numeric(data_matrix[, bend_idx]))
  
  all_data[[file_path]] <- data.frame(
    cilium_length_um = cilium_length_um,
    bending_index = bending_index,
    image_name = image_name,
    stringsAsFactors = FALSE
  )
}

Combined_data <- bind_rows(all_data)
Combined_data$image_name <- str_replace(Combined_data$image_name, "--.*", "")

# ---------- Summarize per image ----------
CiliumCountData <- Combined_data %>%
  group_by(image_name) %>%
  summarise(
    CiliumCount = n(),
    AverageCiliumLength = mean(cilium_length_um, na.rm = TRUE),
    .groups = "drop"
  )

# ---------- Load nuclei summary ----------
nuclei_file <- file.path(MainDirectory, "Combined_Summary_Nuclei.csv")
if (!file.exists(nuclei_file)) {
  stop("❌ Missing nuclei summary file: ", nuclei_file)
}
NucleusCountData <- read.csv(nuclei_file)
NucleusCountData$fileName <- str_replace(NucleusCountData$fileName, "--.*", "")

# ---------- Join nuclei and cilia ----------
CombinedNucCilData <- left_join(NucleusCountData, CiliumCountData,
                                by = c("fileName" = "image_name")) %>%
  mutate(
    CiliumCount = ifelse(is.na(CiliumCount), 0, CiliumCount),
    PercCiliation = ifelse(CiliumCount == 0, 0, CiliumCount / number_of_nuclei * 100)
  )

# ---------- Write outputs ----------
write.table(CombinedNucCilData,
            file = file.path(MainDirectory, "CiliationData.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(Combined_data,
            file = file.path(MainDirectory, "CiliumLength.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")

cat(" ✅ Wrote outputs:\n")
cat("  -", file.path(MainDirectory, "CiliationData.txt"), "\n")
cat("  -", file.path(MainDirectory, "CiliumLength.txt"), "\n")
