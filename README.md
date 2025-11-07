# CiliaTGD

**CiliaTGD** is an R package wrapper for **CiliaQ** that automates high-throughput image segmentation, interactive editing, and downstream analysis of cilia using FIJI’s CiliaQ suite.  
It bridges image analysis in R with ImageJ/FIJI macros, providing a reproducible end-to-end workflow for quantitative ciliary measurements.  
It is currently optimized for **2D images** and has **not been tested on 3D z-stacks**.

---

## ✳️ Features

- **Automated segmentation** of nuclei and cilia in multi-channel images, using either CiliaQ Preparator or custom segmentation methods.
- **Automated SegmentationTester** for testing and optimizing cilium segmentation parameter combinations.
- **Interactive editing** using the CiliaQ Editor within FIJI for optional manual curation of segmented ciliary masks.
- **Direct integration with CiliaQ** for morphometric analysis of cilia (e.g., cilium length, bending index).
- **Automated nucleus counting** for calculation of percent ciliation per image.

---

## 🧰 Installation

You can install the development version directly from GitHub:

```r
# install.packages("remotes")
remotes::install_github("drivas-lab/CiliaTGD")
```

---

## ⚙️ Dependencies

### System-Level Dependencies

The package requires a few system tools to run correctly.

You will need to install **tifffile**, a Python package used for writing ImageJ-compatible TIFF metadata.  

In Terminal, run the following commands:

```bash
python3 -m venv ~/tiffenv
source ~/tiffenv/bin/activate
pip install tifffile
```

### R Package Dependencies

These R packages are installed automatically when you install CiliaTGD:

- `EBImage`
- `parallel`
- `future.apply`
- `progressr`
- `dplyr`
- `stringr`
- `abind`
- `readr`
- `BiocManager`

> **Note:** Installing `EBImage` requires Bioconductor.
> 
> ```r
> install.packages("BiocManager")
> BiocManager::install("EBImage")
> ```

### FIJI and CiliaQ

CiliaTGD requires **FIJI** installed with the **CiliaQ plugin**.

#### Install FIJI:
Download from [https://imagej.net/software/fiji/downloads](https://imagej.net/software/fiji/downloads)

#### Install CiliaQ:
1. Open FIJI  
2. Click **Help → Update...**  
3. Wait for the update check to complete  
4. In the ImageJ Updater window, click **Manage Update Sites**  
5. Search for **CiliaQ** and check the box next to it  
6. Click **Apply and Close**, then restart FIJI  

---

## 📁 Expected Data Directory Structure

CiliaTGD expects image pairs (one nucleus channel, ending in --C00.tif, one cilium channel ending in --C01.tif) arranged as follows:  

```
MainDirectory/
 ├── SubDirectory1/
 │    ├── Image1--C00.tif   # nucleus channel
 │    └── Image1--C01.tif   # cilium channel
 ├── SubDirectory2/
 │    ├── Image2--C00.tif
 │    └── Image2--C01.tif
 └── SubDirectory3/
      ├── ...
```

---

## 🚀 Quick Start Guide

### Available Parameters

| Parameter | Type | Description |
|------------|------|-------------|
| `MainDirectory` | string | Path to directory containing subfolders with images. |
| `pixel_size` | numeric | Microns per pixel. Usually found in microscope metadata or via **Image → Properties** in FIJI. |
| `nuc_length` | numeric | Expected nucleus diameter in µm. Default = 12. |
| `watershed` | numeric | Watershed splitting parameter. Larger values split nuclei less aggressively. Default = 3. |
| `NucleusBorderRemoval` | logical | When counting nuclei, should nuclei on the right and bottom image border be removed? Default is TRUE (remove) |
| `min_cilium_area` | numeric | Minimum area (in pixels) per cilium. Used for filtering in valid cilia by size (excluding small artifacts). |
| `fiji_app_path` | string | Full path to FIJI (e.g., `path.expand("~/Desktop/Fiji.app")`). |
| `SegmentationTester` | logical | Enables parameter grid testing for cilium segmentation to identify ideal segmentation approaches for your data (TRUE/FALSE). Default = FALSE |

---

### When `SegmentationTester = FALSE`

You must choose whether to use **CiliaQ Preparator** (`CiliaQPreparator = TRUE`) or the **custom segmentation** mode (`CiliaQPreparator = FALSE`) to segment cilia.

#### If `CiliaQPreparator = FALSE`:
| Parameter | Description |
|------------|-------------|
| `cilium_threshold_value` | Intensity threshold for initial segmentation. (normalized 0–255). Higher numbers will threshold more and be more selective about what they call cilia. Values between 4 and 12 seem to be best, but test different values for your image set. |
| `background_subtraction` | Rolling ball radius for background subtraction (in pixels). Use higher numbers if you have less background. For IF images with low background, ~20 works best. For live cell images, ~8-10 works best.  |
| `cilum_aspectratio` | Minimum aspect ratio to classify as a cilium. |

#### If `CiliaQPreparator = TRUE`:
| Parameter | Description |
|------------|-------------|
| `subtract_radius_val` | Subtraction brush radius, in pixels, for segmenting cilia from background in CiliaQ Preparator. |
| `smooth_radius_val` | Gaussian blur radius in pixels for applying a blur to cilium images to effectively connect/fill in disjointed cilium objects. |
| `seg_method_val` | Segmentation method to apply in CiliaQ Preparator. Options in clude "RenyiEntropy","Otsu","Triangle","Default","Huang","Li","Mean","Minimum","Moments","Percentile", and "Yen" |

---

### Other Parameters

| Parameter | Type | Description |
|------------|------|-------------|
| `CiliaQEditor` | logical | Run the FIJI CiliaQ Editor for manual cleanup (`TRUE`) or skip directly to analysis (`FALSE`). Default is TRUE. |
| `brightnessfactor` | numeric | Multiplier for cilium channel brightness to aid mask visualization. Default = 2.5. |

---

### Example: Run the SegmentationTester

Automatically iterate over segmentation parameters to identify optimal thresholds for an image:

```r
CiliaTGD(
  MainDirectory       = "/path/to/experiment",
  pixel_size          = 0.0913,
  min_cilium_area     = 50,
  SegmentationTester  = TRUE,
  fiji_app_path       = "/Applications/Fiji.app"
)
```

Expected output will be within your MainDirectory, under a new sub directory called `SegmentationTests` containing:

| Output     | Description |
|------------|-------------|
| `cilia_counts_barplot_png` | A barplot displaying the number of cilia identified per segmentation method. |
| `cilia_counts.csv` and `cilia_counts.txt` | Data files containing the number of cilia identified per segmentation method. |
| *.tif | All the segmented images for manual review |

---

### Example: Run the Pipeline

```r
library(CiliaTGD)

CiliaTGD(
  MainDirectory          = "/path/to/experiment",
  pixel_size             = 0.0913,
  min_cilium_area        = 50,
  background_subtraction = 10,
  cilum_aspectratio      = 3,
  cilium_threshold_value = 15,
  fiji_app_path          = "/Applications/Fiji.app",
  CiliaQPreparator       = FALSE,
  CiliaQEditor           = TRUE
)
```

or 

```r
library(CiliaTGD)

CiliaTGD(
  MainDirectory          = "/path/to/experiment",
  pixel_size             = 0.0913,
  min_cilium_area        = 50,
  subtract_radius_val    = 20,
  smooth_radius_val      = 2,
  seg_method_val         = "RenyiEntropy",
  fiji_app_path          = "/Applications/Fiji.app",
  CiliaQPreparator       = TRUE,
  CiliaQEditor           = FALSE
)
```


Expected output will be within your MainDirectory, and will include:

| Output    | Description |
|------------|-------------|
| `CiliationData.txt` | A tab separated file containing the nucleus count, cilium, count, average cilium length, and percent ciliation per image. |
| `CiliumLength.txt` | A tab separated file containing the cilium length and bending index for each cilium in each image. |
| `CiliaQ_Files` | A directory containing all the output files from CiliaQ analysis. |
| `CiliaQ-Editor_Files` | A directory containing all the output files from the CiliaQ-Editor analysis. |
| `CiliaQ-Preparator_Files` | A directory containing all the output files from the CiliaQ-Preparator analysis. | 


---

## 🧠 Notes

- CiliaTGD outputs segmented masks, CiliaQ files, and summary statistics to subdirectories under your main directory.
- The pipeline creates a full log (`CiliaTGD_RunLog.txt`) documenting all parameters used.
- TIFF metadata and ImageJ compatibility are preserved automatically.
- Currently optimized for 2D data; not tested with z-stacks.

---

## 📜 License

MIT License © 2025 Theodore G. Drivas
