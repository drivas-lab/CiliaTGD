\name{CiliaTGD}
\alias{CiliaTGD}
\title{Automated Pipeline for Cilia and Nuclei Analysis}
\description{Automates the segmentation of nuclei and cilia, interactive editing via Fiji, and CiliaQ processing.}
\usage{
Expects a single directory, with subdirectories each containing an image set, ending *--C00.tif and *--C01.tif, where *--C00.tif is the image of the nuclei and *--C01.tif is the image of the cilia.

#   Example directory structure:
#
#     |MainDirectory
#         |
#         |-> SubDirectory1
#                  |-> Image1--C00.tif
#                  |-> Image1--C01.tif
#         |-> SubDirectory2
#                  |-> Image2--C00.tif
#                  |-> Image2--C01.tif
#         |-> SubDirectory3
#         |-> ...

CiliaTGD(MainDirectory,
         pixel_size,
         nuc_length,
         watershed,
         cilium_threshold_value,
         min_cilium_area,
         background_subtraction,
         cilum_aspectratio,
         CiliaQPreparator,
         subtract_radius_val,
         smooth_radius_val,
         seg_method_val,
         CiliaQEditor,
         brightnessfactor,
         SegmentationTester,
         fiji_app_path)
}
\arguments{
  \item{MainDirectory}{\code{character}. Path to the main directory containing image subdirectories.}
  \item{pixel_size}{\code{numeric}. Microns per pixel. You can get this from ImageJ if needed}
  \item{nuc_length}{\code{numeric}. Expected diameter of a nucleus in microns (12uM works well for hTERT-RPE1).}
  \item{watershed}{\code{numeric}. Watershed splitting parameter. Larger values split less aggressively. Default is 3, but higher numbers will split less, and lower will split more}
  \item{cilium_threshold_value}{\code{numeric}. Threshold for initial cilium segmentation (normalized 0–255). Higher numbers will threshold more and be more selective about what they call cilia. Values between 4 and 12 seem to be best, but test different values for your image set. Only needed if CiliaQ Preparator is set to FALSE (the default)}
  \item{min_cilium_area}{\code{numeric}. Minimum area (in pixels) to keep objects in cilium mask. 50 works well. Only needed if CiliaQ Preparator is set to FALSE (the default)}
  \item{background_subtraction}{\code{numeric}. Brush size for background subtraction (in pixels). Use higher numbers if you have less background. For IF images with low background, ~20 works best. For live cell images, ~8-10 works best. Only needed if CiliaQ Preparator is set to FALSE (the default)}
  \item{cilum_aspectratio}{\code{numeric}. Aspect ratio threshold for filtering cilia. For IF images with long cilia, higher numbers like 3 work well. For fixed images with shorter cilia, lower numebrs like 1.5 work well. Only needed if CiliaQ Preparator is set to FALSE (the default)}
  \item{CiliaQPreparator}{\code{TRUE/FALSE}. Should the background subtraction be done with CiliaQ Preparator? Default is FALSE. If set to true, need to define subtract_radius_val, smooth_radius_val, and seg_method_val}
  \item{subtract_radius_val}{\code{numeric}. Radius in pixels with which to subtract background if CiliaQPreparator is TRUE}
  \item{smooth_radius_val}{\code{numeric}. Radius in pixels with which to blur image to smooth and connect cilia CiliaQPreparator is TRUE}
  \item{seg_method_val}{\code{character}. Background Subtraction method if CiliaQPreparator is TRUE. Options include "RenyiEntropy", "Triangel", "Otsu", etc}
  \item{CiliaQEditor}{\code{TRUE/FALSE}. Should CiliaQ Editor be launched after segmentation for manual clean up (TRUE), or should the segmented image be processed directly (FALSE)?}
  \item{brightnessfactor}{\code{numeric}. If using CiliaQ Editor, this is the factor by which to increase the brightness of the original cilium image when generating the cilium mask overlay for ease of visualization -- this does not change the actual segmentation. Useful if the cilium image is realtively dim so you can see the cilia clearly during editing. The default is 3}
  \item{SegmentationTester}{\code{TRUE/FALSE}. Set to TRUE if you would like CiliaTGD to iterate across a number of segmentation approaches so you can chose the best one. This will not run the whole pipeline, just segmentation and produce a list and graph of results.}
  \item{fiji_app_path}{\code{character}. Path to the \code{Fiji.app} bundle for running ImageJ macros.}
}
\details{
This function wraps three pipeline stages:
\enumerate{
  \item Nucleus and cilium segmentation using \code{EBImage}.
  \item Interactive editing of cilia masks in Fiji (\code{CiliaQ Editor}).
  \item Automated CiliaQ processing and result compilation.
}
All output is written under \code{MainDirectory}, with logging in \code{CiliaTGD_RunLog.txt}.
}
\value{
\code{NULL}. Primary side-effects are writing output files and logs.
}
\author{
Theodore Drivas [aut, cre]
}
\keyword{utilities}
