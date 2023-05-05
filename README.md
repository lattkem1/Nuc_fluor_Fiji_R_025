# Nuc_fluor_Fiji_R_01
Simple nuclear fluorescence intensity image analysis pipeline in Fiji and R

Fiji/R nuclear fluorescence intensity analysis pipeline (v02.5, 5/5/23)

Pipeline usage overview (steps):

1) Basic image read, processing and ROI measurements
2a) (optional): Extract images for presentation
2) Read ROI data into R
3a) Analysis of individual channels
3b) Analysis of intensities of 2 channels (Ch1 vs Ch2)

Required: 
-	Fiji
-	R (including libraries tidyverse, colorRamps, lmerTest)

Input:
-	Leica confocal lif file (Z-stack) 16-bit
o	Other input needs customised input scripts (replacing ‘01_MaxProj greyscale batch_v02_17-12-20_from lif.ijm’, see subfolder ‘alternative input and ROI definitition scripts’)
o	8-bit images may require alternative R scripts for proper visualisation (see subfolder ‘R scripts optimised for 8-bit images’)


Description of pipeline usage steps

1) Basic image read, processing and ROI measurements

-	Copy script directory in your main analysis directory
o	For alternative data formats, replace scripts in main script directory with relevant versions from subdirectories (see “Input” above) 
-	Open images with Fiji

-	Run Fiji script “01_MaxProj greyscale batch_v02_17-12-20_from lif.ijm” 
o	select main analysis directory when asked
o	Will create maximum projection of all opened images in “01_MaxProj greyscale” subfolder 
-	Run Fiji script “02_Nuclei ROI measure batch_v021_18-01-15.ijm” 
o	change DAPI channel if required (script line 2) 
o	select main analysis directory when asked
o	set DAPI threshold for each image to outline nuclei (background has to be red)
o	will create for each image (folder “02_ImageJ analysis”)
	ROI sets 
	csv files with intensity for each nucleus for each channel


2a) (optional): Extract images for presentation

script “02a_Greyscale to RGB sel chan adj brightness selected area batch_v022_20-12-15.ijm”)
o	Can create RGB images of selected channels in selected colors 
o	Define output folder in script
o	Change parameters in script (channels, colors, intensity range)
o	Select input folder when running script


2) Read ROI data into R

-	Run R script “03_R batch input_v02_18-01-15.R”
o	Define channel names in script line 24 (in order of channels)
o	Reads data from “02_ImageJ analysis” for downstream R analyses in “03_R input” folder
o	Creates “image_groups.csv” where you have to define the image groups for the downstream analysis 
	“group” can be e.g. your experimental conditions
	“repl” can be e.g. a number or letter to identify your coverslip, if you have multiple technical replicates, to know which image is from which coverslip


3a) Analysis of individual channels

-	Run R script “04_intens distrib 1 Ch 16bit_v023_230301_all channels new plots.R”
o	Define image groups in “03_R input” folder (see above)
o	Creates QC plots in “04_intens distrib_1Ch” folder, showing intensity distributions for selected channel
-	Run R script “05_threshold testing 1 Ch_v025_230302 with t-test stats.R”
o	Set the channel you want to analyse in line 27
o	From the QC plots and/or the original max proj images get some thresholds that could make sense
o	Set the thresholds for which you want to get summary statistics in line 29 (you can set as many thresholds as you like to test, separated with commas: “thCh1.test.thresh = c(20, 30, 40)” means it will analyse with 20, 30 and 40 as thresholds )
o	Output:
	 plots and summary csv files for the different thresholds
	Basic statistic analysis of first group vs other groups (t-test or GLMM stats with lmertest package)


3b) Analysis of intensities of 2 channels (Ch1 vs Ch2)

-	Run R script “04b_intens distrib 2 Ch 16bit_v026_230301_new plots.R”
o	Define the 2 channels you want to analyse in line 30
o	Gives you folder “04b_intens distrib_2Ch”, with FACS-style scatterplots intensity Ch1 vs Ch2
-	Run R script “05b_threshold analysis 2 Ch_v025_230302 lmer stats included.R”
o	Define dataset that you want to analyse in line 21
o	Set thresholds for both channels in line 28 (only each one threshold possible: c(thresh_Ch1, thresh_Ch2))
o	Output: folder “05b_treshold analysis_2Ch”, with plots of fraction of Ch1/Ch2 -/-, -/+, … of all cells; and basic stats
