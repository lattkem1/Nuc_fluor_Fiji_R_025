//define DAPI channel
Ch_DAPI = 3

//define I/O directories
maxp_dir = getDirectory("Choose Image Directory ('01_MaxProj greyscale' subfolder)");
print("Image directory: ", maxp_dir);
dir = File.getParent(maxp_dir);
print("Main directory: ", dir);
analysis_dir = dir+"/02_ImageJ analysis/";
File.makeDirectory(analysis_dir); 
mask_dir = analysis_dir+"/masks/";
File.makeDirectory(mask_dir); 
res_dir =analysis_dir+"/results/";
File.makeDirectory(res_dir); 

//get image file list
list = getFileList(maxp_dir);

//open threshold tool
run("Threshold...");

//loop for processing of all images in input folder

for(i = 0; i<list.length; i++) {
open(maxp_dir+list[i]);
print(list[i]);

//Set DAPI channel
Stack.setChannel(Ch_DAPI);
//gaussian blur and enhanced contrast to de-noise low intensity dapi
run("Gaussian Blur...", "sigma=2");
run("Enhance Contrast", "saturated=0.35");
run("Apply LUT","slice");

//background subtraction and median filter
run("Subtract Background...", "rolling=50 slice");
run("Median...", "radius=1 slice");

//let user set threshold
waitForUser("Set mask","Threshold set?");
getThreshold(lower, upper);
setThreshold(lower, upper);
setOption("BlackBackground", false);

//convert to mask, watershed, identify nuclei
run("Convert to Mask", "method=Default background=Light calculate only");
run("Watershed", "slice");
run("Analyze Particles...", "size=20-Infinity clear include summarize add slice");

//save ROI-set with nuclei as mask in "masks" folder, close modified image
roiManager("Save", mask_dir+list[i]+"_RoiSet.zip");
close();

//open original image, measure all slices (multi-measure), save in results directory
open(maxp_dir+list[i]);
roiManager("multi-measure measure_all");
saveAs("Results", res_dir+list[i]+"_results.csv");
roiManager("Delete");
close();

};

