
dir = getDirectory("Choose Main Analysis Directory (Images need to be in '00_input' subdirectory)");
in_dir = dir+"/00_input/";
maxp_dir = dir+"/01_MaxProj greyscale/";
File.makeDirectory(maxp_dir); 
// get image IDs of all open images
ids = getFileList(in_dir);

for(i = 0; i<ids.length; i++) {
        run("Bio-Formats Importer", "open=["+in_dir+ids[i]+"] color_mode=Default display_metadata rois_import=[ROI manager] view=[Standard ImageJ] stack_order=Default");
		run("8-bit");
		run("Make Composite", "display=Grayscale");
        title = getTitle;
        print(title);
        saveAs("tiff", maxp_dir+title);
        close();
}

