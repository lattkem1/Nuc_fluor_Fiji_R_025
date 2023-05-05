//specify input/output directory
in_dir = getDirectory("Choose Input Image (MaxP) Directory");
out_dir = File.getParent(in_dir)+"/MaxP_RGB_DAPI HuCD OCT4/";

//specify channel order, colors and intensity ranges
channels = newArray(1,2,3);
channel_colors = newArray("Blue","Green","Magenta");
min_intens = newArray(0,0,0);
max_intens = newArray(200,127,127);



// main

input_file_list = getFileList(in_dir);

File.makeDirectory(out_dir); 

// get image IDs of all open images
ids=newArray(nImages);

//loop through all images

for(i = 0; i<input_file_list.length; i++) {
        
		open(in_dir+input_file_list[i]);

		//allow cropping if required (with standard centered 512x512 cropping window)
        makeRectangle(255, 255, 256, 256);
        waitForUser("Set mask","Crop image as required and continue");
        run("Crop");

        //split image in individual channels for processing
		run("Split Channels");

		//adjust and save channel images
		for (j=0;j<channels.length;j++) {
			selectWindow("C"+channels[j]+"-"+input_file_list[i]);
			setMinAndMax(min_intens[j], max_intens[j]);
			run("Apply LUT");
			run(channel_colors[j]);
			run("RGB Color");
			saveAs("tiff", out_dir+input_file_list[i]+"_C"+channels[j]+".tif");
			close();
		}
		run("Close All");
		
		//reopen channels to create composite image
		for (j=0;j<channels.length;j++) {
			open(out_dir+input_file_list[i]+"_C"+channels[j]+".tif");
		}
		run("Images to Stack", "name=Stack title=[] use");
		run("Z Project...", "projection=[Sum Slices]");
		saveAs("tiff", out_dir+input_file_list[i]+"_merged.tif");
		run("Close All");
		
}


