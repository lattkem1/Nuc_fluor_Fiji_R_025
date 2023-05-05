
// get image IDs of all open images
dir = getDirectory("Choose Main Directory (Images need to be opened)");

rgb_dir = dir+"/02a_MaxP_RGB_Sox9 Ki67 BLBP/";

File.makeDirectory(rgb_dir); 

ids=newArray(nImages);
for (i=1;i<=nImages;i++) {
        selectImage(i);

		Stack.setChannel(1);
		setMinAndMax(0, 255);
		run("Blue");
		
		Stack.setChannel(3);
		setMinAndMax(0, 255);
		run("Green");
		
		Stack.setChannel(2);
		setMinAndMax(0, 255);
		run("Magenta");
		
		Stack.setDisplayMode("composite");
		Stack.setActiveChannels("1111");
		run("Stack to RGB");

        title = getTitle;
        print(title);
        saveAs("tiff", rgb_dir+title);
        close();
}


