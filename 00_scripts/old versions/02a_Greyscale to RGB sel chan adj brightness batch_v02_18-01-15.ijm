
// get image IDs of all open images
dir = getDirectory("Choose Main Directory (Images need to be opened)");

rgb_dir = dir+"/02a_MaxProj_RGB_DAPI GFP Sox9/";

File.makeDirectory(rgb_dir); 

ids=newArray(nImages);
for (i=1;i<=nImages;i++) {
        selectImage(i);

		setMinAndMax(0, 255);
		Stack.setChannel(3);
		run("Blue");
		Stack.setChannel(1);
		run("Green");
		Stack.setChannel(2);
		run("Magenta");
		Stack.setDisplayMode("composite");
		Stack.setActiveChannels("111");
		run("Stack to RGB");

        title = getTitle;
        print(title);
        saveAs("tiff", rgb_dir+title);
        close();
}


