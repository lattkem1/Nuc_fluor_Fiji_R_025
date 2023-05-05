
dir = getDirectory("Choose Main Directory (Images need to be opened)");
maxp_dir = dir+"/01_MaxProj greyscale/";
File.makeDirectory(maxp_dir); 
// get image IDs of all open images
ids=newArray(nImages);
for (i=1;i<=nImages;i++) {
        selectImage(i);
		Stack.setDisplayMode("grayscale");
        title = getTitle;
        print(title);
        saveAs("tiff", maxp_dir+title);
        close();
}

