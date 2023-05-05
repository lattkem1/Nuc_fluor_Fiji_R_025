# read out current directory and set parent directory of "R scripts" folder(scr_dir) as 
# main working directory; warn if R Scripts is not current working directory
getwd()
basename(getwd())
if (basename(getwd()) == "00_scripts"){
  scr_dir = getwd()
  setwd("./..")
  main_dir = getwd()
} else {readline("Check current working directory")}

########define working directories##########################
main_dir = getwd()
exp.name = basename(getwd())
in_dir = "./02_ImageJ analysis/results"
scr_dir = "./00_scripts"
out_dir = "./03_R input"
if (!dir.exists(out_dir)) {dir.create(out_dir)}
setwd(main_dir)

###################    Load packages and define data   ############################
library("tidyverse")

#define channel names
Ch.names = c("DAPI", "HuCD","PAX6")

#input filelist
flist = list.files(in_dir)



#################  main   #######################

########### read in data

int.tab = NULL

#perform for each file (row in group.tab)
i=1
  for (i in 1:length(flist)){
    
    #read file and define ordering parameters
    fname = paste(in_dir,"/",flist[i], sep ="")
    rawdata = read_csv(fname)
    
    #generate table of intensities for each cell for each channel, add group, replicate, file, color in first columns
    int.tab.fl = filter(rawdata, rawdata[,"Ch"] == 1)[,"Mean"]
    if (length(Ch.names) > 1){
      for (h in 2:length(Ch.names)){
        int.tab.fl = cbind(int.tab.fl, filter(rawdata, rawdata[,"Ch"] == h)[,"Mean"])
      }
    }
    int.tab.fl = cbind(group = NA, repl = NA, file = flist[i], color = NA, int.tab.fl)
    
    #add to table with all cells of all files
    int.tab = rbind(int.tab, int.tab.fl)
  }


names(int.tab)[-c(1:4)] = Ch.names

#convert file names from factor to character (else problems with downstream analysis) 
int.tab$file = as.character(int.tab$file) 

#list of files with group assignment, and list of groups
gr.tab = tibble(file = flist, group ="", repl="")

#create ROI_set object (raw intensity dataset)
set0 = list(channels = Ch.names, 
            groups = NULL, 
            replicates = NULL,
            files = flist,
            group_tab = gr.tab,
            intens_tab = int.tab
            )
class(set0) = "ROI_set"


#save group table (keep copy of old version, if existing)
if (file.exists(paste(out_dir,"/image_groups.csv", sep=""))){
  gr.tab.old = read_csv(paste(out_dir,"/image_groups.csv", sep=""))
  write_csv(gr.tab.old, path = paste(out_dir,"/image_groups_", Sys.time(),".csv", sep=""))
}
write_csv(gr.tab, path = paste(out_dir,"/image_groups.csv", sep=""))


#save processed input data for further analysis
save(set0, 
     file = paste(out_dir,"/R rawdata set.RData",sep=""))

