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
in_dir = "./03_R input"
scr_dir = "./00_scripts"
out_dir = "./04b_intens distrib_2Ch"
if (!dir.exists(out_dir)) {dir.create(out_dir)}
setwd(main_dir)

###################    Load packages and define data   ############################
library("tidyverse")
library(ggplot2)

load(file = paste(in_dir,"/R rawdata set.RData",sep=""))

#define groups in csv sheet
group_tab0 = read_csv(paste(in_dir,"/image_groups.csv", sep =""))
group_tab = group_tab0[!is.na(group_tab0$group),]


####define analysis table (Channels, Group order) from original table (int.tab)

#define channel for analysis
set0$channels
Ch1 = c("BLBP", "Ki67")

#select groups and group order and color scheme
glist0 = unique(group_tab$group)
glist0
glist1 = glist0[!is.na(glist0)]
glist = glist1

gr_colors = glist
glist
gr_colors[] = "green"
gr_colors[c()] = c("grey")

gr_colors




################ functions ###############





#################  main   #######################

### generate working dataset set1, define groups, replicates, extract relevant channel and files

set1 = set0

set1$channels = Ch1

set1$group_tab = group_tab
set1$groups = glist
repl = paste(group_tab$group,group_tab$repl, sep="_")
set1$group_tab$repl = repl
set1$replicates = unique(repl)

t1 = set0$intens_tab[set0$intens_tab$file %in% group_tab$file,]
set1$intens_tab = as_tibble(t1[,c(names(t1)[1:4], Ch1)])

set1$intens_tab


#reorder group table by group order
t1=set1$group_tab
t2=NULL
i=glist[1]
for (i in glist){
  files_in_group =  t1$file[t1$group == i]  
  t2 = rbind(t2, t1[t1$file %in% files_in_group,])
}
set1$group_tab = t2


#reorder ROI intensity table and add group information
gr.t = set1$group_tab
int.tab = set1$intens_tab
files = gr.t$file

t2 = NULL

i=2
for (i in 1:length(files)){
  t1 = int.tab[int.tab$file == files[i],]
  t1$group =  gr.t[gr.t$file == files[i],]$group
  t1$repl =  gr.t[gr.t$file == files[i],]$repl
  t1$color = gr_colors[glist == t1$group[1]]
  t2=as_tibble(rbind(t2,t1))
}

set1$intens_tab = t2


#save processed data for channel for further analysis
save(set1,
     file = paste(out_dir,"/dataset ",Ch1[1]," ", Ch1[2],".RData",sep=""))



### QC scatterplots by image/replicate/group


file.name = paste(out_dir, "/QC scatterplots ", Ch1[1]," ", Ch1[2]," single cell intens by file.pdf", sep = "")
pdf(file = file.name, width = length(glist)*1.5, height = length(files)%/%length(glist)*7.5)

par(mfrow=c(length(glist), length(files)%/%length(glist)),lwd = 0.3, cex = 0.5) 
for (i in files){
  
  y.val = as_vector(set1$intens_tab[set1$intens_tab$file == i,Ch1[1]])
  x.val = as_vector(set1$intens_tab[set1$intens_tab$file == i,Ch1[2]])
  
  v1 = strsplit(i, split=".lif - ", fixed=TRUE)[[1]][2]
  plot.title = sub("_results.csv", "", v1)
  
  plot(x.val, y.val,
       ylab = paste(Ch1[1], "intensity (A.U.)"), 
       xlab = paste(Ch1[2], "intensity (A.U.)"), 
       ylim = c(0, 260), xlim=c(0, 260),
       pch = 21, bg = "green", #pch type of blotting point (see "points"), bg (background see "par")
       panel.first=grid(), cex=0.7, 
       main = plot.title) 
  par(axis(1, at = seq(0, 260, by = 10), tick = TRUE, labels = FALSE, lwd = 0.3, cex = 0.5),
      axis(2, at = seq(0, 260, by = 10), tick = TRUE, labels = FALSE, lwd = 0.3, cex = 0.5))
  
}

dev.off()





file.name = paste(out_dir, "/QC scatterplots ", Ch1[1]," ", Ch1[2]," single cell intens by repl.pdf", sep = "")
pdf(file = file.name, width = length(glist)*1.5, height = length(files)%/%length(glist)*7.5)

par(mfrow=c(length(glist), length(files)%/%length(glist)),lwd = 0.3, cex = 0.5) 

for (i in set1$replicates){
  
  y.val = as_vector(set1$intens_tab[set1$intens_tab$repl == i,Ch1[1]])
  x.val = as_vector(set1$intens_tab[set1$intens_tab$repl == i,Ch1[2]])
  
  plot.title = i
  
  plot(x.val, y.val,
       ylab = paste(Ch1[1], "intensity (A.U.)"), 
       xlab = paste(Ch1[2], "intensity (A.U.)"), 
       ylim = c(0, 260), xlim=c(0, 260),
       pch = 21, bg = "green", #pch type of blotting point (see "points"), bg (background see "par")
       panel.first=grid(), cex=0.7, 
       main = plot.title) 
  par(axis(1, at = seq(0, 260, by = 10), tick = TRUE, labels = FALSE, lwd = 0.3, cex = 0.5),
      axis(2, at = seq(0, 260, by = 10), tick = TRUE, labels = FALSE, lwd = 0.3, cex = 0.5))
  
}

dev.off()



file.name = paste(out_dir, "/QC scatterplots ", Ch1[1]," ", Ch1[2]," single cell intens by groups.pdf", sep = "")
pdf(file = file.name, width = length(glist)* 1.5, height = length(files)%/%length(glist)*7.5)

par(mfrow=c(length(glist), length(files)%/%length(glist)),lwd = 0.3, cex = 0.5) 

for (i in glist){
  
  y.val = as_vector(set1$intens_tab[set1$intens_tab$group == i,Ch1[1]])
  x.val = as_vector(set1$intens_tab[set1$intens_tab$group == i,Ch1[2]])
  
  plot.title = i
  
  plot(x.val, y.val,
       ylab = paste(Ch1[1], "intensity (A.U.)"), 
       xlab = paste(Ch1[2], "intensity (A.U.)"), 
       ylim = c(0, 260), xlim=c(0, 260),
       pch = 21, bg = "green", #pch type of blotting point (see "points"), bg (background see "par")
       panel.first=grid(), cex=0.7, 
       main = plot.title) 
  par(axis(1, at = seq(0, 260, by = 10), tick = TRUE, labels = FALSE, lwd = 0.3, cex = 0.5),
      axis(2, at = seq(0, 260, by = 10), tick = TRUE, labels = FALSE, lwd = 0.3, cex = 0.5))
  
}

dev.off()











