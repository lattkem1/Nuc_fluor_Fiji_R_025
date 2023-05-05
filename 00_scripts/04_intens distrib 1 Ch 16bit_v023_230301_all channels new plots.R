### analysis single channel nuclear intensities by cell

# go to main directory (parent directory of scripts)
if (basename(getwd())== "00_scripts"){setwd("../.")} 

#load packages

library("tidyverse")
library(colorRamps)

#define working directories

in_dir = "./03_R input/"
out_dir = "./04_intens_distrib_1Ch/"
if (!dir.exists(out_dir)) {dir.create(out_dir)}


#### load data and set parameters 

load(file = paste0(in_dir,"/R rawdata set.RData"))

#load groups defined in csv sheet

group_tab = read_csv(paste0(in_dir,"/image_groups.csv"))

#define analysis table (Channels, Group order) from original table (int.tab)

#define channel for analysis
Ch = set0$channels

#select groups and group order and color scheme
gr = unique(group_tab$group)

gr_colors = matlab.like(6)
if (length(gr)>6) {gr_colors = matlab.like(length(gr))}




################ functions ###############



#####################################################################
#  update and organise dataset   
#####################################################################


### define groups, replicates, extract relevant channel and files

group_tab$repl = paste(group_tab$group,group_tab$repl, sep="_")
replics = unique(group_tab$repl)

t1 = set0$intens_tab
t1 = t1[t1$file %in% group_tab$file,c(names(t1)[1:4], Ch)]

# simplify filenames (if format "[lif archive].lif - [image]"), keep only image name

if (all(grepl(".lif - ", group_tab$file)) ){
  
  l2 = lapply(group_tab$file, function(x){
    l1 = str_split(x, ".lif - ")
    return(l1[[1]][2])
  })
  group_tab$file = unlist(l2)
  
  l2 = lapply(t1$file, function(x){
    l1 = str_split(x, ".lif - ")
    return(l1[[1]][2])
  })
  t1$file = unlist(l2)
}

#assign group, replicate, color

t1$group = group_tab$group[match(t1$file, group_tab$file)]
t1$repl = group_tab$repl[match(t1$file, group_tab$file)]
t1$color = gr_colors[match(t1$group, gr)]

#reorder cells by group table order

t2 = t1[order(match(t1$file, group_tab$file)),]

#update ROI set, add group information
set1 = set0
set1$groups = gr
set1$replicates = replics
set1$files = group_tab$file
set1$group_tab = group_tab
set1$intens_tab = t2

#save processed data for channel for further analysis
save(set1,
     file = paste(out_dir,"/dataset processed.RData",sep=""))




################################################################
### Intensity histogrammes by image/replicate/group for all channels
################################################################

#plots by group

t1 = set1$intens_tab

pl_file = lapply(Ch, function(Ch1){
  p1 = ggplot(t1, aes(x = as_vector(t1[,Ch1]), fill = group))+geom_histogram(binwidth = 1000)+
    scale_x_continuous(limits = c(-1, 66000), breaks = seq(0, 66000, by = 10000), 
                       minor_breaks = seq(0, 66000, by = 1000))+theme_bw()+
    scale_fill_manual(limits = gr, values = gr_colors)+
    facet_wrap(~factor(group, levels = gr), ncol = 1)+
    labs(title = paste0("Intensity distribution by group - ", Ch1), x = "Intensity (a.u.)")
  return(p1)
})

file.name = paste0(out_dir, "/Intensity_histogrammes_by_group.pdf")

pdf(file = file.name, width = 5, height = length(groups)*1.5+2)
lapply(pl_file, function(x){x})
dev.off()


#plots by replicate

t1 = set1$intens_tab

pl_file = lapply(Ch, function(Ch1){
  p1 = ggplot(t1, aes(x = as_vector(t1[,Ch1]), fill = group))+geom_histogram(binwidth = 1000)+
    scale_x_continuous(limits = c(-1, 66000), breaks = seq(0, 66000, by = 10000), 
                       minor_breaks = seq(0, 66000, by = 1000))+theme_bw()+
    scale_fill_manual(limits = gr, values = gr_colors)+
    facet_wrap(~factor(repl, levels = replics), ncol = 1)+
    labs(title = paste0("Intensity distribution by replicate - ", Ch1), x = "Intensity (a.u.)")
  return(p1)
})

file.name = paste0(out_dir, "/Intensity_histogrammes_by_repl.pdf")

pdf(file = file.name, width = 5, height = length(replics)*1.5+2)
lapply(pl_file, function(x){x})
dev.off()


#plots by file

t1 = set1$intens_tab

files = unique(t1$file)

pl_file = lapply(Ch, function(Ch1){
  p1 = ggplot(t1, aes(x = as_vector(t1[,Ch1]), fill = group))+geom_histogram(binwidth = 1000)+
    scale_x_continuous(limits = c(-1, 66000), breaks = seq(0, 66000, by = 10000), 
                       minor_breaks = seq(0, 66000, by = 1000))+theme_bw()+
    scale_fill_manual(limits = gr, values = gr_colors)+
    facet_wrap(~factor(file, levels = files), ncol = 1)+
    labs(title = paste0("Intensity distribution by group - ", Ch1), x = "Intensity (a.u.)")
  return(p1)
})

file.name = paste0(out_dir, "/Intensity_histogrammes_by_file.pdf")

pdf(file = file.name, width = 5, height = length(files)*1.5+2)
lapply(pl_file, function(x){x})
dev.off()







