### analysis 2 channels nuclear intensities by cell

# go to main directory (parent directory of scripts)
if (basename(getwd())== "00_scripts"){setwd("../.")} 

#load packages

library("tidyverse")
library(colorRamps)

#define working directories

in_dir = "./03_R input/"
out_dir = "./04b_intens_distrib_2Ch/"
if (!dir.exists(out_dir)) {dir.create(out_dir)}


#### load data and set parameters 

load(file = paste0(in_dir,"/R rawdata set.RData"))

#load groups defined in csv sheet

group_tab = read_csv(paste0(in_dir,"/image_groups.csv"))

#define analysis table (Channels, Group order) from original table (int.tab)

#define channel for analysis
set0$channels
Ch = c("HuCD", "PAX6")

#select groups and group order and color scheme
gr = unique(group_tab$group)

gr_colors = matlab.like(6)
if (length(gr)>6) {gr_colors = matlab.like(length(gr))}




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
set1$channels = Ch
set1$groups = gr
set1$replicates = replics
set1$files = group_tab$file
set1$group_tab = group_tab
set1$intens_tab = t2

#save processed data for channel for further analysis
save(set1,
     file = paste(out_dir,"/dataset ",Ch[1]," ", Ch[2],".RData",sep=""))




#####################################################################
#  Scatterplots nuclear intensity by image/replicate/group  
#####################################################################

# basic plot

t1 = set1$intens_tab

p1 = ggplot(data = t1, aes(x=as_vector(t1[Ch[1]]), 
                           y = as_vector(t1[Ch[2]]),
                           color = group))+geom_point(size = 0.5, alpha = 0.5)+
  scale_color_manual(limits = gr, values = gr_colors)+
  theme_bw()


# save plot by group

p2 = p1 + facet_wrap(~factor(group, levels = gr))+
  labs(title = paste0("Nuclear intensities ", Ch[1], " - ", Ch[2], " by group"), 
       x = Ch[1], y = Ch[2])

file.name = paste0(out_dir, "/Intensity_scatterplots_", Ch[1],"_", Ch[2],"_by_group.pdf")

pdf(file = file.name, width = length(gr)*2+1, height = 3)
p2
dev.off()


# save plot by replicate

p2 = p1 + facet_wrap(~factor(repl, levels = replics), ncol = 4, nrow = ceiling(length(replics)/ 4))+
  labs(title = paste0("Nuclear intensities ", Ch[1], " - ", Ch[2], " by replicate"), 
       x = Ch[1], y = Ch[2])

file.name = paste0(out_dir, "/Intensity_scatterplots_", Ch[1],"_", Ch[2],"_by_repl.pdf")

pdf(file = file.name, width = 9, height = ceiling(length(replics)/ 4)*2+1)
p2
dev.off()


# save plot by file

files = unique(t1$file)

p2 = p1 + facet_wrap(~factor(file, levels = files), ncol = 4, nrow = ceiling(length(files)/ 4))+
  labs(title = paste0("Nuclear intensities ", Ch[1], " - ", Ch[2], " by file"), 
       x = Ch[1], y = Ch[2])

file.name = paste0(out_dir, "/Intensity_scatterplots_", Ch[1],"_", Ch[2],"_by_file.pdf")

pdf(file = file.name, width = 9, height = ceiling(length(files)/ 4)*2+1)
p2
dev.off()












