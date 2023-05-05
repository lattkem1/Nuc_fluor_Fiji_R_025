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
out_dir = "./04_intens distrib_1Ch"
if (!dir.exists(out_dir)) {dir.create(out_dir)}
setwd(main_dir)

###################    Load packages and define data   ############################
library("tidyverse")
library(ggplot2)

load(file = paste(in_dir,"/R rawdata set.RData",sep=""))

#define groups in csv sheet
group_tab = read_csv(paste(in_dir,"/image_groups.csv", sep =""))



####define analysis table (Channels, Group order) from original table (int.tab)

#define channel for analysis
set0$channels
Ch1 = "GFAP"



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

QCplot = function(t1, type = "scatter", plot_by = "file", ...) {
  t2 = t1$intens_tab 
  Ch1 = t1$channels
  #define plotting table labels
  if (plot_by == "file"){
    v1 = sapply(strsplit(as_vector(t2[,plot_by]), split=".lif - ", fixed=TRUE), function(x) (x[2]))
    x.labels = sapply(sub("_results.csv", "", v1), function(x) (x))
  } else {x.labels = as_vector(t2[,plot_by])}
  
  y.val = as.numeric(as_vector(t2[,Ch1]))
  
  if (type == "scatter"){
    gg = ggplot() + xlim(unique(x.labels)) + ylim(c(0,255)) +
      geom_dotplot(data = t2, aes(x = x.labels, 
                                  y = y.val,
                                  color = t2$color, 
                                  fill = t2$color), 
                   position = "dodge",
                   binaxis = "y", stackdir = "center", dotsize = 2, binwidth = 0.3) + 
      scale_color_manual(limits = unique(t2$color), values = unique(t2$color))+
      scale_fill_manual(limits = unique(t2$color), values = unique(t2$color))+
      guides(fill = "none", color = "none")+
      labs(title=paste(Ch1," intensities (scatterplot individual cells by ", plot_by ,")", sep = ""), x="", 
           y=paste(Ch1,"intensity (A.U.)"), color = NULL, fill = NULL) +
      theme(axis.line.y = element_line(colour="black"),
            axis.line.x = element_line(colour="black"),
            text = element_text(colour="black", size=14, face="bold"),
            plot.title = element_text(colour="black", size=20, face="bold"),
            axis.title = element_text(colour="black", size=10, face="bold"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0))  ## angeled labels: angle = 45, vjust = 1, hjust = 1
  }
  
  if(type == "violin"){
    gg0 = ggplot() + xlim(unique(x.labels)) + ylim(c(0,255)) + 
      labs(title=paste(Ch1,"distribution","by",plot_by), x="", y=paste(Ch1, "Intensity (A.U)")) 
    gg = gg0 +
      geom_violin(data = t2, aes(x = x.labels, y = y.val)) + 
      geom_boxplot(data = t2, aes(x = x.labels, y = y.val, fill = t2$color),
                   outlier.color = "grey", outlier.size = 0, 
                   size = 0.5, width = 0.1)+
      scale_color_manual(limits = unique(t2$color), values = unique(t2$color))+
      scale_fill_manual(limits = unique(t2$color), values = unique(t2$color))+
      guides(fill = "none", color = "none")+
      theme(axis.line.y = element_line(colour="black"),
            axis.line.x = element_line(colour="black"),
            text = element_text(colour="black", size=14, face="bold"),
            plot.title = element_text(colour="black", size=20, face="bold"),
            axis.title = element_text(colour="black", size=16, face="bold"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0))  ## angeled labels: angle = 45, vjust = 1, hjust = 1
  }
  
  
  
  
  return(gg)
  
}





#################  main   #######################

### generate working dataset set1, define groups, replicates, extract relevant channel

set1 = set0
set1$channels = Ch1
set1$group_tab = group_tab
set1$groups = glist
repl = paste(group_tab$group,group_tab$repl, sep="_")
set1$group_tab$repl = repl
set1$replicates = unique(repl)
set1$intens_tab = as_tibble(set0$intens_tab[,c(names(set0$intens_tab)[1:4], Ch1)])

set1$intens_tab

#reorder group table by group order
t1=set1$group_tab
t2=NULL
t.gr = t1[!is.na(t1$group),]
i=1
for (i in 1:length(glist)){
  files_in_group =  t.gr$file[t.gr$group == glist[i]]  
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
     file = paste(out_dir,"/dataset ",Ch1,".RData",sep=""))



### QC plots by image/replicate/group

pl1 = QCplot(set1, type = "scatter",plot_by = "file")
pl2 = QCplot(set1, type = "violin",plot_by = "file")

pl3 = QCplot(set1, type = "scatter",plot_by = "repl")
pl4 = QCplot(set1, type = "violin",plot_by = "repl")

pl5 = QCplot(set1, type = "scatter",plot_by = "group")
pl6 = QCplot(set1, type = "violin",plot_by = "group")


#save scatterplot, violinplot and boxplot by image/replicate/group

file.name = paste(out_dir, "/QC plots ", Ch1, " single cell intens by file.pdf", sep = "")
pdf(file = file.name, width = 8, height = 6)

plot(pl1)
plot(pl2)

# plot as histogram
i = files[3]

par(mfrow=c(4, length(files)%/%4+1),lwd = 0.3, cex = 0.5) 
for (i in files){

  y.val = as_vector(set1$intens_tab[set1$intens_tab$file == i,Ch1])
  
  v1 = strsplit(i, split=".lif - ", fixed=TRUE)[[1]][2]
  plot.title = sub("_results.csv", "", v1)
  
  hist(y.val, breaks = seq(0, 260, by = 2), col="grey", 
       main = plot.title,
       freq=TRUE, xlim=c(0,260), xlab=paste(Ch1, "Intesity (A.U.)"))
  par(axis(1, at = seq(0, 260, by = 10), tick = TRUE, labels = FALSE, lwd = 0.3, cex = 0.5))
}

dev.off()



file.name = paste(out_dir, "/QC plots ", Ch1, " single cell intens by repl.pdf", sep = "")
pdf(file = file.name, width = 6, height = 6)

plot(pl3)
plot(pl4)

par(mfrow=c(4, length(set1$replicates)%/%4+1),lwd = 0.3, cex = 0.5) 
for (i in set1$replicates){
  y.val = as_vector(set1$intens_tab[set1$intens_tab$repl == i,Ch1])
  hist(y.val, breaks = seq(0, 260, by = 2), col="grey", 
       main = i,
       freq=TRUE, xlim=c(0,260), xlab=paste(Ch1, "Intesity (A.U.)"))
  par(axis(1, at = seq(0, 260, by = 10), tick = TRUE, labels = FALSE, lwd = 0.3, cex = 0.5))
}

dev.off()



file.name = paste(out_dir, "/QC plots ", Ch1, " single cell intens by group.pdf", sep = "")
pdf(file = file.name, width = 5, height = 6)

plot(pl5)
plot(pl6)

par(mfrow=c(4, length(glist)%/%4+1),lwd = 0.3, cex = 0.5) 
for (i in glist){
  y.val = as_vector(set1$intens_tab[set1$intens_tab$group == i,Ch1])
  hist(y.val, breaks = seq(0, 260, by = 2), col="grey", 
       main = i,
       freq=TRUE, xlim=c(0,260), xlab=paste(Ch1, "Intesity (A.U.)"))
  par(axis(1, at = seq(0, 260, by = 10), tick = TRUE, labels = FALSE, lwd = 0.3, cex = 0.5))
}

dev.off()





