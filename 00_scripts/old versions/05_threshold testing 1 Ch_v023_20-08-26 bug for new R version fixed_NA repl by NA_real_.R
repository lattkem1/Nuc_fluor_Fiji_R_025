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
in_dir = "./04_intens distrib_1Ch"
scr_dir = "./00_scripts"
out_dir = "./05_treshold tests_1Ch"
if (!dir.exists(out_dir)) {dir.create(out_dir)}
setwd(main_dir)

###################    Load packages and define data   ############################
library("tidyverse")
library(ggplot2)

Ch1 = "V5"

load(file = paste(in_dir,"/dataset ",Ch1,".RData",sep=""))


#################  parameter settings   #######################

Ch1.test.thresh = c(50, 60)


glist0 = unique(set1$group_tab$group)
glist0
glist = glist0

gr_colors = glist
glist
gr_colors[] = "green"
gr_colors[c()] = c("grey")

gr_colors



############# functions ####################

########function: count cells (in tab1) all/above specified threshold, calculate mean intensity 

count.tab.1Ch = function(s1, Ch1.thresh){
  
  t1 = s1$intens_tab
  gr.tab = s1$group_tab
  flist1 = gr.tab$file
  
  #define vectors for mean intensities and count numbers 
  N.all = vector(length = length(flist1)) 
  N.1p = N.all
  mean.int = N.all
  
  #subset and count all/positive/negative cells
  for (i in 1:length(flist1)){
    file1 = t1[t1$file==flist1[i],]
    mean.int[i] = mean(as_vector(file1[,Ch1]))
    file1.1p = file1[as_vector(file1[,Ch1])>=Ch1.thresh,]
    N.all[i]  = nrow(file1)
    N.1p[i]  = nrow(file1.1p)
  }
  
  file.sum = as_tibble(cbind(gr.tab, mean.int, N.all, N.1p)) 
  
  return(file.sum)
  
}



### function: collapse (summarize) to replicates

collapse.repl = function(tab1){
  #tab1: count table listing counts of single images/files
  
  r.count.tab = NULL
  
  rlist = unique(tab1$repl)
  
  N.files = vector(length = length(rlist))
  for (i in 1:length(rlist)){
    sub.tab = tab1[tab1$repl == rlist[i],]
    sum.sub.tab = sub.tab[1,]
    sum.sub.tab[1 ,"N.all"] = sum(sub.tab[,"N.all"])
    sum.sub.tab[1 ,"mean.int"] = mean(as_vector(sub.tab[,"mean.int"]))
    sum.sub.tab[1 ,"N.1p"] = sum(sub.tab[,"N.1p"])
    r.count.tab = as_tibble(rbind(r.count.tab,sum.sub.tab))
    N.files[i] = nrow(sub.tab)
  }
  
  r.count.tab2 = as_tibble(cbind(r.count.tab, N.files))
  
  return(r.count.tab2)
  
}


bar_scatterplot_dyn = function(x,y,mean_y, sd_y, color = "grey", title = "", y_title = "", ...){
  
  #determine y range
  y.min = min(na.omit(mean_y-sd_y))
  y.max = max(na.omit(mean_y+sd_y))
  #if single values are beyond mean +/- sd range, or no sd is given, adjust accordingly
  if (min(na.omit(y))<y.min | is.infinite(y.min)){y.min = min(na.omit(y))}
  if (max(na.omit(y))>y.max | is.infinite(y.max)){y.max = max(na.omit(y))}
  if (y.min > 0){y.min = 0}
  
  gg.bar.scatter = ggplot() + xlim(unique(x)) + ylim(c(y.min, y.max)) +
    geom_bar(aes(x = x, y = mean_y), 
             stat="identity", width = 0.6, colour = "black", fill = gr_colors) + 
    geom_dotplot(aes(x = x, y = y), position = position_dodge(0.6),
                 binaxis = "y", stackdir = "center", dotsize = 0.015*y.max, binwidth = 2) + 
    geom_errorbar(aes(x = x, ymax = mean_y+sd_y, ymin = mean_y-sd_y), width = 0.3)+
    labs(title=title, x="", y=y_title, fill = NULL)   +
    geom_hline(yintercept = 0, colour = "darkgrey")+
    theme_classic()+
    theme(axis.line.y = element_line(colour="black"),
          axis.line.x = element_line(colour="black"),
          text = element_text(colour="black", size=14, face="bold"),
          plot.title = element_text(colour="black", size=20, face="bold"),
          axis.title = element_text(colour="black", size=16, face="bold")) 
 
  return(gg.bar.scatter) 
}







#################  main   #######################


### generate list of statistics tables for each test-threshold

#list of stats for individual test thresholds
test.thresh.stat.list = vector("list", length = length(Ch1.test.thresh))

#run analysis for each test threshold
i=5
for (i in 1:length(Ch1.test.thresh)){
  
  #generate table with cell counts per image/replicate (total and above threshold)
  Ch1.counts.f = count.tab.1Ch(set1, Ch1.thresh = Ch1.test.thresh[i])
  Ch1.counts.r = collapse.repl(Ch1.counts.f)

  #fraction calculations, add columns for group statistics 
  count.stat = mutate(Ch1.counts.r, 
                             N_field = N.all/N.files, 
                             frac_1p = N.1p/N.all, 
                             mean_int_group = NA_real_, sd_int_group = NA_real_, 
                             mean_N_field = NA_real_, sd_N_field = NA_real_, 
                             mean_frac_1p = NA_real_, sd_frac_1p = NA_real_)
  
  #calculate group statistics
  for (j in count.stat$group){
    t1 = count.stat[count.stat$group == j,]
    t1$mean_int_group[1] = mean(t1$mean.int)
    t1$sd_int_group[1] = sd(t1$mean.int)
    t1$mean_N_field[1] = mean(t1$N_field)
    t1$sd_N_field[1] = sd(t1$N_field)
    t1$mean_frac_1p[1] = mean(t1$frac_1p)
    t1$sd_frac_1p[1] = sd(t1$frac_1p)
    count.stat[count.stat$group == j,] = t1
  }
  
  test.thresh.stat.list[[i]] = count.stat
  write_csv(count.stat, path = 
              paste(out_dir,"/",exp.name, "_count stat ",Ch1,"_thresh ", Ch1.test.thresh[i],".csv", sep="")
  )
  
}




##### generate plots


file.name = paste(out_dir, "/",Ch1," threshholding test plots.pdf", sep = "")

#plot environment: width groups number, height image number/group
pdf(file = file.name, width = 16, height = 8)

for (a in 1:length(Ch1.test.thresh)){
  
  #load specific plotting table for threshold a
  t1 = test.thresh.stat.list[[a]]
  
  pl.title = paste("% ",Ch1," pos, Treshold = ",Ch1.test.thresh[a], sep = "")
  y.title = paste("%",Ch1,"pos")
  
  gg = bar_scatterplot_dyn(x = t1$group, y = t1$frac_1p*100, 
                           mean_y = t1$mean_frac_1p*100, sd_y = t1$sd_frac_1p*100,
                           color = gr.colors, 
                           title = pl.title,
                           y_title = y.title)
  
  plot(gg)
  
}

dev.off()


