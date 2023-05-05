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
in_dir = "./04b_intens distrib_2Ch"
scr_dir = "./00_scripts"
out_dir = "./05c_treshold analysis_2Ch_by Ch1 n_p"
if (!dir.exists(out_dir)) {dir.create(out_dir)}
setwd(main_dir)

###################    Load packages and define data   ############################
library("tidyverse")
library(ggplot2)

load(file = paste(in_dir,"/dataset V5_GFP GS.RData",sep=""))


#################  parameter settings   #######################

Ch = set1$channels
Ch
Ch.thresh = c(30, 50)

#select groups and group order and color scheme
glist0 = unique(set1$group_tab$group)
glist0
glist = glist0

gr_colors = glist
glist
gr_colors[] = "grey"
gr_colors[c()] = c()

gr_colors


############# functions ####################

########function: count cells all/above/below specified thresholds 

count.tab.2Ch = function(s1, Ch.thresh){
  
  t1 = s1$intens_tab
  gr.tab = s1$group_tab
  flist1 = gr.tab$file
  
  #define vectors for count numbers 
  N.all = vector(length = length(flist1)) 
  N.nn = N.all
  N.np = N.all
  N.pn = N.all
  N.pp = N.all
  
  #subset and count all/positive/negative cells
  for (i in 1:length(flist1)){
    t2 = t1[t1$file==flist1[i],]
    N.all[i]  = nrow(t2)
    Ch1.int = as.numeric(as_vector(t2[,Ch[1]]))
    Ch2.int = as.numeric(as_vector(t2[,Ch[2]]))
    N.nn[i] = nrow(t2[(Ch1.int < Ch.thresh[1] & Ch2.int < Ch.thresh[2]),])
    N.np[i] = nrow(t2[(Ch1.int < Ch.thresh[1] & Ch2.int >= Ch.thresh[2]),])
    N.pn[i] = nrow(t2[(Ch1.int >= Ch.thresh[1] & Ch2.int < Ch.thresh[2]),])
    N.pp[i] = nrow(t2[(Ch1.int >= Ch.thresh[1] & Ch2.int >= Ch.thresh[2]),])
                              
  }
  
  file.sum = as_tibble(cbind(gr.tab, N.all, N.nn, N.np, N.pn, N.pp)) 
  
  return(file.sum)
  
}



### function: collapse (summarize) to replicates

merge_files_to_repl = function(tab1){
  #tab1: count table listing counts of single images/files
  
  r.count.tab = NULL
  
  rlist = unique(tab1$repl)
  
  N.files = vector(length = length(rlist))
  for (i in 1:length(rlist)){
    sub.tab = tab1[tab1$repl == rlist[i],]
    sum.sub.tab = sub.tab[1,]
    sum.sub.tab[1 ,-c(1:3)] = apply(sub.tab[,-c(1:3)],2,sum)
    sum.sub.tab
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
    geom_bar(aes(x = x, y = mean_y, fill = color), 
             stat="identity", width = 0.6, colour = "black") + 
    geom_dotplot(aes(x = x, y = y), position = position_dodge(0.6),
                 binaxis = "y", stackdir = "center", dotsize = 0.005*y.max, binwidth = 2) + 
    geom_errorbar(aes(x = x, ymax = mean_y+sd_y, ymin = mean_y-sd_y), width = 0.3)+
    labs(title=title, x="", y=y_title, fill = NULL)   +
    geom_hline(yintercept = 0, colour = "darkgrey")+
    scale_fill_manual(limits = unique(color), values = unique(color))+
    guides(fill = "none", color = "none")+
    theme_classic()+
    theme(axis.line.y = element_line(colour="black"),
          axis.line.x = element_line(colour="black"),
          text = element_text(colour="black", size=14, face="bold"),
          plot.title = element_text(colour="black", size=20, face="bold"),
          axis.title = element_text(colour="black", size=16, face="bold"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0))
 
  return(gg.bar.scatter) 
}







#################  main   #######################

  
### generate table with cell counts per image/replicate (total and above threshold)
  Ch1.counts.f = count.tab.2Ch(set1, Ch.thresh = Ch.thresh)
  Ch1.counts.r = merge_files_to_repl(Ch1.counts.f)
  
  write_csv(Ch1.counts.f, path = 
              paste(out_dir,"/",exp.name, "_counts per file_",Ch[1]," ", Ch.thresh[1],
                    "_",Ch[2]," ", Ch.thresh[2],".csv", sep=""))
  write_csv(Ch1.counts.r, path = 
              paste(out_dir,"/",exp.name, "_counts per repl_",Ch[1]," ", Ch.thresh[1],
                    "_",Ch[2]," ", Ch.thresh[2],".csv", sep=""))
  
  
### calculate fraction of cells for pos/neg for both channels combined, add columns for group statistics 
  count.stat = mutate(Ch1.counts.r, 
                      N_field = N.all/N.files, 
                      frac_1n = (N.nn+N.np)/N.all, 
                      frac_1p = (N.pn+N.pp)/N.all,
                      frac_2p_of_1n = N.np/(N.nn+N.np),
                      frac_2p_of_1p = N.pp/(N.pn+N.pp),
                      mean_N_field = NA, sd_N_field = NA, 
                      mean_frac_1n = NA, sd_frac_1n = NA,
                      mean_frac_1p = NA, sd_frac_1p = NA,
                      mean_frac_2p_of_1n = NA, sd_frac_2p_of_1n = NA,
                      mean_frac_2p_of_1p = NA, sd_frac_2p_of_1p = NA
                      )
  
  #calculate group statistics
  for (j in count.stat$group){
    t1 = count.stat[count.stat$group == j,]
    t1$mean_N_field[1] = mean(t1$N_field)
    t1$sd_N_field[1] = sd(t1$N_field)
    t1$mean_frac_1n[1] = mean(t1$frac_1n)
    t1$sd_frac_1n[1] = sd(t1$frac_1n)
    t1$mean_frac_1p[1] = mean(t1$frac_1p)
    t1$sd_frac_1p[1] = sd(t1$frac_1p)
    
    t1$mean_frac_2p_of_1n[1] = mean(t1$frac_2p_of_1n)
    t1$sd_frac_2p_of_1n[1] = sd(t1$frac_2p_of_1n)
    t1$mean_frac_2p_of_1p[1] = mean(t1$frac_2p_of_1p)
    t1$sd_frac_2p_of_1p[1] = sd(t1$frac_2p_of_1p)

    count.stat[count.stat$group == j,]= t1
  }

  write_csv(count.stat, path = 
              paste(out_dir,"/",exp.name, "_count stat thresh_",Ch[1]," ", Ch.thresh[1],
                    "_",Ch[2]," ", Ch.thresh[2],".csv", sep="")
            )
  


  
  
  
  


##### generate plots

### generate table for bar+scatterplot by groups by n/p populations (gr1_nn/np/pn/pp-gr2_nn...)
# labels: gr_nn/np/pn/pp, columns: label/single replicate values/mean/sd/color
  
t1 = count.stat
#remove absolute numbers
t2 = t1[,-c(4:10)]

t4 = NULL

for (i in 1:length(glist)){
    t3 = t2[t2$group == glist[i],]
    frac_2p_of_1n = as_tibble(cbind(label = paste(glist[i],"_",Ch[1],"_neg",sep=""), repl = t3$frac_2p_of_1n,
                         mean = t3$mean_frac_2p_of_1n, sd = t3$sd_frac_2p_of_1n))
    frac_2p_of_1p = as_tibble(cbind(label = paste(glist[i],"_",Ch[1],"_pos",sep=""), repl = t3$frac_2p_of_1p,
                                    mean = t3$mean_frac_2p_of_1p, sd = t3$sd_frac_2p_of_1p))
    t3_1 = as_tibble(rbind(frac_2p_of_1n,frac_2p_of_1p))
    t3_2 = as_tibble(cbind(t3_1, color = gr_colors[i]))
    t4 = as_tibble(rbind(t4,t3_2))
}

t4[,-c(1,5)] = apply(t4[,-c(1,5)],2,as.numeric)
t4[,5] = as.character(t4$color)


#plot scatter/barplot

file.name = paste(out_dir, "/fract ",Ch[2]," pos of ",Ch[1],
                  " pos_neg_thresh ", Ch.thresh[2],"_", Ch.thresh[1],".pdf", sep = "")

#plot environment: width groups number, height image number/group
pdf(file = file.name, width = 12, height = 8)
  
  pl.title = paste("% ",Ch[2]," pos (thresh = ",Ch.thresh[2],") of ",Ch[1],
                  " pos/neg cells (thresh = ", Ch.thresh[1],")", sep = "")
  y.title = paste("% " ,Ch[2]," pos")
  
  gg = bar_scatterplot_dyn(x = t4$label, y = t4$repl*100, 
                           mean_y = t4$mean*100, sd_y = t4$sd*100,
                           color = t4$color, 
                           title = pl.title,
                           y_title = y.title)
  
  plot(gg)

dev.off()


