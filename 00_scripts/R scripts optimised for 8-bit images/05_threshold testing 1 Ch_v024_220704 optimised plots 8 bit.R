### analysis nuclear intensities single channel with different thresholds

# go to main directory (parent directory of scripts)
if (basename(getwd())== "00_scripts"){setwd("../.")} 

#load packages

library("tidyverse")
library(colorRamps)
library(lmerTest)

#define working directories

in_dir = "./04_intens_distrib_1Ch/"
out_dir = "./05_intens_thresh_analysis_1Ch/"
if (!dir.exists(out_dir)) {dir.create(out_dir)}


#### load data and set parameters 

load(file = paste0(in_dir,"/dataset processed.RData"))

#define channel and thresholds for analysis
set1$channels

Ch = "PAX6"

thresh = c(50, 70, 100)




#############################################
#create summary statistics for each threshold
#############################################

t1 = set1$intens_tab
t1$Ch = t1[,Ch]

replics = unique(t1$repl)
gr = unique(t1$group)


#stats by file (with mean/sd by replicate and group)

stats_list_by_file = lapply(thresh, function(thresh1){
  
  t2 = t1 %>% group_by(group, repl, file) %>% summarise(N_cells = n(), 
                                                        N_pos = length(Ch[Ch>thresh1]))
  
  t2$fract_pos = t2$N_pos/t2$N_cells
  
  t3 = t2 %>% group_by(group, repl) %>% summarise(N_cells_repl_mean = mean(N_cells), 
                                                  N_cells_repl_sd = sd(N_cells), 
                                                  fract_pos_repl_mean = mean(fract_pos),
                                                  fract_pos_repl_sd = sd(fract_pos))
  
  t4 = t2 %>% group_by(group) %>% summarise(N_cells_group_mean = mean(N_cells), 
                                            N_cells_group_sd = sd(N_cells), 
                                            fract_pos_group_mean = mean(fract_pos),
                                            fract_pos_group_sd = sd(fract_pos))
  
  t5 = cbind(t2, t3[match(t2$repl, t3$repl),-c(1,2)], t4[match(t2$group, t4$group),-1])
  
  t5 = t5[order(match(t5$repl, replics)),]
  
  write_csv(t5, file = paste0(out_dir, "/Stats_by_image_", Ch, "_thresh_", thresh1, ".csv"))
  
  return(t5)
  
})

names(stats_list_by_file) = paste0(Ch, "_thresh_", thresh)



#stats by replicate (with mean/sd by group)

stats_list_by_repl = lapply(thresh, function(thresh1){
  
  t2 = t1 %>% group_by(group, repl, file) %>% summarise(N_cells = n(), 
                                                        N_pos = length(Ch[Ch>thresh1]))
  
  t3 = t2 %>% group_by(group, repl) %>% summarise(N_cells = mean(N_cells), N_pos = mean(N_pos))
  t3$fract_pos = t3$N_pos/t3$N_cells
  
  t4 = t3 %>% group_by(group) %>% summarise(N_cells_mean = mean(N_cells), 
                                            N_cells_sd = sd(N_cells), 
                                            fract_pos_mean = mean(fract_pos),
                                            fract_pos_sd = sd(fract_pos) )
  
  t5 = cbind(t3, t4[match(t3$group, t4$group),-1])
  
  write_csv(t5, file = paste0(out_dir, "/Stats_by_repl_", Ch, "_thresh_", thresh1, ".csv"))
  
  return(t5)
  
})

names(stats_list_by_repl) = paste0(Ch, "_thresh_", thresh)



#############################################
#plot summary statistics for each threshold
#############################################

#define plot colors (by group)

gr_colors = matlab.like(6)
if (length(gr)>6) {gr_colors = matlab.like(length(gr))}



### plot by replicate (individual images)

pl_repl = lapply(thresh, function(thresh1){
  
  t1 = stats_list_by_file[[paste0(Ch, "_thresh_", thresh1)]]
  
  p1 = ggplot(data = t1)+
    geom_col(aes(x = repl, y = fract_pos_repl_mean, color = group), fill = "grey90", 
             position = position_dodge(), width = 0.5, lwd = 0.4)+
    geom_errorbar(aes(x = repl, ymin = fract_pos_repl_mean-fract_pos_repl_sd,
                      y = fract_pos_repl_mean,
                      ymax =  fract_pos_repl_mean+fract_pos_repl_sd,
                      color = group),
                  position = position_dodge(), width = 0.2, lwd = 0.4)+
    geom_point(data = t1, aes(x = repl, y = fract_pos, color = group), 
               position = position_dodge(width = 0.5), size = 1, stroke = 0.5)+
    scale_color_manual(limits = gr, values = gr_colors)+
    scale_x_discrete(limits = replics)+
    geom_hline(yintercept = 0)+
    labs(title = paste0("Fraction ", Ch, " pos (thresh = ", thresh1,"; individual images)"), 
         x = "", y = "fraction positive")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
})

names(pl_repl) = paste0(Ch, "_thresh_", thresh)


# add plot with cell numbers by replicate (individual images)

t1 = stats_list_by_file[[1]]

pl_repl[["N_cells"]] = ggplot(data = t1)+
  geom_col(aes(x = repl, y = N_cells_repl_mean, color = group), fill = "grey90", 
           position = position_dodge(), width = 0.5, lwd = 0.4)+
  geom_errorbar(aes(x = repl, ymin = N_cells_repl_mean-N_cells_repl_sd,
                    y = N_cells_repl_mean,
                    ymax =  N_cells_repl_mean+N_cells_repl_sd,
                    color = group),
                position = position_dodge(), width = 0.2, lwd = 0.4)+
  geom_point(data = t1, aes(x = repl, y = N_cells, color = group), 
             position = position_dodge(width = 0.5), size = 1, stroke = 0.5)+
  scale_color_manual(limits = gr, values = gr_colors)+
  scale_x_discrete(limits = replics)+
  geom_hline(yintercept = 0)+
  labs(title = "Number of cells/field (individual images)", 
       x = "", y = "N cells")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


#save plots by replicate

file.name = paste0(out_dir, "/Plots_cell_stats_",Ch,"_by_replicate.pdf")

pdf(file = file.name, width = length(replics)/4+2, height = 3)
lapply(pl_repl, function(x){x})
dev.off()




### plots by group (individual images/replicates)

# threshold plots individual images by group

pl_gr = lapply(thresh, function(thresh1){
  
  t1 = stats_list_by_file[[paste0(Ch, "_thresh_", thresh1)]]
  
  p1 = ggplot(data = t1)+
    geom_col(aes(x = group, y = fract_pos_group_mean, color = group), fill = "grey90", 
             position = position_dodge(), width = 0.5, lwd = 0.4)+
    geom_errorbar(aes(x = group, ymin = fract_pos_group_mean-fract_pos_group_sd,
                      y = fract_pos_group_mean,
                      ymax =  fract_pos_group_mean+fract_pos_group_sd,
                      color = group),
                  position = position_dodge(), width = 0.2, lwd = 0.4)+
    geom_point(data = t1, aes(x = group, y = fract_pos, color = group), 
               position = position_dodge(width = 0.5), size = 1, stroke = 0.5)+
    scale_color_manual(limits = gr, values = gr_colors)+
    scale_x_discrete(limits = gr)+
    geom_hline(yintercept = 0)+
    labs(title = paste0( Ch, " (thresh = ", thresh1,"; indiv images)"), 
         x = "", y = "fraction positive")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
})

names(pl_gr) = paste0(Ch, "_thresh_", thresh, "_indiv_images")


# add plot with cell numbers individual images by group

t1 = stats_list_by_file[[1]]

pl_gr[["N_cells_indiv_images"]] = ggplot(data = t1)+
  geom_col(aes(x = group, y = N_cells_group_mean, color = group), fill = "grey90", 
           position = position_dodge(), width = 0.5, lwd = 0.4)+
  geom_errorbar(aes(x = group, ymin = N_cells_group_mean-N_cells_group_sd,
                    y = N_cells_group_mean,
                    ymax =  N_cells_group_mean+N_cells_group_sd,
                    color = group),
                position = position_dodge(), width = 0.2, lwd = 0.4)+
  geom_point(data = t1, aes(x = group, y = N_cells, color = group), 
             position = position_dodge(width = 0.5), size = 1, stroke = 0.5)+
  scale_color_manual(limits = gr, values = gr_colors)+
  scale_x_discrete(limits = gr)+
  geom_hline(yintercept = 0)+
  labs(title = "Number of cells/field (individual images)", 
       x = "", y = "N cells")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))



# create threshold plots individual replicates by group

l1 = lapply(thresh, function(thresh1){
  
  t1 = stats_list_by_repl[[paste0(Ch, "_thresh_", thresh1)]]
  
  p1 = ggplot(data = t1)+
    geom_col(aes(x = group, y = fract_pos_mean, color = group), fill = "grey90", 
             position = position_dodge(), width = 0.5, lwd = 0.4)+
    geom_errorbar(aes(x = group, ymin = fract_pos_mean-fract_pos_sd,
                      y = fract_pos_mean,
                      ymax =  fract_pos_mean+fract_pos_sd,
                      color = group),
                  position = position_dodge(), width = 0.2, lwd = 0.4)+
    geom_point(data = t1, aes(x = group, y = fract_pos, color = group), 
               position = position_dodge(width = 0.5), size = 1, stroke = 0.5)+
    scale_color_manual(limits = gr, values = gr_colors)+
    scale_x_discrete(limits = gr)+
    geom_hline(yintercept = 0)+
    labs(title = paste0(Ch, " (thresh = ", thresh1,"; indiv replicates)"), 
         x = "", y = "fraction positive")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
})

names(l1) = paste0(Ch, "_thresh_", thresh, "_indiv_replicates")


# add plot with cell numbers individual replicates by group

t1 = stats_list_by_repl[[1]]

l1[["N_cells_indiv_repl"]] = ggplot(data = t1)+
  geom_col(aes(x = group, y = N_cells_mean, color = group), fill = "grey90", 
           position = position_dodge(), width = 0.5, lwd = 0.4)+
  geom_errorbar(aes(x = group, ymin = N_cells_mean-N_cells_sd,
                    y = N_cells_mean,
                    ymax =  N_cells_mean+N_cells_sd,
                    color = group),
                position = position_dodge(), width = 0.2, lwd = 0.4)+
  geom_point(data = t1, aes(x = group, y = N_cells, color = group), 
             position = position_dodge(width = 0.5), size = 1, stroke = 0.5)+
  scale_color_manual(limits = gr, values = gr_colors)+
  scale_x_discrete(limits = gr)+
  geom_hline(yintercept = 0)+
  labs(title = "Number of cells/field (individual replicates)", 
       x = "", y = "N cells")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


#merge plots of individual images and replicates and save plots by group
pl_gr = c(pl_gr, l1)

file.name = paste0(out_dir, "/Plots_cell_stats_",Ch,"_by_group.pdf")

pdf(file = file.name, width = length(gr)/4+2, height = 3)
lapply(pl_gr, function(x){x})
dev.off()





### plot fraction positive vs N_cells for all images (by group)

pl_repl = lapply(thresh, function(thresh1){
  
  t1 = stats_list_by_file[[paste0(Ch, "_thresh_", thresh1)]]
  
  p1 = ggplot(data = t1)+
    geom_point(data = t1, aes(x = N_cells, y = fract_pos, color = group), 
               position = position_dodge(width = 0.5), size = 1, stroke = 0.5)+
    scale_color_manual(limits = gr, values = gr_colors)+
    geom_hline(yintercept = 0)+
    geom_vline(xintercept = 0)+
    labs(title = paste0("N_cells vs ", Ch, " pos (thresh = ", thresh1,"; individual images)"), 
         x = "N cells per field", y = "fraction positive")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
})

names(pl_repl) = paste0(Ch, "_thresh_", thresh)

#save plots

file.name = paste0(out_dir, "/Plots_cell_stats_",Ch,"_vs_N_cells.pdf")

pdf(file = file.name, width = 4, height = 4)
lapply(pl_repl, function(x){x})
dev.off()



#############################################################################################
# statistical analysis with GLMM (lmerTest algorithm) approach from by Yu et al., 2022(Neuron) 
#############################################################################################

# create list to collect lmertest results

l1 = list()

# GLMM model for number of cells per field vs group 

t1 = stats_list_by_file[[1]]

l1[["N_cells"]] = lmerTest::lmer(N_cells ~ group+(1|repl), data=t1)


# GLMM intensity by cell vs group

t1 = set1$intens_tab
t1$intensity = t1[,Ch] 

l1[["Intensity"]] = lmerTest::lmer(intensity ~ group+(1|repl) + (1|file), data=t1)


# GLMM fract pos by field vs group for each treshhold

l2 = lapply(thresh, function(thresh1){
  
  t1 = stats_list_by_file[[paste0(Ch, "_thresh_", thresh1)]]
  t2 = lmerTest::lmer(fract_pos ~ group+(1|repl), data=t1)
  return(t2)
  
})
names(l2) = paste0(Ch, "_thresh_", thresh)


#merge GLMM results and extract summary

l3 = c(l1, l2)

l4 = lapply(names(l3), function(x){
  
  t3 = summary(l3[[x]])
  t4 = t3$coefficients
  t5 = as_tibble(cbind(comp = x, group = rownames(t4), t4))
  return(t5)
  
})
names(l4) = names(l3)

t6 = NULL
for (i in names(l4)){
  t6 = rbind(t6, l4[[i]])
}

t6$`Pr(>|t|)` = as.numeric(t6$`Pr(>|t|)`)
t6$sign = ""
t6$sign[t6$`Pr(>|t|)`<=0.05] = "*"
t6$sign[t6$`Pr(>|t|)`<=0.01] = "**"
t6$sign[t6$`Pr(>|t|)`<=0.001] = "***"

write_csv(t6, file = paste0(out_dir, "lmerTest_results_",Ch, ".csv"))




