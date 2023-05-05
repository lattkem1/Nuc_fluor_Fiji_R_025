### analysis nuclear intensities single channel with different thresholds

# go to main directory (parent directory of scripts)
if (basename(getwd())== "00_scripts"){setwd("../.")} 

#load packages

library("tidyverse")
library(colorRamps)
library(lmerTest)

#define working directories

in_dir = "./04b_intens_distrib_2Ch/"
out_dir = "./05b_intens_thresh_analysis_2Ch/"
if (!dir.exists(out_dir)) {dir.create(out_dir)}


#### load data and set parameters 

load(file = paste0(in_dir,"/dataset PAX6 KI67.RData"))

#define channel and thresholds for analysis
set1$channels

Ch = set1$channels

thresh = c(3000, 1500)




#############################################
#create summary statistics by file(image)/replicate
#############################################

#extract groups and replicates

replics = set1$replicates
gr = set1$groups

# set channels for comparison as generic Ch1/Ch2 variable

t1 = set1$intens_tab
t1$Ch1 = t1[,Ch[1]]
t1$Ch2 = t1[,Ch[2]]

int_tab = t1


### create each one stats table by image for cell number and fraction of specific subsets 
#       (Ch1 neg/pos Ch2 neg/pos (nn, np, pn, pp)) and merge

t1 = int_tab

t2 = t1 %>% group_by(group, repl, file) %>% 
  summarise(N_cells = n(),  subset = paste0(Ch[1],"/", Ch[2], " -/-"), 
            N_subset = length(Ch1[Ch1<thresh[1]& Ch2<thresh[2]]))

t3 = t1 %>% group_by(group, repl, file) %>% 
  summarise(N_cells = n(),  subset = paste0(Ch[1],"/", Ch[2], " -/+"), 
            N_subset = length(Ch1[Ch1<thresh[1]& Ch2>=thresh[2]]))

t4 = t1 %>% group_by(group, repl, file) %>% 
  summarise(N_cells = n(),  subset = paste0(Ch[1],"/", Ch[2], " +/-"), 
            N_subset = length(Ch1[Ch1>=thresh[1]& Ch2<thresh[2]]))

t5 = t1 %>% group_by(group, repl, file) %>% 
  summarise(N_cells = n(),  subset = paste0(Ch[1],"/", Ch[2], " +/+"), 
            N_subset = length(Ch1[Ch1>=thresh[1]& Ch2>=thresh[2]]))

t6 = rbind(t2,t3,t4,t5)

t6$fract_subset = t6$N_subset/t6$N_cells

stats_by_file = t6

# calculate mean and for each subset by replicate, add to stats table, then same for group

t1 = stats_by_file

t2 = t1 %>% group_by(group, repl, subset) %>% 
  summarize(fract_subset_repl_mean = mean(fract_subset),
            fract_subset_repl_sd = sd(fract_subset) )

t1$fract_subset_repl_mean = t2$fract_subset_repl_mean[match(
  paste0(t1$repl, t1$subset), paste0(t2$repl, t2$subset))]
t1$fract_subset_repl_sd = t2$fract_subset_repl_sd[match(
  paste0(t1$repl, t1$subset), paste0(t2$repl, t2$subset))]


t2 = t1 %>% group_by(group, subset) %>% 
  summarize(fract_subset_group_mean = mean(fract_subset),
            fract_subset_group_sd = sd(fract_subset) )

t1$fract_subset_group_mean = t2$fract_subset_group_mean[match(
  paste0(t1$group, t1$subset), paste0(t2$group, t2$subset))]
t1$fract_subset_group_sd = t2$fract_subset_group_sd[match(
  paste0(t1$group, t1$subset), paste0(t2$group, t2$subset))]

stats_by_file = t1

write_csv(t1, file = paste0(out_dir, "stats_by_file_",Ch[1],"_", Ch[2],
                            "_thresh_", thresh[1],"_", thresh[2], ".csv"))



### create each one stats table by replicate for cell number and fraction of specific subsets 
#       (Ch1 neg/pos Ch2 neg/pos (nn, np, pn, pp)) and merge

t1 = int_tab

t2 = t1 %>% group_by(group, repl) %>% 
  summarise(N_cells = n(),  subset = paste0(Ch[1],"/", Ch[2], " -/-"), 
            N_subset = length(Ch1[Ch1<thresh[1]& Ch2<thresh[2]]))

t3 = t1 %>% group_by(group, repl) %>% 
  summarise(N_cells = n(),  subset = paste0(Ch[1],"/", Ch[2], " -/+"), 
            N_subset = length(Ch1[Ch1<thresh[1]& Ch2>=thresh[2]]))

t4 = t1 %>% group_by(group, repl) %>% 
  summarise(N_cells = n(),  subset = paste0(Ch[1],"/", Ch[2], " +/-"), 
            N_subset = length(Ch1[Ch1>=thresh[1]& Ch2<thresh[2]]))

t5 = t1 %>% group_by(group, repl) %>% 
  summarise(N_cells = n(),  subset = paste0(Ch[1],"/", Ch[2], " +/+"), 
            N_subset = length(Ch1[Ch1>=thresh[1]& Ch2>=thresh[2]]))

t6 = rbind(t2,t3,t4,t5)

t6$fract_subset = t6$N_subset/t6$N_cells

stats_by_repl = t6

# calculate mean and for each subset by replicate, add to stats table, then same for group

t1 = stats_by_repl

t2 = t1 %>% group_by(group, subset) %>% 
  summarize(fract_subset_group_mean = mean(fract_subset),
            fract_subset_group_sd = sd(fract_subset) )

t1$fract_subset_group_mean = t2$fract_subset_group_mean[match(
  paste0(t1$group, t1$subset), paste0(t2$group, t2$subset))]
t1$fract_subset_group_sd = t2$fract_subset_group_sd[match(
  paste0(t1$group, t1$subset), paste0(t2$group, t2$subset))]

stats_by_repl = t1

write_csv(t1, file = paste0(out_dir, "stats_by_repl_",Ch[1],"_", Ch[2],
                            "_thresh_", thresh[1],"_", thresh[2], ".csv"))




#############################################
#plot summary statistics 
#############################################

#define plot colors (by group)

gr_colors = matlab.like(6)
if (length(gr)>6) {gr_colors = matlab.like(length(gr))}



### plot by replicate (individual images)

t1 = stats_by_file

p1 = ggplot(data = t1)+
  geom_col(aes(x = subset, y = fract_subset_repl_mean, group = repl, color = group), 
           fill = "grey90", position = position_dodge(), width = 0.8, lwd = 0.4)+
  geom_errorbar(aes(x = subset, ymin = fract_subset_repl_mean-fract_subset_repl_sd,
                    y = fract_subset_repl_mean,
                    ymax =  fract_subset_repl_mean+fract_subset_repl_sd,
                    color = group,
                    group = repl),
                position = position_dodge(width = 0.8), width = 0.2, lwd = 0.4)+
  geom_point(data = t1, aes(x = subset, y = fract_subset, color = group, group = repl), 
             position = position_dodge(width = 0.8), size = 0.5, stroke = 0.5)+
  scale_color_manual(limits = gr, values = gr_colors)+
  scale_x_discrete(limits = unique(t1$subset))+
  geom_hline(yintercept = 0)+
  labs(title = paste0("Fraction by replicate (individual images)"), 
       x = "", y = "fraction of cells")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


file.name = paste0(out_dir, "/Plot_cell_stats_",Ch[1],"_",Ch[2],
                   "_tresh_", thresh[1],"_",thresh[2],"_by_repl.pdf")

pdf(file = file.name, width = length(replics)/4+2, height = 3)
p1
dev.off()



### plot by group (individual images/replicates)

#plot individual images

t1 = stats_by_file

p1 = ggplot(data = t1)+
  geom_col(aes(x = subset, y = fract_subset_group_mean, group = group, color = group), 
           fill = "grey90", position = position_dodge(), width = 0.8, lwd = 0.4)+
  geom_errorbar(aes(x = subset, ymin = fract_subset_group_mean-fract_subset_group_sd,
                    y = fract_subset_group_mean,
                    ymax =  fract_subset_group_mean+fract_subset_group_sd,
                    color = group,
                    group = group),
                position = position_dodge(width = 0.8), width = 0.2, lwd = 0.4)+
  geom_point(data = t1, aes(x = subset, y = fract_subset, color = group, group = group), 
             position = position_dodge(width = 0.8), size = 0.5, stroke = 0.5)+
  scale_color_manual(limits = gr, values = gr_colors)+
  scale_x_discrete(limits = unique(t1$subset))+
  geom_hline(yintercept = 0)+
  labs(title = paste0("Fraction of cells (individual images)"), 
       x = "", y = "fraction of cells")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))



#plot individual replicates

t1 = stats_by_repl

p2 = ggplot(data = t1)+
  geom_col(aes(x = subset, y = fract_subset_group_mean, group = group, color = group), 
           fill = "grey90", position = position_dodge(), width = 0.8, lwd = 0.4)+
  geom_errorbar(aes(x = subset, ymin = fract_subset_group_mean-fract_subset_group_sd,
                    y = fract_subset_group_mean,
                    ymax =  fract_subset_group_mean+fract_subset_group_sd,
                    color = group,
                    group = group),
                position = position_dodge(width = 0.8), width = 0.2, lwd = 0.4)+
  geom_point(data = t1, aes(x = subset, y = fract_subset, color = group, group = group), 
             position = position_dodge(width = 0.8), size = 0.5, stroke = 0.5)+
  scale_color_manual(limits = gr, values = gr_colors)+
  scale_x_discrete(limits = unique(t1$subset))+
  geom_hline(yintercept = 0)+
  labs(title = paste0("Fraction of cells (individual replicates)"), 
       x = "", y = "fraction of cells")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


#save plots

file.name = paste0(out_dir, "/Plot_cell_stats_",Ch[1],"_",Ch[2],
                   "_tresh_", thresh[1],"_",thresh[2],"_by_group.pdf")

pdf(file = file.name, width = length(gr)/4+3, height = 3)
p1
p2
dev.off()





#############################################################################################
# statistical analysis multiple t-tests by subset (for each group vs first group)
#############################################################################################

t1 = stats_by_repl

subsets = unique(t1$subset)

t6 = NULL

for (subs in subsets){
  
  t2 = t1[t1$subset == subs,]
  t3 = pairwise.t.test(t2$fract_subset, t2$group, p.adjust.method = "none")
  t4 = t3$p.value
  t5 = as_tibble(cbind(subset = subs, comp1 = colnames(t4)[1], comp2 = rownames(t4), p = t4[,1]))
  t6 = rbind(t6, t5)
  
}

t6$padj = p.adjust(t6$p, method = "BH")

write_csv(t6, file = paste0(out_dir, "/stats_ttest_",Ch[1],"_",Ch[2],
                            "_tresh_", thresh[1],"_",thresh[2],".csv"))



#############################################################################################
# statistical analysis with GLMM (lmerTest algorithm) approach from by Yu et al., 2022(Neuron)
# (works only for multiple fields/replicate)
#############################################################################################


# GLMM model for number of cells per field vs group (if more files than replicates, else not meaningful)

t1 = stats_by_file

subsets = unique(t1$subset)


if (length(unique(t1$repl))<length(unique(t1$file))){
  
  t6 = NULL
  
  # GLMM fract pos by field vs group for each subset
  
  t6 = NULL
  
  for (subs in subsets){
    
    t2 = t1[t1$subset == subs,]
    t3 = lmerTest::lmer(fract_subset ~ group+(1|repl), data=t2)
    t3 = summary(t3)
    t4 = t3$coefficients
    t5 = as_tibble(cbind(subset = subs, group = rownames(t4), t4))
    t6 = rbind(t6, t5)
    
  }
  
  write_csv(t6, file = paste0(out_dir, "lmerTest_results_by_subset_",Ch[1],"_",Ch[2],
                              "_tresh_", thresh[1],"_",thresh[2],".csv"))
  
}





