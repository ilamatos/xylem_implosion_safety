# Load packages # -----------------
library(ggpubr)
library(viridis)
library(BIEN)
library(maps)
library(mapdata)
library(raster)
library(sp)
library(smatr) 
library(vegan)
library(tidyverse)
library(readxl)
library(ggtree)
library(smatr)

R.Version() 
citation()

# Import data # -------------------

# "spp_list.csv" = list of 122 plant species
tx<-read_csv("data/spp_list.csv")
glimpse(tx)

# "Implosion_safety_v3.xlsx" = anatomical and functional data
dt<- read_xlsx("data/implosion_safety_v3.xlsx", sheet = "data")
glimpse(dt)

unique(dt$spp_code) # 122 species

# Order factor/categorical variables
unique(dt$clade)
dt$clade <- factor(dt$clade, levels = c("ferns", "basal angiosperms", "monocots", "basal eudicots", "rosids", "asterids"))

dt$species_names <- as.factor(dt$species_names)

unique(dt$growth_form)
dt$growth_form <- factor(dt$growth_form, levels = c("aquatic", "climbing", "herb", "shrub", "tree"))

unique(dt$vein_order)
dt$vein_order <- factor(dt$vein_order, levels = c("minor", "medium", "major", NA))

glimpse(dt)

# converting dataset into data.frame format to be able to run multcomp analysis on SMATR
dt <- as.data.frame(dt)

# converting chr in dbl
dt$ppt_annual_mean <- as.double(dt$ppt_annual_mean)
dt$min_Diameter <- as.double(dt$min_Diameter)
dt$e <- as.double(dt$e)

glimpse(dt)# final dataset ready for analysis

# Building phylogenetic tree # ----
# using V.PhyloMaker 2 Jin and Qian 2022 - https://www.sciencedirect.com/science/article/pii/S2468265922000580

# installing and loading PhyloMaker 2 
#library(devtools)
#devtools::install_github("jinyizju/V.PhyloMaker2")
library(V.PhyloMaker2) 

# creating the phylogenetic tree
tree122 <- phylo.maker(sp.list = tx%>%select (Species_WP, Genus_WP, Family_WP), tree = GBOTB.extended.WP, nodes = nodes.info.1.WP, scenarios = "S3")
tree122
str(tree122)

# checking if all 122 species were included
tree122$scenario.3$tip.label

# saving the phylogenetic tree
write.tree(tree122$scenario.3, "data/UCBG_122spp.tre") 

# convert the tree from tre to txt file (just need to change the file extension)

# Figure S1 # ---------------------

# installing and loading ggtree
#install.packages("BiocManager", repos = "https://cloud.r-project.org")
#library(BiocManager)
#BiocManager::install("ggtree")
library(ggtree)

# reading tree txt
tree <- read.tree("data/UCBG_122spp.txt")

glimpse(tree)
# plotting tree with ggtree
ggtree(tree, layout="circular")+
geom_tiplab(size=3)


# Figure S2 # ---------------------
spp_occ<-read_csv("data/spp_occurrences_v3.csv")

# loading packages
library(hrbrthemes)
library(mapdata)
library(hexbin)

glimpse(spp_occ)

a<-ggplot(spp_occ, aes(x=longitude, y=latitude)) +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="grey60", alpha=0.3)+
  geom_hex(bins=100)+
  theme_void()+
  scale_fill_viridis(
    option="D",
    trans = "log", 
    breaks = c(1,7,54,403,3000),
    name=" Species occurrences", 
    #guide = guide_legend( keyheight = unit(2.5, units = "mm"), keywidth=unit(10, units = "mm"), label.position = "bottom", title.position = 'top', nrow=1) 
  )+ 
  ggtitle( "" ) +
  theme(
    legend.position = "bottom",
    legend.title=element_text(color="black", size=12),
    text = element_text(color = "#22211d"),
    plot.title = element_text(size= 16, hjust=0.1, color = "#4e4d47", margin = margin(b = -0.1, t = 0.4, l = 2, unit = "cm")),
  );a
a
ggsave("figures/Figure_S2.png", width = 18, height = 12, units = "cm")

# Calculating ovality and Pcri # --------------

dt<- dt%>%
  mutate(Pcri1_low = (2*100)/(1-0.28^2) *((Thickness/2)/max_Diameter)^3)%>% # calculate lower range Pcri1, with E = 100
  mutate(Pcri1_up = (2*300)/(1-0.28^2) *((Thickness/2)/max_Diameter)^3)%>% # calculate upper range Pcri1, with E = 300
  mutate(Pcri2_low = (10/0.25)*(Thickness/max_Diameter)^2)%>% # calculate lower range Pcri2, with W = 10
  mutate(Pcri2_up = (80/0.25)*(Thickness/max_Diameter)^2)%>% # calculate upper range Pcri2, with W = 80
  mutate (O = (max_Diameter - min_Diameter )/ (max_Diameter + min_Diameter)) # calculate conduit ovality
glimpse(dt)

# Figure S4 # ---------------------
# Relationship between Pcr1 and Pcri 2
S4a<-ggplot(dt, aes(x = Pcri1_low, y = Pcri2_low))+
  geom_point(size = 4, shape =1, alpha = 0.5)+
  geom_line(size = 1, col = "red" )+
  theme_classic()+
  labs(x = "Pcri1_low (-MPa)", y = "Pcri2_low (-MPa)"); S4a

S4b<-ggplot(dt, aes(x = Pcri1_up, y = Pcri2_up))+
  geom_point(size = 4, shape =1, alpha = 0.5)+
  geom_line(size = 1, col = "red" )+
  theme_classic()+
  labs(x = "Pcri1_up (-MPa)", y = "Pcri2_up (-MPa)"); S4b

# countor plot of Pcri values 
# Pcri1 
S4c<-expand_grid(Diameter = seq(min(dt$max_Diameter),max(dt$max_Diameter),length.out = 100), 
                 Thickness = seq(min(dt$Thickness),max(dt$Thickness), length.out= 100)) %>% 
  mutate(Pcri1_low = (2*100)/(1-0.28^2) *((Thickness/2)/Diameter)^3) %>% 
  ggplot(data=., aes(Diameter, Thickness,fill=Pcri1_low))+
  geom_tile()+
  geom_point(data =dt,
             size = 3,
             alpha = 0.5,
             shape =21,
             color = 'grey20',
             inherit.aes = F,
             aes(max_Diameter,Thickness,color=Pcri1_low))+
  scale_fill_viridis_c(
    trans = "log10",
    name=" Pcri1_low (-MPa)",
    limits = c(0.0001, 10000),
    oob = scales::squish)+
  coord_cartesian(expand = F)+
  labs( x = "Diameter (um)", y = "Thickness (um)"); S4c


S4d<-expand_grid(Diameter = seq(min(dt$max_Diameter),max(dt$max_Diameter),length.out = 100), 
                 Thickness = seq(min(dt$Thickness),max(dt$Thickness), length.out= 100)) %>% 
  mutate(Pcri1_up = (2*300)/(1-0.28^2) *((Thickness/2)/Diameter)^3) %>% 
  ggplot(data=., aes(Diameter, Thickness,fill=Pcri1_up))+
  geom_tile()+
  geom_point(data =dt,
             size = 3,
             alpha = 0.5,
             shape =21,
             color = 'grey20',
             inherit.aes = F,
             aes(max_Diameter,Thickness,color=Pcri1_up))+
  scale_fill_viridis_c(
    trans = "log10",
    name=" Pcri1_up (-MPa)",
    limits = c(0.0001, 10000),
    oob = scales::squish)+
  coord_cartesian(expand = F)+
  labs( x = "Diameter (um)", y = "Thickness (um)"); S4d

# Pcri2

S4e<-expand_grid(Diameter = seq(min(dt$max_Diameter),max(dt$max_Diameter),length.out = 100), 
                 Thickness = seq(min(dt$Thickness),max(dt$Thickness), length.out= 100)) %>% 
  mutate(Pcri2_low = (10/0.25)*(Thickness/Diameter)^2) %>% 
  ggplot(data=., aes(Diameter, Thickness,fill=Pcri2_low))+
  geom_tile()+
  geom_point(data =dt,
             size = 3,
             alpha = 0.5,
             shape =21,
             color = 'grey20',
             inherit.aes = F,
             aes(max_Diameter,Thickness,color=Pcri2_low))+
  scale_fill_viridis_c(
    trans = "log10",
    name=" Pcri2_low (-MPa)",
    limits = c(0.0001, 10000),
    oob = scales::squish)+
  coord_cartesian(expand = F)+
  labs( x = "Diameter (um)", y = "Thickness (um)"); S4e

S4f<-expand_grid(Diameter = seq(min(dt$max_Diameter),max(dt$max_Diameter),length.out = 100), 
                 Thickness = seq(min(dt$Thickness),max(dt$Thickness), length.out= 100)) %>% 
  mutate(Pcri2_up = (80/0.25)*(Thickness/Diameter)^2) %>% 
  ggplot(data=., aes(Diameter, Thickness,fill=Pcri2_up))+
  geom_tile()+
  geom_point(data =dt,
             size = 3,
             alpha = 0.5,
             shape =21,
             color = 'grey20',
             inherit.aes = F,
             aes(max_Diameter,Thickness,color=Pcri2_up))+
  scale_fill_viridis_c(
    trans = "log10",
    name=" Pcri2_up (-MPa)",
    limits = c(0.0001, 10000),
    oob = scales::squish)+
  coord_cartesian(expand = F)+
  labs( x = "Diameter (um)", y = "Thickness (um)"); S4f

FigS4<- ggarrange(S4a, S4b, S4c, S4d, S4e, S4f,  
                  ncol =2, nrow=3, common.legend = T, legend = "right",
                  labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"))
FigS4
ggsave(FigS4,filename="figures/Figure_S4.png", width = 18, height = 20, limitsize = FALSE, units = "cm")

# Descriptive analysis # ----------
glimpse(dt)

# number of vessels and species per clade
dt%>% group_by(clade)%>%
  summarize(n_spp = n_distinct(spp_code), n_vessels = n()) 

# Ovality range
hist(dt$O)
min(dt$O, na.rm = T)
max(dt$O, na.rm = T)
median(dt$O, na.rm = T)

# variation in conduits diameters
hist(dt$max_Diameter)
min(dt$max_Diameter)
max(dt$max_Diameter)
median(dt$max_Diameter)

# variation in conduits double wall thickness
hist(dt$Thickness)
min(dt$Thickness)
max(dt$Thickness)
median(dt$Thickness)

# variation in conduits implosion safety
hist(dt$TD)
min(dt$TD)
max(dt$TD)
median(dt$TD)

#Pcri1 
min(dt$Pcri1_low)
max(dt$Pcri1_low)
median(dt$Pcri1_low)
min(dt$Pcri1_up)
max(dt$Pcri1_up)
median(dt$Pcri1_up)

#Pcri2 
min(dt$Pcri2_low)
max(dt$Pcri2_low)
median(dt$Pcri2_low)
min(dt$Pcri2_up)
max(dt$Pcri2_up)
median(dt$Pcri2_up)

# variation in Kleaf max
hist(dt$Kmax)
min(dt$Kmax)
max(dt$Kmax)
median(dt$TD)

# variation in e
hist(dt$e)
min(dt$e, na.rm = TRUE)
max(dt$e, na.rm = TRUE)
median(dt$TD)

# variation in LMA
hist(dt$LMA)
min(dt$LMA)
max(dt$LMA)
median(dt$LMA)

# sma1 - all species # ------------

# Fit SMA model to test if overall log10 T x log10 D relationship is isometric (slope =1)
m1<-sma( log_Thickness ~ log_Diameter, data = dt, slope.test = 1, robust = T) # robust = T to reduce effect of outliers

# model results
summary(m1)

# residual analysis
plot(m1, which = "residual")  
abline(a = 0, b = 0)
# QQ plot - testing normal distribution of residuals
plot(m1, which = "qq")

# Figure 2 # ----------------------
# create new data set of log_Diameter at a high resolution (200 points from min to max)
set.seed(1) # for reproducibility
preds <- data.frame(expand.grid(log_Diameter = seq(min(dt$log_Diameter, na.rm = T), 
                                                   max(dt$log_Diameter, na.rm = T), 
                                                   length.out = 100), stringsAsFactors = FALSE))

# mean coefficients
preds_mean<-mutate(preds, log_Thickness = coef(m1)[1] + coef(m1)[2]*log_Diameter)%>%
  mutate(Diameter = 10**log_Diameter, Thickness = 10**log_Thickness) # back-transform values

# upper confidence interval
preds_up<-mutate(preds, log_Thickness = m1$coef[[1]][3]$`upper limit`[1] + m1$coef[[1]][3]$`upper limit`[2]*log_Diameter)%>%
  mutate(Diameter = 10**log_Diameter, Thickness = 10**log_Thickness) # back-transform values

# lower confidence interval
preds_low<-mutate(preds, log_Thickness = m1$coef[[1]][2]$`lower limit`[1] +m1$coef[[1]][2]$`lower limit`[2]*log_Diameter)%>%
  mutate(Diameter = 10**log_Diameter, Thickness = 10**log_Thickness)

# plot with ggplot
dt %>% 
  mutate(Diameter = 10**log_Diameter,
         Thickness = 10**log_Thickness) %>% 
  ggplot(data=., aes(Diameter, Thickness)) +
  geom_point(shape =1, alpha =0.5, size = 2, col = "grey40")+ 
  geom_abline (slope =1,intercept = 0, linetype = "dashed", linewidth = 0.8, color = "grey20")+
  geom_abline (slope =1.5, intercept = 0, linetype = "dashed", linewidth = 0.8, color = "grey20")+
  geom_line(data = preds_mean,
            #inherit.aes = F, 
            aes(Diameter,Thickness),
            col = "black", size =1.2)+
  geom_line(data = preds_up, col = "black", size =0.5) +
  geom_line(data = preds_low, col = "black", size =0.5) +
  scale_x_log10(breaks = c(2, 5, 10, 20, 50))+
  scale_y_log10(breaks = c(0.5, 1, 2, 5, 10))+
  #xlim(2,50)+
  #ylim(0.5,10)+
  #annotate("text", label = "y = -0.12 + 0.55x ", x = 2.5, y =10, size = 4, col = "darkred")+
  #annotate("text", label = "R2 = 0.19, p < 2.22e-16 ", x = 2.5, y =8.5, size = 4, col = "darkred")+
  theme_classic()+
  labs(x = "Diameter (µm)", y = "Thickness (µm)")+
  theme(axis.text = element_text(size = 14), text = element_text(size = 16),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent'), #transparent legend panel
        panel.border = element_rect(colour = "grey20", fill=NA, size=1)
  )

# save figure
ggsave(filename="figures/Figure_2_raw.png", width = 8, height = 8, bg='transparent')

# sma2 - each species # -----------
# test for differences in slopes across species 
m2<-sma(log_Thickness ~ log_Diameter*spp_code, data = dt, robust = T, slope.test = 1)

# summary results
summary(m2)
m2
coef(m2)

spp_test_results <-m2$groupsummary
glimpse(spp_test_results)
median(spp_test_results$Slope)
max(spp_test_results$Slope)
min(spp_test_results$Slope)
# save results
write_csv(spp_test_results, "data/Table_S2.csv")

# plot residuals
plot(m2, which = "residual")  
abline(a = 0, b = 0)
# QQ plot
plot(m2, which = "qq") 

# test for differences in intercepts across species
m2y<-sma(log_Thickness ~ log_Diameter+spp_code, data = dt, robust = T)

# Figure S5 # ---------------------

# plot results  for each species with ggplot
vec_spp <- sort(unique(dt$spp_code)) # create a vector with species names in alphabetic order
vec_spp

p<-list() # create an empty list to store plots

# create a loop to make and save a plot for each species
for (i in 1:length(vec_spp)){
  # filter the data
  dt_spp<-dt%>%filter(spp_code == vec_spp[i])
  # fit the sma model for each species
  m<-sma(log_Thickness ~ log_Diameter, data = dt_spp, robust = T, slope.test = 1)
  
  set.seed(1) # for reproducibility
  preds <- data.frame(expand.grid(log_Diameter = seq(min(dt_spp$log_Diameter, na.rm = T), 
                                                     max(dt_spp$log_Diameter, na.rm = T), 
                                                     length.out = 100), stringsAsFactors = FALSE))
  
  # mean coefficients
  preds_mean<-mutate(preds, log_Thickness = coef(m)[1] + coef(m)[2]*log_Diameter)%>%
    mutate(Diameter = 10**log_Diameter, Thickness = 10**log_Thickness) # back-transform values
  
  # plot with ggplot
  dt_spp %>% 
    mutate(Diameter = 10**log_Diameter,
           Thickness = 10**log_Thickness)
  g<-ggplot(data=dt_spp, aes(max_Diameter, Thickness)) +
    geom_point(shape =1, alpha =0.8, size = 6, col = "grey40", stroke =1)+ 
    geom_abline (slope =1, linetype = "dashed", linewidth = 1.5, col = "grey20")+
    geom_line(data = preds_mean,
              aes(Diameter,Thickness),
              col = "red", size =2)+
    scale_x_log10(breaks = c(2, 5, 10, 20, 50))+
    scale_y_log10(breaks = c(0.5, 1, 2, 5, 10))+
    theme_classic()+
    labs(x = "Diameter (µm)", y = "Thickness (µm)", title = unique(dt_spp$spp_code))+
    theme(axis.text = element_text(size = 24), text = element_text(size = 32))
  #save plot for each species
 # ggsave(filename=paste0(vec_spp[i],".png"), width = 6, height = 6)
  # add plot to the p list
  p[[i]]<-g
  print(vec_spp[i])
}
p

# make a figure with all species plots 
gg_spp<-ggarrange(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],
                  p[[11]], p[[12]], p[[13]], p[[14]], p[[15]],p[[16]],p[[17]],p[[18]],p[[19]],p[[20]],
                  p[[21]], p[[22]], p[[23]], p[[24]], p[[25]],p[[26]],p[[27]],p[[28]],p[[29]],p[[30]],
                  p[[31]], p[[32]], p[[33]], p[[34]], p[[35]],p[[36]],p[[37]],p[[38]],p[[39]],p[[40]],
                  p[[41]], p[[42]], p[[43]], p[[44]], p[[45]],p[[46]],p[[47]],p[[48]],p[[49]],p[[50]],
                  p[[51]], p[[52]], p[[53]], p[[54]], p[[55]],p[[56]],p[[57]],p[[58]],p[[59]],p[[60]],
                  p[[61]], p[[62]], p[[63]], p[[64]], p[[65]],p[[66]],p[[67]],p[[68]],p[[69]],p[[70]],
                  p[[71]], p[[72]], p[[73]], p[[74]], p[[75]],p[[76]],p[[77]],p[[78]],p[[79]],p[[80]],
                  p[[81]], p[[82]], p[[83]], p[[84]], p[[85]],p[[86]],p[[87]],p[[88]],p[[89]],p[[90]],
                  p[[91]], p[[92]], p[[93]], p[[94]], p[[95]],p[[96]],p[[97]],p[[98]],p[[99]],p[[100]],
                  p[[101]], p[[102]], p[[103]], p[[104]], p[[105]],p[[106]],p[[107]],p[[108]],p[[109]],p[[110]],
                  p[[111]], p[[112]], p[[113]], p[[114]], p[[115]],p[[116]],p[[117]],p[[118]],p[[119]],p[[120]],
                  p[[121]],p[[122]],
                  ncol =11, nrow=12, common.legend = T)
gg_spp
# save figure
ggsave(gg_spp,filename="figures/Figure_S5_raw.png", width = 60, height = 45, limitsize = FALSE)

# Figure 3 #----------------------
# Dap_him plot 
dt_spp<-dt%>%filter(spp_code == "Dap_him")
glimpse(dt_spp)
# fit the sma model 
m<-sma(log_Thickness ~ log_Diameter, data = dt_spp, robust = T, slope.test = 1)

summary(m)
set.seed(1) # for reproducibility
preds <- data.frame(expand.grid(log_Diameter = seq(min(dt_spp$log_Diameter, na.rm = T), 
                                                   max(dt_spp$log_Diameter, na.rm = T), 
                                                   length.out = 100), stringsAsFactors = FALSE))

# mean coefficients
preds_mean<-mutate(preds, log_Thickness = coef(m)[1] + coef(m)[2]*log_Diameter)%>%
  mutate(Diameter = 10**log_Diameter, Thickness = 10**log_Thickness) # back-transform values

# plot with ggplot
g1<-ggplot(data=dt_spp, aes(log_Diameter, log_Thickness))+
  geom_point(shape =1, alpha =0.8, size = 4, col = "grey40", stroke =1)+ 
  geom_abline (slope =1, intercept = 0.05, linetype = "dashed", linewidth = 0.8, col = "grey20")+
  geom_line(data = preds_mean,
            aes(log_Diameter,log_Thickness),
            col = "black", size =1.5)+
  xlim(0,1.5)+
  ylim(0,1)+
  theme_classic()+
  labs(x = "log10 Diameter (µm)", y = "Log10 Thickness (µm)", title = unique(dt_spp$Species_WP))+
  theme(axis.text = element_text(size = 12), text = element_text(size = 14),plot.title = element_text(face = "italic"),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent'), #transparent legend panel
        panel.border = element_rect(colour = "grey20", fill=NA, size=1)); g1


# Ari_bae plot 
dt_spp<-dt%>%filter(spp_code == "Ari_bae")
glimpse(dt_spp)
# fit the sma model 
m<-sma(log_Thickness ~ log_Diameter, data = dt_spp, robust = T, slope.test = 1)

summary(m)
set.seed(1) # for reproducibility
preds <- data.frame(expand.grid(log_Diameter = seq(min(dt_spp$log_Diameter, na.rm = T), 
                                                   max(dt_spp$log_Diameter, na.rm = T), 
                                                   length.out = 100), stringsAsFactors = FALSE))

# mean coefficients
preds_mean<-mutate(preds, log_Thickness = coef(m)[1] + coef(m)[2]*log_Diameter)%>%
  mutate(Diameter = 10**log_Diameter, Thickness = 10**log_Thickness) # back-transform values

# plot with ggplot
g2<-ggplot(data=dt_spp, aes(log_Diameter, log_Thickness))+
  geom_point(shape =1, alpha =0.8, size = 4, col = "grey40", stroke =1)+ 
  geom_abline (slope =1, intercept = -0.34, linetype = "dashed", linewidth = 0.8, col = "grey20")+
  geom_line(data = preds_mean,
            #inherit.aes = F, 
            aes(log_Diameter,log_Thickness),
            col = "black", size =1.5)+
  xlim(0,1.5)+
  ylim(0,1)+
  theme_classic()+
  labs(x = "log10 Diameter (µm)", y = "Log10 Thickness (µm)", title = unique(dt_spp$Species_WP))+
  theme(axis.text = element_text(size = 12), text = element_text(size = 14),plot.title = element_text(face = "italic"),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent'), #transparent legend panel
        panel.border = element_rect(colour = "grey20", fill=NA, size=1)); g2


# Equ_tel plot 
dt_spp<-dt%>%filter(spp_code == "Equ_tel")
glimpse(dt_spp)
# fit the sma model 
m<-sma(log_Thickness ~ log_Diameter, data = dt_spp, robust = T, slope.test = 1)

summary(m)
set.seed(1) # for reproducibility
preds <- data.frame(expand.grid(log_Diameter = seq(min(dt_spp$log_Diameter, na.rm = T), 
                                                   max(dt_spp$log_Diameter, na.rm = T), 
                                                   length.out = 100), stringsAsFactors = FALSE))

# mean coefficients
preds_mean<-mutate(preds, log_Thickness = coef(m)[1] + coef(m)[2]*log_Diameter)%>%
  mutate(Diameter = 10**log_Diameter, Thickness = 10**log_Thickness) # back-transform values

# plot with ggplot
g3<-ggplot(data=dt_spp, aes(log_Diameter, log_Thickness))+
  geom_point(shape =1, alpha =0.8, size = 4, col = "grey40", stroke =1)+ 
  geom_abline (slope =1, intercept = -0.74, linetype = "dashed", linewidth = 0.8, col = "grey20")+
  geom_line(data = preds_mean,
            #inherit.aes = F, 
            aes(log_Diameter,log_Thickness),
            col = "black", size =1.5)+
  xlim(0,1.5)+
  ylim(0,1)+
  theme_classic()+
  labs(x = "log10 Diameter (µm)", y = "Log10 Thickness (µm)", title = unique(dt_spp$Species_WP))+
  theme(axis.text = element_text(size = 12), text = element_text(size = 14),plot.title = element_text(face = "italic"),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent'), #transparent legend panel
        panel.border = element_rect(colour = "grey20", fill=NA, size=1)); g3


# calculations 
# use the line equation of each species y-intercept + slope * log10(Diameter) to find log10(Thickness)
Thicklog = -0.3361035 + 0.8684820* 1.4
# back-transform thickness 
Thick = 10^Thicklog
Thick
# back-transform diameter
Diam = 10^ 1.4
Diam
# calculate implosion safety
safety = Thick/Diam
safety
# calculate critical pressures (low)
Pcri1 = (2*100)/(1-0.28^2)*((Thick/2)/Diam)^3
Pcri1
Pcri2 = (10/0.25)* (Thick/Diam)^2
Pcri2

# make figure 3 
fig3<-ggarrange( g1, g2, g3,ncol =3, nrow=1, common.legend = T)
fig3
# save figure
ggsave(fig3,filename="figures/Figure_3_raw.png", width = 12, height = 4, limitsize = FALSE, bg='transparent')

# sma3 - clades # -----------------

# test for differences in slopes across clades 
m3<-sma(log_Thickness ~ log_Diameter*clade, data = dt, slope.test = 1, robust = T, multcomp = T, multcompmethod = "adjusted")
summary(m3)

clade_test_results <-m3$groupsummary
glimpse(clade_test_results)

# plot residuals
plot(m3, which = "residual")  
abline(a = 0, b = 0)
# QQ plot
plot(m3, which = "qq") 

# comparing slopes across clades 
multcompmatrix(m3)

# test for differences in intercepts across species
m3y<-sma(log_Thickness ~ log_Diameter+clade, data = dt, robust = T)

# sma4 - habitas #----------------
# test for differences across species habitats
# divide climate data into categorical bins: arid, mesic and hydric
hist(dt$ppt_annual_mean)
dt_ppt1 <- dt %>%
  mutate(ppt_binned = case_when(ppt_annual_mean <= 500 ~ "arid",
                                ppt_annual_mean > 500 & ppt_annual_mean <= 2000 ~ "mesic",
                                ppt_annual_mean > 2000 ~ "hydric"))

glimpse(dt_ppt1)

m4<-sma(log_Thickness ~ log_Diameter*ppt_binned, data = dt_ppt1, robust = T, slope.test = 1, multcomp = T, multcompmethod = "adjusted")

# summary results
summary(m4)
m4
coef(m4)

# plot results
plot(m4)
abline(a = 0, b = 1, lwd =2, lty=2)

# plot residuals
plot(m4, which = "residual")  
abline(a = 0, b = 0)
# QQ plot
plot(m4, which = "qq") 

# comparing across species habitats 
multcompmatrix(m4)

# test for differences in intercepts across species
m4y<-sma(log_Thickness ~ log_Diameter+ppt_binned, data = dt_ppt1, robust = T)

# sma5 - growth forms # -----------
# test for differences across growth forms
glimpse(dt)
unique(dt$growth_form) 
m5<-sma(log_Thickness ~ log_Diameter*growth_form, data = dt, robust = T, slope.test = 1, multcomp = T, multcompmethod = "adjusted")

# summary results
summary(m5)
m5
coef(m5)

# plot results
plot(m5)
abline(a = 0, b = 1, lwd =2, lty=2)

# plot residuals
plot(m5, which = "residual")  
abline(a = 0, b = 0)
# QQ plot
plot(m5, which = "qq") 

# comparing across growth forms
multcompmatrix(m5)

# test for differences in intercepts across species
m5y<-sma(log_Thickness ~ log_Diameter+growth_form, data = dt, robust = T)

# sma6 - vein orders # ------------
# test for differences across vein orders 
# filtering out samples with NA values in vein orders
unique(dt$vein_order)
dt_order <- dt %>% drop_na (vein_order) # filtering out NAs in vein_order
class(dt_order)  

m6<-sma(log_Thickness ~ log_Diameter*vein_order, data = dt_order, robust = T, slope.test = 1, multcomp = T, multcompmethod = "adjusted")
# summary results
summary(m6)
m6

# plot residuals
plot(m6, which = "residual")  
abline(a = 0, b = 0)
# QQ plot
plot(m6, which = "qq") 

# comparing across vein orders
multcompmatrix(m6)

# test for differences in intercepts across species
m6y<-sma(log_Thickness ~ log_Diameter+vein_order, data = dt, robust = T)

# Figure 4 # ----------------------
# plot T x D scaling for each clade 
vec_clade <- unique(dt$clade) # create a vector with clades names
vec_clade
p<-list() # create an empty list to store plots

# create a loop to make and save a plot for each clade
for (i in 1:length(vec_clade)){
  # filter the data
  dt_clade<-dt%>%filter(clade == vec_clade[i])
  # fit the sma model for each clade
  m<-sma(log_Thickness ~ log_Diameter, data = dt_clade, robust = T, slope.test = 1)
  
  set.seed(1) # for reproducibility
  preds <- data.frame(expand.grid(log_Diameter = seq(min(dt_clade$log_Diameter, na.rm = T), 
                                                     max(dt_clade$log_Diameter, na.rm = T), 
                                                     length.out = 100), stringsAsFactors = FALSE))
  
  # mean coefficients
  preds_mean<-mutate(preds, log_Thickness = coef(m)[1] + coef(m)[2]*log_Diameter)%>%
    mutate(Diameter = 10**log_Diameter, Thickness = 10**log_Thickness) # back-transform values
  
  # upper confidence interval
  preds_up<-mutate(preds, log_Thickness = m$coef[[1]][3]$`upper limit`[1] + m$coef[[1]][3]$`upper limit`[2]*log_Diameter)%>%
    mutate(Diameter = 10**log_Diameter, Thickness = 10**log_Thickness) # back-transform values
  
  # lower confidence interval
  preds_low<-mutate(preds, log_Thickness = m$coef[[1]][2]$`lower limit`[1] +m$coef[[1]][2]$`lower limit`[2]*log_Diameter)%>%
    mutate(Diameter = 10**log_Diameter, Thickness = 10**log_Thickness)
  
  # plot with ggplot
  dt_clade %>% 
    mutate(Diameter = 10**log_Diameter,
           Thickness = 10**log_Thickness)
  
  g<-ggplot(data=dt_clade, aes(max_Diameter, Thickness)) +
    geom_point(shape =1, alpha =0.5, size = 2, col = "grey50")+ 
    geom_abline (slope =1, linetype = "dashed", linewidth = 1)+
    geom_line(data = preds_mean,
              aes(Diameter,Thickness),
              col = "black", size =1)+
    geom_line(data = preds_up, aes(Diameter,Thickness), col = "black", size =0.4) +
    geom_line(data = preds_low, aes(Diameter,Thickness), col = "black", size =0.4) +
    scale_x_log10(breaks = c(2, 5, 10, 20, 50))+ 
    scale_y_log10(breaks = c(0.5, 1, 2, 5, 10))+
    theme_classic()+
    labs(x = "Diameter (µm)", y = "Thickness (µm)", title = vec_clade[i])+
    theme(axis.text = element_text(size = 10), text = element_text(size = 12))
  #save plot
 # ggsave(filename=paste0(vec_clade[i],".png"), width = 6, height = 6)
  # add plot to the p list
  p[[i]]<-g
  print(vec_clade[i])
}

fig4a<- p[[2]]
fig4b<- p[[4]]
fig4c<- p[[1]]
fig4d <-p[[5]]
fig4e <- p[[3]]
fig4f<- p[[6]]


# plot T x D scaling for each species habitat 
vec_ppt1 <- sort(unique(dt_ppt1$ppt_binned)) # create a vector with species names in alphabetic order
vec_ppt1
p<-list() # create an empty list to store plots

# create a loop to make and save a plot for each vein order
for (i in 1:length(vec_ppt1)){
  # filter the data
  dt2<-dt_ppt1%>%filter(ppt_binned == vec_ppt1[i])
  
  # fit the sma model for each species
  m<-sma(log_Thickness ~ log_Diameter, data = dt2, robust = T, slope.test = 1)
  
  set.seed(1) # for reproducibility
  preds <- data.frame(expand.grid(log_Diameter = seq(min(dt2$log_Diameter, na.rm = T), 
                                                     max(dt2$log_Diameter, na.rm = T), 
                                                     length.out = 100), stringsAsFactors = FALSE))
  
  # mean coefficients
  preds_mean<-mutate(preds, log_Thickness = coef(m)[1] + coef(m)[2]*log_Diameter)%>%
    mutate(Diameter = 10**log_Diameter, Thickness = 10**log_Thickness) # back-transform values
  
  # upper confidence interval
  preds_up<-mutate(preds, log_Thickness = m$coef[[1]][3]$`upper limit`[1] + m$coef[[1]][3]$`upper limit`[2]*log_Diameter)%>%
    mutate(Diameter = 10**log_Diameter, Thickness = 10**log_Thickness) # back-transform values
  
  # lower confidence interval
  preds_low<-mutate(preds, log_Thickness = m$coef[[1]][2]$`lower limit`[1] +m$coef[[1]][2]$`lower limit`[2]*log_Diameter)%>%
    mutate(Diameter = 10**log_Diameter, Thickness = 10**log_Thickness)
  
  # plot with ggplot
  dt2 %>% 
    mutate(Diameter = 10**log_Diameter,
           Thickness = 10**log_Thickness)
  g<-ggplot(data=dt2, aes(max_Diameter, Thickness)) +
    geom_point(shape =1, alpha =0.5, size = 2, col = "grey50")+ 
    geom_abline (slope =1, linetype = "dashed", linewidth = 1)+
    geom_line(data = preds_mean,
              #inherit.aes = F, 
              aes(Diameter,Thickness),
              col = "black", size =1)+
    geom_line(data = preds_up, aes(Diameter,Thickness), col = "black", size =0.4) +
    geom_line(data = preds_low, aes(Diameter,Thickness), col = "black", size =0.4) +
    scale_x_log10(breaks = c(2, 5, 10, 20, 50))+ 
    scale_y_log10(breaks = c(0.5, 1, 2, 5, 10))+
    theme_classic()+
    labs(x = "Diameter (µm)", y = "Thickness (µm)", title =vec_ppt1[i])+
    theme(axis.text = element_text(size = 10), text = element_text(size = 12))
  g
  #save plot
 # ggsave(filename=paste0(vec_ppt1[i],".png"), width = 6, height = 6)
  # add plot to the p list
  p[[i]]<-g
  print(vec_ppt1[i])
}


fig4g <- p[[1]]
fig4h <- p[[3]]
fig4i <- p[[2]]

# plot T x D scaling for each vein order 
vec_order <- sort(unique(dt_order$vein_order)) # create a vector with species names in alphabetic order
vec_order
p<-list() # create an empty list to store plots

# create a loop to make and save a plot for each vein order
for (i in 1:length(vec_order)){
  # filter the data
  dt2<-dt_order%>%filter(vein_order == vec_order[i])
  # fit the sma model for each species
  m<-sma(log_Thickness ~ log_Diameter, data = dt2, robust = T, slope.test = 1)
  
  set.seed(1) # for reproducibility
  preds <- data.frame(expand.grid(log_Diameter = seq(min(dt2$log_Diameter, na.rm = T), 
                                                     max(dt2$log_Diameter, na.rm = T), 
                                                     length.out = 100), stringsAsFactors = FALSE))
  
  # mean coefficients
  preds_mean<-mutate(preds, log_Thickness = coef(m)[1] + coef(m)[2]*log_Diameter)%>%
    mutate(Diameter = 10**log_Diameter, Thickness = 10**log_Thickness) # back-transform values
  
  # upper confidence interval
  preds_up<-mutate(preds, log_Thickness = m$coef[[1]][3]$`upper limit`[1] + m$coef[[1]][3]$`upper limit`[2]*log_Diameter)%>%
    mutate(Diameter = 10**log_Diameter, Thickness = 10**log_Thickness) # back-transform values
  
  # lower confidence interval
  preds_low<-mutate(preds, log_Thickness = m$coef[[1]][2]$`lower limit`[1] +m$coef[[1]][2]$`lower limit`[2]*log_Diameter)%>%
    mutate(Diameter = 10**log_Diameter, Thickness = 10**log_Thickness)
  
  # plot with ggplot
  dt2 %>% 
    mutate(Diameter = 10**log_Diameter,
           Thickness = 10**log_Thickness)
  g<-ggplot(data=dt2, aes(max_Diameter, Thickness)) +
    geom_point(shape =1, alpha =0.5, size = 2, col = "grey50")+ 
    geom_abline (slope =1, linetype = "dashed", linewidth = 1)+
    geom_line(data = preds_mean,
              aes(Diameter,Thickness),
              col = "black", size =1)+
    geom_line(data = preds_up, aes(Diameter,Thickness), col = "black", size =0.4) +
    geom_line(data = preds_low, aes(Diameter,Thickness), col = "black", size =0.4) +
    scale_x_log10(breaks = c(2, 5, 10, 20, 50))+ 
    scale_y_log10(breaks = c(0.5, 1, 2, 5, 10))+
    theme_classic()+
    labs(x = "Diameter (µm)", y = "Thickness (µm)", title =vec_order[i])+
    theme(axis.text = element_text(size = 10), text = element_text(size = 12))
  #save plot
  #ggsave(filename=paste0(vec_order[i],".png"), width = 6, height = 6)
  # add plot to the p list
  p[[i]]<-g
  print(vec_order[i])
}

fig4j <- p[[3]]
fig4k <- p[[2]]
fig4l <- p[[1]]

# prepare Figure 4
fig4<-ggarrange(fig4a, fig4b, fig4c, fig4d, fig4e, fig4f,
                fig4g, fig4h, fig4i,
                fig4j, fig4k, fig4l,
                ncol =3, nrow=4, common.legend = T)
fig4

# save Figure 3
ggsave(fig4,filename="figures/Figure_4_raw.png", width = 17.5, height = 24.8, units = "cm", dpi = 300)

# Figure S6 # ---------------------
#plot T x D scaling for each growth form 
# create a vector with species names in alphabetic order
vec_form <- sort(unique(dt$growth_form)) 
vec_form
p<-list() # create an empty list to store plots

# create a loop to make and save a plot for each vein order
for (i in 1:length(vec_form)){
  # filter the data
  dt2<-dt%>%filter(growth_form == vec_form[i])
  
  # fit the sma model for each species
  m<-sma(log_Thickness ~ log_Diameter, data = dt2, robust = T, slope.test = 1)
  
  set.seed(1) # for reproducibility
  preds <- data.frame(expand.grid(log_Diameter = seq(min(dt2$log_Diameter, na.rm = T), 
                                                     max(dt2$log_Diameter, na.rm = T), 
                                                     length.out = 100), stringsAsFactors = FALSE))
  
  # mean coefficients
  preds_mean<-mutate(preds, log_Thickness = coef(m)[1] + coef(m)[2]*log_Diameter)%>%
    mutate(Diameter = 10**log_Diameter, Thickness = 10**log_Thickness) # back-transform values
  
  # upper confidence interval
  preds_up<-mutate(preds, log_Thickness = m$coef[[1]][3]$`upper limit`[1] + m$coef[[1]][3]$`upper limit`[2]*log_Diameter)%>%
    mutate(Diameter = 10**log_Diameter, Thickness = 10**log_Thickness) # back-transform values
  
  # lower confidence interval
  preds_low<-mutate(preds, log_Thickness = m$coef[[1]][2]$`lower limit`[1] +m$coef[[1]][2]$`lower limit`[2]*log_Diameter)%>%
    mutate(Diameter = 10**log_Diameter, Thickness = 10**log_Thickness)
  
  # plot with ggplot
  dt %>% 
    mutate(Diameter = 10**log_Diameter,
           Thickness = 10**log_Thickness)
  g<-ggplot(data=dt, aes(max_Diameter, Thickness)) +
    geom_point(shape =1, alpha =0.5, size = 2, col = "grey50")+ 
    geom_abline (slope =1, linetype = "dashed", linewidth = 1)+
    geom_line(data = preds_mean,
              aes(Diameter,Thickness),
              col = "black", size =1)+
    geom_line(data = preds_up,aes(Diameter,Thickness), col = "black", size =0.4) +
    geom_line(data = preds_low, aes(Diameter,Thickness), col = "black", size =0.4) +
    scale_x_log10(breaks = c(2, 5, 10, 20, 50))+ 
    scale_y_log10(breaks = c(0.5, 1, 2, 5, 10))+
    theme_classic()+
    labs(x = "Diameter (µm)", y = "Thickness (µm)", title =vec_form[i])+
    theme(axis.text = element_text(size = 10), text = element_text(size = 12))
  g
  #save plot
  #ggsave(filename=paste0(vec_form[i],".png"), width = 6, height = 6)
  # add plot to the p list
  p[[i]]<-g
  print(vec_form[i])
}

# make a figure with all species plots 
gg_form<-ggarrange(p[[1]], p[[3]], p[[2]], p[[5]],p[[4]],
                   ncol =3, nrow=2, common.legend = T)
gg_form
# save figure
ggsave(gg_form,filename="figures/Figure_S6_raw.png", width = 17.5, height = 12.4, units = "cm", dpi = 300)

# Kruskal Wallis tests # --------
# Differences across clades 

dt%>% group_by(clade)%>%
  summarise (n = n(),medD = median (max_Diameter, na.rm = T), 
             medT = median (Thickness, na.rm = T),
             med_TD = median (TD, na.rm = T), 
             med_Pcri1_low = median (Pcri1_low, na.rm = T),
             med_Pcri1_up = median (Pcri1_up, na.rm = T),
             med_Pcri2_low = median (Pcri2_low, na.rm = T),
             med_Pcri2_up = median (Pcri2_up, na.rm = T))

kruskal.test(TD ~ clade, data = dt)
kruskal.test(max_Diameter ~ clade, data = dt)
kruskal.test(Thickness ~ clade, data = dt)
kruskal.test(Pcri1_low ~ clade, data = dt)
kruskal.test(Pcri1_up ~ clade, data = dt)
kruskal.test(Pcri2_low ~ clade, data = dt)
kruskal.test(Pcri2_up ~ clade, data = dt)

pairwise.wilcox.test(dt$TD, dt$clade, p.adjust.method = "BH")
pairwise.wilcox.test(dt$max_Diameter, dt$clade, p.adjust.method = "BH")
pairwise.wilcox.test(dt$Thickness, dt$clade, p.adjust.method = "BH")
pairwise.wilcox.test(dt$Pcri1_low, dt$clade, p.adjust.method = "BH")
pairwise.wilcox.test(dt$Pcri1_up, dt$clade, p.adjust.method = "BH")
pairwise.wilcox.test(dt$Pcri2_low, dt$clade, p.adjust.method = "BH")
pairwise.wilcox.test(dt$Pcri2_up, dt$clade, p.adjust.method = "BH")


# Differences species habitats
glimpse(dt_ppt1)

dt_ppt1%>% group_by(ppt_binned)%>%
  summarise (n = n(),med_TD = median (TD, na.rm = T), 
             medD = median (max_Diameter, na.rm = T),
             medT = median (Thickness, na.rm = T),
             med_Pcri1_low = median (Pcri1_low, na.rm = T),
             med_Pcri1_up = median (Pcri1_up, na.rm = T),
             med_Pcri2_low = median (Pcri2_low, na.rm = T),
             med_Pcri2_up = median (Pcri2_up, na.rm = T))

# NAs is for Nym spp

kruskal.test(TD ~ ppt_binned, data = dt_ppt1)
kruskal.test(max_Diameter ~ ppt_binned, data = dt_ppt1)
kruskal.test(Thickness ~ ppt_binned, data = dt_ppt1)
kruskal.test(Pcri1_low ~ ppt_binned, data = dt_ppt1)
kruskal.test(Pcri1_up ~ ppt_binned, data = dt_ppt1)
kruskal.test(Pcri2_low ~ ppt_binned, data = dt_ppt1)
kruskal.test(Pcri2_up ~ ppt_binned, data = dt_ppt1)

pairwise.wilcox.test(dt_ppt1$TD, dt_ppt1$ppt_binned,p.adjust.method = "BH")
pairwise.wilcox.test(dt_ppt1$Thickness, dt_ppt1$ppt_binned,p.adjust.method = "BH")
pairwise.wilcox.test(dt_ppt1$max_Diameter, dt_ppt1$ppt_binned,p.adjust.method = "BH")
pairwise.wilcox.test(dt_ppt1$Pcri1_low, dt_ppt1$ppt_binned,p.adjust.method = "BH")
pairwise.wilcox.test(dt_ppt1$Pcri1_up, dt_ppt1$ppt_binned,p.adjust.method = "BH")
pairwise.wilcox.test(dt_ppt1$Pcri2_low, dt_ppt1$ppt_binned,p.adjust.method = "BH")
pairwise.wilcox.test(dt_ppt1$Pcri2_up, dt_ppt1$ppt_binned,p.adjust.method = "BH")

# Differences across growth forms 

dt%>% group_by(growth_form)%>%
  summarise (n = n(),med_TD = median (TD, na.rm = T), 
             medD = median (max_Diameter, na.rm = T),
             medT = median (Thickness, na.rm = T),
             mean_Pcri1_low = mean(Pcri1_low, na.rm = T),
             mean_Pcri1_up = mean(Pcri1_up, na.rm = T),
             mean_Pcri2_low = mean(Pcri2_low, na.rm = T),
             mean_Pcri2_up = mean(Pcri2_up, na.rm = T))



kruskal.test(TD ~ growth_form, data = dt)
kruskal.test(max_Diameter ~ growth_form, data = dt)
kruskal.test(Thickness ~ growth_form, data = dt)
kruskal.test(Pcri1_low ~ growth_form, data = dt)
kruskal.test(Pcri1_up ~ growth_form, data = dt)
kruskal.test(Pcri2_low ~ growth_form, data = dt)
kruskal.test(Pcri2_up ~ growth_form, data = dt)

pairwise.wilcox.test(dt$TD, dt$growth_form, p.adjust.method = "BH")
pairwise.wilcox.test(dt$max_Diameter, dt$growth_form, p.adjust.method = "BH")
pairwise.wilcox.test(dt$Thickness, dt$growth_form, p.adjust.method = "BH")
pairwise.wilcox.test(dt$Pcri1_low, dt$growth_form, p.adjust.method = "BH")
pairwise.wilcox.test(dt$Pcri1_up, dt$growth_form, p.adjust.method = "BH")
pairwise.wilcox.test(dt$Pcri2_low, dt$growth_form, p.adjust.method = "BH")
pairwise.wilcox.test(dt$Pcri2_up, dt$growth_form, p.adjust.method = "BH")

# Differences across vein orders 

dt_order%>% group_by(vein_order)%>%
  summarise (n = n(), med_TD = median (TD, na.rm = T), 
             medD = median (max_Diameter, na.rm = T),
             medT = median (Thickness, na.rm = T),
             med_Pcri1_low = median(Pcri1_low, na.rm = T),
             med_Pcri1_up = median(Pcri1_up, na.rm = T),
             med_Pcri2_low = median(Pcri2_low, na.rm = T),
             med_Pcri2_up = median(Pcri2_up, na.rm = T))

kruskal.test(TD ~ vein_order, data = dt_order)
kruskal.test(max_Diameter ~ vein_order, data = dt_order)
kruskal.test(Thickness ~ vein_order, data = dt_order)
kruskal.test(Pcri1_low ~ vein_order, data = dt_order)
kruskal.test(Pcri1_up ~ vein_order, data = dt_order)
kruskal.test(Pcri2_low ~ vein_order, data = dt_order)
kruskal.test(Pcri2_up ~ vein_order, data = dt_order)

pairwise.wilcox.test(dt_order$TD, dt_order$vein_order,p.adjust.method = "BH")
pairwise.wilcox.test(dt_order$max_Diameter, dt_order$vein_order, p.adjust.method = "BH")
pairwise.wilcox.test(dt_order$Thickness, dt_order$vein_order, p.adjust.method = "BH")
pairwise.wilcox.test(dt_order$Pcri1_low, dt_order$vein_order, p.adjust.method = "BH")
pairwise.wilcox.test(dt_order$Pcri1_up, dt_order$vein_order, p.adjust.method = "BH")
pairwise.wilcox.test(dt_order$Pcri2_low, dt_order$vein_order, p.adjust.method = "BH")
pairwise.wilcox.test(dt_order$Pcri2_up, dt_order$vein_order, p.adjust.method = "BH")

# Trade-offs - PCA # --------------
# Prepare data for PCA
dt_m<- dt%>% group_by (Species_name)%>%
  summarise(TD_m = median (TD, na.rm=T),
            LMA_m = median(LMA, na.rm = T),
            Kmax_m = median(Kmax, na.rm = T),
            e_m = median(e, na.rm = T),
            VTotV_m = median(VTotV, na.rm = T),
            clade = unique(Clade))%>%
  na.omit(e_m)%>% # removing species with missing elasticity data
  na.omit(TotV_m)%>% # removing species with missing total vlume of veins
  ungroup()%>%
  dplyr::select(TD_m, LMA_m, Kmax_m, e_m, VTotV_m, clade)

glimpse(dt_m)

# Run PCA
pca1 <- prcomp(dt_m[,-6],center=TRUE,scale=TRUE) # z-transforming variables prior to the PCA
pca1

# variance explained
varexp = 100*pca1$sdev^2/sum(pca1$sdev^2)
varexp


# get loadings
pca_rotations = as.data.frame(pca1$rotation[,1:2])
pca_rotations$var = row.names(pca_rotations)
pca_rotations$x0 = 0
pca_rotations$y0 = 0
row.names(pca_rotations) = NULL

scale_factor = 3

# overall PCA space 

dt_pca <-data.frame(pca1$x[,c(1,2,3)], Clades=as.factor(dt_m$clade)) # getting scores
glimpse(dt_pca)
dt_pca$Clades <- factor(dt_pca$Clades, levels = c("ferns", "basal angiosperms", "monocots", "basal eudicots", "rosids", "asterids"))


# Figure S7 # ---------------------
#  Broken Stick criterion method
library(vegan)
bstick(pca1)

png(file="figures/Figure_S7.png")
screeplot(pca1, bstick = T)
dev.off()

# Figure 5 # ----------------------
pca_plot<-ggplot(dt_pca,
                 aes(x = PC1, y = PC2, group= Clades, col = Clades))+
  geom_point(alpha=0.7, size = 3)+
  stat_ellipse() + 
  geom_segment(data=pca_rotations,aes(x=x0,y=y0,xend=scale_factor*PC1,yend=scale_factor*PC2),inherit.aes = FALSE,size=0.8,color='black') +
  geom_text(data=pca_rotations,aes(x=scale_factor*PC1-0.2,y=scale_factor*PC2-0.1,label=c("implosion safety", "cost (LMA)","efficiency", "support", "cost (VTotV)")),inherit.aes = FALSE,size=5,color='black') +  
  theme_bw() +
  xlab(sprintf("PC1 (%.2f%%)",varexp[1])) +
  ylab(sprintf("PC2 (%.2f%%)",varexp[2])) +
  scale_color_manual(values = c("asterids"="#D55E00",# dark orange
                                "rosids"="#CC79A7", # pink
                                "basal angiosperms" = "#56B4E9",# blue
                                "monocots" = "#009E73",# green
                                "basal eudicots" = "#F0E442",#yellow
                                "ferns" = "#999999"))+
  theme( legend.title=element_text(size=18),legend.text=element_text(size=14)); pca_plot
pca_plot
# save plot
ggsave(pca_plot, file='figures/Figure_5.png',width=8,height=6.5)

# Trade-offs - lm # ---------------
lm1<-lm(TD_m ~ Kmax_m, data = dt_m)
summary(lm1)
lm2<- lm(TD_m ~ e_m, data = dt_m)
summary(lm2)
lm3<- lm(TD_m ~ LMA_m, data = dt_m)
summary(lm3)
lm4<- lm(TD_m ~ VTotV_m, data = dt_m)
summary(lm4)

# Figure 6 ----------------------
glimpse(dt_m)
dt_m$clade <- factor(dt_m$clade, levels = c("ferns", "basal angiosperms", "monocots", "basal eudicots", "rosids", "asterids"))

fig6a<-ggplot(data = dt_m, aes(x = Kmax_m, y = TD_m))+
  geom_point(aes(col = clade), size = 3, alpha =0.7)+
  geom_smooth (method = "lm", col = "black")+
  stat_cor(label.x = 0, label.y = 0.9)+
  scale_color_manual(values = c("asterids"="#D55E00",# dark orange
                                "rosids"="#CC79A7", # pink
                                "basal angiosperms" = "#56B4E9",# blue
                                "monocots" = "#009E73",# green
                                "basal eudicots" = "#F0E442",#yellow
                                "ferns" = "#999999"))+
  labs(x = "Efficiency \n Kleaf max (mmol m-2 s-1 MPa-1)",
       y = "Implosion safety \n T/D", col = "Clades")+
  theme_classic(); fig6a

fig6b<-ggplot(data = dt_m, aes(x = e_m, y = TD_m))+
  geom_point(aes(col = clade), size = 3, alpha =0.7)+
  geom_smooth (method = "lm", col = "black")+
  stat_cor(label.x = 0, label.y = 0.9)+
  scale_color_manual(values = c("asterids"="#D55E00",# dark orange
                                "rosids"="#CC79A7", # pink
                                "basal angiosperms" = "#56B4E9",# blue
                                "monocots" = "#009E73",# green
                                "basal eudicots" = "#F0E442",#yellow
                                "ferns" = "#999999"))+
  labs(x = "Support \n ε (MN m-2)",
       y = "Implosion safety \n T/D", col = "Clades")+
  theme_classic(); fig6b

fig6c<-ggplot(data = dt_m, aes(x = LMA_m, y = TD_m))+
  geom_point(aes(col = clade), size = 3, alpha =0.7)+
  geom_smooth (method = "lm", col = "black")+
  stat_cor(label.x = 0, label.y = 0.9)+
  scale_color_manual(values = c("asterids"="#D55E00",# dark orange
                                "rosids"="#CC79A7", # pink
                                "basal angiosperms" = "#56B4E9",# blue
                                "monocots" = "#009E73",# green
                                "basal eudicots" = "#F0E442",#yellow
                                "ferns" = "#999999"))+
  labs(x = "Cost \n LMA (g m-2)",
       y = "Implosion safety \n T/D", col = "Clades")+
  theme_classic(); fig6c

fig6d<-ggplot(data = dt_m, aes(x = LMA_m, y = VTotV_m))+
  geom_point(aes(col = clade), size = 3, alpha =0.7)+
  geom_smooth (method = "lm", col = "black")+
  stat_cor(label.x = 0, label.y = 0.9)+
  scale_color_manual(values = c("asterids"="#D55E00",# dark orange
                                "rosids"="#CC79A7", # pink
                                "basal angiosperms" = "#56B4E9",# blue
                                "monocots" = "#009E73",# green
                                "basal eudicots" = "#F0E442",#yellow
                                "ferns" = "#999999"))+
  labs(x = "Cost \n VTotV (mm3 mm-2)",
       y = "Implosion safety \n T/D", col = "Clades")+
  theme_classic(); fig6d

Fig6<- ggarrange(fig6a, fig6b, fig6c, fig6d,
                 ncol =2, nrow=2, common.legend = T,
                 labels = c("(a)", "(b)", "(c)", "(d)",
                 legend = "bottom")
Fig6
ggsave(Fig6,filename="figures/Figure_6.png", width = 23, height = 9, limitsize = FALSE, units = "cm")
