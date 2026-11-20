#### Libraries ####
library(adespatial) # yes
library(ecotraj) # yes
library(ggplot2) # yes
library(ggpubr) # yes
library(ade4) # yes
library(ape) # yes 
library(cluster)
library(paletteer) # yes
library(multcompView) # yes
library(RColorBrewer) # yes
library(vegan) # yes
library(dplyr) # yes
library(stringr) # yes
library(tidyr) # yes
library(tibble) # yes

rm(list=ls()) # Clear environment
cat("\014")   # Clear console

Flore_tot <- read.csv(file ="data/Flore_totale_fin.CSV")
Flore_tot <- Flore_tot[, -1]

# Select data
Flore_cut <- Flore_tot %>% 
  filter(grepl("An", Sample_ID)) %>% 
  filter(grepl("Aut", Sample_ID)) 
# %>%  
#   select(-Hydropunctaria_maura, -Lichina_pygmaea, 
#          -Caloplaca_marina, -Wahlenbergiella_mucosa)


Flore <- Flore_cut[,-1]
rownames(Flore) <- Flore_cut[,1]

##############################################################################
#### Distance matrix creation ####
# Distance matrix creation
D_flore <- dist.ldc(Flore, "hellinger") # adespatial

is.euclid(D_flore)
is.metric(D_flore)


#############################################################################
# Create metadata ####
meta <- Flore %>%
  rownames_to_column("Sample_ID") %>%
  distinct(Sample_ID) %>%
  separate(Sample_ID, c("Site", "Season", "Year", "Point"),
           sep ="_",remove = FALSE) %>%
  unite(Site_point, Site, Point, sep = "_", remove = FALSE) %>%
  mutate(Year = as.numeric(Year)) %>%
  mutate_if(is.character, as.factor)


meta$Survey <- ifelse(meta$Season =="Aut" & meta$Year%in%c(2004:2007),1, 
                    ifelse(meta$Season =="Aut" & meta$Year%in%c(2008:2010),2.75, 
                     ifelse(meta$Season =="Aut" & meta$Year%in%c(2011:2013),3.75,
                      ifelse(meta$Season =="Aut" & meta$Year%in%c(2014:2016),4.75,
                       ifelse(meta$Season =="Aut" & meta$Year%in%c(2017:2019),5.75,
                        ifelse(meta$Season =="Aut" & meta$Year%in%c(2020:2022),6.75,
                                 ""))))))
meta$Survey <- ifelse(meta$Season =="Hiv" & meta$Year%in%c(2004:2007),1.25, 
                    ifelse(meta$Season =="Ete" & meta$Year%in%c(2004:2007),1.5,
                           meta$Survey))
meta$Survey <- ifelse(meta$Season =="Ptps" & meta$Year%in%c(2004:2007),1.75, 
                    ifelse(meta$Season =="Ptps" & meta$Year%in%c(2008:2010),2.25, 
                     ifelse(meta$Season =="Ptps" & meta$Year%in%c(2011:2013),3.25,
                      ifelse(meta$Season =="Ptps" & meta$Year%in%c(2014:2016),4.25,
                       ifelse(meta$Season =="Ptps" & meta$Year%in%c(2017:2019),5.25,
                        ifelse(meta$Season =="Ptps" & meta$Year%in%c(2020:2022),6.25,
                               meta$Survey))))))

# Create objects with sites and surveys
Sites <- meta$Site_point
Surveys <- meta$Survey

##############################################################################
############################# Trajectories ##################################
############################################################################
#### Display PCOA and trajectories ####  

DomCover <- read.csv(file ="Cover_AnFves.CSV")

DomCover2 <- DomCover %>% select(., -X) %>%
  spread(Espece, Cover, fill = 0) %>%
  mutate(., DOM = ifelse(Ascophyllum_nodosum + Fucus_vesiculosus <= 10, "NC",
                         ifelse(Ascophyllum_nodosum > Fucus_vesiculosus, "An",
                                ifelse(Ascophyllum_nodosum < Fucus_vesiculosus,
                                       "Fves", 
                                       ifelse(Ascophyllum_nodosum == Fucus_vesiculosus,
                                              "ND", NA )))))

meta2 <- meta %>%
  mutate(Dominance = DomCover2$DOM)%>%
  mutate_if(is.character, as.factor)

col.vec3_dark <- c('#E89242','#ABCD72', '#2D7DB4', 'black')

# PCOA computation
Flore_PCoA <- pcoa(D_flore, correction="none", rn=NULL)
# Calculate weighted average scores using the species data and PCoA result
weighted_scores <- wascores(Flore_PCoA$vectors, Flore)

# Extract the scores for PCoA1 and PCoA2 dimensions
weighted_scores_pcoa1 <- weighted_scores[, 1]
weighted_scores_pcoa2 <- weighted_scores[, 2]

# Select the top contributing species (e.g., top 10)
top_species <- names(weighted_scores_pcoa1)[order(weighted_scores_pcoa1, decreasing = TRUE)][1:9]
# create clusters
Df <- cbind(Flore_PCoA$vectors[, 1], Flore_PCoA$vectors[, 2])
Df <- Df %>% as.data.frame(.) %>%
  rownames_to_column("Sample_ID")
Df <- merge(Df, meta2, by = "Sample_ID")
centroids <- aggregate(cbind(V1,V2)~Dominance,Df,mean)

# Create a scatter plot of PCoA1 and PCoA2 dimensions
plot(Flore_PCoA$vectors[, 1], Flore_PCoA$vectors[, 2], 
     xlab = "PCoA1", ylab = "PCoA2", 
     col = col.vec3_dark, pch = 16, cex = 1, xlim=c(-0.8, 0.8), 
     ylim=c(-0.8, 0.8))

axis(side = 1, at = seq(-0.8, 0.8, by = 0.2))  # x-axis breaks
axis(side = 2, at = seq(-0.8, 0.8, by = 0.2))  # y-axis breaks

points(centroids$V1, centroids$V2, pch = 15, cex = 2, col = col.vec3_dark[centroids$Dominance])
points(centroids$V1, centroids$V2, pch = 22, cex = 2, col = "black")
# Add the weighted average projection as arrows onto the plot
arrows(0, 0, weighted_scores_pcoa1[top_species], weighted_scores_pcoa2[top_species], col = "red", length = 0.1)
# Add species labels for selected species
text(weighted_scores_pcoa1[top_species], weighted_scores_pcoa2[top_species], 
     labels = top_species, col = "red", pos = 4)
legend("bottomleft", legend = levels(meta2$Dominance), fill = col.vec3_dark, title = "Dominance")


# Display trajectory with trajectoryPCoA
par(mfrow = c(1,1))

col.vec <- c('#5ABCEB','#5ABCEB','#5ABCEB','#050708','#050708','#050708',
                      '#D33C32','#D33C32','#D33C32','#97A1A7','#97A1A7','#97A1A7',
                      '#ABCD72','#ABCD72','#ABCD72', '#FFCD00','#FFCD00','#FFCD00',
                      '#E89242','#E89242','#E89242', '#8A4198','#8A4198','#8A4198',
                      '#FB6467','#FB6467','#FB6467','#E7C9C6','#E7C9C6','#E7C9C6')
                      
name.vec.rac <- c('Be1','Be2','Be3','Br1','Br2','Br3','Cr1','Cr2','Cr3',
                  'De1','De2','De3','Mo1','Mo2','Mo3','Po1','Po2','Po3',
                  'Qb1','Qb2','Qb3','Si1','Si2','Si3','Sb1','Sb2','Sb3',
                  'Tg1','Tg2','Tg3')
traj <- defineTrajectories(D_flore, Sites, Surveys)
Obj <- trajectoryPCoA(traj, Sites,  traj.colors = col.vec, lwd = 2, 
                      survey.labels = F)
legend("bottomleft", col = col.vec, legend=name.vec.rac, bty="n", lty=1, lwd = 2)

# Display trajectory with ggplot
Col.cont2 <- paletteer_c("grDevices::Purple-Orange", 7)
Objet2 <- Obj$points %>% as.data.frame(Obj$points)%>%
  bind_cols(meta)

#col.vec2 <- c('#5ABCEB','#8A4198','#D33C32','grey40','#CF4E9C','#97A1A7', '#E89242','#FFCD00','#ABCD72','#E7C9C6')

desired_order <- c("StBriac", "Brehat", "SeptIles", "Portsall", 
                   "Molene", "Dellec", "Tregunc", "Quiberon", 
                   "Berchis", "Croisic")
Objet2$Site <- factor(Objet2$Site, levels = desired_order)
                      
col.vec2 <- c('#E89242','#97A1A7','#FFCD00','#D33C32','#ABCD72','grey40',
              '#CF4E9C','#5ABCEB','#8A4198','#E7C9C6')

ggplot(Objet2, aes(x = V1, y = V2, color = Site, group = Site_point)) +
  geom_point(size = 0.5) +
  geom_path(aes(x = V1, y = V2, color = Site), linewidth = 1,
            arrow = arrow(length = unit(0.3, "cm"), type = "closed")) +
  scale_color_manual(values = col.vec2)  +
  ylim(c(-0.75, 0.75))+
  xlim(c(-0.75, 0.75) )+
  xlab(paste("PCoA1", round(100*(Flore_PCoA$values$Eigenvalues[1]/sum(Flore_PCoA$values$Eigenvalues)), 1), "%", sep = " "))+
  ylab(paste("PCoA2", round(100*(Flore_PCoA$values$Eigenvalues[2]/sum(Flore_PCoA$values$Eigenvalues)), 1), "%", sep = " "))+
  theme_classic() 
#+  facet_wrap(as.factor(Objet2$Site), scales = "free")

################################################################
color_palette <- brewer.pal(9, "YlGn")
PALLAITE <- paletteer::paletteer_c("grDevices::ag_Sunset", 9, direction = -1)
PALLAITE <- PALLAITE[-1]
Objet2$Survey <- as.numeric(Objet2$Survey)
# Darken the arrow with time spent
# Define breaks and colors for arrows
time_breaks <- seq(min(Objet2$Survey), max(Objet2$Survey), length.out = length(color_palette) + 1)
time_colors <- PALLAITE[findInterval(Objet2$Survey, time_breaks, all.inside = TRUE)]

Objet2$Survey <- as.factor(Objet2$Survey)
TrajSep <- ggplot(Objet2, aes(x = V1, y = V2, color = time_colors, 
                              group = Site_point)) +
  geom_point(size = 1) +
  geom_path(aes(x = V1, y = V2, color = time_colors), linewidth = 1.1,
            arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
  scale_color_identity()+
  xlab(paste("PCoA1", round(100*(Flore_PCoA$values$Eigenvalues[1]/sum(Flore_PCoA$values$Eigenvalues)), 1), "%", sep = " "))+
  ylab(paste("PCoA2", round(100*(Flore_PCoA$values$Eigenvalues[2]/sum(Flore_PCoA$values$Eigenvalues)), 1), "%", sep = " "))+  
  ylim(c(-0.75, 0.75))+
  xlim(c(-0.75, 0.75))+
  facet_wrap(as.factor(Objet2$Site))+
  theme_classic()+
  theme(axis.line = element_line(colour = "black", linewidth = 1, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", size = 13)) + 
  theme(axis.title = element_text(face="bold", size = 13))+ 
  theme(axis.text.y = element_text(face="bold", size = 13))+
  theme(strip.text = element_text(face="bold",size=13))+ 
  theme(panel.grid.major = element_line(colour = "grey90", linewidth = 0.5))

TrajSep

##############################################################################
#############################################################################
############################################################################
Objet3 <- Objet2

# Add a new plot with the additional points
TrajSep <- ggplot(Objet2, aes(x = V1, y = V2)) +
  geom_point(size = 1, aes(color = time_colors, group = Site_point)) +
  geom_path(aes(x = V1, y = V2, color = time_colors, group = Site_point), 
            size = 1.1, arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
  scale_color_identity() +
  xlab(paste("PCoA1", round(100 * (Flore_PCoA$values$Eigenvalues[1] / 
                                     sum(Flore_PCoA$values$Eigenvalues)), 1), "%", sep = " ")) +
  ylab(paste("PCoA2", round(100 * (Flore_PCoA$values$Eigenvalues[2] / 
                                     sum(Flore_PCoA$values$Eigenvalues)), 1), "%", sep = " ")) +  
  ylim(c(-0.75, 0.75)) +
  xlim(c(-0.75, 0.75)) +
  facet_wrap(~ Site) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid")) +
  theme(axis.text.x = element_text(face = "bold", size = 13)) + 
  theme(axis.title = element_text(face = "bold", size = 13)) + 
  theme(axis.text.y = element_text(face = "bold", size = 13)) +
  theme(strip.text = element_text(face = "bold", size = 13)) + 
  theme(panel.grid.major = element_line(colour = "grey90", size = 0.5))

# Add the second geom_point with subset data
TrajSep + 
  geom_point(data = subset(Objet2, select = c(-Site, -Site_point)), aes(x = V1, y = V2), 
             color = "grey", pch = 20, size = 1)

##############################################################################
########################## Metric computation ###############################
############################################################################
#### Trajectory metrics ####
# Length metrics
library(plotrix)

Net_changes <- trajectoryLengths(D_flore, Sites, Surveys, relativeToInitial = T) 
head(Net_changes)

segment_lengths <- trajectoryLengths(D_flore, Sites, Surveys, relativeToInitial = F) 
head(segment_lengths)

mean(c(1.2239988, 1.3840103, 1.2567784))
std.error(c(1.2239988, 1.3840103, 1.2567784))

# Direction metrics
Angles <- trajectoryAngles(D_flore, Sites, Surveys)
head(Angles)
Direction <- trajectoryDirectionality(D_flore, Sites, Surveys)
head(Direction)
DD <- as.data.frame(Direction)

# distances between trajectories
Ds_flore <- trajectoryDistances(D_flore, Sites, surveys = NULL, 
                                distance.type = "DSPD", symmetrization = 'mean', 
                                add = TRUE, verbose = FALSE)
Ds_flore_m <- as.matrix(Ds_flore)
heatmap(Ds_flore_m, symm = TRUE)

is.euclid(Ds_flore)
# distances between directed segments
Ds <- segmentDistances(D_flore, Sites, Surveys)$Dseg
Ds <- as.matrix(Ds)
Ds

################################################################################
########################## Box plots and tests ################################
##############################################################################
# Segment Lentgh ####
Net_changesb <- tibble::rownames_to_column(Net_changes, "ID") 
Net_changesb <- Net_changesb %>% separate(ID, c("Site","Point"), sep ="_")

finalNC <- Net_changesb[, c(1,7)]
finalNC$Lt1_t6[1:3] <- Net_changesb$Lt1_t5[1:3]

# Define replacement patterns and their corresponding substitutions
patterns <- c("Berchis", "Brehat", "Croisic", "Dellec", "Molene", "Portsall", 
              "Quiberon", "StBriac", "SeptIles", "Tregunc")
replacements <- c("Be", "Br", "Cr", "De", "Mo", "Po", "Qb", "SB", "SI", "Tg")

# Use str_replace_all to replace all occurrences of the patterns with their corresponding replacements
finalNC$Site <- str_replace_all(finalNC$Site, 
                                         setNames(replacements, patterns))
finalNC$Site <- as.factor(finalNC$Site)

mod <- aov(Lt1_t6 ~ Site, data = finalNC)
# normality:
hist(mod$residuals)
shapiro.test(mod$residuals) 
# variance homogeneity:
bartlett.test(Lt1_t6 ~ Site, data = finalNC)

anova(mod)

TUKEY <- TukeyHSD(mod)
TUKEY

#######################################################################
generate_label_df <- function(TUKEY, variable){
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  Tukey.labels$Site <- rownames(Tukey.labels)
  Tukey.labels <- Tukey.labels[order(Tukey.labels$Site) , ]
  return(Tukey.labels)
}
LABELSNC <- generate_label_df(TUKEY , "Site")

# A panel of colors to draw each group with the same color :
my_colors=c("#001889", "#9900CC", "#FF7F00", "#FDE333")

a <- boxplot(finalNC$Lt1_t6 ~ finalNC$Site , ylim=c(min(finalNC$Lt1_t6) , 1.1*max(finalNC$Lt1_t6)) ,col=my_colors[as.factor(LABELSNC[,1])], ylab="final net change", xlab = "Sites" , main="")

over <- 0.05*max( a$stats[nrow(a$stats),] )
text(c(1:nlevels(finalNC$Site)) , 
     a$stats[nrow(a$stats),]+over , LABELSNC[,1] , 
     col=my_colors[as.factor(LABELSNC[,1])], font = 2)

################################################################################
# Segment Lentgh ####
segment_lengthsb <- tibble::rownames_to_column(segment_lengths, "ID") 
segment_lengthsb <- segment_lengthsb %>% separate(ID, c("Site","Point"), sep ="_")

# Define replacement patterns and their corresponding substitutions
patterns <- c("Berchis", "Brehat", "Croisic", "Dellec", "Molene", "Portsall", 
              "Quiberon", "StBriac", "SeptIles", "Tregunc")
replacements <- c("Be", "Br", "Cr", "De", "Mo", "Po", "Qb", "SB", "SI", "Tg")

# Use str_replace_all to replace all occurrences of the patterns with their corresponding replacements
segment_lengthsb$Site <- str_replace_all(segment_lengthsb$Site, 
                                         setNames(replacements, patterns))
segment_lengthsb$Site <- as.factor(segment_lengthsb$Site)

mod <- aov(Trajectory ~ Site, data = segment_lengthsb)
# normality:
shapiro.test(mod$residuals) 
# variance homogeneity:
bartlett.test(Trajectory ~ Site, data = segment_lengthsb)

anova(mod)

TUKEY <- TukeyHSD(mod)
TUKEY

###############################################################################
generate_label_df <- function(TUKEY, variable){
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  Tukey.labels$Site <- rownames(Tukey.labels)
  Tukey.labels <- Tukey.labels[order(Tukey.labels$Site) , ]
  return(Tukey.labels)
}
LABELSSL <- generate_label_df(TUKEY , "Site")

a <- boxplot(segment_lengthsb$Trajectory ~ segment_lengthsb$Site , 
             ylim=c(min(segment_lengthsb$Trajectory) , 1.1*max(segment_lengthsb$Trajectory)) ,col=my_colors[as.factor(LABELSSL[,1])], 
             ylab="Total segment length", xlab = "Sites", main="")

over <- 0.05*max( a$stats[nrow(a$stats),] )
text(c(1:nlevels(segment_lengthsb$Site)) , 
     a$stats[nrow(a$stats),]+over , LABELSSL[,1] , 
     col=my_colors[as.factor(LABELSSL[,1])], font = 2)


###############################################################################
#################### NC and SL temporal  selected ############################
#############################################################################
# SL temporal selected ####
# Put the NA's as the first segment for Berchis
SL <- segment_lengthsb[,-8]
SL <- SL %>% unite(Site_point, Site, Point, sep = "_", remove = FALSE)
Be <- SL[1:3,]
SL <- SL %>% filter(Site != "Be")
colnames(Be) <- c("ID","Site","Point", 2:5, 1)
Be <- Be %>% relocate("1", .before="2")
colnames(SL) <- c("ID","Site","Point", 1:5)
SLtot <- rbind(Be, SL)


# prepare the subset
SLline <- tidyr::pivot_longer(SLtot, 4:8, names_to ="Time_since_start", values_to="Segment_lengths", names_transform = function(x) {as.numeric(x)+1}) 

SLline <- SLline[order(SLline$Time_since_start, decreasing=FALSE),]
SLline <- as.data.frame(SLline)

SLline <- SLline %>% mutate(Site_ceinture  = ID)
SLline$Site_ceinture <- str_replace(SLline$Site_ceinture, "1", "")
SLline$Site_ceinture <- str_replace(SLline$Site_ceinture, "2", "")
SLline$Site_ceinture <- str_replace(SLline$Site_ceinture, "3", "")

SLline_means <- SLline %>% group_by(Site_ceinture, Time_since_start) %>%
  summarise(MeanSL = mean(Segment_lengths), .groups = "keep")
SLline_means <- as.data.frame(SLline_means)

SLline_select <- SLline %>% 
  filter(Site %in% c("Be", "De", "Qb"))
SLline_means_select <- SLline_means %>% 
  filter(Site_ceinture %in% c("Be_An", "De_An", "Qb_An"))

# create the plot to display net changes time series for all individuals
SLplot <- ggplot(data=SLline_select, aes(x=Time_since_start, y=Segment_lengths))+
  geom_path(data=SLline_select, aes(x=Time_since_start, y=Segment_lengths, group=ID), size =1, color = "#DFDEFF")+
  geom_path(data=SLline_means_select,aes(x=Time_since_start, y=MeanSL, group=Site_ceinture), size =1.3, color = "#9900CC")+
  facet_wrap(~Site_ceinture)+
  labs(x = "Survey Period", y = "Segment Length")+
  theme_classic()+
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", size = 16)) + 
  theme(axis.title = element_text(face="bold", size = 15))+ 
  theme(axis.text.y = element_text(face="bold", size = 16))+
  theme(strip.text = element_text(face="bold",size=15))
SLplot

###############################################################################
# NC temporal selected ####
# Use str_replace_all to replace all occurrences of the patterns with their corresponding replacements
Net_changesb$Site <- str_replace_all(Net_changesb$Site, 
                                setNames(replacements, patterns))

# Put the NA's as the first segment for Berchis
NC <- Net_changesb[,-8]
NC <- NC %>% unite(Site_point, Site, Point, sep = "_", remove = FALSE)
Be <- NC[1:3,]
NC <- NC %>% filter(Site != "Be")
colnames(Be) <- c("ID","Site","Point", 2:5, 1)
Be <- Be %>% relocate("1", .before="2")
colnames(NC) <- c("ID","Site","Point", 1:5)
NCtot <- rbind(Be, NC)


# prepare the subset
NCline <- tidyr::pivot_longer(NCtot, 4:8, names_to ="Time_since_start", values_to="Net_changes", names_transform = function(x) {as.numeric(x)+1}) 

NCline <- NCline[order(NCline$Time_since_start, decreasing=FALSE),]
NCline <- as.data.frame(NCline)

NCline <- NCline %>% mutate(Site_ceinture  = ID)
NCline$Site_ceinture <- str_replace(NCline$Site_ceinture, "1", "")
NCline$Site_ceinture <- str_replace(NCline$Site_ceinture, "2", "")
NCline$Site_ceinture <- str_replace(NCline$Site_ceinture, "3", "")

NCline_means <- NCline %>% group_by(Site_ceinture, Time_since_start) %>%
  summarise(MeanNC = mean(Net_changes), .groups = "keep")
NCline_means <- as.data.frame(NCline_means)
NCline_means_bis_for_datas <- NCline %>% group_by(Site_ceinture, Time_since_start) %>%
  summarise(MeanNC = mean(Net_changes),
            seNC = std.error(Net_changes), .groups = "keep")
NCline_means_bis_for_datas <- as.data.frame(NCline_means_bis_for_datas)
  
NCline_select <- NCline %>% 
  filter(Site %in% c("Be", "De", "Qb"))
NCline_means_select <- NCline_means %>% 
  filter(Site_ceinture %in% c("Be_An", "De_An", "Qb_An"))

# create the plot to display net changes time series for all individuals
NCplot <- ggplot(data=NCline_select, aes(x=Time_since_start, y=Net_changes))+
  geom_path(data=NCline_select,aes(x=Time_since_start, y=Net_changes, group=ID), 
            size =1, color = "#FFDCD4")+
  geom_path(data=NCline_means_select,aes(x=Time_since_start, y=MeanNC, 
                                  group= Site_ceinture), size =1.3, 
            color = "#E14BA8")+
  facet_wrap(~Site_ceinture)+
  labs(x = "Survey Period", y = "Net change")+
  theme_classic()+
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", size = 16)) + 
  theme(axis.title = element_text(face="bold", size = 15))+ 
  theme(axis.text.y = element_text(face="bold", size = 16))+
  theme(strip.text = element_text(face="bold",size=15))
NCplot

ggarrange(SLplot, NCplot, labels = c("A", "B"), ncol = 1, nrow = 2)

###############################################################################
#################### NC and SL temporal cluster ##############################
#############################################################################

# Center data ####
D_center <- centerTrajectories(D_flore, Sites, verbose = FALSE)
Ds_center <- trajectoryDistances(D_center, Sites, surveys = NULL, 
                                distance.type = "DSPD", symmetrization = 'mean', 
                                add = TRUE, verbose = FALSE)
#Cluster on centered data
# hierarchical clustering on the symmetric matrix
par(mfrow = c(1,1))
colstd <- c("black","yellow","green","blue","grey","red")

hsxy <- hclust(Ds_center, "ward.D2")

plot(hsxy,hang = -1, main = "Distance entre trajectoires", cex=.6)
Hst=0.95 # Cutting height
x <- rect.hclust(hsxy, h = Hst, border = colstd)

# Merge the groups and flore cut data
# cut the dendrogram at height to obtain a vector of cluster membership and copy it in furseals as a factor
Groups <- cutree(hsxy, h=Hst)
Groups <- as.data.frame(Groups)
Clusters <- tibble::rownames_to_column(Groups, "Site_point")
Clusters <- Clusters %>% separate(Site_point, c("Site", "Point"))
Clusters$Site <- str_replace_all(Clusters$Site, setNames(replacements, patterns))
Clusters <- Clusters %>% unite(Site_point, Site, Point, sep = "_")

# Merge Groups and SL data
SLline2 <- rename(SLline, Site_point = ID)
SLline2 <- merge(SLline2, Clusters, by = "Site_point")
SLline2$Groups <- as.factor(SLline2$Groups)

SLline2$Site_ceinture = factor(SLline2$Site_ceinture, levels=c("SB_An", "Br_An", "SI_An", "Po_An", "Mo_An", "De_An", "Tg_An", "Qb_An", "Be_An", "Cr_An"))
SLline_means$Site_ceinture = factor(SLline_means$Site_ceinture, levels=c("SB_An", "Br_An", "SI_An", "Po_An", "Mo_An", "De_An", "Tg_An", "Qb_An", "Be_An", "Cr_An"))


SLline2 <- SLline2[order(SLline2$Time_since_start, decreasing=FALSE),]
SLline2 <- as.data.frame(SLline2)
# Plot SL for all sites with clusters ####
# create the plot to display net changes time series for all individuals
ggplot(data=SLline2, aes(x=Time_since_start, y=Segment_lengths, color=Groups))+
  geom_path(data=SLline2,aes(x=Time_since_start, y=Segment_lengths, 
                             group=Site_point, color=Groups), 
            arrow = arrow(length = unit(0.10, "cm")), size =1 )+
  geom_path(data=SLline_means,aes(x=Time_since_start, y=MeanSL, 
                                  group= Site_ceinture), size =1.3, 
            color = "grey", linetype = "dashed")+
  scale_color_manual(values = c("#FFA347", "violetred",  "#5E0CCA"))+
  facet_wrap(~Site_ceinture, 
             labeller = as_labeller(c(`SB_An` = "Saint-Briac", `Br_An` = "Brehat",
                                      `SI_An` = "Sept-Iles",`Po_An` = "Portsall",
                                      `Mo_An` = "Molene", `De_An` = "Dellec",
                                      `Tg_An` = "Tregunc", `Qb_An` = "Quiberon",
                                      `Be_An` = "Berchis",`Cr_An` = "Croisic" )))+
  labs(x = "Survey Period", y = "Segment Length")+
  theme_classic()+
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", size = 16)) + 
  theme(axis.title = element_text(face="bold", size = 15))+ 
  theme(axis.text.y = element_text(face="bold", size = 16))+
  theme(strip.text = element_text(face="bold",size=15))



NCline2 <- rename(NCline, Site_point = ID)
NCline2 <- merge(NCline2, Clusters, by = "Site_point")
NCline2$Groups <- as.factor(NCline2$Groups)

NCline2 <- NCline2[order(NCline2$Time_since_start, decreasing=FALSE),]
NCline2 <- as.data.frame(NCline2)

NCline2$Site_ceinture = factor(NCline2$Site_ceinture, levels=c("SB_An", "Br_An", "SI_An", "Po_An", "Mo_An", "De_An", "Tg_An", "Qb_An", "Be_An", "Cr_An"))
NCline_means$Site_ceinture = factor(NCline_means$Site_ceinture, levels=c("SB_An", "Br_An", "SI_An", "Po_An", "Mo_An", "De_An", "Tg_An", "Qb_An", "Be_An", "Cr_An"))

# Plot NC for all sites with clusters ####
# create the plot to display net changes time series for all individuals
ggplot(data=NCline2, aes(x=Time_since_start, y=Net_changes, color=Groups))+
  geom_path(data=NCline2,aes(x=Time_since_start, y=Net_changes, 
                             group=Site_point, color=Groups), 
            arrow = arrow(length = unit(0.10, "cm")), size =1 )+
  geom_path(data=NCline_means,aes(x=Time_since_start, y=MeanNC, 
                                         group= Site_ceinture), size =1.3, 
            color = "grey", linetype = "dashed")+
  scale_color_manual(values = c("#FFA347", "violetred",  "#5E0CCA"))+
  facet_wrap(~Site_ceinture, 
             labeller = as_labeller(c(`SB_An` = "Saint-Briac", `Br_An` = "Brehat",
                                      `SI_An` = "Sept-Iles",`Po_An` = "Portsall",
                                      `Mo_An` = "Molene", `De_An` = "Dellec",
                                      `Tg_An` = "Tregunc", `Qb_An` = "Quiberon",
                                      `Be_An` = "Berchis",`Cr_An` = "Croisic")))+
  labs(x = "Survey Period", y = "Net change")+
  theme_classic()+
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid")) +
  theme(axis.text.x = element_text(face="bold", size = 16)) + 
  theme(axis.title = element_text(face="bold", size = 15))+ 
  theme(axis.text.y = element_text(face="bold", size = 16))+
  theme(strip.text = element_text(face="bold",size=15))

###############################################################################
#################### NC and SL as maps #######################################
#############################################################################
##MAP
library(tidyterra)
library(terra)
library(ggplot2)
library(ggspatial)

setwd("~/FAC BIOLOGIE/Stage IUEM 2023/Figures/Figs/QGIS")

maptry <- vect("Europe_coastline.shp")
maptry2 <- vect("DEPARTEMENT.shp")

lon_min <- -5 
lon_max <- -1 
lat_min <- 47 
lat_max <- 49 
extent <- ext(lon_min, lon_max, lat_min, lat_max)

mapTRY4326 <- project(maptry2, "epsg:4326")
mapHR_WestEurope <- crop(mapTRY4326, extent)

SitesGPS <- read.csv('Points_GPS_2.CSV')
SitesGPS <- SitesGPS %>% 
  filter(Site != "Ile Callot") %>% 
  filter(Site != "Le Caro")

# Define replacement patterns and their corresponding substitutions
patterns <- c("Berchis", "Brehat", "Croisic", "Dellec", "Molene", "Portsall", "Quiberon", "Saint-Briac", "Sept-Iles", "Tregunc")
replacements <- c("Be", "Br", "Cr", "De", "Mo", "Po", "Qb", "SB", "SI", "Tg")

###############################################################################
#### Map x2 with significance as color ####
# Use str_replace_all to replace all occurrences of the patterns with their corresponding replacements
SitesGPS$Site <- str_replace_all(SitesGPS$Site, setNames(replacements, patterns))

finalNC_mean <- finalNC %>% group_by(Site) %>% summarise(MeanNC = mean(Lt1_t6), se = std.error(Lt1_t6))
GroupsNC <- merge(finalNC_mean, LABELSNC, by = "Site")

OverallSL_mean <- segment_lengthsb %>% select(., Site, Trajectory)%>% 
  group_by(Site) %>% 
  summarise(MeanSL = mean(Trajectory), seSL = std.error(Trajectory))
GroupsSL <- merge(OverallSL_mean, LABELSSL, by = "Site")


SitesGPS_SL <- merge(SitesGPS, GroupsSL, by = "Site")
SitesGPS_NC <- merge(SitesGPS, GroupsNC, by = "Site")

SitesGPS_SL$Letters <- as.factor(SitesGPS_SL$Letters)
SitesGPS_NC$Letters <- as.factor(SitesGPS_NC$Letters)

size_breaks1 <- seq(3, 4.6, by = 0.4)
size_breaks2 <- seq(0.5, 1.5, by = 0.1)

my_colors2=c("#FF7F00", "#9900CC", "#001889")
#Carte europe france + brest + echelle + NOES + EEA 
P <- ggplot() +
  geom_spatvector(data=mapHR_WestEurope, color="grey30", fill="gray95") +
  coord_sf(xlim=c(lon_min,lon_max), ylim=c(lat_min,lat_max)) +
  geom_point(data = SitesGPS_SL, mapping = aes(x = Longitude, y = Latitude, size = MeanSL, color = Letters), pch = 20)+
  geom_text(data = SitesGPS_SL,
            aes(x = Longitude, y = Latitude, label = Site, color = Letters), hjust = +0.9, vjust = -1, fontface = "bold", size = 5)+
  scale_size_continuous(range = c(5, 20), breaks = size_breaks1) +
  scale_color_manual(values = my_colors) +
  annotation_scale(location="bl",
                   width_hint=0.40,
                   line_width=0.25,
                   height=unit(0.15,"cm")) +
  annotation_north_arrow(location="tl",
                         style=north_arrow_nautical,
                         pad_x=unit(1.10,"cm"),
                         pad_y=unit(0.75,"cm")) +
  guides(color = "none")+
  labs(size = "Trajectory Path Length")+
  theme_bw()+
  theme(panel.background = element_rect(fill = "white"))+
  theme(axis.text.x = element_text(face="bold", size = 13)) + 
  theme(axis.title = element_text(face="bold", size = 15))+ 
  theme(axis.text.y = element_text(face="bold", size = 13))+
  theme(legend.title = element_text(size = 15),
        legend.position = "top", 
        legend.justification = "center")



Q <- ggplot() +
  geom_spatvector(data=mapHR_WestEurope, color="grey30", fill="gray95") +
  coord_sf(xlim=c(lon_min,lon_max), ylim=c(lat_min,lat_max)) +
  geom_point(data = SitesGPS_NC, 
             mapping = aes(x = Longitude, y = Latitude, size = MeanNC,
                           color = Letters), pch = 20)+
  geom_text(data = SitesGPS_NC,
            aes(x = Longitude, y = Latitude, label = Site, color = Letters), hjust = +0.9, vjust = -1, fontface = "bold", size = 5)+
  scale_size_continuous(range = c(5, 20), breaks = size_breaks2) +
  scale_color_manual(values = my_colors2) +
  annotation_scale(location="bl",
                   width_hint=0.40,
                   line_width=0.25,
                   height=unit(0.15,"cm")) +
  annotation_north_arrow(location="tl",
                         style=north_arrow_nautical,
                         pad_x=unit(1.10,"cm"),
                         pad_y=unit(0.75,"cm")) +
  guides(color = "none")+
  labs(size = "Total Net Change")+
  theme_bw()+
  theme(panel.background = element_rect(fill = "white"))+
  theme(axis.text.x = element_text(face="bold", size = 13)) + 
  theme(axis.title = element_text(face="bold", size = 15))+ 
  theme(axis.text.y = element_text(face="bold", size = 13))+
  theme(legend.title = element_text(size = 15),
        legend.position = "top", 
        legend.justification = "center")


ggarrange(P, Q, labels = c("A", "B"), ncol = 2)


color_palette <- paletteer_c("grDevices::Zissou 1", 70)
color_breaks <- seq(0.8, 1.4, by = 0.2)
color_breaks2 <- seq(3, 4.5, by = 0.5)

size_breaks1 <- seq(3.2, 4.4, by = 0.4)
size_breaks2 <- seq(0.8, 1.2, by = 0.1)

###############################################################################
#### Map x1 with color and size as scales ####
library(ggnewscale)

ggplot() +
  geom_spatvector(data=mapHR_WestEurope, color="grey30", fill="gray95") +
  coord_sf(xlim=c(lon_min,lon_max), ylim=c(lat_min,lat_max)) +
  geom_point(data = SitesGPS_NC, 
             mapping = aes(x = Longitude, y = Latitude, 
                           color = MeanNC, size = MeanNC), 
             pch = 17)+
  scale_color_gradientn(limits = c(min(color_breaks), max(color_breaks)),
                        breaks = color_breaks,
                        colors = color_palette)+
  scale_size_continuous(range = c(4, 10), breaks = size_breaks2) +
  new_scale_colour() +
  new_scale("size")+
  geom_point(data = SitesGPS_SL, 
             mapping = aes(x = (Longitude-0.18), y = Latitude, 
                           color = MeanSL, size = MeanSL), 
             pch = 20)+
  scale_color_gradientn(limits = c(min(color_breaks2), max(color_breaks2)),
                        breaks = color_breaks2,
                        colors = color_palette)+
  scale_size_continuous(range = c(5, 15), 
                        breaks = size_breaks1) +
  geom_text(data = SitesGPS_NC,
            aes(x = Longitude, y = Latitude, label = Site), 
            hjust = +0.9, vjust = -1, 
            fontface = "bold", size = 5)+
  annotation_scale(location="bl",
                   width_hint=0.40,
                   line_width=0.25,
                   height=unit(0.15,"cm")) +
  annotation_north_arrow(location="tl",
                         style=north_arrow_nautical,
                         pad_x=unit(1.10,"cm"),
                         pad_y=unit(0.75,"cm")) +
  guides(color = "none")+
  theme_bw()+
  theme(panel.background = element_rect(fill = "white"))+
  theme(axis.text.x = element_text(face="bold", size = 13)) + 
  theme(axis.title = element_text(face="bold", size = 15))+ 
  theme(axis.text.y = element_text(face="bold", size = 13))+
  theme(legend.title = element_text(size = 15),
        legend.position = "top", 
        legend.justification = "center")

###############################################################################
#### Map x2 with color and size as scales ####
#Carte europe france + brest + echelle + NOES + EEA 
TPL <- ggplot() +
  geom_spatvector(data=mapHR_WestEurope, color="grey30", fill="gray95") +
  coord_sf(xlim=c(lon_min,lon_max), ylim=c(lat_min,lat_max)) +
  geom_point(data = SitesGPS_SL, mapping = aes(x = Longitude, 
                                               y = Latitude, size = MeanSL), 
             color = "#001889", pch = 20)+
  geom_text(data = SitesGPS_SL,
            aes(x = Longitude, y = Latitude, label = Site), 
            hjust = +0.9, vjust = -2, fontface = "bold", size = 5)+
  scale_size_continuous(range = c(4, 28), breaks = size_breaks1) +
  annotation_scale(location="bl",
                   width_hint=0.40,
                   line_width=0.25,
                   height=unit(0.15,"cm")) +
  annotation_north_arrow(location="tl",
                         style=north_arrow_nautical,
                         pad_x=unit(0.6,"cm"),
                         pad_y=unit(0.6,"cm")) +
  labs(size = "Trajectory Path Length", color = "")+
  theme_bw()+
  theme(panel.background = element_rect(fill = "white"))+
  theme(axis.text.x = element_text(face="bold", size = 18)) + 
  theme(axis.title = element_text(face="bold", size = 20))+ 
  theme(axis.text.y = element_text(face="bold", size = 18))+
  theme(legend.title = element_text(size = 18),
        legend.position = "top", 
        legend.justification = "center", 
        legend.text = element_text(size = 18))



TNC <- ggplot() +
  geom_spatvector(data=mapHR_WestEurope, color="grey30", fill="gray95") +
  coord_sf(xlim=c(lon_min,lon_max), ylim=c(lat_min,lat_max)) +
  geom_point(data = SitesGPS_NC, 
             mapping = aes(x = Longitude, y = Latitude, size = MeanNC),
             color = "#FFA347", pch = 20)+
  geom_text(data = SitesGPS_NC,
            aes(x = Longitude, y = Latitude, label = Site), hjust = +0.9, 
            vjust = -2, fontface = "bold", size = 5)+
  scale_size_continuous(range = c(4, 28), breaks = size_breaks2) +
  annotation_scale(location="bl",
                   width_hint=0.40,
                   line_width=0.25,
                   height=unit(0.15,"cm")) +
  annotation_north_arrow(location="tl",
                         style=north_arrow_nautical,
                         pad_x=unit(0.6,"cm"),
                         pad_y=unit(0.6,"cm")) +
  labs(size = "Total Net Change", color = "")+
  theme_bw()+
  theme(panel.background = element_rect(fill = "white"))+
  theme(axis.text.x = element_text(face="bold", size = 18)) + 
  theme(axis.title = element_text(face="bold", size = 20))+ 
  theme(axis.text.y = element_text(face="bold", size = 18))+
  theme(legend.title = element_text(size = 18),
        legend.position = "top", 
        legend.justification = "center", 
        legend.text = element_text(size = 18))


ggarrange(TPL, TNC, labels = c("A", "B"), ncol = 2)

###############################################################################
#### Map x1 with only size as scale ####
library(ggnewscale)

ggplot() +
  geom_spatvector(data=mapHR_WestEurope, color="grey30", fill="gray95") +
  geom_point(data = SitesGPS_NC, 
             mapping = aes(x = (Longitude-0.08), y = Latitude, size = MeanNC), 
             pch = 20, color =  "#FFA347")+
  scale_size_continuous(range = c(4, 16), breaks = size_breaks2) +
  new_scale("size")+
  geom_point(data = SitesGPS_SL, 
             mapping = aes(x = (Longitude+0.08), y = Latitude, 
                           color = MeanSL, size = MeanSL), 
             pch = 20, color = "#001889")+
  scale_size_continuous(range = c(4, 16), breaks = size_breaks1) +
  geom_text(data = SitesGPS_NC,
            aes(x = Longitude, y = Latitude, label = Site), 
             vjust = -1.5, fontface = "bold", size = 5)+
  annotation_scale(location="bl",
                   width_hint=0.40,
                   line_width=0.25,
                   height=unit(0.15,"cm")) +
  annotation_north_arrow(location="tl",
                         style=north_arrow_nautical,
                         pad_x=unit(0.9,"cm"),
                         pad_y=unit(0.6,"cm")) +
  labs(x = "Longitude")+
  coord_sf(xlim=c(lon_min,lon_max), ylim=c(lat_min,lat_max)) +
  theme_bw()+
  theme(axis.text.x = element_text(face="bold", size = 13)) + 
  theme(axis.title = element_text(face="bold", size = 15))+ 
  theme(axis.text.y = element_text(face="bold", size = 13))+
  theme(legend.title = element_text(size = 15, face = "bold"), 
        legend.text = element_text(size = 13, face = "bold"),
        legend.position = "top", 
        legend.justification = "center")
# +  theme(panel.border = element_blank()) # to remove the frame
