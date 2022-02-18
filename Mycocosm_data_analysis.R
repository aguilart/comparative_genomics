

####Mycocosms data analysis: comparative genomics of functional groups

rm(list=ls())

library(tidyverse)
library(vegan)

#Loading files assembled in code Mycocosm_data.R
geneFreq_dat_wFunct <- read.csv("gene_freq.csv", header = T, stringsAsFactors = F)# Species x gene frequency table
funct_groups <- read.csv("funct_groups.csv", header = T, stringsAsFactors = F)#Functional group info associated to a given species
taxonomy <- read.csv("taxonomy.csv", header = T, stringsAsFactors = F)# Standardized taxonomy for each species
gene_descriptions <- readRDS("gene_descriptions.RDS")#Functional description associated to each 

#1.For some reason there are 350 gene families for which all entries are 0
c0 <- sapply(geneFreq_dat_wFunct[, -c(1:5)], function(x){all(x == 0)})
length(which(c0 ==T))

gene_descriptions$description[gene_descriptions$family%in%names(c0[which(c0==T)])]

#However almost half of the families have an average lower of 0.1 and a variance
#lower of 0.1. 
varianza<-sapply(geneFreq_dat_wFunct[,-c(1:5)],var)
promedio<-sapply(geneFreq_dat_wFunct[,-c(1:5)],mean)
s_deviation<-sapply(geneFreq_dat_wFunct[,-c(1:5)],sd)
c_variation<-s_deviation/promedio
el_max<-sapply(geneFreq_dat_wFunct[,-c(1:5)],max)
el_min<-sapply(geneFreq_dat_wFunct[,-c(1:5)],min)
rango<-el_max-el_min

#I would assume it makes sense to exclude those because they driver almost no variation
t_e<-names(promedio)[which(promedio<0.1&varianza<0.1)]

t_e<-which(names(geneFreq_dat_wFunct)%in%t_e)


#PERMANOVA 
model_perm<-adonis(geneFreq_dat_wFunct[,-c(1:5,t_e)]~group,data = funct_groups,
                   method = "bray")


#To visualize the patterns:

#1. Principal component analysis (based on euclidean distances)

model_pca<-rda(geneFreq_dat_wFunct[,-c(1:5,t_e)],scale=T)

#2. I am using Principal Coordinate analysis (because euclidean distances
#do not hold well 0 values (which happen a lot here))
model_pcoa <- vegdist(geneFreq_dat_wFunct[,-c(1:5,t_e)])#the default are bray curtis distances
model_pcoa <- wcmdscale(model_pcoa, eig = TRUE)

#or using dbRDA (because it allows to obtain species scores)
model_pcoa2 <- vegdist(geneFreq_dat_wFunct[,-c(1:5,t_e)])#the default are bray curtis distances
model_pcoa2 <- dbrda(model_pcoa2~1)
sppscores(model_pcoa2)<-geneFreq_dat_wFunct[,-c(1:5,t_e)]


#3. NMDS. This because PCoA does not give me species scores, thus I have no way
#to know which variables drive most of each axis.
# model_nmds<-metaMDS(geneFreq_dat_wFunct[,-c(1:5,t_e)],
#                     distance = "bray",
#                     k=3,maxit=999,trymax = 1000)
# 
# model_nmds1<-metaMDS(geneFreq_dat_wFunct[,-c(1:5,t_e)],
#                      distance = "bray",
#                      k=3,maxit=999,trymax = 1000,previous.best = model_nmds)
# 
# model_nmds2<-metaMDS(geneFreq_dat_wFunct[,-c(1:5,t_e)],
#                      distance = "bray",
#                      k=3,maxit=999,trymax = 1000,previous.best = model_nmds1)
# 
# model_nmds3<-metaMDS(geneFreq_dat_wFunct[,-c(1:5,t_e)],
#                      distance = "bray",
#                      k=3,maxit=999,trymax = 1000,previous.best = model_nmds2)
# 
#So, NMDS does not work! It nevers find a stable solution



#Some aesthetic parameters
my_theme<-
  theme(#title = element_text(size = 18),
    axis.title=element_text(size=20),
    panel.background = element_blank(),
    panel.grid.major.x = element_line(size=0.25,linetype = "longdash",colour = "gray"),
    panel.grid.major.y = element_line(size=0.25,linetype = "longdash",colour = "gray"),
    #axis.text.x = element_text(size = 20,angle = 45,hjust = 1),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    strip.text.x = element_text(size = 20),
    strip.text.y = element_text(size = 25),
    legend.position = "bottom"
  )

#Because there are so many groups, here I am creating this vector to subset
#only the taxa with exclusive functional groups

grupos<-
  names(
    table(funct_groups$group[-grep("-",funct_groups$group)])[
      table(funct_groups$group[-grep("-",funct_groups$group)])>3])


#PCA
mod<-model_pca
e1<-round(summary(mod)[["cont"]][["importance"]][2,1],2)*100
e2<-round(summary(mod)[["cont"]][["importance"]][2,2],2)*100

#Extracting the scores for the species
vectors<-as.data.frame(scores(mod,display = "sites",scaling = "sites",choices=c(1,2)))
vectors$group<-funct_groups$group;vectors<-vectors[,c(3,1,2)]
vectors$taxon<-funct_groups$taxon
vectors$phylum<-taxonomy$phylum

#Extracting the scores for the gene families
e_vectors_traits_scaled_a<-as.data.frame(
  scores(mod,display = "species",choices = c(1,2),scaling = "sites"));
e_vectors_traits_scaled_a$traits<-rownames(e_vectors_traits_scaled_a)

names(e_vectors_traits_scaled_a)[1]<-"RDA1_species"
names(e_vectors_traits_scaled_a)[2]<-"RDA2_species"

#further organizing the information on the gene family scores
e_vectors_traits_scaled_b<-left_join(e_vectors_traits_scaled_a,gene_descriptions %>% 
                                       rename(traits=family))

e_vectors_traits_scaled_b$RDA1_species<-abs(e_vectors_traits_scaled_b$RDA1_species)
e_vectors_traits_scaled_b$RDA2_species<-abs(e_vectors_traits_scaled_b$RDA2_species)

e_vectors_traits_scaled_a<-left_join(e_vectors_traits_scaled_a,gene_descriptions %>% 
                                       rename(traits=family))

#Visualizing the data: 

#The top 20 gene families and their function driving the separation in PCA

gene_fam_PCscores_a<-e_vectors_traits_scaled_a
gene_fam_PCscores_b<-e_vectors_traits_scaled_b

#First axis
gene_fam_PCscores_b[order(gene_fam_PCscores_b$RDA1_species,
                          decreasing = T),][1:20,c(3,4)]

#Second axis
gene_fam_PCscores_b[order(gene_fam_PCscores_b$RDA2_species,
                          decreasing = T),][1:20,c(3,4)]

#Plotting:
pca_plot<-
  ggplot()+
  geom_point(aes(x=PC1,y=PC2,col=phylum),size=2,#colour="red",
             data=vectors %>% 
               filter(group%in%grupos))+
  scale_color_brewer(palette = "Set1")+
  
  ggrepel::geom_text_repel(size=5,aes(x=RDA1_species,y=RDA2_species,label=description,fontface="bold"),
                           #hjust = 1.1, vjust =-0.5,#segment.size = 0,segment.color = 'transparent',
                           data = e_vectors_traits_scaled_a[e_vectors_traits_scaled_b$RDA1_species>1.65,])+
  geom_segment(aes(x=0, y=0, xend=RDA1_species, yend=RDA2_species), 
               arrow = arrow(),color="gray",
               data = e_vectors_traits_scaled_a[e_vectors_traits_scaled_b$RDA1_species>1.65,])+
  
  labs(y=paste("PC2",e2,"%",sep = " "),
       x=paste("PC1",e1,"%",sep = " "),title = "PCA")+
  my_theme      



#PCoA
mod<-model_pcoa2
e1<-round(summary(mod)[["cont"]][["importance"]][2,1],2)*100
e2<-round(summary(mod)[["cont"]][["importance"]][2,2],2)*100

#Extracting the scores for the species
vectors<-as.data.frame(scores(mod,display = "sites",scaling = "sites",choices=c(1,2)))
vectors$group<-funct_groups$group;vectors<-vectors[,c(3,1,2)]
vectors$taxon<-funct_groups$taxon
vectors$phylum<-taxonomy$phylum

#Extracting the scores for the gene families
e_vectors_traits_scaled_a<-as.data.frame(
  scores(mod,display = "species",choices = c(1,2),scaling = "sites"));
e_vectors_traits_scaled_a$traits<-rownames(e_vectors_traits_scaled_a)

names(e_vectors_traits_scaled_a)[1]<-"RDA1_species"
names(e_vectors_traits_scaled_a)[2]<-"RDA2_species"

#further organizing the information on the gene family scores
e_vectors_traits_scaled_b<-left_join(e_vectors_traits_scaled_a,gene_descriptions %>% 
                                       rename(traits=family))

e_vectors_traits_scaled_b$RDA1_species<-abs(e_vectors_traits_scaled_b$RDA1_species)
e_vectors_traits_scaled_b$RDA2_species<-abs(e_vectors_traits_scaled_b$RDA2_species)

e_vectors_traits_scaled_a<-left_join(e_vectors_traits_scaled_a,gene_descriptions %>% 
                                       rename(traits=family))

#Visualizing the data: 

#The top 20 gene families and their function driving the separation in PCA

gene_fam_PCoAscores_a<-e_vectors_traits_scaled_a
gene_fam_PCoAscores_b<-e_vectors_traits_scaled_b

#First axis
gene_fam_PCoAscores_b[order(gene_fam_PCoAscores_b$RDA1_species,
                            decreasing = T),][1:20,c(3,4)]

#Second axis
gene_fam_PCoAscores_b[order(gene_fam_PCoAscores_b$RDA2_species,
                            decreasing = T),][1:20,c(3,4)]

#Plotting:
pcoa_plot2<-
  ggplot()+
  geom_point(aes(x=MDS1,y=MDS2,col=phylum),size=2,#colour="red",
             data=vectors %>% 
               filter(group%in%grupos))+
  scale_color_brewer(palette = "Set1")+
  
  ggrepel::geom_text_repel(size=4,aes(x=RDA1_species/7,y=RDA2_species/7,label=description,fontface="bold"),
                           #hjust = 1.1, vjust =-0.5,#segment.size = 0,segment.color = 'transparent',
                           data = e_vectors_traits_scaled_a[e_vectors_traits_scaled_b$RDA1_species>1.8,])+
  geom_segment(aes(x=0, y=0, xend=RDA1_species/7, yend=RDA2_species/7), 
               arrow = arrow(),color="gray",
               data = e_vectors_traits_scaled_a[e_vectors_traits_scaled_b$RDA1_species>1.8,])+
  
  labs(y=paste("MDS2",e2,"%",sep = " "),
       x=paste("MDS1",e1,"%",sep = " "),title = "PCoA 2")+
  my_theme      

#Based on the coversation with Jeff the strong differences in groups may be an artifact caused
#by different methods in the assembly of genome. A ways to subset the genomes that have been 
#exclusively assembled by JGI.

#NMDS did not work



#Here I getting parameters to get polygons that represent each group
grupos<-
  names(
    table(funct_groups$group[-grep("-",funct_groups$group)])[
      table(funct_groups$group[-grep("-",funct_groups$group)])>3])

mod<-model_pcoa2
vectors<-as.data.frame(scores(mod,display = "sites",scaling = "sites",choices=c(1,2)))
vectors$group<-funct_groups$group;vectors<-vectors[,c(3,1,2)]
centroids<-aggregate(cbind(PC1,PC2)~group,data=vectors,mean)
# to create the lines from the centroids to each point we will put it in a format that ggplot can handle
seg.data<-cbind(vectors[,1:3],centroids[rep(1:nrow(centroids), as.data.frame(table(vectors$group))$Freq),2:3])
names(seg.data)<-c("group","v.PC1","v.PC2","PC1","PC2")
seg.data

trial<-lapply(
  split(seg.data,seg.data$group),function(x){
    x[,c(1:3)][chull(x[,c(2:3)]),]
  })

all.hull<-do.call("rbind",trial)
names(all.hull)[2:3]<-c("p.PC1","p.PC2")
vectors$taxon<-funct_groups$taxon
vectors$phylum<-taxonomy$phylum


viendo_grupos<-
  ggplot()+
  geom_polygon(aes(x=p.PC1,y=p.PC2,fill=group),
               alpha=0.5,linetype="dashed",
               data=all.hull %>% 
                 filter(group%in%grupos)
  )+#filter(!grepl("-",group))) +
  #scale_fill_manual(values = c25)+
  geom_point(aes(x=MDS1,y=MDS2),size=2,#colour="red",
             data=vectors %>% 
               filter(group%in%grupos))+
  scale_color_brewer(palette = "Set1")+
  
  # geom_text(aes(x=PC1,y=PC2,label=taxon,fontface="bold.italic"),size=5,
  #           data=vectors %>% filter(grepl("Arbuscular",group)) %>% filter(!grepl("-",group))) +
  labs(y=paste("MDS2",e2,"%",sep = " "),
       x=paste("MDS1",e1,"%",sep = " "))+
  my_theme



############################## Older plots/ Outdated analysis ########################

m<-"Arbuscular_Mycorrhizal"
uno<-
  ggplot()+
  geom_polygon(aes(x=p.PC1,y=p.PC2,fill=group),
               alpha=0.5,linetype="dashed",
               data=all.hull %>% 
                 filter(group==m))+
  geom_point(aes(x=PC1,y=PC2,col=group),size=2,#colour="red",
             data=vectors %>% 
               filter(group==m))+
  labs(y=paste("PC2",e2,"%",sep = " "),
       x=paste("PC1",e1,"%",sep = " "))+
  my_theme

m<-"E_Mycorrhiza"
dos<-
  ggplot()+
  geom_polygon(aes(x=p.PC1,y=p.PC2,fill=group),
               alpha=0.5,linetype="dashed",
               data=all.hull %>% 
                 filter(group==m))+
  geom_point(aes(x=PC1,y=PC2,col=group),size=2,#colour="red",
             data=vectors %>% 
               filter(group==m))+
  labs(y=paste("PC2",e2,"%",sep = " "),
       x=paste("PC1",e1,"%",sep = " "))+
  my_theme

m<-"Endophyte"
tres<-
  ggplot()+
  geom_polygon(aes(x=p.PC1,y=p.PC2,fill=group),
               alpha=0.5,linetype="dashed",
               data=all.hull %>% 
                 filter(group==m))+
  geom_point(aes(x=PC1,y=PC2,col=group),size=2,#colour="red",
             data=vectors %>% 
               filter(group==m))+
  labs(y=paste("PC2",e2,"%",sep = " "),
       x=paste("PC1",e1,"%",sep = " "))+
  my_theme

m<-"Fungal_Parasite"
cuatro<-
  ggplot()+
  geom_polygon(aes(x=p.PC1,y=p.PC2,fill=group),
               alpha=0.5,linetype="dashed",
               data=all.hull %>% 
                 filter(group==m))+
  geom_point(aes(x=PC1,y=PC2,col=group),size=2,#colour="red",
             data=vectors %>% 
               filter(group==m))+
  labs(y=paste("PC2",e2,"%",sep = " "),
       x=paste("PC1",e1,"%",sep = " "))+
  my_theme

m<-"Human_Pathogen"
cinco<-
  ggplot()+
  geom_polygon(aes(x=p.PC1,y=p.PC2,fill=group),
               alpha=0.5,linetype="dashed",
               data=all.hull %>% 
                 filter(group==m))+
  geom_point(aes(x=PC1,y=PC2,col=group),size=2,#colour="red",
             data=vectors %>% 
               filter(group==m))+
  labs(y=paste("PC2",e2,"%",sep = " "),
       x=paste("PC1",e1,"%",sep = " "))+
  my_theme

m<-"Insect_Pathogen"
seis<-
  ggplot()+
  geom_polygon(aes(x=p.PC1,y=p.PC2,fill=group),
               alpha=0.5,linetype="dashed",
               data=all.hull %>% 
                 filter(group==m))+
  geom_point(aes(x=PC1,y=PC2,col=group),size=2,#colour="red",
             data=vectors %>% 
               filter(group==m))+
  labs(y=paste("PC2",e2,"%",sep = " "),
       x=paste("PC1",e1,"%",sep = " "))+
  my_theme

m<-"Insect_Symbiont"
siete<-
  ggplot()+
  geom_polygon(aes(x=p.PC1,y=p.PC2,fill=group),
               alpha=0.5,linetype="dashed",
               data=all.hull %>% 
                 filter(group==m))+
  geom_point(aes(x=PC1,y=PC2,col=group),size=2,#colour="red",
             data=vectors %>% 
               filter(group==m))+
  labs(y=paste("PC2",e2,"%",sep = " "),
       x=paste("PC1",e1,"%",sep = " "))+
  my_theme

m<-"Invertebrate_pathogen"
ocho<-
  ggplot()+
  geom_polygon(aes(x=p.PC1,y=p.PC2,fill=group),
               alpha=0.5,linetype="dashed",
               data=all.hull %>% 
                 filter(group==m))+
  geom_point(aes(x=PC1,y=PC2,col=group),size=2,#colour="red",
             data=vectors %>% 
               filter(group==m))+
  labs(y=paste("PC2",e2,"%",sep = " "),
       x=paste("PC1",e1,"%",sep = " "))+
  my_theme

m<-"Lichenized"
nueve<-
  ggplot()+
  geom_polygon(aes(x=p.PC1,y=p.PC2,fill=group),
               alpha=0.5,linetype="dashed",
               data=all.hull %>% 
                 filter(group==m))+
  geom_point(aes(x=PC1,y=PC2,col=group),size=2,#colour="red",
             data=vectors %>% 
               filter(group==m))+
  labs(y=paste("PC2",e2,"%",sep = " "),
       x=paste("PC1",e1,"%",sep = " "))+
  my_theme

m<-"Plant_Pathogen"
diez<-
  ggplot()+
  geom_polygon(aes(x=p.PC1,y=p.PC2,fill=group),
               alpha=0.5,linetype="dashed",
               data=all.hull %>% 
                 filter(group==m))+
  geom_point(aes(x=PC1,y=PC2,col=group),size=2,#colour="red",
             data=vectors %>% 
               filter(group==m))+
  labs(y=paste("PC2",e2,"%",sep = " "),
       x=paste("PC1",e1,"%",sep = " "))+
  my_theme

m<-"Saprotroph"
once<-
  ggplot()+
  geom_polygon(aes(x=p.PC1,y=p.PC2,fill=group),
               alpha=0.5,linetype="dashed",
               data=all.hull %>% 
                 filter(group==m))+
  geom_point(aes(x=PC1,y=PC2,col=group),size=2,#colour="red",
             data=vectors %>% 
               filter(group==m))+
  labs(y=paste("PC2",e2,"%",sep = " "),
       x=paste("PC1",e1,"%",sep = " "))+
  my_theme

m<-"warm_vertebrate_symbiont"
doce<-
  ggplot()+
  geom_polygon(aes(x=p.PC1,y=p.PC2,fill=group),
               alpha=0.5,linetype="dashed",
               data=all.hull %>% 
                 filter(group==m))+
  geom_point(aes(x=PC1,y=PC2,col=group),size=2,#colour="red",
             data=vectors %>% 
               filter(group==m))+
  labs(y=paste("PC2",e2,"%",sep = " "),
       x=paste("PC1",e1,"%",sep = " "))+
  my_theme


library(gridExtra)



#checking with families are the most important

e_vectors_traits_scaled_a<-as.data.frame(
  scores(model3,display = "species",choices = c(1,2),scaling = "sites"));
e_vectors_traits_scaled_a$traits<-rownames(e_vectors_traits_scaled_a)

names(e_vectors_traits_scaled_a)[1]<-"RDA1_species"
names(e_vectors_traits_scaled_a)[2]<-"RDA2_species"

e_vectors_traits_scaled_b<-left_join(e_vectors_traits_scaled_a,gene_descriptions %>% 
                                       rename(traits=family))

e_vectors_traits_scaled_b$RDA1_species<-abs(e_vectors_traits_scaled_b$RDA1_species)
e_vectors_traits_scaled_b$RDA2_species<-abs(e_vectors_traits_scaled_b$RDA2_species)

library(tidytext)

e_vectors_traits_scaled_a %>% 
  #mutate(traits=new_names) %>% 
  pivot_longer(RDA1_species:RDA2_species,names_to = "Axis",values_to="loading") %>% 
  ggplot()+
  aes(x=reorder_within(traits,-abs(loading),Axis),y=abs(loading))+
  geom_bar(stat="identity",position = "dodge")+
  facet_wrap(.~Axis,scales = "free_x")+
  labs(y="Loadings")+
  #scale_fill_manual(values = c("#7FCDBB","#0C2C84"))+
  scale_x_reordered()+
  theme(title = element_text(size = 18),
        axis.title.x=element_blank(),
        panel.background = element_blank(),
        panel.grid.major.y = element_line(size=0.25,linetype = "longdash",colour = "gray"),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 15,angle = 45,hjust = 1),
        #axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 25),
        legend.position = "none"
  )




plot(model2,display="sites")
ordihull(model2,funct_groups$group)


#Adding the two RDA axis to the dataframe
temporal<-as.data.frame(scores(model2,display = "sites",scaling = "sites",choices=c(1,2)))
funct_groups0<-cbind(funct_groups,temporal);rm(temporal)

e1<-round(summary(model2)[["cont"]][["importance"]][2,1],2)*100
e2<-round(summary(model2)[["cont"]][["importance"]][2,2],2)*100
#m<-2.5

funct_groups0%>% 
  ggplot()+
  aes(x=RDA1,y=RDA2,col=group) +
  geom_point(size=2)+
  #geom_text(aes(label=easy_code,fontface="bold.italic",color=phylum),size=5) +
  labs(y=paste("RDA2",e2,"%",sep = " "),
       x=paste("RDA1",e1,"%",sep = " "))+
  my_theme


c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  #"black", "gold1",
  "skyblue2", #"#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70",# "khaki2",
  #"maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", #"yellow4", "yellow3",
  "darkorange4", "brown"
)











#Adding the two RDA axis to the dataframe
temporal<-as.data.frame(scores(model3,display = "sites",scaling = "sites",choices=c(1,2)))

funct_groups1<-cbind(funct_groups,temporal);rm(temporal)


e1<-round(summary(model3)[["cont"]][["importance"]][2,1],2)*100
e2<-round(summary(model3)[["cont"]][["importance"]][2,2],2)*100
#m<-2.5


funct_groups1%>% 
  ggplot()+
  aes(x=PC1,y=PC2,col=group) +
  geom_point(size=2)+
  #geom_text(aes(label=easy_code,fontface="bold.italic",color=phylum),size=5) +
  labs(y=paste("PC2",e2,"%",sep = " "),
       x=paste("PC1",e1,"%",sep = " "))+
  my_theme



mod<-model2
centroids<-data.frame(grps=rownames(as.data.frame(scores(mod,display = "cn",scaling = "sites",choices=c(1,2)))),
                      as.data.frame(scores(mod,display = "cn",scaling = "sites",choices=c(1,2))))
vectors<-as.data.frame(scores(mod,display = "sites",scaling = "sites",choices=c(1,2)))
vectors$group<-funct_groups$group;vectors<-vectors[,c(3,1,2)]
# to create the lines from the centroids to each point we will put it in a format that ggplot can handle
seg.data<-cbind(vectors[,1:3],centroids[rep(1:nrow(centroids), as.data.frame(table(vectors$group))$Freq),2:3])
names(seg.data)<-c("group","v.RDA1","v.RDA2","RDA1","RDA2")
seg.data






#The idea is now to collapse all the as a mean for each functional group 
#(also when they are reported in more than one function). Then I will do 
#again, PCA´s or networks with collapsed data.


#To do that I need to use the format from functional_dat. First I need to subset it 
#to include only the fungi from which we have gene frequencies

functional_dat2<-left_join(geneFreq_dat_wFunct[,c("X","taxon")],functional_dat)


trial<-sapply(functional_dat2[,-c(1,2)],function(x){which(x==1)})

s<-sapply(trial,function(x){length(x)==0})
trial<-trial[!s]



functional_summary<-
  bind_rows(
    sapply(geneFreq_dat_wFunct[trial$Saprotroph,-c(1:5)],mean),
    sapply(geneFreq_dat_wFunct[trial$Plant_Pathogen,-c(1:5)],mean),
    sapply(geneFreq_dat_wFunct[trial$Insect_Pathogen,-c(1:5)],mean),
    sapply(geneFreq_dat_wFunct[trial$Human_Pathogen,-c(1:5)],mean),
    sapply(geneFreq_dat_wFunct[trial$Invertebrate_pathogen,-c(1:5)],mean),
    sapply(geneFreq_dat_wFunct[trial$cold_vertebrate_pathogen,-c(1:5)],mean),
    sapply(geneFreq_dat_wFunct[trial$warm_vertebrate_pathogen,-c(1:5)],mean),
    sapply(geneFreq_dat_wFunct[trial$Endophyte,-c(1:5)],mean),
    sapply(geneFreq_dat_wFunct[trial$Lichenized,-c(1:5)],mean),
    sapply(geneFreq_dat_wFunct[trial$Arbuscular_Mycorrhizal,-c(1:5)],mean),
    sapply(geneFreq_dat_wFunct[trial$E_Mycorrhiza,-c(1:5)],mean),
    sapply(geneFreq_dat_wFunct[trial$warm_vertebrate_symbiont,-c(1:5)],mean),
    sapply(geneFreq_dat_wFunct[trial$Insect_Symbiont,-c(1:5)],mean),
    sapply(geneFreq_dat_wFunct[trial$Fungal_Parasite,-c(1:5)],mean),
    sapply(geneFreq_dat_wFunct[trial$Ectomycorrhizal,-c(1:5)],mean),
    sapply(geneFreq_dat_wFunct[trial$Wood_Saprotroph,-c(1:5)],mean),
    sapply(geneFreq_dat_wFunct[trial$Dung_Saprotroph,-c(1:5)],mean),
    sapply(geneFreq_dat_wFunct[trial$Epiphyte,-c(1:5)],mean),
    sapply(geneFreq_dat_wFunct[trial$Lichen_Parasite,-c(1:5)],mean),
    sapply(geneFreq_dat_wFunct[trial$Ericoid_Mycorrhizal,-c(1:5)],mean),
    sapply(geneFreq_dat_wFunct[trial$Orchid_Mycorrhizal,-c(1:5)],mean)
  )

functional_summary$functional_group<-names(trial)

#THis confirms everything is right
all(
  sapply(geneFreq_dat_wFunct[trial$Orchid_Mycorrhizal,-c(1:5)],mean)==functional_summary[21,-8295])

#Now the RDAs

model4<-rda(functional_summary[,-8295],scale=T)

#Adding the two RDA axis to the dataframe
temporal<-as.data.frame(scores(model4,display = "sites",scaling = "sites",choices=c(1,2)))

funct_groups2<-cbind(functional_summary,temporal);rm(temporal)


e1<-round(summary(model4)[["cont"]][["importance"]][2,1],2)*100
e2<-round(summary(model4)[["cont"]][["importance"]][2,2],2)*100
#m<-2.5


funct_groups2%>% 
  ggplot()+
  aes(x=PC1,y=PC2) +
  geom_point(size=2)+
  geom_text(aes(label=functional_group,fontface="bold.italic"),size=3) +
  labs(y=paste("PC2",e2,"%",sep = " "),
       x=paste("PC1",e1,"%",sep = " "))+
  my_theme




#Heatmap
library(pheatmap)
library(RColorBrewer)

#https://stackoverflow.com/questions/41608294/r-heatmap-2-manual-grouping-of-rows-and-columns


funct_groups<-funct_groups[order(funct_groups$group),]
geneFreq_dat_wFunct<-geneFreq_dat[match(funct_groups$taxon,geneFreq_dat$taxon),]
rownames(funct_groups)<-NULL
rownames(geneFreq_dat_wFunct)<-NULL

#rownames(funct_groups)<-funct_groups$taxon
# funct_groups$X<-NULL
# funct_groups$taxon<-NULL


p<-data.frame(row.names = rownames(funct_groups)[1:20],group=funct_groups[c(1:20),])

o<-p
o$id<-c(1:length(o$group))
o<-split(o,o$group)
o<-sapply(o, function(x){x[length(x$id),2]})

pheatmap(geneFreq_dat[c(1:20),c(6:16)],cluster_cols = F, cluster_rows = F,
         scale = "column",annotation_row = p,gaps_row = o)

o<-funct_groups
o$id<-c(1:length(o$group))
o<-split(o,o$group)
o<-sapply(o, function(x){x[length(x$id),2]})

pheatmap(geneFreq_dat[,-c(1:5)],cluster_cols = F, cluster_rows = F,
         scale = "column",annotation_row = funct_groups,
         #gaps_row = o,
         show_rownames =F,
         show_colnames =F,
         legend = F)


pheatmap(geneFreq_dat_wFunct2[c(1:1000),c(6:1000)],cluster_cols = F, cluster_rows = F,
         scale = "column",annotation_row = as.data.frame(funct_groups$group[1:1000]),
         #gaps_row = o,
         show_rownames =F,
         show_colnames =F,
         legend = F)


pheatmap(functional_summary[,-8295],cluster_cols = F, cluster_rows = F,
         scale = "column",annotation_row = as.data.frame(functional_summary$functional_group),
         #gaps_row = o,
         show_rownames =F,
         show_colnames =F,
         legend = F)









#Network approach
dat<-geneFreq_dat[,-c(1,3:5)]
# dat<-geneFreq_dat[c(1:5),-c(1,3:5)]
# dat<-dat[,c(1:7)]

dat<-as.data.frame(t(dat))

names(dat)<-dat[1,]
dat<-dat[-1,]

g_names<-rownames(dat)

dat<-lapply(dat,as.numeric)
dat<-as.data.frame(dat)
dat<-as.matrix(dat)
dat[dat!=0]<-1

#There are 95 genes that are common to all fungi
length(which(rowSums(dat)==1294))

d1<-dat[-which(rowSums(dat)==1294),]

#d1<-sqrt(d1)

d1<-crossprod(d1)

#Cap of 3 to 5 species per genus
#PF Protein Family Database
#Doing ANOVA´s per single gene ()
#

#####

library(igraph)

functional_network<-graph_from_adjacency_matrix(d1,mode = "undirected",weighted = T)
functional_network<-igraph::simplify(functional_network,remove.loops = TRUE)

#This one is a copy to make calculations of centrality
#functional_network_calc<-functional_network

#E(functional_network)$pesos<-E(functional_network)$weight

#functional_network<-delete_edge_attr(functional_network,"weight")

V(functional_network)$size<-diag(dat)


###
#edge.color=adjustcolor(1, alpha.f = 0.15)
par(mar=c(0,0,0,0))
plot(functional_network,edge.curved=0,arrow.mode=0,
     edge.width=log10(E(functional_network)$weight),
     vertex.size=log10(V(functional_network)$size),
     vertex.label=NA,
     display.isolates=F,
     edge.color=adjustcolor(1, alpha.f = 0.15),
     #rescale=FALSE,
     #vertex.label.dist=1,
     #layout=layout_nicely
     #layout=layout_in_circle
     #layout=layout.fruchterman.reingold(functional_network)*10000
     #layout=layout.sphere
     layout=layout.fruchterman.reingold
     #layout=test.layout
     #layout=l
     #layout=layout_with_dh
     #layout = layout_with_graphopt
     #layout=layout_with_lgl
     #layout=layout_with_kk
)








