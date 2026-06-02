#REQUIREMENTS
library("NetHypGeom") 
library("igraph")
library("dplyr")
library("ggplot2")
library("ggpubr")
library("enrichR")
library("lattice") 
library("RColorBrewer")
library("tidyr")
library("writexl")
library("enrichR")
library("plyr")
library("tibble")
library("ggrepel")
library("colorspace")
library("org.Hs.eg.db")
library("stringr")

#set pathway
setwd("~/R-4.2.1/FU1_FU2_new_analysis")
options(scipen = 999)

#import file 
coor_multi_omics<-readRDS("coordinates_protein_lipids_new.RData")
hippie_v2.3_mapped<-read.csv("hippie_0.71_mapped_uniprot.csv", fileEncoding = "UTF-8-BOM", header=TRUE)
prot_lipids_metabolic<-read.csv("lipids_proteins_metabolic_uniprot.csv", fileEncoding = "UTF-8-BOM", header=TRUE)
prot_lipids_genetic<-read.csv("lipids_proteins_genetic_uniprot.csv", fileEncoding = "UTF-8-BOM", header=TRUE)
prot_metabol_interactions<-read.csv("protein_metabolites_reaction3.csv", fileEncoding = "UTF-8-BOM", header=TRUE)
nodes_category<-read.csv("nodes_multi_omics_category.csv", fileEncoding = "UTF-8-BOM", header=TRUE)
nodes_cluster<-read.csv("multi_omics_clusters.csv", fileEncoding = "UTF-8-BOM", header=TRUE)
edges_multi_omics_category<-read.csv("edges_multi_omics.csv", fileEncoding = "UTF-8-BOM", header=TRUE)
#FU1_FU2_signatures<-read.csv("FU1_FU2_protein_signature.csv", fileEncoding = "UTF-8-BOM", header=TRUE)
#nodes_all_features<-read.csv("nodes_multi_omics_gene_symbol.csv", fileEncoding = "UTF-8-BOM", header=TRUE)

setwd("~/R-4.2.1/tissue_multi_omics")


#import files
lipids_tissues<-read.csv("lipid_tissues_swisslipid3.csv", fileEncoding = "UTF-8-BOM", header=TRUE)
metabolite_tissues<-read.csv("metabolite_tissue_HMDB.csv", fileEncoding = "UTF-8-BOM", header=TRUE)
protein_tissue<-read.delim("filtered_rna_tissue_with_uniprot.tsv",sep="\t")
Supp_tab_1<-read.csv("Suppl_Table_S1_nodes.csv", fileEncoding = "UTF-8-BOM", header=TRUE)

setwd("C:/Users/avagiona/Desktop/multi_omics/multi_omics_paper/Supplementary_files_january_2025")
Supp_tab_1_last<-read.csv("Suppl_Table_S1.csv", fileEncoding = "UTF-8-BOM", header=TRUE)
Supp_tab_1_last_coord<-Supp_tab_1_last%>%
  dplyr::left_join(coor_multi_omics$polar,by="id")
library("data.table")
#fwrite(Supp_tab_1_last_coord,"Supp_tab_1_last_coord.csv")

#Lipids
lipids_tissues$Network = ifelse(lipids_tissues$Lipid.ID %in% nodes_category$id,'1','0')

lipids_tissues_counts <- lipids_tissues %>%
  dplyr::group_by(Tissue.Cell.name, Network) %>%
  dplyr::summarise(count = n(), .groups = 'drop') %>%
  ungroup()

#lipids_tissues_counts$Network <- factor(lipids_tissues_counts$Network, levels = c("1", "0"))
lipids_tissues_network<-ggplot(lipids_tissues_counts, aes(x = reorder(Tissue.Cell.name,-count), 
                                                          y=count, fill= Network)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(x = "Tissue", y = "Counts",title = "Lipids") +
  scale_fill_manual(values = c("1" = "#44AA99", "0" = "#F8766D")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
lipids_tissues_network

#Metabolites
metabolite_tissues$Network = ifelse(metabolite_tissues$PC_CID %in% nodes_category$id,'1','0')

metabolite_tissues_counts <- metabolite_tissues %>%
  dplyr::group_by(tissue_location,Network) %>%
  dplyr::summarise(count = n(), .groups = 'drop') %>%
  ungroup()

tissue_order <- metabolite_tissues_counts %>%
  dplyr::group_by(tissue_location) %>%
  dplyr::summarise(total_count = sum(count)) %>%
  dplyr::arrange(-total_count) %>%
  dplyr::pull(tissue_location)

metabolite_tissues_counts$tissue_location <- factor(metabolite_tissues_counts$tissue_location, levels = tissue_order)

#metabolite_tissues_counts$Network <- factor(metabolite_tissues_counts$Network, levels = c("1", "0"))
metabolite_tissues_network<-ggplot(metabolite_tissues_counts, 
                               aes(x =tissue_location, y=count, fill= Network)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(x = "Tissue", y = "Counts",title = "Metabolites") +
  scale_fill_manual(values = c("1" = "#44AA99", "0" = "#F8766D")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
metabolite_tissues_network

#SOS SOS SOS#
#analyze again proteins annotated with brain tissues
tissue_of_interest<-c("brain", "substantia nigra",
                      "caudate","amygdala",
                      "hypothalamus","nucleus accumbens")

protein_tissue_brain <- protein_tissue %>%
  dplyr::select(3,4)%>%
  dplyr::group_by(UniProt) %>%
  dplyr::summarise(all_in_tissues = all(tissue %in% tissue_of_interest)) %>%
  dplyr::filter(all_in_tissues) %>%
  dplyr::pull(UniProt)

proteins_all_brain <- protein_tissue %>%
  dplyr::filter(UniProt %in% protein_tissue_brain & tissue %in% tissue_of_interest) %>%
  dplyr::select(UniProt, tissue)


proteins_all_brain$Network = ifelse(proteins_all_brain$UniProt %in% nodes_category$id,'1','0')

protein_tissue_unique_counts <- proteins_all_brain %>%
  dplyr::group_by(tissue, Network) %>%
  dplyr::summarise(count = n(), .groups = 'drop') %>%
  ungroup()

#protein_tissue_unique_counts$Network <- factor(protein_tissue_unique_counts$Network, levels = c("1", "0"))
#protein_tissue_unique_counts_network<-ggplot(protein_tissue_unique_counts, 
#                                   aes(x = reorder(tissue, -count), y=count, fill= Network)) +
#  geom_bar(position = "stack", stat = "identity") +
#  labs(x = "Tissue", y = "Counts",title = "Proteins") +
#  scale_fill_manual(values = c("1" = "#44AA99", "0" = "#F8766D")) +
#  theme_bw() +
#  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
#protein_tissue_unique_counts_network


#final<-ggarrange(lipids_tissues_network, metabolite_tissues_network,
#                 protein_tissue_unique_counts_network,
#                 ncol = 2, nrow = 2,
#                labels = c("a", "b","c"))
#final

#ggexport(final,filename = "final.png",
#         width = 3500,height = 2000,
#         pointsize = 100,res=200)


#visualize network 
nodes_euclidean_coordinates<-coor_multi_omics$cartesian
nodes_polar_coordinates<-coor_multi_omics$polar

#lipids brain
lipids_brain<-lipids_tissues%>%
  dplyr::filter(Tissue.Cell.name == "brain")
lipids_brain_network<-lipids_brain%>%
  dplyr::left_join(nodes_euclidean_coordinates, by=c("Lipid.ID"="id"))%>%
  dplyr::filter(Network == 1)%>%
  dplyr::select(1,4,5)%>%
  dplyr::rename(id=1)%>%
  dplyr::mutate(category= "lipid")

#metabolites brain
metabolites_brain<-metabolite_tissues%>%
  dplyr::filter(tissue_location == "brain")
metabolites_brain$PC_CID <- as.character(metabolites_brain$PC_CID)
metabolites_brain_network<-metabolites_brain%>%
  dplyr::left_join(nodes_euclidean_coordinates, by=c("PC_CID"="id"))%>%
  dplyr::filter(Network == 1)%>%
  dplyr::select(2,5,6)%>%
  dplyr::rename(id=1)%>%
  dplyr::mutate(category="metabolite")

#proteins brain
#proteins_brain<-protein_tissue_unique%>%
#  dplyr::filter(tissue == "brain")
proteins_brain_network<-proteins_all_brain%>%
  dplyr::left_join(nodes_euclidean_coordinates, by=c("UniProt"="id"))%>%
  dplyr::filter(Network == 1)%>%
  dplyr::select(1,4,5)%>%
  dplyr::rename(id=1)%>%
  dplyr::mutate(category="protein")


nodes_brain<-rbind(lipids_brain_network,metabolites_brain_network,proteins_brain_network)
labels_brian<-rbind(lipids_brain_network,metabolites_brain_network,proteins_brain_network)
nodes_brain <- nodes_brain[order(nodes_brain$category,decreasing = TRUE), ]
labels_brian$id <- ifelse(grepl("^[0-9]", labels_brian$id), paste0("pubchem:", labels_brian$id), labels_brian$id)

labels_brian_2<-labels_brian%>%
  dplyr::distinct()
nodes_brain_cluster<-nodes_brain%>%
  dplyr::left_join(nodes_cluster,by="id")%>%
  dplyr::select(1,2,3,4,6)


#setwd("~/R-4.2.1/tissue_multi_omics")
#library(data.table)
#fwrite(nodes_brain_cluster,"nodes_all_brain_cluster.csv")

#there are not interactions between the elements of the network related to brain 
network_brain_<-edges_multi_omics_category[with(edges_multi_omics_category,id1 %in% nodes_brain$id & 
                                     id2 %in% nodes_brain$id ),]


library(data.table)
#fwrite(nodes_brain,"nodes_brain.csv")
brain_multi_omics_2<-ggplot(nodes_brain,aes(x=x,y=y))+
  theme_classic() + 
  theme(legend.position = "none"  )+
  geom_point(aes(colour = as.factor(category)), size = 0.5)+
  scale_color_manual(values = c("lipid" = "cyan3","protein" = "gray15", "metabolite"="red"))+
  geom_text_repel(data=labels_brian_2, aes(label= id,colour = as.factor(category)),
                  max.overlaps = 100,size=1.5, segment.size = 0.05,fontface = "bold")
brain_multi_omics_2
#setwd("/Users/avagiona/Desktop/Figure_1_multiomics")


setwd("~/R-4.2.1/tissue_multi_omics")
#nodes_pubchem<-read.csv("multi_omics_clusters_pubchem.csv", fileEncoding = "UTF-8-BOM", header=TRUE)
#add to the supplementary table 1 the clusters and the brain annotation 
Supp_tab_1_cluster<-Supp_tab_1 %>%
  dplyr::left_join(nodes_pubchem, by="id")%>%
  dplyr::select(1,2,4)%>%
  dplyr::rename("cluster"=3)

#nodes_brain_pubchem<-read.csv("nodes_brain_cluster.csv", fileEncoding = "UTF-8-BOM", header=TRUE)
nodes_brain_cluster$id <- ifelse(grepl("^[0-9]", nodes_brain_cluster$id), paste0("pubchem:", nodes_brain_cluster$id), nodes_brain_cluster$id)
Supp_tab_1_cluster$brain = ifelse(Supp_tab_1_cluster$id %in% nodes_brain_cluster$id,'1','0')
#fwrite(Supp_tab_1_cluster,"Supp_tab_1_cluster_brain_13_12.csv")

nodes_polar_coordinates$id <- ifelse(grepl("^[0-9]", nodes_polar_coordinates$id), 
                                     paste0("pubchem:", nodes_polar_coordinates$id), nodes_polar_coordinates$id)

setwd("C:/Users/avagiona/Desktop/alex_paper_2025_version")
Suppl_Table_S1_nodes<-Supp_tab_1_cluster%>%
  dplyr::left_join(nodes_polar_coordinates, by="id")
#fwrite(Suppl_Table_S1_nodes,"Suppl_Table_S1_nodes.csv")

Supp_tab_1_cluster_brain_only<-Supp_tab_1_cluster%>%
  dplyr::filter(brain == 1)

#calculate the number of proteins per cluster 
proteins_per_cluster<-Supp_tab_1_cluster%>%
  dplyr::filter(category=="protein")

# visualize distributions of polar coordinates 
nodes_brain_polar<-nodes_brain%>%
  dplyr::left_join(nodes_polar_coordinates)

#lipids
lipids_all_network<-lipids_tissues%>%
  dplyr::filter(Network ==1)%>%
  dplyr::mutate(category = "all")%>%
  dplyr::select(1,4)%>%
  dplyr::left_join(nodes_polar_coordinates, by=c("Lipid.ID"="id"))%>%
  dplyr::rename("id"=1)

lipids_brain_network_polar<-lipids_brain_network%>%
  dplyr::left_join(nodes_polar_coordinates)%>%
  dplyr::select(1,5,6)%>%
  dplyr::mutate(category="brain")%>%
  dplyr::select(1,4,2,3)

lipids_hist<-rbind(lipids_all_network,lipids_brain_network_polar)

hist_lipids_theta<-ggplot(lipids_hist, aes(x=theta,colour = category,fill=category)) + 
  geom_histogram(aes(y = ..density..),alpha=0.7,bins=10,position = "dodge")+ 
  theme_bw()+
  xlab("theta")+ylab("Density")+
  theme(plot.tag=element_text(face="bold"))+
  scale_color_manual(values=c("all"="grey", "brain"="#F8766D"))+
  scale_fill_manual(values=c("all"="grey", "brain"="#F8766D"))+
  scale_x_continuous(breaks=seq(0, 7, 0.5))+
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12))+
  ggtitle("theta distribution of lipids")
hist_lipids_theta

hist_lipids_r<-ggplot(lipids_hist, aes(x=r,colour = category,fill=category)) + 
  geom_histogram(aes(y = ..density..),alpha=0.7,bins=10,position = "dodge")+ 
  theme_bw()+
  xlab("r")+ylab("Density")+
  theme(plot.tag=element_text(face="bold"))+
  scale_color_manual(values=c("all"="grey", "brain"="#F8766D"))+
  scale_fill_manual(values=c("all"="grey", "brain"="#F8766D"))+
  #scale_x_continuous(breaks=seq(0.5))+
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12))+
  ggtitle("r distribution of lipids")
hist_lipids_r


#metabolites
metabolite_tissues$PC_CID<-as.character(metabolite_tissues$PC_CID)
metabolites_all_network<-metabolite_tissues%>%
  dplyr::filter(Network ==1)%>%
  dplyr::mutate(category = "all")%>%
  dplyr::select(2,5)%>%
  dplyr::left_join(nodes_polar_coordinates, by=c("PC_CID"="id"))%>%
  dplyr::rename("id"=1)

metabolites_brain_network_polar<-metabolites_brain_network%>%
  dplyr::left_join(nodes_polar_coordinates)%>%
  dplyr::select(1,5,6)%>%
  dplyr::mutate(category="brain")%>%
  dplyr::select(1,4,2,3)

metabolites_hist<-rbind(metabolites_all_network,metabolites_brain_network_polar)

hist_metabolites_theta<-ggplot(metabolites_hist, aes(x=theta,colour = category,fill=category)) + 
  geom_histogram(aes(y = ..density..),alpha=0.8,bins=10,position = "dodge")+ 
  theme_bw()+
  xlab("theta")+ylab("Density")+
  theme(plot.tag=element_text(face="bold"))+
  scale_color_manual(values=c("all"="grey", "brain"="#44AA99"))+
  scale_fill_manual(values=c("all"="grey", "brain"="#44AA99"))+
  scale_x_continuous(breaks=seq(0, 7, 0.5))+
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12))+
  ggtitle("theta distribution of metabolites")
hist_metabolites_theta

hist_metabolites_r<-ggplot(metabolites_hist, aes(x=r,colour = category,fill=category)) + 
  geom_histogram(aes(y = ..density..),alpha=0.8,bins=10,position = "dodge")+ 
  theme_bw()+
  xlab("r")+ylab("Density")+
  theme(plot.tag=element_text(face="bold"))+
  scale_color_manual(values=c("all"="grey", "brain"="#44AA99"))+
  scale_fill_manual(values=c("all"="grey", "brain"="#44AA99"))+
  #scale_x_continuous(breaks=seq(0, 7, 0.5))+
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12))+
  ggtitle("r distribution of metabolites")
hist_metabolites_r

#proteins
proteins_all_network<-nodes_category%>%
  dplyr::filter(category=="protein")%>%
  dplyr::left_join(nodes_polar_coordinates)%>%
  dplyr::select(1,3,4)%>%
  dplyr::mutate(category="all")%>%
  dplyr::select(1,4,2,3)

proteins_brain_network_polar<-proteins_brain_network%>%
  dplyr::left_join(nodes_polar_coordinates)%>%
  dplyr::select(1,5,6)%>%
  dplyr::mutate(category="brain")%>%
  dplyr::select(1,4,2,3)
  
proteins_hist<-rbind(proteins_all_network,proteins_brain_network_polar)

hist_proteins_theta<-ggplot(proteins_hist, aes(x=theta,colour = category,fill=category)) + 
  geom_histogram(aes(y = ..density..),alpha=0.8,bins=10,position = "dodge")+ 
  theme_bw()+
  xlab("theta")+ylab("Density")+
  theme(plot.tag=element_text(face="bold"))+
  scale_color_manual(values=c("all"="grey", "brain"="#DDCC77"))+
  scale_fill_manual(values=c("all"="grey", "brain"="#DDCC77"))+
  scale_x_continuous(breaks=seq(0, 7, 0.5))+
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12))+
  ggtitle("theta distribution of proteins")
hist_lipids_theta

hist_proteins_r<-ggplot(proteins_hist, aes(x=r,colour = category,fill=category)) + 
  geom_histogram(aes(y = ..density..),alpha=0.8,bins=10,position = "dodge")+ 
  theme_bw()+
  xlab("r")+ylab("Density")+
  theme(plot.tag=element_text(face="bold"))+
  scale_color_manual(values=c("all"="grey", "brain"="#DDCC77"))+
  scale_fill_manual(values=c("all"="grey", "brain"="#DDCC77"))+
  #scale_x_continuous(breaks=seq(0, 7, 0.5))+
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12))+
  ggtitle("r distribution of proteins")
hist_lipids_r


distibutions<-ggarrange(hist_lipids_r, hist_lipids_theta,
                        hist_metabolites_r, hist_metabolites_theta,
                        hist_proteins_r, hist_proteins_theta,
                 ncol = 2, nrow = 3,
                labels = c("a", "b","c","d", "e","f"))
distibutions

ggexport(distibutions,filename = "distibutions_multi_omics.png",
         width = 3500,height = 2500,
         pointsize = 100,res=250)

#try to change the plot accordingly to the reviews requests



