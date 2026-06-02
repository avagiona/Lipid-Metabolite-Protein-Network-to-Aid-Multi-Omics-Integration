library("dplyr")
library("igraph")
library("ggplot2")
library("NetHypGeom")


#set pathway
#setwd("~/R-4.2.1/multi_omics_embedding")
options(scipen = 999)

#import files 
hippie_v2.3_mapped<-read.csv("hippie_0.71_mapped_uniprot.csv", fileEncoding = "UTF-8-BOM", header=TRUE)
prot_lipids_metabolic<-read.csv("lipids_proteins_metabolic_uniprot.csv", fileEncoding = "UTF-8-BOM", header=TRUE)
prot_lipids_genetic<-read.csv("lipids_proteins_genetic_uniprot.csv", fileEncoding = "UTF-8-BOM", header=TRUE)
prot_metabol_interactions<-read.csv("protein_metabolites_reaction3.csv", fileEncoding = "UTF-8-BOM", header=TRUE)

#hippie keep the uniprot ids only
prot_prot_inter<-hippie_v2.3_mapped%>%
  dplyr::select(3,6)%>%
  dplyr::rename("V1"=1, "V2"=2)%>%  
  dplyr::filter(!(V1 == ""))%>%
  dplyr::filter(!(V2 == ""))

#convert the numeric into character 
prot_metabol_interactions$V2<-as.character(prot_metabol_interactions$V2)

multi_omics_data<-prot_prot_inter%>%
  dplyr::bind_rows(prot_lipids_metabolic)%>%
  dplyr::bind_rows(prot_lipids_genetic)%>%
  dplyr::bind_rows(prot_metabol_interactions)

#graphs from edgelists
net_multi_omics <- graph_from_data_frame(multi_omics_data, directed = FALSE)

#simplify to remove self and parallel interactions 
net_multi_omics_simplify<-simplify(net_multi_omics)

#parameter estimation netork 0.71
nodes_net_multi_omics<-data.frame(name = igraph::V(net_multi_omics_simplify)$name) #Nodes
edges_net_multi_omics<-igraph::E(net_multi_omics_simplify) #Edges
clusters_net_multi_omics<-decompose(net_multi_omics_simplify)#get the clusters
giant_net_multi_omics<-clusters_net_multi_omics[[ which.max(sapply(clusters_net_multi_omics, vcount)) ]] #get the LCC
#properties of giant cluster
nodes_giant_net_multi_omics<-data.frame(name = igraph::V(giant_net_multi_omics)$name) #Nodes LCC
edges_giant_net_multi_omics<-igraph::E(giant_net_multi_omics) #Edges LCC
edgeslist_giant_giant_net_multi_omics<-as_long_data_frame(giant_net_multi_omics)
#network properties of the network 
gma_net_multi_omics<-fit_power_law(degree(giant_net_multi_omics))$alpha #estimate_gma
c_net_multi_omics<-transitivity(giant_net_multi_omics,type = "localaverage") #estimate clustering coefficient
mean_degree_net_multi_omics<-mean(degree(giant_net_multi_omics))#estimate average node degree for artificial networks

#estimate Temperature
#psmodels_multi_omics=list()
#clustering_psmodels_multi_omics=list()
#mean_clustering_models_multi_omics=NULL

#for (i in 1:10){
#  psmodels_multi_omics[[i]]<-ps_model(N = 16625, avg.k = 22.72, gma = 2.98, Temp = 0)
#  clustering_psmodels_multi_omics[[i]]<-transitivity(psmodels_multi_omics[[i]][["network"]],type = "localaverage")
#  mean_clustering_models_multi_omics<-(mean(clustering_psmodels_multi_omics[[i]]))
#}

# Some maths by hand here, based on SF3 (LaBNE+HM paper) 
# c=a*T+b
# when T=0, c=0.77
# when T=1, c=0 so the equation is: c=-0.77*T + 0.77
# when c=0.12 (clustering coefficient of the LCC), T=0.84

#apply LaBNE+HM 
coordinates_multi_omics<-labne_hm(giant_net_multi_omics,gma = 2.98, Temp = 0.84, w = 2*pi)
saveRDS(coordinates_multi_omics,file="coordinates_protein_lipids_new.RData")

#coor_multi_omics<-readRDS("coordinates_protein_lipids_new.RData")



