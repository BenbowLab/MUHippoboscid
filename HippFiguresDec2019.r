library(vegan)
library(MASS)
library(ggplot2)
library(plyr)
library(dplyr)
library(magrittr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(randomForest)
library(knitr)
library(ape)
library(ggpubr)
library(agricolae)
library(Rmisc)
library(multcompView)
set.seed(1256)

parse_taxonomy_silva_128 <- function(char.vec){
  # Use default to assign names to elements in case problem with greengenes prefix
  char.vec = parse_taxonomy_default(char.vec)
  # Check for unassigned taxa
  if (char.vec["Rank1"] == "Unassigned") {
    char.vec <- c(Rank1="D_0__Unassigned", Rank2="D_1__Unassigned", Rank3="D_2__Unassigned", Rank4="D_3__Unassigned",
                  Rank5="D_4__Unassigned", Rank6="D_5__Unassigned", Rank7="D_6__Unassigned")
  }
  # Define the meaning of each prefix according to SILVA taxonomy
  Tranks = c(D_0="Kingdom", D_1="Phylum", D_2="Class", D_3="Order", D_4="Family", D_5="Genus", D_6="Species")
  # Check for prefix using regexp, warn if there were none. trim indices, ti
  ti = grep("[[:alpha:]]\\_[[:digit:]]{1}\\_\\_", char.vec)
  if( length(ti) == 0L ){
    warning(
      "No silva prefixes were found. \n",
      "Consider using parse_taxonomy_delfault() instead if true for all OTUs. \n",
      "Dummy ranks may be included among taxonomic ranks now."
    )
    # Will want to return without further modifying char.vec
    taxvec = char.vec
    # Replace names of taxvec according to prefix, if any present...
  } else {
    # Format character vectors for Ambiguous taxa
    if( length(ti) < 7 ){
      for (key in names(char.vec)) {
        if ( char.vec[key] == "Ambiguous_taxa" ) {
          tax_no <- (as.numeric(substr(key, 5, 5)) - 1)
          char.vec[key] = sprintf("D_%s__Ambiguous_taxa", tax_no)
        }
      }
      # Reset the trimmed indicies if Ambiguous taxa
      ti = grep("[[:alpha:]]\\_[[:digit:]]{1}\\_\\_", char.vec)
    }
    # Remove prefix using sub-"" regexp, call result taxvec
    taxvec = gsub("[[:alpha:]]\\_[[:digit:]]{1}\\_\\_", "", char.vec)
    # Define the ranks that will be replaced
    repranks = Tranks[substr(char.vec[ti], 1, 3)]
    # Replace, being sure to avoid prefixes notK present in Tranks
    names(taxvec)[ti[!is.na(repranks)]] = repranks[!is.na(repranks)]
  }
  return(taxvec)
}



cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#000000","#CC79A7")
theme_set(theme_bw(base_size = 18)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
biom=import_biom("HippoWTax9.18.19.biom",parseFunction= parse_taxonomy_silva_128)
#tax_table(biom) <- tax_table(biom)[,-c(5:10,14)]#remove dummy ranks

metadata=read.table("HippoMetadataWDiversity9.18.19.tsv",header = TRUE)

tree=read_tree("tree9.18.19.nwk")

sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$id
physeq=merge_phyloseq(biom,sampdat,tree)


#FIgure 1A

GPr  = transform_sample_counts(physeq, function(x) x / sum(x) ) #transform samples based on relative abundance
PhylumAll=tax_glom(GPr, "Phylum")

PhylumLevel = filter_taxa(PhylumAll, function(x) mean(x) > 1e-4, TRUE) #filter out any taxa lower tha 0.1%
df <- psmelt(PhylumLevel)
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Phylum", "BSpecies"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)


Figure1A=ggplot(Trtdata, aes(x=BSpecies,y=mean))+geom_bar(aes(fill = Phylum),colour="black", stat="identity")+xlab("Bird Species")+ylab("Relative Abundance (> 0.01%)")+theme(axis.title.x=element_blank())+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+ scale_fill_manual(values=cbPalette)
Figure1A#+  scale_x_discrete(labels=c("0 hrs", "24 hrs", "48 hrs","72 hrs"))


#Figure 1B
GenusAll=tax_glom(GPr,"Genus")
GenusLevel = filter_taxa(GenusAll, function(x) mean(x) > 1e-3, TRUE) #filter out any taxa lower tha 1%

df <- psmelt(GenusLevel)
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Genus", "BSpecies"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
sapply(Trtdata, class)
summary(Trtdata$Genus)
levels(Trtdata$Genus)<-c("Arsenophonus","Exiguobacterium","Ignatzschineria","Koukoulia","Pectobacterium","Endosymbiont of \n P. canariensis","Tissierella")
(Trtdata)
Figure1B=ggplot(Trtdata, aes(x=BSpecies,y=mean))+geom_bar(aes(fill = Genus),colour="black", stat="identity")+xlab("Bird BSpecies")+ylab("Relative Abundance (> 0.1%)") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(axis.title.x=element_blank())+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+ scale_fill_manual(values=cbPalette)#+geom_text(aes(x=BSpecies, y=mean+se+10,label=vec))
Figure1B#+  scale_x_discrete(labels=c("0 hrs", "24 hrs", "48 hrs","72 hrs"))


#Figure 1C
av=kruskal.test( faith_pd~ BSpecies, data=metadata) #fig.width = 8, fig.height = 8 for line above to change size
posthoc<-compare_means(faith_pd ~ BSpecies, data = metadata, p.adjust.method = "fdr")



Hyphenated<-as.character(paste0(posthoc$group1,"-",posthoc$group2))
difference<-posthoc$p.format
names(difference)<-Hyphenated

Letters<-multcompLetters(difference)
stats<-summarySE(metadata,measurevar="faith_pd",groupvars=c("BSpecies",na.rm=TRUE))
Figure1C<-ggplot(stats,aes(x=BSpecies,y=faith_pd,fill=BSpecies))+geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=faith_pd-se,ymax=faith_pd+se),color="black")+
  geom_text(aes(x=BSpecies, y=faith_pd+se+0.3,label=Letters$Letters), position=position_dodge(width=0.9), size=4,color="black")+
  guides(fill=FALSE)+  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ scale_fill_manual(values=cbPalette)+xlab("")+
  ylab("Faith's PD")#+ 
Figure1C

#Figure1D
df <- psmelt(PhylumLevel)
df$Abundance=df$Abundance*100

PhylumAll3<-ddply(df, c("Phylum","Month","BSpecies"), summarise, #For Rel abu tabs
                  N    = length(Abundance),
                  mean = mean(Abundance),
                  sd   = sd(Abundance),
                  se   = sd / sqrt(N)
)
PhylumAll3$Month = factor(PhylumAll3$Month, levels = c("Sept","Oct","Nov")) #fixes x-axis labels
Figure1D=ggplot(PhylumAll3, aes(x=Month,y=mean))+geom_bar(aes(fill = Phylum),colour="black", stat="identity")+xlab("Month")+ylab("Relative Abundance (> 0.01%)")+theme(axis.title.x=element_blank())+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+facet_wrap(~BSpecies)+ scale_fill_manual(values=cbPalette)
Figure1D

ggarrange(Figure1A,Figure1B,Figure1C,Figure1D,
          labels = c("A", "B","C","D"),
          ncol = 2, nrow = 2)
#ggsave("name.eps", width = 20, height = 20, units = "cm")

theme_set(theme_bw(base_size = 7.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))


ggsave(filename = "HippFig1EPS.eps",
       (plot = ggarrange(Figure1A,Figure1B,Figure1C,Figure1D,
                        labels = c("A", "B  ","C","D"),
                        ncol = 2, nrow = 2)),
       
       device = "eps",dpi = 600,width = 5.2, height = 5.2,unit="in") #MAx width 5.2, height = 8.75 min width = 2.63 height = 

dev.off()
tiff("HippFigure1Tiff.tiff", width = 5.2, height = 5.2, units = 'in', res = 600)
ggarrange(Figure1A,Figure1B,Figure1C,Figure1D,
          labels = c("", "","",""),
          ncol = 2, nrow = 2)
dev.off()


#Figure2
ord=ordinate(physeq,"PCoA", "jaccard")
ordplot=plot_ordination(physeq, ord,"samples", color="BSpecies")+geom_point(size=2)+scale_colour_manual(values=cbPalette)+scale_fill_manual(values=cbPalette)+ theme(legend.title = element_blank())
Figure2<-ordplot+ stat_ellipse(type= "norm",geom = "polygon", alpha = 1/4, aes(fill = BSpecies))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())#+ 
Figure2


dev.off()
tiff("HippFigure2Tiff.tiff", width = 2.63, height = 2.63, units = 'in', res = 600)
Figure2
dev.off()

