## R script to plot phylo tree (Fig.2)
##### phylogenetic tree from the study: "Growth–reproduction trade-offs are common but changing in woody plants: a meta-analysis" 

## open packages

#### 
ipak <- function(pkg){
  # Identify packages not currently installed on your system
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  # If there are any new packages, install them (including dependencies)
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE) #return T or F if loaded and installed
}

packages <- c('tidyverse', 'cowplot',"patchwork", "confintr", "betareg", "broom", "mgcv", "ape", 'here',
              "phytools", "phylolm", 'sjPlot', "sjmisc", 'sjlabelled', "RColorBrewer", 'scatterpie', 'ggtree',
              "scales", "climwin", "readxl" , "gghalves", "rstatix", "ggpubr", "ggtree","ggtreeExtra","ggnewscale")

ipak(packages)

### should print all true

## load packages
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)

################################################################################
## final plots with the data prepared in R script (Growth-reproduction_trade-offs_meta_MKB2026.Rmd)

## updated phylogeny from Zanne et al. 2016
new_phylo_updated

## open tree of life phylogeny
my_phylo_branched

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

## copy the file with effect_sizes and species names
ring = meta_Z

ring$label <- sub("_", " ", ring$tip_label) #remove underscore
## annotate summary of data 
#ring$tip_label_an <- paste0(ring$label," (",ring$pop,",",ring$meanY,")")

#tip_label_an <- sub("_", " ", tip_label_an) #remove underscore
#ring$tip_label_an <- sub("_", " ", tip_label_an) #remove underscore

### sum of years
ring$tip_label_an <- paste0(meta_Z$tip_label," (",meta_Z$pop,",",meta_Z$sumY,")")

### mean number of years
ring$tip_label_an2 <- paste0(meta_Z$tip_label," (",meta_Z$pop,",",as.numeric(meta_Z$meanY),")") 

## annotate summary of data 
#tip_label_an <- paste0(meta_Z$tip_label," (",meta_Z$pop,",",meta_Z$sumY,")")
#meta_Z$tip_label_an2 <- paste0(meta_Z$tip_label," (",meta_Z$pop,",",as.numeric(meta_Z$meanY),")") 

str(ring)

## shorten the dataset for analysis (phylo function)

meta_Z <- meta_Z %>%
  dplyr::select(tip_label,meanZ)


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

## function for phylogenetical analysis (see description in the bottom)

phyloAnalysisFun<-function(phylogenetictree = new_phylo_updated, 
                           subsetphyloanalysisdf = meta_Z, 
                           variableEvo = meanZ){
  dfPhylo<- subsetphyloanalysisdf[!is.na(subsetphyloanalysisdf[variableEvo]),]
  Drop<- phylogenetictree$tip.label[-c(which(phylogenetictree$tip.label %in% dfPhylo$tip_label))]
  phylo_up<- drop.tip(phylogenetictree,Drop)
  dfPhylo<- dfPhylo[match(phylo_up$tip.label, dfPhylo$tip_label),] %>% as.data.frame()
  
  rownames(dfPhylo) <- dfPhylo$tip_label
  
  dfPhylo$tip_label<- NULL
  X <- as.matrix(dfPhylo)
  phyAnalyse <- phylosig(phylo_up, X[,1], method="lambda", test=TRUE, nsim=10000, se=NULL, start=NULL,
                         control=list())
  
  return(phyAnalyse)
}

#check lambda phylo signal 

result_phylo <- phyloAnalysisFun(
  phylogenetictree     = new_phylo_updated,
  subsetphyloanalysisdf= meta_Z,
  variableEvo          = "meanZ"
)

print(result_phylo)

## no signal as lambda is very close to 0

## check Blomberg's K phylo signal

signal_data <- meta_Z %>%
  dplyr::select(tip_label, meanZ) %>%
  filter(
    !is.na(tip_label),
    !is.na(meanZ)
  ) %>%
  distinct(tip_label, .keep_all = TRUE)

#Check whether the species names agree (if so they equal 0)
  
setdiff(signal_data$tip_label, new_phylo_updated$tip.label)


#Now retain only species present in both objects:
common_species <- intersect(
  new_phylo_updated$tip.label,
  signal_data$tip_label
)

tree_K <- keep.tip(
  new_phylo_updated,
  common_species
)


## 
# General checks of the tree 
is.rooted(tree_K)
is.binary(tree_K)
is.ultrametric(tree_K)
is.null(tree_K$edge.length)

# Problematic labels
anyDuplicated(tree_K$tip.label)
tree_K$tip.label[duplicated(tree_K$tip.label)]

# Problematic branch lengths
which(tree_K$edge.length == 0)
which(tree_K$edge.length < 0)
range(tree_K$edge.length)

# Tree covariance matrix
C <- vcv.phylo(tree_K)

dim(C)
qr(C)$rank
nrow(C)
kappa(C)


signal_data_K <- signal_data %>%
  filter(tip_label %in% common_species)

#Create the named vector:
  
  Z_vector <- signal_data_K$meanZ
names(Z_vector) <- signal_data_K$tip_label

#Reorder it to match the tree exactly (they return TRUE if correct)
  
Z_vector <- Z_vector[tree_K$tip.label]

identical(names(Z_vector), tree_K$tip.label)


# now run Blomberg’s K:
  
 set.seed(123)

K_result <- phylosig(
  tree = tree_K,
  x = Z_vector,
  method = "K",
  test = TRUE,
  nsim = 9999
)

K_result

K_summary <- data.frame(
  trait = "Z",
  n_species = length(Z_vector),
  Blomberg_K = K_result$K,
  P_value = K_result$P
)

K_summary

# Interpretation: K≈1: similarity is approximately consistent with Brownian-motion evolution;
# K<1: close relatives are less similar than expected under Brownian motion;
# K>1: close relatives are more similar than expected under Brownian motion;
# a small permutation P-value indicates evidence that the trait is non-randomly distributed across the phylogeny.
 
## check also for open tree of life phylogeny

setdiff(signal_data$tip_label, my_phylo_branched$tip.label)

common_species <- intersect(
  my_phylo_branched$tip.label,
  signal_data$tip_label
)

tree_K <- keep.tip(
  my_phylo_branched,
  common_species
)

###  plot tree

plotphyloGGtree <- function(phylogenetictree = new_phylo_updated, 
                            subsetphyloanalysisdf = meta_Z, 
                            variableEvo = 'meanZ', 
                            offset = 50){
  #basic merging dataset phylo and mastree+
  dfPhylo <- subsetphyloanalysisdf[!is.na(subsetphyloanalysisdf[variableEvo]),]
  Drop<- new_phylo_updated$tip.label[-c(which(new_phylo_updated$tip.label %in% dfPhylo$tip_label))]
  phylo_up<- drop.tip(new_phylo_updated,Drop)
  dfPhylo<- dfPhylo[match(phylo_up$tip.label, dfPhylo$tip_label),] %>% as.data.frame()
  
  rownames(dfPhylo) <- dfPhylo$tip_label
  dfPhylo$tip_label<- NULL
  svl <- as.matrix(dfPhylo)[,1]
  
  fit <- phytools::fastAnc(phylo_up,svl,vars=TRUE,CI=TRUE)
  td <- data.frame(node = nodeid(phylo_up, names(svl)),
                   var = svl)
  nd <- data.frame(node = names(fit$ace), var = fit$ace)
  
  d <- rbind(td, nd)
  d$node <- as.numeric(d$node)
  tree <- full_join(phylo_up, d, by = 'node')
  tree@phylo$tip.label <- sub("_", " ", tree@phylo$tip.label) #remove underscore
  

#color gradient 
  #  pal <- rev(wesanderson::wes_palette(name="Rushmore", type = "continuous", n = 42) )

    pal <- RColorBrewer::brewer.pal(n = 11, name = "RdBu")
  variableEvolegend <- "Mean Z"
  mid <- median(tree@extraInfo$var)
  palo <- scales::gradient_n_pal(scales::brewer_pal(palette = "RdBu", direction = -1)(11))
  # palo <- scales::gradient_n_pal(wesanderson::wes_palette(name="Rushmore",direction = -1) (11) )
    
  library(RColorBrewer)
  library(scales) # needed for rescale
  
cols <- brewer.pal(n = 11, name = "RdBu") 
#   cols <- rev(wesanderson::wes_palette(n = 42, type = "continuous", name = "Rushmore") )

  offset = offset
  phyloAutoTemp <- ggtree(tree, layout="fan", 
                          open.angle=10, 
                          aes(col=var)) +
    geom_tree(size = 1) +   geom_fruit(geom = geom_tile,width = 8,
                               aes(fill=var),
                              width = 0.02 ) +
    scale_fill_gradientn(colours = cols, 
                         values = scales::rescale(c(-1.5, -0.5, 0, 0.5, 1)),
                         guide = "colorbar", 
                         limits=c(-1.5, 1),
                         name = paste0(variableEvolegend))  +
    scale_colour_gradientn(colours = cols, 
                           values = scales::rescale(c(-1.5, -0.5, 0, 0.5, 1)),
                            guide = "colorbar", 
                            limits=c(-1.5, 1),
                           name = paste0(variableEvolegend)) +
    theme(legend.position="right", legend.key.height = unit(8, "mm"))   
   # +
  ### geom_cladelab(node=ggtree::MRCA(tree, c('Pinus sylvestris')), label="Pinus", barsize=1.2,
  ###              col="black", offset.text=5, angle = "auto",hjust=0.5, horizontal = F, offset = offset)+
  ### geom_cladelab(node=ggtree::MRCA(tree, c('Fagus sylvatica')), label="Fagus", barsize=1.2)
  ### col="black", offset.text=5, angle = "auto",hjust=0.5, horizontal = F, offset = offset)
  
  phyloAutoTemp
  
  return(phyloAutoTemp)
}


#------------------
#### plot the tree

my_plot <- plotphyloGGtree(
  phylogenetictree       = new_phylo_updated,
  subsetphyloanalysisdf  = meta_Z,
  variableEvo            = "meanZ",
  offset                 = 50
)

print(my_plot)

### simple tree
my_plot + geom_tiplab(size= 3, aes(angle=angle), colour="black", hjust = -0.3) + theme(
  legend.position = "none")  


my_plot + geom_tiplab(size= 3, aes(angle=angle), colour="black", hjust = -0.3) 

### 

my_plot    +  geom_tiplab(hjust = -.1, size = 5)   + 
  theme(legend.position = c(.02, .85))     


plot_phylo  <-  my_plot    +  geom_tiplab(hjust = -.1, size = 5, aes(angle=angle), colour="black", fontface = 1, offset = 6) +   
  theme(legend.position = c(.02, .85),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14))   +  #xlim(0, 1.3)  
  xlim(0,450)

## base plot produced
plot_phylo 
  
### produce more advanced plot 

## prepare pallettes 
pal <- RColorBrewer::brewer.pal(n = 11, name = "RdBu")
pal2 <- RColorBrewer::brewer.pal(n=9,name="Greys") 
#pal2 <- RColorBrewer::brewer.pal(n=9,name="Blues") 


p1leg <- my_plot + geom_tree(linewidth = 1) +
  new_scale_fill() +
  geom_fruit(
    data = ring,
    geom = geom_tile,
    mapping = aes(y = label, fill = pop, col = as.numeric(meanY)),
    width = 6,      # thickness of the ring
    offset = 0.05 ,)  +  # gap between tips and ring ,
   scale_fill_gradientn(colours = pal2, name = "N, n") +
  geom_fruit(
    data = ring,
    geom = geom_tile,
    mapping = aes(y = label, fill = as.numeric(meanY), col = as.numeric(meanY)),
    offset = 0.04 ,   # gap between tips and ring 
    width = 6,      # thickness of the ring   # gap between tips and ring
  ) +
  geom_tiplab(aes(label = label), align = FALSE, linetype = NA,
              hjust = -.1, size = 6.5, aes(angle=angle),  colour = "black", alpha = 0.8,
              fontface = 1, offset = 40) +
  #scale_fill_viridis_c(name = "n", option = "C", na.value = "grey90") +
  scale_fill_gradientn(colours = pal2, name = "N, n") +
  theme(legend.position ="right" )  
 
#geom_treescale(x = 250, y =0, width = 100, offset = 0.1,linesize = 0.3, fontsize = 5) ## adding time scale

p1leg  

### set margins

p1 <- p1leg + theme(legend.position ="none" ) + theme(plot.margin = margin(175, 155, 155, 155))

p1

### save the plot 

ggsave(
  filename = "phylo_plot_ring.pdf",
  plot     = p1,
  #  device   = cairo_pdf,   # needs Cairo installed; on Windows it ships with R, on macOS use XQuartz if needed
  width    = 13,
  height   = 13,
  units    = "in"
)
 
################################################################################
## save the plot's legend

leg_phylo <- ggpubr::get_legend(p1leg)

ggsave(
  filename = "leg_phylo.pdf",
  plot     = leg_phylo,
  #  device   = cairo_pdf,   # needs Cairo installed; on Windows it ships with R, on macOS use XQuartz if needed
  width    = 14,
  height   = 12,
  units    = "in"
)

################################################################################ 
################################################################################
###### DENDROGRAM - supplementary plot
## only secondary growth ## only secondary growth ## only secondary growth

ring2 = meta_Z_sec
ring2$label <- sub("_", " ", ring2$tip_label) #remove underscore
## annotate summary of data 
ring2$tip_label_an <- paste0(ring2$label," (",ring2$pop,",",as.numeric(ring2$meanY),")") 

##

meta_Z_sup <- ring2 %>%
  dplyr::select(tip_label_an,meanZ)

## only secondary growth ## only secondary growth ## only secondary growth
meta_Z_sec

##### plot for supplement - dendrogram with annotations and secondary growth data only


plotphyloGGtree3 <- function(phylogenetictree = new_phylo_updated, 
                             subsetphyloanalysisdf = meta_Z_sec, 
                             variableEvo = 'meanZ', 
                             offset = 70){
  #basic merging dataset phylo and mastree+
  dfPhylo <- subsetphyloanalysisdf[!is.na(subsetphyloanalysisdf[variableEvo]),]
  Drop<- new_phylo_updated$tip.label[-c(which(new_phylo_updated$tip.label %in% dfPhylo$tip_label))]
  phylo_up<- drop.tip(new_phylo_updated,Drop)
  dfPhylo<- dfPhylo[match(phylo_up$tip.label, dfPhylo$tip_label),] %>% as.data.frame()
  
  rownames(dfPhylo) <- dfPhylo$tip_label
  dfPhylo$tip_label<- NULL
  svl <- as.matrix(dfPhylo)[,1]
  
  fit <- phytools::fastAnc(phylo_up,svl,vars=TRUE,CI=TRUE)
  td <- data.frame(node = nodeid(phylo_up, names(svl)),
                   var = svl)
  nd <- data.frame(node = names(fit$ace), var = fit$ace)
  
  d <- rbind(td, nd)
  d$node <- as.numeric(d$node)
  tree <- full_join(phylo_up, d, by = 'node')
  tree@phylo$tip.label <- sub("_", " ", tree@phylo$tip.label) #remove underscore
  
  
  #color gradient 
  #pal <- rev(wesanderson::wes_palette(n = 42, type = "continuous", name = "Cavalcanti1") )
  pal <- RColorBrewer::brewer.pal(n = 9, name = "RdBu")
  variableEvolegend <- "Mean Z"
  mid <- median(tree@extraInfo$var)
  palo <- scales::gradient_n_pal(scales::brewer_pal(palette = "RdBu", direction = -1)(11))
 # palo <- scales::gradient_n_pal(rev(wesanderson::wes_palette(n = 42, type = "continuous", name = "Cavalcanti1") ))
  
  library(RColorBrewer)
  library(ggtreeExtra)
  library(scales) # needed for rescale
  
  cols <- brewer.pal(n = 9, name = "RdBu") 
 # cols <- rev(wesanderson::wes_palette(n = 42, type = "continuous", name = "Cavalcanti1") )
  
  offset = offset
  phyloAutoTemp <- ggtree(tree, layout="rectangular", 
                          open.angle=10 ) +
    geom_tree(size = 1) + # geom_fruit(geom = geom_tile,
                           #           aes(fill=var),
                           #           width = 0.02 ) +
    scale_fill_gradientn(colours = cols, 
                         values = scales::rescale(c(-1.5, -0.5, 0, 0.5, 1)),
                         guide = "colorbar", 
                         limits=c(-1.5, 1),
                         name = paste0(variableEvolegend)
    ) +
    theme(legend.position="right", legend.key.height = unit(8, "mm"))   
  # +
  ### geom_cladelab(node=ggtree::MRCA(tree, c('Pinus sylvestris')), label="Pinus", barsize=1.2,
  ###              col="black", offset.text=5, angle = "auto",hjust=0.5, horizontal = F, offset = offset)+
  ### geom_cladelab(node=ggtree::MRCA(tree, c('Fagus sylvatica')), label="Fagus", barsize=1.2)
  ### col="black", offset.text=5, angle = "auto",hjust=0.5, horizontal = F, offset = offset)
  
  phyloAutoTemp
  
  return(phyloAutoTemp)
  
}

## phylogeny with only secondary growth 

my_new_plot_dend <- plotphyloGGtree3(
  phylogenetictree       = new_phylo_updated,
  subsetphyloanalysisdf  = meta_Z_sec,
  variableEvo            = "meanZ",
  offset                 = 50
)

print(my_new_plot_dend)

my_new_plot_dend + geom_tree(linewidth = 1.2, colour = "black") +
   geom_fruit(
    data = ring2,
    geom = geom_tile,
    mapping = aes(y = label, fill = meanZ),
    width = 3, 
    offset = 0      # thickness of the ring
  ) +
 geom_tiplab(aes(label = label), align = TRUE, linetype = NA,
              hjust = -0.1, size = 4, aes(angle=angle), colour="black", fontface = 1, offset = 3) +
theme(legend.position ="left")

cols <- brewer.pal(n = 11, name = "RdBu") 
cols <- brewer.pal(n = 11, name = "YlGn")  

######
old = ring2$tip_label
new = ring2$tip_label_an

syn_lookup <- as.tibble(old)
syn_lookup$new <- new
syn_lookup$value <- sub("_", " ", old) #remove underscore

dat <- my_new_plot_dend$data %>%
  filter(isTip) %>%
  dplyr::select(label) %>%
  left_join(syn_lookup, by = c("label" = "value")) %>%
  mutate(label_new = coalesce(new, label))  # fallback to original if no match

##################

suppl_phylo <- my_new_plot_dend %<+% dat +
  geom_tree() +
  geom_fruit(
    data = ring2,
    geom = geom_tile,
    mapping = aes(y = label, fill = meanZ),
    width = 6, 
    offset = 0.017      # thickness of the ring
  ) +
  geom_tiplab(aes(label = label_new), align = FALSE, linetype = NA,
              hjust = -0.1, size = 5, aes(angle=angle), fontface = 3 , offset = 4)   +    xlim(0, 460)  +
 theme(aspect.ratio = (1.2)) + theme_tree2() + theme(
   legend.position = c(.25, .75),
   legend.justification = c("right", "top"),
   legend.box.just = "right",
   legend.margin = margin(6, 6, 6, 6),
   legend.title = element_text(size = 18),
   legend.text = element_text(size = 16))    

suppl_phylo

### SAVE THE FIG.S7

ggsave(
  filename = "plot_suppl_phylogeny_annot.pdf",
  plot     = suppl_phylo,
  #device   = cairo_pdf,   # needs Cairo installed; on Windows it ships with R, on macOS use XQuartz if needed
  width    = 14,
  height   = 15,
  units    = "in"
)


suppl_phylo2 <- my_new_plot_dend + geom_tree(linewidth = 1.2, colour = "black") +
  geom_fruit(
    data = ring2,
    geom = geom_tile,
    mapping = aes(y = label, fill = meanZ),
    width = 5, 
    offset = 0      # thickness of the ring
  ) +
  geom_tiplab(aes(label = label), align = TRUE, linetype = NA,
              hjust = -0.1, size = 6, aes(angle=angle), colour="black", fontface = 1, offset = 5) +
  xlim(0,600)  +
  theme(legend.position = "left",
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16)) + theme(aspect.ratio = (1.2))

suppl_phylo2


###### DENDROGRAM - supplementary plot but with ROTL tree
## only secondary growth ## only secondary growth ## only secondary growth
meta_Z_sec

##### plot for supplement - dendrogram with annotations and secondary growth data only


plotphyloGGtree_rotl <- function(phylogenetictree = my_phylo_branched, 
                             subsetphyloanalysisdf = meta_Z_sec, 
                             variableEvo = 'meanZ', 
                             offset = 70){
  #basic merging dataset phylo and mastree+
  dfPhylo <- subsetphyloanalysisdf[!is.na(subsetphyloanalysisdf[variableEvo]),]
  Drop<- my_phylo_branched$tip.label[-c(which(my_phylo_branched$tip.label %in% dfPhylo$tip_label))]
  phylo_up<- drop.tip(my_phylo_branched,Drop)
  dfPhylo<- dfPhylo[match(phylo_up$tip.label, dfPhylo$tip_label),] %>% as.data.frame()
  
  rownames(dfPhylo) <- dfPhylo$tip_label
  dfPhylo$tip_label<- NULL
  svl <- as.matrix(dfPhylo)[,1]
  
  fit <- phytools::fastAnc(phylo_up,svl,vars=TRUE,CI=TRUE)
  td <- data.frame(node = nodeid(phylo_up, names(svl)),
                   var = svl)
  nd <- data.frame(node = names(fit$ace), var = fit$ace)
  
  d <- rbind(td, nd)
  d$node <- as.numeric(d$node)
  tree <- full_join(phylo_up, d, by = 'node')
  tree@phylo$tip.label <- sub("_", " ", tree@phylo$tip.label) #remove underscore
  
  
  #color gradient 
  #pal <- rev(wesanderson::wes_palette(n = 42, type = "continuous", name = "Cavalcanti1") )
  pal <- RColorBrewer::brewer.pal(n = 9, name = "RdBu")
  variableEvolegend <- "Mean Z"
  mid <- median(tree@extraInfo$var)
  palo <- scales::gradient_n_pal(scales::brewer_pal(palette = "RdBu", direction = -1)(11))
  # palo <- scales::gradient_n_pal(rev(wesanderson::wes_palette(n = 42, type = "continuous", name = "Cavalcanti1") ))
  
  library(RColorBrewer)
  library(ggtreeExtra)
  library(scales) # needed for rescale
  
  cols <- brewer.pal(n = 9, name = "RdBu") 
  # cols <- rev(wesanderson::wes_palette(n = 42, type = "continuous", name = "Cavalcanti1") )
  
  offset = offset
  phyloAutoTemp <- ggtree(tree, layout="rectangular", 
                          open.angle=10 ) +
    geom_tree(size = 1) + # geom_fruit(geom = geom_tile,
    #           aes(fill=var),
    #           width = 0.02 ) +
    scale_fill_gradientn(colours = cols, 
                         values = scales::rescale(c(-1.5, -0.5, 0, 0.5, 1)),
                         guide = "colorbar", 
                         limits=c(-1.5, 1),
                         name = paste0(variableEvolegend)
    ) +
    theme(legend.position="right", legend.key.height = unit(8, "mm"))   
  # +
  ### geom_cladelab(node=ggtree::MRCA(tree, c('Pinus sylvestris')), label="Pinus", barsize=1.2,
  ###              col="black", offset.text=5, angle = "auto",hjust=0.5, horizontal = F, offset = offset)+
  ### geom_cladelab(node=ggtree::MRCA(tree, c('Fagus sylvatica')), label="Fagus", barsize=1.2)
  ### col="black", offset.text=5, angle = "auto",hjust=0.5, horizontal = F, offset = offset)
  
  phyloAutoTemp
  
  return(phyloAutoTemp)
  
}

## phylogeny with only secondary growth 

rotl_dend <- plotphyloGGtree_rotl(
  phylogenetictree       = my_phylo_branched,
  subsetphyloanalysisdf  = meta_Z_sec,
  variableEvo            = "meanZ",
  offset                 = 50
)

print(rotl_dend)

rotl_dend + geom_tree(linewidth = 1.2, colour = "black") +
  geom_fruit(
    data = ring2,
    geom = geom_tile,
    mapping = aes(y = label, fill = meanZ),
    width = 0.01, 
    offset = 0.01      # thickness of the ring
  ) +
  geom_tiplab(aes(label = label), align = TRUE, linetype = NA,
              hjust = -0.1, size = 4, aes(angle=angle), colour="black", fontface = 3, offset = 0.01) +
  theme(legend.position ="left")

cols <- brewer.pal(n = 11, name = "RdBu") 
cols <- brewer.pal(n = 11, name = "YlGn")  

######
old = ring2$tip_label
new = ring2$tip_label_an

syn_lookup <- as.tibble(old)
syn_lookup$new <- new
syn_lookup$value <- sub("_", " ", old) #remove underscore

dat <- rotl_dend$data %>%
  filter(isTip) %>%
  dplyr::select(label) %>%
  left_join(syn_lookup, by = c("label" = "value")) %>%
  mutate(label_new = coalesce(new, label))  # fallback to original if no match

##################

suppl_phylo_rotl <- rotl_dend %<+% dat +
  geom_tree(linewidth = 1.2, colour = "black") +
  geom_fruit(
    data = ring2,
    geom = geom_tile,
    mapping = aes(y = label, fill = meanZ),
    width = 0.01, 
    offset = 0.01      # thickness of the ring
  ) +
  geom_tiplab(aes(label = label_new), align = TRUE, linetype = NA,
              hjust = -0.1, size = 5, aes(angle=angle), colour="black", fontface = 3, offset = 0.01) + xlim(0,1.3) +
  theme(aspect.ratio = (1.2)) + theme_tree2() + theme(
    legend.position = c(.25, .75),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)) 

suppl_phylo_rotl

### SAVE THE FIG.S8

ggsave(
  filename = "plot_suppl_phylogeny_annot_rotl.pdf",
  plot     = suppl_phylo_rotl,
  #device   = cairo_pdf,   # needs Cairo installed; on Windows it ships with R, on macOS use XQuartz if needed
  width    = 14,
  height   = 15,
  units    = "in"
)
##############################################################################

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## supplementary material - TABLE S3

## checking the phylogenetic signal across functional traits and CVp of seed production

##### ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

## check the phylogenetic signal in other traits

meta_trait_seed <- meta_FINAL_1y %>%
  dplyr::select(tip_label,seed_mass)

meta_trait_wood <- meta_FINAL_1y %>%
  dplyr::select(tip_label,wood_density_imp)

meta_trait_SLA <- meta_FINAL_1y %>%
  dplyr::select(tip_label,SLA)

meta_trait_height <- meta_FINAL_1y %>%
  dplyr::select(tip_label,plant_height)

meta_trait_CVp <- meta_FINAL_1y %>%
  dplyr::select(tip_label,meanCVp)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

## function for phylogenetical analysis for seed mass

phyloAnalysisFun<-function(phylogenetictree = new_phylo_updated, 
                           subsetphyloanalysisdf = meta_trait_seed, 
                           variableEvo = seed_mass){
  dfPhylo<- subsetphyloanalysisdf[!is.na(subsetphyloanalysisdf[variableEvo]),]
  Drop<- phylogenetictree$tip.label[-c(which(phylogenetictree$tip.label %in% dfPhylo$tip_label))]
  phylo_up<- drop.tip(phylogenetictree,Drop)
  dfPhylo<- dfPhylo[match(phylo_up$tip.label, dfPhylo$tip_label),] %>% as.data.frame()
  
  rownames(dfPhylo) <- dfPhylo$tip_label
  
  dfPhylo$tip_label<- NULL
  X <- as.matrix(dfPhylo)
  phyAnalyse <- phylosig(phylo_up, X[,1], method="lambda", test=TRUE, nsim=10000, se=NULL, start=NULL,
                         control=list())
  
  return(phyAnalyse)
}


#check lambda phylo signal for seed mass

result_phylo_seed <- phyloAnalysisFun(
  phylogenetictree     = new_phylo_updated,
  subsetphyloanalysisdf= meta_trait_seed,
  variableEvo          = "seed_mass"
)

print(result_phylo_seed)

## function for phylogenetical analysis for wood density

phyloAnalysisFun<-function(phylogenetictree = new_phylo_updated, 
                           subsetphyloanalysisdf = meta_trait_wood, 
                           variableEvo = wood_density_imp){
  dfPhylo<- subsetphyloanalysisdf[!is.na(subsetphyloanalysisdf[variableEvo]),]
  Drop<- phylogenetictree$tip.label[-c(which(phylogenetictree$tip.label %in% dfPhylo$tip_label))]
  phylo_up<- drop.tip(phylogenetictree,Drop)
  dfPhylo<- dfPhylo[match(phylo_up$tip.label, dfPhylo$tip_label),] %>% as.data.frame()
  
  rownames(dfPhylo) <- dfPhylo$tip_label
  
  dfPhylo$tip_label<- NULL
  X <- as.matrix(dfPhylo)
  phyAnalyse <- phylosig(phylo_up, X[,1], method="lambda", test=TRUE, nsim=10000, se=NULL, start=NULL,
                         control=list())
  
  return(phyAnalyse)
}

#check lambda phylo signal for wood density

result_phylo_wood <- phyloAnalysisFun(
  phylogenetictree     = new_phylo_updated,
  subsetphyloanalysisdf= meta_trait_wood,
  variableEvo          = "wood_density_imp"
)

print(result_phylo_wood)


## function for phylogenetical analysis for SLA

phyloAnalysisFun<-function(phylogenetictree = new_phylo_updated, 
                           subsetphyloanalysisdf = meta_trait_SLA, 
                           variableEvo = SLA){
  dfPhylo<- subsetphyloanalysisdf[!is.na(subsetphyloanalysisdf[variableEvo]),]
  Drop<- phylogenetictree$tip.label[-c(which(phylogenetictree$tip.label %in% dfPhylo$tip_label))]
  phylo_up<- drop.tip(phylogenetictree,Drop)
  dfPhylo<- dfPhylo[match(phylo_up$tip.label, dfPhylo$tip_label),] %>% as.data.frame()
  
  rownames(dfPhylo) <- dfPhylo$tip_label
  
  dfPhylo$tip_label<- NULL
  X <- as.matrix(dfPhylo)
  phyAnalyse <- phylosig(phylo_up, X[,1], method="lambda", test=TRUE, nsim=10000, se=NULL, start=NULL,
                         control=list())
  
  return(phyAnalyse)
}


#check lambda phylo signal for SLA 

result_phylo_SLA <- phyloAnalysisFun(
  phylogenetictree     = new_phylo_updated,
  subsetphyloanalysisdf= meta_trait_SLA,
  variableEvo          = "SLA"
)

print(result_phylo_SLA)

## function for phylogenetical analysis for plant height


phyloAnalysisFun<-function(phylogenetictree = new_phylo_updated, 
                           subsetphyloanalysisdf = meta_trait_height, 
                           variableEvo = plant_height){
  dfPhylo<- subsetphyloanalysisdf[!is.na(subsetphyloanalysisdf[variableEvo]),]
  Drop<- phylogenetictree$tip.label[-c(which(phylogenetictree$tip.label %in% dfPhylo$tip_label))]
  phylo_up<- drop.tip(phylogenetictree,Drop)
  dfPhylo<- dfPhylo[match(phylo_up$tip.label, dfPhylo$tip_label),] %>% as.data.frame()
  
  rownames(dfPhylo) <- dfPhylo$tip_label
  
  dfPhylo$tip_label<- NULL
  X <- as.matrix(dfPhylo)
  phyAnalyse <- phylosig(phylo_up, X[,1], method="lambda", test=TRUE, nsim=10000, se=NULL, start=NULL,
                         control=list())
  
  return(phyAnalyse)
}

###
#check lambda phylo signal for height

result_phylo_height <- phyloAnalysisFun(
  phylogenetictree     = new_phylo_updated,
  subsetphyloanalysisdf= meta_trait_height,
  variableEvo          = "plant_height"
)

print(result_phylo_height)

## function for phylogenetical analysis for mean CVp

phyloAnalysisFun<-function(phylogenetictree = new_phylo_updated, 
                           subsetphyloanalysisdf = meta_trait_CVp, 
                           variableEvo = meanCVp){
  dfPhylo<- subsetphyloanalysisdf[!is.na(subsetphyloanalysisdf[variableEvo]),]
  Drop<- phylogenetictree$tip.label[-c(which(phylogenetictree$tip.label %in% dfPhylo$tip_label))]
  phylo_up<- drop.tip(phylogenetictree,Drop)
  dfPhylo<- dfPhylo[match(phylo_up$tip.label, dfPhylo$tip_label),] %>% as.data.frame()
  
  rownames(dfPhylo) <- dfPhylo$tip_label
  
  dfPhylo$tip_label<- NULL
  X <- as.matrix(dfPhylo)
  phyAnalyse <- phylosig(phylo_up, X[,1], method="lambda", test=TRUE, nsim=10000, se=NULL, start=NULL,
                         control=list())
  
  return(phyAnalyse)
}

#check lambda phylo signal for CVp

result_phylo_CVp <- phyloAnalysisFun(
  phylogenetictree     = new_phylo_updated,
  subsetphyloanalysisdf= meta_trait_CVp, 
  variableEvo          = "meanCVp"
)

print(result_phylo_CVp)

### Blomberg's K for traits

signal_data_seed <- meta_trait_seed %>%
  dplyr::select(tip_label, seed_mass) %>%
  filter(
    !is.na(tip_label),
    !is.na(seed_mass)
  ) %>%
  distinct(tip_label, .keep_all = TRUE)

#Check whether the species names agree (if so they equal 0)

setdiff(signal_data_seed$tip_label, new_phylo_updated$tip.label)

#Now retain only species present in both objects:

 
Z_vector_seed <- signal_data_seed$seed_mass
names(Z_vector_seed) <- signal_data_seed$tip_label
 
#Reorder it to match the tree exactly (they return TRUE if correct)

Z_vector_seed <- Z_vector_seed[tree_K$tip.label]

identical(names(Z_vector_seed), tree_K$tip.label)


# now run Blomberg’s K:

set.seed(123)

K_result_seed <- phylosig(
  tree = tree_K,
  x = Z_vector_seed,
  method = "K",
  test = TRUE,
  nsim = 9999
)

K_result_seed
 
## wood density 

signal_data_wood <- meta_trait_wood %>%
  dplyr::select(tip_label, wood_density_imp) %>%
  filter(
    !is.na(tip_label),
    !is.na(wood_density_imp)
  ) %>%
  distinct(tip_label, .keep_all = TRUE)

#Check whether the species names agree (if so they equal 0)

setdiff(signal_data_wood$tip_label, new_phylo_updated$tip.label)

#Now retain only species present in both objects:

Z_vector_wood <- signal_data_wood$wood_density_imp
names(Z_vector_wood) <- signal_data_wood$tip_label

#Reorder it to match the tree exactly (they return TRUE if correct)

Z_vector_wood <- Z_vector_wood[tree_K$tip.label]

identical(names(Z_vector_wood), tree_K$tip.label)


# now run Blomberg’s K:

set.seed(123)

K_result_wood <- phylosig(
  tree = tree_K,
  x = Z_vector_wood,
  method = "K",
  test = TRUE,
  nsim = 9999
)

K_result_wood


## SLA

signal_data_SLA <- meta_trait_SLA %>%
  dplyr::select(tip_label, SLA) %>%
  filter(
    !is.na(tip_label),
    !is.na(SLA)
  ) %>%
  distinct(tip_label, .keep_all = TRUE)

#Check whether the species names agree (if so they equal 0)

setdiff(signal_data_SLA$tip_label, new_phylo_updated$tip.label)

#Now retain only species present in both objects:

Z_vector_SLA<- signal_data_SLA$SLA
names(Z_vector_SLA) <- signal_data_SLA$tip_label

#Reorder it to match the tree exactly (they return TRUE if correct)

Z_vector_SLA <- Z_vector_SLA[tree_K$tip.label]

identical(names(Z_vector_SLA), tree_K$tip.label)


# now run Blomberg’s K:

set.seed(123)

K_result_SLA<- phylosig(
  tree = tree_K,
  x = Z_vector_SLA,
  method = "K",
  test = TRUE,
  nsim = 9999
)

K_result_SLA

## plant height

signal_data_height <- meta_trait_height %>%
  dplyr::select(tip_label, plant_height) %>%
  filter(
    !is.na(tip_label),
    !is.na(plant_height)
  ) %>%
  distinct(tip_label, .keep_all = TRUE)

#Check whether the species names agree (if so they equal 0)

setdiff(signal_data_height$tip_label, new_phylo_updated$tip.label)

#Now retain only species present in both objects:

Z_vector_height<- signal_data_height$plant_height
names(Z_vector_height) <- signal_data_height$tip_label

#Reorder it to match the tree exactly (they return TRUE if correct)

Z_vector_height <- Z_vector_height[tree_K$tip.label]

identical(names(Z_vector_height), tree_K$tip.label)


# now run Blomberg’s K:

set.seed(123)

K_result_height <- phylosig(
  tree = tree_K,
  x = Z_vector_height,
  method = "K",
  test = TRUE,
  nsim = 9999
)

K_result_height


## CVp

signal_data_CVp <- meta_trait_CVp %>%
  dplyr::select(tip_label, meanCVp) %>%
  filter(
    !is.na(tip_label),
    !is.na(meanCVp)
  ) %>%
  distinct(tip_label, .keep_all = TRUE)

#Check whether the species names agree (if so they equal 0)

setdiff(signal_data_CVp$tip_label, new_phylo_updated$tip.label)

#Now retain only species present in both objects:
common_species2 <- intersect(
  new_phylo_updated$tip.label,
  signal_data_CVp$tip_label
)

tree_K2 <- keep.tip(
  new_phylo_updated,
  common_species2
)

# General checks of the tree 
is.rooted(tree_K2)
is.binary(tree_K2)
is.ultrametric(tree_K2)
is.null(tree_K$edge.length)

# Problematic labels
anyDuplicated(tree_K2$tip.label)
tree_K2$tip.label[duplicated(tree_K2$tip.label)]

# Problematic branch lengths
which(tree_K2$edge.length == 0)
which(tree_K2$edge.length < 0)
range(tree_K2$edge.length)

# Tree covariance matrix
C2 <- vcv.phylo(tree_K2)

dim(C2)
qr(C2)$rank
nrow(C2)
kappa(C2)


signal_data_K2 <- signal_data_CVp %>%
  filter(tip_label %in% common_species2)

#Create the named vector:

Z_vector_CVp <- signal_data_CVp$meanCVp
names(Z_vector_CVp) <- signal_data_CVp$tip_label

#Reorder it to match the tree exactly (they return TRUE if correct)

Z_vector_CVp <- Z_vector_CVp[tree_K2$tip.label]

identical(names(Z_vector_CVp), tree_K2$tip.label)


# now run Blomberg’s K:

set.seed(123)

K_result_CVp <- phylosig(
  tree = tree_K2,
  x = Z_vector_CVp,
  method = "K",
  test = TRUE,
  nsim = 9999
)

K_result_CVp


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## function for phylogenetical analysis on subset of data (secondary growth only)

phyloAnalysisFun <-function(phylogenetictree = new_phylo_updated, 
                           subsetphyloanalysisdf = meta_Z_sec, 
                           variableEvo = meanZ){
  dfPhylo<- subsetphyloanalysisdf[!is.na(subsetphyloanalysisdf[variableEvo]),]
  Drop<- phylogenetictree$tip.label[-c(which(phylogenetictree$tip.label %in% dfPhylo$tip_label))]
  phylo_up<- drop.tip(phylogenetictree,Drop)
  dfPhylo<- dfPhylo[match(phylo_up$tip.label, dfPhylo$tip_label),] %>% as.data.frame()
  
  rownames(dfPhylo) <- dfPhylo$tip_label
  
  dfPhylo$tip_label<- NULL
  X <- as.matrix(dfPhylo)
  phyAnalyse <- phylosig(phylo_up, X[,1], method="lambda", test=TRUE, nsim=10000, se=NULL, start=NULL,
                         control=list())
  
  return(phyAnalyse)
}


#check lambda phylo signal 

result_phylo_secondary_gr <- phyloAnalysisFun(
  phylogenetictree     = new_phylo_updated,
  subsetphyloanalysisdf= meta_Z_sec,
  variableEvo          = "meanZ"
)

print(result_phylo_secondary_gr)


signal_data_sec <- meta_Z_sec %>%
  dplyr::select(tip_label, meanZ) %>%
  filter(
    !is.na(tip_label),
    !is.na(meanZ)
  ) %>%
  distinct(tip_label, .keep_all = TRUE)

#Check whether the species names agree (if so they equal 0)

setdiff(signal_data_sec$tip_label, new_phylo_updated$tip.label)
setdiff(signal_data_sec$tip_label, my_phylo_branched$tip.label)

#Now retain only species present in both objects:

common_species_sec <- intersect(
  my_phylo_branched$tip.label,
  signal_data_sec$tip_label
)

tree_K_sec <- keep.tip(
  my_phylo_branched,
  common_species_sec
)

##

common_species_sec <- intersect(
  new_phylo_updated$tip.label,
  signal_data_sec$tip_label
)

tree_K_sec <- keep.tip(
  new_phylo_updated,
  common_species_sec
)
# General checks of the tree 
is.rooted(tree_K_sec)
is.binary(tree_K_sec)
is.ultrametric(tree_K_sec)
is.null(tree_K_sec$edge.length)

# Problematic labels
anyDuplicated(tree_K_sec$tip.label)
tree_K_sec$tip.label[duplicated(tree_K_sec$tip.label)]

# Problematic branch lengths
which(tree_K_sec$edge.length == 0)
which(tree_K_sec$edge.length < 0)
range(tree_K_sec$edge.length)

# Tree covariance matrix
C_sec <- vcv.phylo(tree_K_sec)

dim(C_sec)
qr(C_sec)$rank
nrow(C_sec)
kappa(C_sec)


signal_data_K_sec <- signal_data_sec %>%
  filter(tip_label %in% common_species_sec)

#Create the named vector:

Z_vector_sec <- signal_data_K_sec$meanZ
names(Z_vector_sec) <- signal_data_K_sec$tip_label

#Reorder it to match the tree exactly (they return TRUE if correct)

Z_vector_sec <- Z_vector_sec[tree_K_sec$tip.label]

identical(names(Z_vector_sec), tree_K_sec$tip.label)


# now run Blomberg’s K:

set.seed(123)

K_result_sec <- phylosig(
  tree = tree_K_sec,
  x = Z_vector_sec,
  method = "K",
  test = TRUE,
  nsim = 9999
)

K_result_sec

#---------------------------------------------------------------------------------------
# phyloAnalysisFun (description)
#---------------------------------------------------------------------------------------
# This function extracts a subset of species data for a specified trait (variableEvo),
# prunes the given phylogenetic tree to retain only those species, and then calculates
# the phylogenetic signal (λ) using the phylosig() function from the 'phytools' package.
#
# Arguments:
#   phylogenetictree     : A phylogenetic tree of class 'phylo' (e.g., from read.tree()).
#   subsetphyloanalysisdf: A data frame containing at least two columns:
#                         - 'Species' : Species names matching those in the phylogenetic tree tip labels.
#                         - variableEvo: The trait or variable of interest used for phylogenetic analysis.
#   variableEvo          : A string specifying the name of the trait/column in subsetphyloanalysisdf
#                         to be analyzed (e.g., "SSD", "Height", etc.).
#
# Workflow:
#   1) Filter out rows in 'subsetphyloanalysisdf' where 'variableEvo' is NA.
#   2) Identify tree tips not present in the filtered data frame, and drop them from 'phylogenetictree'
#      using drop.tip(), creating 'phylo_up'.
#   3) Reorder the data frame rows to match the tip labels in 'phylo_up'.
#   4) Convert the specified trait column to a matrix 'X'.
#   5) Run 'phylosig()' on 'phylo_up' and the trait column 'X[,1]', using method="lambda" to estimate λ
#      (Pagel's lambda), and perform a significance test (nsim=10,000).
#   6) Return the result of 'phylosig()'.
#
# Return Value:
#   An object of class 'phylosig', typically containing:
#     - The estimate of Pagel's lambda (λ).
#     - A p-value for the test of λ vs. no phylogenetic signal.
#     - Additional information on model fitting and the simulation test.
#
# Note:
#   - Ensure that 'phylogenetictree$tip.label' matches the 'Species' column in your data.
#   - If there are mismatches, those tips are dropped, and the corresponding data rows
#     are removed or reordered for consistent alignment.
#   - The 'phylosig()' function is from the 'phytools' package.
#
# Usage Example:
#   result_phylo <- phyloAnalysisFun(
#       phylogenetictree     = my_tree,
#       subsetphyloanalysisdf= my_data,
#       variableEvo          = "SSD"
#   )
#   print(result_phylo)
#
#------------------

