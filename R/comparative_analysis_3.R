# Load library
library(ggplot2) # use ggplot2 method of plotting
library(cowplot) # organise figures
library(tidyverse) 
library(colorspace)
library(brms)
library(rstan)
library(ggeffects)
library(ape)
library(phytools)
library(diversitree)
library(corHMM)
library(hisse)

# Set directory
setwd('C:/Users/nicho/Dropbox (Personal)/Banded behaviour study')

# Functions
mytheme <- function() {
  theme_bw() + 
    theme(panel.border          = element_rect(fill = NA, colour = "black"), # set border around plot.
          panel.grid.major      = element_blank(), # remove major grid lines
          panel.grid.minor      = element_blank(), # remove minor grid lines
          axis.line             = element_blank(), # remove axis lines
          axis.ticks            = element_line(colour = "black"),
          axis.text             = element_text(size = 10, colour = "black"), # axis text size
          axis.title            = element_text(size = 10), # axis title size
          axis.title.y          = element_text(vjust = 3), # increase distance from the y-axis
          axis.title.x          = element_text(vjust = -1), # increase distance from the x-axis
          panel.background      = element_rect(fill = NA),
          plot.background       = element_rect(fill = NA, color = NA), # remove background colour
          plot.margin           = unit(c(1, 1, 1, 1), units = , "cm"), 
          legend.background     = element_rect(fill = NA, color = NA), # get rid of legend bg
          legend.box.background = element_rect(fill = NA, color = NA), # get rid of legend panel bg
          strip.text.x          = element_text(size = 10, color = "black", face = "bold"), # for facet plots
          strip.background      = element_rect(fill = NA, color = NA)
    )
} # set up plot theme

# Load data
diet_dat   <- read.csv("snake_diet_phylogeny.csv")
#phylo_tree <- ape::read.tree("squam_shl_new_Consensus_9755.tre") # Fully-sampled phylogeny
phylo_tree <- ape::read.tree("squam_shl_new.tre") # Tree with molecular supermatrix 

diet_dat$species_tree[duplicated(diet_dat$species_tree)] # check if there are duplicate species names

diet_dat <- diet_dat %>%
  dplyr::mutate(amphibian   = ifelse(amphibian == "Yes",1,0),
                snake       = ifelse(snake == "Yes",1,0),
                lizard      = ifelse(lizard == "Yes",1,0),
                fish        = ifelse(fish == "Yes",1,0),
                insect      = ifelse(insect == "Yes",1,0),
                annelid     = ifelse(annelid == "Yes",1,0),
                arthropod   = ifelse(arthropod == "Yes",1,0),
                crustacean  = ifelse(crustacean == "Yes",1,0),
                mollusc     = ifelse(mollusc == "Yes",1,0),
                bird        = ifelse(bird == "Yes",1,0),
                bird_egg    = ifelse(bird_egg == "Yes",1,0),
                reptile_egg = ifelse(reptile_egg == "Yes",1,0),
                mammal      = ifelse(mammal == "Yes",1,0),
                unknown     = ifelse(unknown == "Yes",1,0))

diet_dat <- diet_dat %>%
  dplyr::mutate(diet_total     = rowSums(diet_dat[c(12:25)] == 1), # count rows with number = 1
                ophio_freq     = round(ifelse(diet_total != 0, snake / diet_total, 0), digits = 2), # change NaN to zero, when zero divided by zero
                rept_special   = rowSums(diet_dat[c("snake", "lizard")] == 1),
                ophiophagy     = ifelse(ophio_freq >= 0.50, "Yes", "No"),
                ophiophagy_bin = ifelse(ophiophagy == "Yes",1,0),
                state          = case_when(
                  band_bin %in% 0 & ophiophagy_bin %in% 0 ~ 1, # state 1 = no bands and no ophiophagy
                  band_bin %in% 1 & ophiophagy_bin %in% 0 ~ 2, # state 2 = bands and no ophiophagy
                  band_bin %in% 0 & ophiophagy_bin %in% 1 ~ 3, # state 3 = no bands and ophiophagy
                  band_bin %in% 1 & ophiophagy_bin %in% 1 ~ 4)) # state 4 = bands and ophiophagy

## DATA SUMMARY ##--------------------------------------------------------------------------------------------------------
diet_dat %>%
  dplyr::group_by(ophiophagy) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n) * 100)

pattern_sum <- as.data.frame(diet_dat %>%
  dplyr::group_by(ophiophagy, pattern_type) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n) * 100) %>%
  dplyr::ungroup() %>%
  tidyr::complete(pattern_type, ophiophagy, fill = list(n = 0, freq = 0))) # Turns implicit missing values into explicit missing values

## SNAKE PHYLOGENY ##-----------------------------------------------------------------------------------------------------
tree_tip_label  <- phylo_tree$tip.label # extract tree tip names
species_list    <- diet_dat$species_tree # extract species name from diet analysis

pruned_tree     <- ape::drop.tip(phylo_tree, setdiff(phylo_tree$tip.label, species_list)) # prune phylo_tree to keep species from the diet analysis
diet_phylo_dat  <- subset(diet_dat, species_tree %in% pruned_tree$tip.label) # remove species from dataset not available in the pruned_tree
diet_phylo_dat$species_tree <- factor(diet_phylo_dat$species_tree, levels = c(tree_tip_label))

# check if lengths match for both data and tree
length(unique(pruned_tree$tip.label)) # 1543
length(unique(diet_phylo_dat$species_tree)) # 1543

# reorder factor levels and add species names to rowname
diet_phylo_dat$species_tree <- factor(diet_phylo_dat$species_tree, levels = pruned_tree$tip.label)
rownames(diet_phylo_dat)    <- NULL
trait_dat                   <- tibble::column_to_rownames(diet_phylo_dat, var = 'species_tree')
trait_dat$pattern_type      <- factor(trait_dat$pattern_type, 
                                      levels = c(
                                        "ringed",
                                        "stripe line",
                                        "uniform",
                                        "blotchy",
                                        "uniform with spots",
                                        "zig-zag",
                                        "bold ventral",
                                        "diamond"))

#plot.phylo(pruned_tree, show.tip.label = FALSE)
ape::is.binary(pruned_tree) # check if polytomies exist

# Time-calibrated tree via Penalised Likelihood (http://phylobotanist.blogspot.com/2018/04/time-calibrated-or-at-least-ultrametric.html)
tree_calibration <- ape::makeChronosCalib(pruned_tree, node = "root", age.min = 81.1, age.max = 137.0) # set 95% HPDI estimated root age (myr) from Tonini et al 2016

# Run all 3 clock models
timetree_corr <- ape::chronos(pruned_tree, lambda = 1, model = "correlated", calibration = tree_calibration, control = chronos.control()) # adjacent parts of the phylogeny are not allowed to evolve at rates that are very different
timetree_disc <- ape::chronos(pruned_tree, lambda = 1, model = "discrete", calibration = tree_calibration, control = chronos.control(nb.rate.cat = 1)) # models different parts of the tree as evolving at different rates
timetree_relax <- ape::chronos(pruned_tree, lambda = 1, model = "relaxed", calibration = tree_calibration, control = chronos.control()) # rates to vary freely

# log-likelihood higher is better, PHIIC lower is better
# timetree_corr = log-Lik = -359.7154, PHIIC = 9971.43 
# timetree_disc = log-Lik = -366.5213 , PHIIC = 3819.04
# timetree_relax = log-Lik = -364.5472 , PHIIC = 10004.29

pruned_tree_cal <- timetree_disc # discrete model best

#saveRDS(pruned_tree_cal, "pruned_tree_cal.rds")
pruned_tree_cal  <- readRDS("pruned_tree_cal.rds")

plot(pruned_tree_cal, cex = 0.3)
axisPhylo()

# check tree is ultrametric
is.ultrametric(pruned_tree_cal) # TRUE

# Create correlation matrix for analaysis
phylo_cor <- vcv(pruned_tree_cal, cor = T)
#saveRDS(phylo_cor, "phylo_cor.rds")
phylo_cor <- readRDS("phylo_cor.rds")

# Pearson correlation of subsetted data
family_n_subset <- diet_phylo_dat %>% dplyr::group_by(family) %>% summarise(n_subset = n())
family_n_full   <- diet_dat %>% dplyr::group_by(family) %>% summarise(n_full = n())
family_n        <- merge(family_n_full, family_n_subset, by = "family", all = TRUE)
family_cor      <- cor.test(family_n$n_full, family_n$n_subset, method = "pearson", use = "complete.obs")
family_cor$estimate
family_cor$p.value

as.data.frame(family_n %>% group_by(family) %>% summarise(n_subset / n_full * 100)) # proportion of species represented

ggplot(family_n, aes(x = log(n_full + 1), y = log(n_subset + 1))) + geom_point()+ mytheme() + xlab("Species per family (ln(n); full)") + ylab("Species per family (ln(n); strict)")

## PHYLOGENETIC SIGNAL AND ANCESTRAL STATE RECONSTRUCTION ## ---------------------------------------------------------------------------------------
## ASR of Ophiophagy (Continuous) ##
ophio_freq     <- setNames(trait_dat$ophio_freq, rownames(trait_dat))
#ophio_asr      <- phytools::fastAnc(pruned_tree_cal, ophio_freq, vars = TRUE, CI = TRUE)

# Plot ASR
ophio_freq_obj <- phytools::contMap(pruned_tree_cal, ophio_freq, plot = FALSE)
ophio_freq_obj <- phytools::setMap(ophio_freq_obj, c("grey", "red"))

# Fig 1a - internal (original pruned_tree - no ultrametric)
#plot(ophio_freq_obj, type = "fan", legend = 0.7 * max(nodeHeights(pruned_tree)), lwd = 0.8, outline = FALSE, ftype = "off")

## ASR of Ophiophagy and Pattern type (Discrete) ##
# Map continuous trait evolution on the tree
#ophio_matrix   <- setNames(trait_dat$ophiophagy_bin, rownames(trait_dat))
pattern_matrix <- setNames(as.numeric(trait_dat$pattern_type), rownames(trait_dat))

col_range <- colorRampPalette(c("grey", "red"))
ophio_col <- col_range(length(unique(trait_dat$ophio_freq)))
names(ophio_col) <- levels(factor(trait_dat$ophio_freq))

# Fig 1a
pattern_col <- viridis::viridis(8) # set 8 discrete colours
#diversitree::trait.plot(pruned_tree_cal, trait_dat, cols = list(ophiophagy_bin = c("grey", "red"), pattern_type = pattern_col), cex.lab = 0.00001, w = 0.06)
diversitree::trait.plot(pruned_tree_cal, trait_dat,
                        cols = list(pattern_type = pattern_col, ophio_freq = ophio_col),
                        cex.lab = 0.00001, w = 0.06)
plot(ophio_freq_obj$tree, 
     colors  = ophio_freq_obj$cols, 
     type    = "fan",
     add     = TRUE, 
     ftype   = "off", 
     lwd     = 0.5,
     outline = FALSE,
     xlim    = get("last_plot.phylo", envir =.PlotPhyloEnv)$x.lim,
     ylim    = get("last_plot.phylo", envir =.PlotPhyloEnv)$y.lim)

# Phylogenetic signal via Pagel's Lambda: 0 = no correlation, 1 = correlation between species equal to Brownian expectation

phytools::phylosig(pruned_tree, ophio_freq, method = "lambda", test = TRUE, nsim = 1000)
#phytools::phylosig(pruned_tree_cal, ophio_matrix, method = "lambda", test = TRUE, nsim = 1000)
phytools::phylosig(pruned_tree, pattern_matrix, method = "lambda", test = TRUE, nsim = 1000)


## ANCESTRAL STATE INFERENCE FOR OPHIOPHAGY STATUS VIA corHMM2.1 - ALL ## --------------------------------------------------------------------------
ophio_dat <- diet_phylo_dat %>% dplyr::select(species_tree, ophiophagy_bin)
geiger::name.check(pruned_tree$tree.label, ophio_dat$species_tree, data.names=NULL) #using geiger
#picante::match.phylo.data(pruned_tree, ophio_dat)  #using picante

## ER: Equal rates model ##
ophio_ER <- corHMM::corHMM(pruned_tree_cal, ophio_dat, model = c("ER"), 
                           node.states = "marginal", lewis.asc.bias = FALSE, root.p = c(1,0), rate.cat = 1, get.tip.states = TRUE, nstarts = 100)
#saveRDS(ophio_ER_ARD, "ophio_ER_ARD.rds")
#ophio_ER <- readRDS("ophio_ER.rds")

# Visualise model
corHMM::plotMKmodel(ophio_ER) 

## ARD: All rates different model ##
ophio_ARD <- corHMM::corHMM(pruned_tree_cal, ophio_dat, model = c("ARD"),
                            node.states = "marginal", lewis.asc.bias = FALSE, root.p = c(1,0), rate.cat = 1, get.tip.states = TRUE, nstarts = 100) 
#ophio_ARD <- readRDS("ophio_ARD.rds")

# Plot tree
ophio_ARD_model       <- ophio_ARD$solution
ophio_ARD_model[is.na(ophio_ARD_model)] <- 0
diag(ophio_ARD_model) <- -rowSums(ophio_ARD_model)
ophio_ARD_simmap      <- corHMM::makeSimmap(tree = pruned_tree_cal, data = ophio_dat, model = ophio_ARD_model, rate.cat = 1, nSim = 1, nCores = 1)
ophio_ARD_cols        <- setNames(c("grey","#DF536B"), c("1","2"))
phytools::plotSimmap(ophio_ARD_simmap[[1]], fsize = 0.001, lwd = 0.6, type = "fan", colors = ophio_ARD_cols)

## Model ARD/ARD: Hidden Markov model using ARD model in both matrices ##
ophio_ARD_ARD <- corHMM::corHMM(pruned_tree_cal, ophio_dat, model = c("ARD"), 
                                node.states = "marginal", lewis.asc.bias = FALSE, root.p = c(1,0), rate.cat = 2, get.tip.states = TRUE, nstarts = 100)
#ophio_ARD_ARD <- readRDS("ophio_ARD_ARD.rds")

## Model ER/ARD: Hidden Markov model using ARD model and ER model as the two rate matrices ##
# Construct two 'within' rate category rate.mat objects (R1 and R) - state-dependent processes
# Used to index the unique parameters to be estimated
ophio_R1_ER <- corHMM::getStateMat4Dat(ophio_dat, model = "ER")$rate.mat # R1 assume a drift-like hypothesis where all transition rates are equal
ophio_R2_ER <- corHMM::getStateMat4Dat(ophio_dat, model = "ER")$rate.mat # R2
ophio_R2_ARD <- corHMM::getStateMat4Dat(ophio_dat, model = "ARD")$rate.mat # R2 assume differences in transition rate

# getRateCatMat function to generate matrix for transitions 'among' the different rate classes - parameter process
RateClassMat  <- corHMM::getRateCatMat(2)

# Create list where first element corresponds to R1 and the second to R2 to create full model using getFullMat
ophio_FullMat <- corHMM::getFullMat(list(ophio_R1_ER, ophio_R2_ARD), RateClassMat)

# Visualize model setup
#corHMM::plotMKmodel(ophio_FullMat, rate.cat = 2, text.scale = 0.7) # check that the model is correct

ophio_ER_ARD <- corHMM::corHMM(pruned_tree_cal, ophio_dat, rate.cat = 2, rate.mat = ophio_FullMat, 
                               node.states = "marginal", lewis.asc.bias = FALSE, root.p = c(1,0), get.tip.states = TRUE, nstarts = 100) 
#ophio_ER_ARD <- readRDS("ophio_ER_ARD.rds")

## Model ER/ER: Hidden Markov model using two ER models as the two rate matrices ##
ophio_FullMat2 <- corHMM::getFullMat(list(ophio_R1_ER, ophio_R2_ER), RateClassMat) # R1 and R2 assume a drift-like hypothesis where all transition rates are equal

# Visualize model setup
#corHMM::plotMKmodel(ophio_FullMat2, rate.cat = 2, text.scale = 0.7)

opio_ER_ER <- corHMM::corHMM(pruned_tree_cal, ophio_dat, rate.cat = 2, rate.mat = ophio_FullMat2, 
                             node.states  = "marginal", lewis.asc.bias = FALSE,  root.p = c(1,0), get.tip.states = TRUE, nstarts = 100) 
#opio_ER_ER <- readRDS("opio_ER_ER.rds")

## ANCESTRAL STATE INFERENCE FOR BANDS STATUS VIA corHMM2.1 - ALL ## --------------------------------------------------------------------------
band_dat <- diet_phylo_dat %>% dplyr::select(species_tree, band_bin)

## ER: Equal rates model ##
band_ER <- corHMM::corHMM(pruned_tree_cal, band_dat, model = c("ER"), 
                           node.states = "marginal", lewis.asc.bias = FALSE, root.p = c(1,0), rate.cat = 1, get.tip.states = TRUE, nstarts = 100)
#band_ER <- readRDS("band_ER.rds")

## ARD: All rates different model ##
band_ARD <- corHMM::corHMM(pruned_tree_cal, band_dat, model = c("ARD"),
                            node.states = "marginal", lewis.asc.bias = FALSE, root.p = c(1,0), rate.cat = 1, get.tip.states = TRUE, nstarts = 100) 
#band_ARD <- readRDS("band_ARD.rds")

## Model ARD/ARD: Hidden Markov model using ARD model in both matrices ##
band_ARD_ARD <- corHMM::corHMM(pruned_tree_cal, band_dat, model = c("ARD"), 
                                node.states = "marginal", lewis.asc.bias = FALSE, root.p = c(1,0), rate.cat = 2, get.tip.states = TRUE, nstarts = 100)
#band_ARD_ARD <- readRDS("band_ARD_ARD.rds")

## Model ER/ARD: Hidden Markov model using ARD model and ER model as the two rate matrices ##
# Construct two 'within' rate category rate.mat objects (R1 and R) - state-dependent processes
# Used to index the unique parameters to be estimated
band_R1_ER  <- corHMM::getStateMat4Dat(band_dat, model = "ER")$rate.mat  # R1
band_R2_ER  <- corHMM::getStateMat4Dat(band_dat, model = "ER")$rate.mat  # R2
band_R2_ARD <- corHMM::getStateMat4Dat(band_dat, model = "ARD")$rate.mat # R2 

# Create list where first element corresponds to R1 and the second to R2 to create full model using getFullMat
band_FullMat <- corHMM::getFullMat(list(band_R1_ER, band_R2_ARD), RateClassMat)
band_ER_ARD <- corHMM::corHMM(pruned_tree_cal, band_dat, rate.cat = 2, rate.mat = band_FullMat, 
                               node.states = "marginal", lewis.asc.bias = FALSE, root.p = c(1,0), get.tip.states = TRUE, nstarts = 100) 
#band_ER_ARD <- readRDS("band_ER_ARD.rds")

## Model ER/ER: Hidden Markov model using two ER models as the two rate matrices ##
band_FullMat2 <- corHMM::getFullMat(list(band_R1_ER, band_R2_ER), RateClassMat) # R1 and R2 assume a drift-like hypothesis where all transition rates are equal
band_ER_ER <- corHMM::corHMM(pruned_tree_cal, band_dat, rate.cat = 2, rate.mat = band_FullMat2, 
                             node.states  = "marginal", lewis.asc.bias = FALSE,  root.p = c(1,0), get.tip.states = TRUE, nstarts = 100) 
#band_ER_ER <- readRDS("band_ER_ER.rds")

##Plotting the tree and model
par(mfrow = c(1,2))
# Ophiophagy
corHMM::plotRECON(pruned_tree_cal, ophio_ER_ARD$states, 
                  pie.cex = 0.3, 
                  piecolors =(c("black","orange", "grey","#DF536B")), 
                  show.tip.label = FALSE,
                  cex = 0.25, 
                  adj = 0.25) 
add.scale.bar(pruned_tree_cal)
cols <- setNames(c("grey", "#DF536B"), levels(ophio_dat$ophiophagy_bin))
ophio_dat2 <- ophio_dat %>% dplyr::arrange(factor(species_tree, levels = levels(ophio_dat$species_tree))) # reorder whole dataset to match phylogeny
ape::tiplabels(pie = to.matrix(ophio_dat2$ophiophagy_bin, sort(unique(ophio_dat$ophiophagy_bin))), piecol = cols, cex = 0.3)

# Banded rings
corHMM::plotRECON(pruned_tree_cal, band_ER_ARD$states, 
                  pie.cex = 0.3, 
                  piecolors =(c("black","#53d1df", "grey","#5384df")),
                  show.tip.label = FALSE,
                  cex = 0.25, 
                  adj = 0.25) 
add.scale.bar(pruned_tree_cal)
cols <- setNames(c("grey", "#5384df"), levels(band_dat$band_bin))
band_dat2 <- band_dat %>% dplyr::arrange(factor(species_tree, levels = levels(band_dat$species_tree))) # reorder whole dataset to match phylogeny
ape::tiplabels(pie = to.matrix(band_dat2$band_bin, sort(unique(band_dat$band_bin))), piecol = cols, cex = 0.3)
par(mfrow = c(1,1))

## ANCESTRAL STATE INFERENCE FOR OPHIOPHAGY AND BANDS STATUS VIA corHMM2.1 - ALL ## --------------------------------------------------------------------------
band_ophio_dat <- diet_phylo_dat %>% dplyr::select(species_tree, band_bin, ophiophagy_bin)

## ER: Equal rates model ##
band_ophio_ER <- corHMM::corHMM(pruned_tree_cal, band_ophio_dat, model = c("ER"), 
                          node.states = "marginal", lewis.asc.bias = FALSE, root.p = c(1,0,0,0), rate.cat = 1, get.tip.states = TRUE, nstarts = 100)
#saveRDS(band_ophio_ER, "band_ophio_ER.rds")
#band_ophio_ER <- readRDS("band_ophio_ER.rds")

#corHMM::plotMKmodel(band_ophio_ER) 

## SYM: All rates different model ##
band_ophio_SYM <- corHMM::corHMM(pruned_tree_cal, band_ophio_dat, model = c("SYM"),
                                 node.states = "marginal", lewis.asc.bias = FALSE, root.p = c(1,0,0,0), rate.cat = 1, get.tip.states = TRUE, nstarts = 100) 
#saveRDS(band_ophio_SYM, "band_ophio_SYM.rds")
#band_ophio_SYM <- readRDS("band_ophio_SYM.rds")

## ARD: All rates different model ##
band_ophio_ARD <- corHMM::corHMM(pruned_tree_cal, band_ophio_dat, model = c("ARD"),
                           node.states = "marginal", lewis.asc.bias = FALSE, root.p = c(1,0,0,0), rate.cat = 1, get.tip.states = TRUE, nstarts = 100) 
#saveRDS(band_ophio_ARD, "band_ophio_ARD.rds")
#band_ophio_ARD <- readRDS("band_ophio_ARD.rds")

# Construct two 'within' rate category rate.mat objects (R1 and R) - state-dependent processes
# Used to index the unique parameters to be estimated
band_ophio_R1_ER  <- corHMM::getStateMat4Dat(band_ophio_dat, model= "ER")$rate.mat # R1 assume a drift-like hypothesis where all transition rates are equal
band_ophio_R2_ER  <- corHMM::getStateMat4Dat(band_ophio_dat, model= "ER")$rate.mat # R2
band_ophio_R1_SYM <- corHMM::getStateMat4Dat(band_ophio_dat, model = "SYM")$rate.mat # R1 
band_ophio_R2_SYM <- corHMM::getStateMat4Dat(band_ophio_dat, model = "SYM")$rate.mat # R2 
band_ophio_R1_ARD <- corHMM::getStateMat4Dat(band_ophio_dat, model = "ARD")$rate.mat # R1 
band_ophio_R2_ARD <- corHMM::getStateMat4Dat(band_ophio_dat, model = "ARD")$rate.mat # R2 

## Model ER/ARD: Hidden Markov model using ARD model and ER model as the two rate matrices ##
# Create list where first element corresponds to R1 and the second to R2 to create full model using getFullMat
RateClassMat           <- corHMM::getRateCatMat(2)
RateClassMat1          <- corHMM::equateStateMatPars(RateClassMat, c(1, 2))
band_ophio_FullMat_mix <- corHMM::getFullMat(list(band_ophio_R1_ER, band_ophio_R2_ARD), RateClassMat1)

# Visualize model setup
#corHMM::plotMKmodel(band_ophio_FullMat, rate.cat = 2, text.scale = 0.7) # check that the model is correct

band_ophio_ER_ARD <- corHMM::corHMM(pruned_tree_cal, band_ophio_dat, rate.cat = 2, rate.mat = band_ophio_FullMat_mix, 
                              node.states = "marginal", lewis.asc.bias = FALSE, root.p = c(1,0,0,0), get.tip.states = TRUE, nstarts = 100) 
#saveRDS(band_ophio_ER_ARD, "band_ophio_ER_ARD.rds")
#band_ophio_ER_ARD <- readRDS("band_ophio_ER_ARD.rds")

## Model ER/ER: Hidden Markov model using two ER models as the two rate matrices ##
band_FullMat_ER  <- corHMM::getFullMat(list(band_ophio_R1_ER, band_ophio_R2_ER), RateClassMat) # R1 and R2 assume a drift-like hypothesis where all transition rates are equal
band_ophio_ER_ER <- corHMM::corHMM(pruned_tree_cal, band_ophio_dat, rate.cat = 2, rate.mat = band_FullMat_ER, 
                             node.states  = "marginal", lewis.asc.bias = FALSE,  root.p = c(1,0,0,0), get.tip.states = TRUE, nstarts = 100) 
#saveRDS(band_ophio_ER_ER, "band_ophio_ER_ER.rds")
#band_ophio_ER_ER <- readRDS("band_ophio_ER_ER.rds")

## Model SYM/SYM: Hidden Markov model using two ARD models as the two rate matrices ##
band_FullMat_SYM   <- corHMM::getFullMat(list(band_ophio_R1_SYM, band_ophio_R2_SYM), RateClassMat) # R1 and R2 assume a drift-like hypothesis where all transition rates are equal
band_ophio_SYM_SYM <- corHMM::corHMM(pruned_tree_cal, band_ophio_dat, rate.cat = 2, rate.mat = band_FullMat_SYM, 
                                     node.states  = "marginal", lewis.asc.bias = FALSE,  root.p = c(1,0,0,0), get.tip.states = TRUE, nstarts = 100) 
#saveRDS(band_ophio_SYM_SYM, "band_ophio_SYM_SYM.rds")
#band_ophio_SYM_SYM <- readRDS("band_ophio_SYM_SYM.rds")

## Model ARD/ARD: Hidden Markov model using two ARD models as the two rate matrices ##
band_FullMat_ARD   <- corHMM::getFullMat(list(band_ophio_R1_ARD, band_ophio_R2_ARD), RateClassMat) # R1 and R2 assume a drift-like hypothesis where all transition rates are equal
band_ophio_ARD_ARD <- corHMM::corHMM(pruned_tree_cal, band_ophio_dat, rate.cat = 2, rate.mat = band_FullMat_ARD, 
                             node.states  = "marginal", lewis.asc.bias = FALSE, root.p = c(1,0,0,0), get.tip.states = TRUE, nstarts = 100) 
#saveRDS(band_ophio_ARD_ARD, "band_ophio_ARD_ARD.rds")
#band_ophio_ARD_ARD <- readRDS("band_ophio_ARD_ARD.rds")

band_ophio_R1_ARD_dual  <- corHMM::getStateMat4Dat(band_ophio_dat, model = "ARD", dual = TRUE)$rate.mat #R1
band_ophio_R2_ARD_dual  <- corHMM::getStateMat4Dat(band_ophio_dat, model = "ARD", dual = TRUE)$rate.mat #R2
band_FullMat_ARD_dual   <- corHMM::getFullMat(list(band_ophio_R1_ARD_dual, band_ophio_R2_ARD_dual), RateClassMat) # R1 and R2 assume a drift-like hypothesis where all transition rates are equal
band_ophio_ARD_ARD_dual <- corHMM::corHMM(pruned_tree_cal, band_ophio_dat, rate.cat = 2, rate.mat = band_FullMat_ARD_dual, 
                                     node.states  = "marginal", lewis.asc.bias = FALSE,  root.p = c(1,0,0,0), get.tip.states = TRUE, nstarts = 100) 
saveRDS(band_ophio_ARD_ARD_dual, "band_ophio_ARD_ARD_dual.rds")
#band_ophio_ARD_ARD_dual <- readRDS("band_ophio_ARD_ARD_dual.rds")

corHMM::plotMKmodel(band_ophio_ARD_ARD_dual) 

band_ophio_ARD_ARD_2 <- corHMM::corHMM(pruned_tree, band_ophio_dat, rate.cat = 2, rate.mat = band_FullMat_ARD, 
                                     node.states  = "marginal", lewis.asc.bias = FALSE,  root.p = c(1,0,0,0), get.tip.states = TRUE, nstarts = 100)
saveRDS(band_ophio_ARD_ARD_2, "band_ophio_ARD_ARD_2.rds")

corHMM::plotMKmodel(band_ophio_ARD_ARD_2) 

##Plotting the tree and model
corHMM::plotRECON(pruned_tree_cal, band_ophio_ARD_ARD$states, 
                  pie.cex = 0.3, 
                  piecolors =(c("#808082", "#b98fc0", "#e6a0b8", "#fecab3",
                                "#000004", "#721F81", "#CD4071", "#FD9567")),  
                  cex = 0.2, 
                  adj = 0.25) 
add.scale.bar(pruned_tree_cal)
cols <- setNames(c("grey", "#5384df"), levels(band_dat$band_bin))
ape::tiplabels(pie = to.matrix(band_dat$band_bin, sort(unique(band_dat$band_bin))), piecol = cols, cex = 0.2)
par(mfrow = c(1,1))

# Plot tree
band_ophio_ARD_model       <- band_ophio_ARD$solution
band_ophio_ARD_model[is.na(band_ophio_ARD_model)] <- 0
diag(band_ophio_ARD_model) <- -rowSums(band_ophio_ARD_model)
band_ophio_ARD_simmap      <- corHMM::makeSimmap(tree = pruned_tree_cal, data = band_ophio_dat, model = band_ophio_ARD_model, rate.cat = 1, nSim = 1, nCores = 1)
band_ophio_ARD_cols        <- setNames(c("#000004", "#721F81", "#CD4071", "#FD9567"), c("1","2","3","4"))

phytools::plotSimmap(band_ophio_ARD_simmap[[1]], fsize = 0.001, lwd = 0.6, type = "fan", colors = band_ophio_ARD_cols)
phytools::add.simmap.legend(colors = band_ophio_ARD_cols)

band_ophio_ARD_ARD_model       <- band_ophio_ARD_ARD$solution
band_ophio_ARD_ARD_model[is.na(band_ophio_ARD_ARD_model)] <- 0
diag(band_ophio_ARD_ARD_model) <- -rowSums(band_ophio_ARD_ARD_model)
band_ophio_ARD_ARD_simmap      <- corHMM::makeSimmap(tree = pruned_tree_cal, data = band_ophio_dat, model = band_ophio_ARD_ARD_model, rate.cat = 2, nSim = 1, nCores = 1)
band_ophio_ARD_ARD_cols        <- setNames(c("#000004", "#721F81", "#CD4071", "#FD9567", "#808082", "#b98fc0", "#e6a0b8", "#fecab3"), 
                                       c("1","2","3","4","5","6","7","8"))

phytools::plotSimmap(band_ophio_ARD_ARD_simmap[[1]], fsize = 0.001, lwd = 0.6, type = " fan", colors = band_ophio_ARD_ARD_cols)
phytools::add.simmap.legend(colors = band_ophio_ARD_ARD_cols)

## EVOLUTIONARY SCENARIO - MuHiSSE ## ---------------------------------------------------------------------------------------
## Sampling fractions
sampling_fraction <- c(0.47, 0.43, 0.47, 0.52)

## MuSSE null 
trans_rate <- hisse::TransMatMakerMuHiSSE(hidden.traits = 0)
null_MuS   <- hisse::MuHiSSE(phy = pruned_tree_cal, data = band_ophio_dat, f = sampling_fraction, turnover = c(1,1,1,1),
                            eps = c(1,1,1,1), root.p = c(1,0,0,0), hidden.states = FALSE,
                            trans.rate = trans_rate)
saveRDS(null_MuS, "null_MuS.rds")
#null_MuS <- readRDS("null_MuS.rds")

## MuSSE true
true_MuS <- hisse::MuHiSSE(phy = pruned_tree_cal, data = band_ophio_dat, f = sampling_fraction, turnover = c(1,2,3,4),
                           eps = c(1,1,1,1), root.p = c(1,0,0,0), hidden.states = FALSE,
                           trans.rate = trans_rate)
saveRDS(true_MuS, "true_MuS.rds")

## Character-dependent MuHiSSE model with two hidden states
trans_rate_CharDep_TwoHidden <- hisse::TransMatMakerMuHiSSE(hidden.traits = 1)
MuH_CD2 <- hisse::MuHiSSE(phy = pruned_tree_cal, data = band_ophio_dat, f = sampling_fraction, turnover = c(1,2,3,4,5,6,7,8),
                          eps = rep(1, 8), root.p = c(1,0,0,0,1,0,0,0), hidden.states = TRUE,
                          trans.rate = trans_rate_CharDep_TwoHidden)
saveRDS(MuH_CD2, "MuH_CD2.rds")

## Character-independent MuHiSSE model with two hidden states 
trans_rate_CharIndep_TwoHidden <- hisse::TransMatMakerMuHiSSE(hidden.traits = 1, make.null = TRUE)
MuH_CID2 <- hisse::MuHiSSE(phy = pruned_tree_cal, data = band_ophio_dat, f = sampling_fraction, turnover = c(1,1,1,1,2,2,2,2),
                           eps = rep(1, 8), root.p = c(1,0,0,0,1,0,0,0), hidden.states = TRUE,
                           trans.rate = trans_rate_CharIndep_TwoHidden)
saveRDS(MuH_CID2, "MuH_CID2.rds")

## Character-dependent MuHiSSE model with three hidden states
trans_rate_CharDep_ThreeHidden <- hisse::TransMatMakerMuHiSSE(hidden.traits = 2)
MuH_CD3 <- hisse::MuHiSSE(phy = pruned_tree_cal, data = band_ophio_dat, f = sampling_fraction, turnover = c(1,2,3,4,5,6,7,8,9,10,11,12),
                          eps = rep(1, 12), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0), hidden.states = TRUE,
                          trans.rate = trans_rate_CharDep_ThreeHidden)
saveRDS(MuH_CD3, "MuH_CD3.rds")

## Character-independent MuHiSSE model with three hidden states 
trans_rate_CharIndep_ThreeHidden <- hisse::TransMatMakerMuHiSSE(hidden.traits = 2, make.null = TRUE)
MuH_CID3 <- hisse::MuHiSSE(phy = pruned_tree_cal, data = band_ophio_dat, f = sampling_fraction, turnover = c(rep(1,4), rep(2,4), rep(3,4)),
                           eps = rep(1, 12), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0), hidden.states = TRUE,
                           trans.rate = trans_rate_CharIndep_ThreeHidden)
saveRDS(MuH_CID3, "MuH_CID3.rds")

## Character-dependent MuHiSSE model with four hidden states
trans_rate_CharDep_FourHidden <- hisse::TransMatMakerMuHiSSE(hidden.traits = 3)
MuH_CD4 <- hisse::MuHiSSE(phy = pruned_tree_cal, data = band_ophio_dat, f = sampling_fraction, turnover = c(1:16),
                          eps = rep(1, 16), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states = TRUE,
                          trans.rate = trans_rate_CharDep_FourHidden)
saveRDS(MuH_CD4, "MuH_CD4.rds")

## Character-independent MuHiSSE model with four hidden states 
trans_rate_CharIndep_FourHidden <- hisse::TransMatMakerMuHiSSE(hidden.traits = 3, make.null = TRUE)
MuH_CID4 <- hisse::MuHiSSE(phy = pruned_tree_cal, data = band_ophio_dat, f = sampling_fraction, turnover = c(rep(1,4), rep(2,4), rep(3,4), rep(4,4)),
                           eps = rep(1, 16), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states = TRUE,
                           trans.rate = trans_rate_CharIndep_FourHidden)
saveRDS(MuH_CID4, "MuH_CID4.rds")

## Character-dependent MuHiSSE model with five hidden states
trans_rate_CharDep_FiveHidden <- hisse::TransMatMakerMuHiSSE(hidden.traits = 4)
MuH_CD5 <- hisse::MuHiSSE(phy = pruned_tree_cal, data = band_ophio_dat, f = sampling_fraction, turnover = c(1:20),
                          eps = rep(1, 20), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states = TRUE,
                          trans.rate = trans_rate_CharDep_FiveHidden)
saveRDS(MuH_CD5, "MuH_CD5.rds")

## Character-independent MuHiSSE model with five hidden states 
trans_rate_CharIndep_FiveHidden <- hisse::TransMatMakerMuHiSSE(hidden.traits = 4, make.null = TRUE)
MuH_CID5 <- hisse::MuHiSSE(phy = pruned_tree_cal, data = band_ophio_dat, f = sampling_fraction, turnover = c(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4)),
                           eps = rep(1, 20), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states = TRUE,
                           trans.rate = trans_rate_CharIndep_FiveHidden)
saveRDS(MuH_CID5, "MuH_CID5.rds")

## Character-dependent MuHiSSE model with six hidden states
trans_rate_CharDep_SixHidden <- hisse::TransMatMakerMuHiSSE(hidden.traits = 5)
MuH_CD6 <- hisse::MuHiSSE(phy = pruned_tree_cal, data = band_ophio_dat, f = sampling_fraction, turnover = c(1:24),
                          eps = rep(1, 24), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states = TRUE,
                          trans.rate = trans_rate_CharDep_SixHidden)
saveRDS(MuH_CD6, "MuH_CD6.rds")

## Character-independent MuHiSSE model with six hidden states 
trans_rate_CharIndep_SixHidden <- hisse::TransMatMakerMuHiSSE(hidden.traits = 5, make.null = TRUE)
MuH_CID6 <- hisse::MuHiSSE(phy = pruned_tree_cal, data = band_ophio_dat, f = sampling_fraction, turnover = c(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4), rep(6,4)),
                           eps = rep(1, 24), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states = TRUE,
                           trans.rate = trans_rate_CharIndep_SixHidden)
saveRDS(MuH_CID6, "MuH_CID6.rds")

## Character-dependent MuHiSSE model with seven hidden states
trans_rate_CharDep_SevenHidden <- hisse::TransMatMakerMuHiSSE(hidden.traits = 6)
MuH_CD7 <-  hisse::MuHiSSE(phy = pruned_tree_cal, data = band_ophio_dat, f = sampling_fraction, turnover = c(1:28),
                           eps = rep(1, 28), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states = TRUE,
                           trans.rate = trans_rate_CharDep_SevenHidden)
saveRDS(MuH_CD7, "MuH_CD7.rds")

## Character-independent MuHiSSE model with seven hidden states 
trans_rate_CharIndep_SevenHidden <- hisse::TransMatMakerMuHiSSE(hidden.traits = 6, make.null = TRUE)
MuH_CID7 <- hisse:: MuHiSSE(phy = pruned_tree_cal, data = band_ophio_dat, f = sampling_fraction, turnover = c(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4), rep(6,4), rep(7, 4)),
                            eps = rep(1, 28), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states = TRUE,
                            trans.rate = trans_rate_CharIndep_SevenHidden)
saveRDS(MuH_CID7, "MuH_CID7.rds")

## Character-dependent MuHiSSE model with eight hidden states
trans_rate_CharDep_EightHidden <- hisse::TransMatMakerMuHiSSE(hidden.traits = 7)
MuH_CD8 <- hisse::MuHiSSE(phy = pruned_tree_cal, data = band_ophio_dat, f = sampling_fraction, turnover = c(1:32),
                          eps = rep(1, 32), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states = TRUE,
                          trans.rate = trans_rate_CharDep_EightHidden)
saveRDS(MuH_CD8, "MuH_CD8.rds")

## Character-independent MuHiSSE model with eight hidden states 
trans_rate_CharIndep_EightHidden <- hisse::TransMatMakerMuHiSSE(hidden.traits = 7, make.null = TRUE)
MuH_CID8 <- hisse::MuHiSSE(phy = pruned_tree_cal, data = band_ophio_dat, f = sampling_fraction, turnover = c(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4), rep(6,4),rep(7,4), rep(8,4)),
                           eps = rep(1, 32), root.p = c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), hidden.states = TRUE,
                           trans.rate = trans_rate_CharIndep_EightHidden)
saveRDS(MuH_CID8, "MuH_CID8.rds")

#Best fit model for Full Strict dataset is MuH_CD4
MuH_CD4  <- readRDS("MuH_CD4.rds")

# Extracting transition rate estimates
Solution            <- MuH_CD4$solution 
Solution_Transposed <- as.data.frame(t(as.matrix(Solution))) # Transpose the dataframe
Solution_Rates      <- Solution_Transposed %>% select(contains("_"))

## 1
Solution_Rates1      <- Solution_Rates %>% select(contains("00") & contains("10"))
Solution_Rates_00_10 <- Solution_Rates1 %>% select(contains("q00")) # q00 > q10
Solution_Rates_10_00 <- Solution_Rates1 %>% select(contains("q10")) # q10 > q00

## 2
Solution_Rates_q01q00 <- Solution_Rates %>% select(contains("01") & contains("00"))
Solution_Rates_00_01  <- Solution_Rates_q01q00 %>% select(contains("q00")) # q00 > q01
Solution_Rates_01_00  <- Solution_Rates_q01q00 %>% select(contains("q01")) # q01 > q00

## 3
Solution_Rates_q01q11 <- Solution_Rates %>% select(contains("01") & contains("11"))
Solution_Rates_01_11  <- Solution_Rates_q01q11 %>% select(contains("q01")) # q01 > q11
Solution_Rates_11_01  <- Solution_Rates_q01q11 %>% select(contains("q11")) # q11 > q01

## 4
Solution_Rates_q10q11 <- Solution_Rates %>% select(contains("10") & contains("11"))
Solution_Rates_10_11  <- Solution_Rates_q10q11 %>% select(contains("q10")) # q10 > q11
Solution_Rates_11_10  <- Solution_Rates_q10q11 %>% select(contains("q11")) # q11 > q10


StrictFull_qrates <- list(Solution_Rates_00_10, 
                          Solution_Rates_10_00, 
                          Solution_Rates_00_01, 
                          Solution_Rates_01_00, 
                          Solution_Rates_01_11, 
                          Solution_Rates_11_01, 
                          Solution_Rates_10_11, 
                          Solution_Rates_11_10) %>%
  setNames(c("00.10",
             "10.00",
             "00.01",
             "01.00",
             "01.11",
             "11.01",
             "10.11",
             "11.10"))

#Just sample one
for(i in 1:length(StrictFull_qrates)){
  MuHiSSEStrictFullRates <- sample(StrictFull_qrates[[i]][2]) * 100
  MuHiSSEStrictFullRatesPrint <- print(round(MuHiSSEStrictFullRates, 4))
}

# Ancestral State Reconstruction
MuH_CD4_recon <- hisse::MarginReconMuHiSSE(phy     = MuH_CD4$phy,
                                           data    = MuH_CD4$data,
                                           f       = MuH_CD4$f,
                                           pars    = MuH_CD4$solution,
                                           hidden.states = 4,
                                           AIC     = MuH_CD4$AICc,
                                           verbose = TRUE)
saveRDS(MuH_CD4_recon, "MuH_CD4_recon.rds")
#MuH_CD4_recon  <- readRDS("MuH_CD4_recon.rds")



devtools::install_github("discindo/gghisse")
library(ggtree)
library(gghisse)

MuH_CD4_phylo <- MuH_CD4_recon
class(MuH_CD4_phylo$phy) <- class(MuH_CD4_phylo$phy)[2] # remove chronos class
MuH_CD4_proc <- gghisse::m_process_recon(MuH_CD4_phylo)

# ASR plot
gghisse::m_trait_recon(processed_recon            = MuH_CD4_proc,
                       states_of_first_character  = c('0','1'),
                       states_of_second_character = c('0','1'),
                       cutoff                     = as.numeric(c('0.5','0.5')),
                       colors                     = c("#000004FF", "#721F81FF", "#CD4071FF", "#FD9567FF"),
                       show_tip_labels            = FALSE,
                       tree_layout                = 'fan',
                       tree_direction             = 'up',
                       time_axis_ticks            = 10,
                       open_angle                 = 5)

# Net diverisification plot
gghisse::m_scatterplot(processed_recon = MuH_CD4_proc,
                         parameter       = 'net.div',
                         states_names    = c('00','01','10','11'), 
                         colors          = c("#000004FF", "#721F81FF", "#CD4071FF", "#FD9567FF"), 
                         plot_as_waiting_time = FALSE) + 
  mytheme() +
  ylab("Net Diversfication")

## EVOLUTIONARY SCENARIO - Pagel's model ## ---------------------------------------------------------------------------------------
ophio_depend <- setNames(trait_dat$ophiophagy, rownames(trait_dat))
band_depend  <- setNames(trait_dat$band_bin, rownames(trait_dat))

# http://phytools.org/mexico2018/ex/8/Pagel94-method.html
independ_fit <- phytools::fitPagel(pruned_tree_cal, band_depend, ophio_depend) # Scenario 1 and 2 - Independent and dependent (correlated) evolution of ophiophagy
band_fit     <- phytools::fitPagel(pruned_tree_cal, band_depend, ophio_depend, dep.var = "x") # Scenario 3 - Evolutionary change in ophiophagy depends upon the state of banded stripe
ophio_fit    <- phytools::fitPagel(pruned_tree_cal, band_depend, ophio_depend, dep.var = "y") # Scenario 4 - Evolutionary change in banded stripe depends upon the state of ophiophagy
snake_aic    <- setNames(c(independ_fit$independent.AIC, band_fit$dependent.AIC, ophio_fit$dependent.AIC, independ_fit$dependent.AIC),
                         c("independent", "dependent band", "dependent ophio", "dependent band_ophio"))

plot(independ_fit, lwd.by.rate = TRUE)
plot(band_fit, lwd.by.rate = TRUE)
plot(ophio_fit, lwd.by.rate = TRUE)


## OPHIOPHAGY ANALYSIS ## ---------------------------------------------------------------------------------------------------------------------------
## PLOT FIG 1b - PERCENTAGE OF OPHIOPHAGY DIET BY FAMILY ##
family_sum_full <- as.data.frame(diet_dat %>%
                              group_by(family, ophiophagy) %>%
                              summarise(n = n()) %>%
                              mutate(freq = n / sum(n) * 100) %>%
                              ungroup() %>%
                              complete(family, ophiophagy, fill = list(n = 0, freq = 0))) # fill-in zeros for the summary variables

family_sum_strict <- as.data.frame(diet_phylo_dat %>%
                              group_by(family, ophiophagy) %>%
                              summarise(n = n()) %>%
                              mutate(freq = n / sum(n) * 100) %>%
                              ungroup() %>%
                              complete(family, ophiophagy, fill = list(n = 0, freq = 0))) # fill-in zeros for the summary variables

family_sum_full <- family_sum_full %>%
  dplyr::filter(ophiophagy == "Yes") %>%
  dplyr::filter(family %in% c("Elapidae", "Atractaspididae", "Cyclocoridae", "Lamprophiidae", "Cyclocoridae", 
                "Colubridae", "Pseudoxyrhophiidae", "Psammophiidae", "Boidae", "Viperidae", "Xenotyphlopidae")) %>%
  dplyr::mutate(family = recode_factor(family, "Xenotyphlopidae" = "Other"))


family_sum_strict <- family_sum_strict %>%
  dplyr::filter(ophiophagy == "Yes") %>%
  dplyr::filter(family %in% c("Elapidae", "Atractaspididae", "Cyclocoridae", "Lamprophiidae", "Cyclocoridae", 
                              "Colubridae", "Pseudoxyrhophiidae", "Psammophiidae", "Boidae", "Viperidae", "Xenotyphlopidae")) %>%
  dplyr::mutate(family = recode_factor(family, "Xenotyphlopidae" = "Other"))

family_plot <- family_sum_full %>%
  ggplot(aes(x = reorder(family, freq), y = freq)) +
  geom_point(position = position_nudge(x = 0.1)) +
  geom_segment(aes(x = family, xend = family, y = 0, yend = freq), colour = "grey", position = position_nudge(x = 0.1)) +
  geom_point(size = 2.5, colour = "#440154FF", position = position_nudge(x = 0.1)) +
  geom_segment(data = family_sum_strict, aes(x = family, xend = family, y = 0, yend = freq), colour = "grey", linetype = "dashed", position = position_nudge(x = -0.1)) +
  geom_point(data = family_sum_strict, aes(x = family, y = freq), size = 2.5, colour = "#277F8EFF", pch = 18, position = position_nudge(x = -0.1)) +
  ylab("Percentage of ophiophagus species (%)") + xlab(NULL) +
  mytheme() + coord_flip()


# PROBABILITY OF OPHIOPHAGY ##
str(diet_dat)
set.seed(1)

# Set options in Rstan
rstan::rstan_options(auto_write = TRUE) # translate to STAN platform for running Bayesian model
options(mc.cores = parallel::detectCores()) # detects how many cores available to use

diet_model <- brms::brm(ophiophagy ~ -1 + pattern_type, 
                        data    = diet_phylo_dat, 
                        family  = 'bernoulli', 
                        prior   = prior(normal(0, 3), "b"),
                        iter    = 5e3, 
                        warmup  = 2.5e3, 
                        chains  = 4, 
                        cores   = 4, 
                        control = list(adapt_delta = 0.99, max_treedepth = 15))
summary(diet_model)

# phylogenetic model
diet_dat$pattern_type <- as.factor(diet_dat$pattern_type)

set.seed(1)

library(Matrix)
tree_vcvc_new <- as.matrix(Matrix::nearPD(phylo_cor)$mat) # use nearest positive definite matrix
diet_phylo_model <- brms::brm(ophiophagy ~ -1 + pattern_type + (1 | gr(species_tree, cov = A)), 
                              data    = diet_phylo_dat, 
                              family  = 'bernoulli',
                              prior  = prior(normal(0, 3), "b"),
                              data2   = list(A = tree_vcvc_new),
                              iter    = 5e3, warmup = 2.5e3, chains = 4, cores = 4,
                              control = list(adapt_delta = 0.999, max_treedepth = 15))

summary(diet_phylo_model)
brms::fixef(diet_model)
brms::fixef(diet_phylo_model)
brms::pp_check(diet_phylo_model)

#saveRDS(diet_phylo_model, "diet_phylo_model.rds")
diet_phylo_model <- readRDS("diet_phylo_model.rds")

# compare non-phylogenetic model and phylogenetic-corrected model
performance::compare_performance(diet_model, diet_phylo_model, rank = TRUE)

# Phylogenetic signal via hypothesis method
#parnames(diet_phylo_model)
(hyp <- hypothesis(diet_phylo_model, "sd_species_tree__Intercept^2 / (sd_species_tree__Intercept^2 + sd_species_tree__Intercept) = 0", class = NULL))

# PLOT FIG 1c
# Extract estimates and 95% CI
predPattern <- as.data.frame(ggeffects::ggpredict(diet_phylo_model, "pattern_type")) %>%
  dplyr::rename(pattern_type = x) %>%
  dplyr::mutate(pattern_type = fct_reorder(pattern_type, predicted))

predPattern_plot <- predPattern %>%
  ggplot(aes(x = reorder(pattern_type, predicted), y = predicted, colour = pattern_type)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), size = 1, width = 0.1) + 
  colorspace::scale_colour_discrete_sequential(palette = "Viridis") +
  xlab(NULL) + ylab("Probability of pattern type with ophiophagus diet") +
  mytheme() +
  coord_flip()

# COMBINE FIG 1b and 1c
plot_grid(family_plot + theme(legend.position = "none"), 
          predPattern_plot + theme(legend.position = "none"),
          align = "v", axis = "lr", nrow = 2, labels = c("b", "c"),
          rel_heights = c(1, 0.7))


## ALTERNATIVE HYPOTHESIS ## ------------------------------------------------------------------------------------ 
# Frequency of banded stripes with ophiophagy, reptiles diet, mammal diet
as.data.frame(diet_dat %>%
                dplyr::group_by(pattern_type) %>%
                dplyr::summarise(ophio = sum(ophiophagy_bin == 1),
                                 reptile = sum(rept_special == 1),
                                 mammal = sum(mammal == 1)) %>%
                dplyr::mutate(ophio_freq = ophio / sum(ophio) * 100,
                              reptile_freq = reptile / sum(reptile) * 100,
                              mammal_freq = mammal / sum(mammal) * 100) %>%
                dplyr::ungroup() %>%
                tidyr::complete(pattern_type, fill = list(n = 0, freq = 0))) # Turns implicit missing values into explicit missing values) 


              
