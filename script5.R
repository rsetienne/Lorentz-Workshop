# fix random seed for reproducibility
seed <- 42
set.seed(seed)
# model parameters, modify as you see fit
lambda <- 2 # speciation rate
mu <- 0.3 # extinction rate
t_max <- 1 # duration of the simulation
tree <- TreeSim::sim.bd.age(age = t_max,
                           numbsim = 1, # only simulate one tree
                           lambda = lambda,
                           mu = mu,
                           mrca = FALSE, # age is time since origin, tree has a root edge
                           complete = TRUE)[[1]] # get the complete tree
# inspect the tree
tree
plot(tree,
     root.edge = TRUE,
     main = "Birth-death tree")
axis(1)
mtext("Time since origin [Myr]", side = 1, line = 2)
psi <- 5 # fossil sampling rate
fossils <- FossilSim::sim.fossils.poisson(rate = psi,
                                         tree = tree,
                                         root.edge = FALSE)
# inspect the fossils
fossils
plot(fossils,
     tree = tree,
     show.tip.label = TRUE,
     rho = 0) # suppress highlighting of extant tips
plot(fossils,
     tree = tree, 
     show.ranges = TRUE,
     rho = 0)
rho = 0.5 # proportion of extant tips that are sampled
fossils2 = FossilSim::sim.extant.samples(fossils = fossils,
                                         tree = tree,
                                         rho = rho)
# plot the output
plot(fossils2,
     tree = tree,
     extant.col = "red",
     show.tip.label = TRUE)
plot(fossils2,
     tree = tree,
     extant.col = "red",
     reconstructed = FALSE,
     main = "Simulated tree and fossils")
plot(fossils2,
     tree = tree,
     extant.col = "red",
     reconstructed = TRUE,
     main = "Sampled tree")
# generate sampled ancestor tree
SAt = FossilSim::SAtree.from.fossils(tree = tree,
                                     fossils = fossils2)
plot(SAt$tree,
     main = "Sampled ancestor tree\n simulated")
# select extant tips based on the fossil object
sampled_extant_tips = SAt$fossils$tip.label[SAt$fossils$hmin == 0]
# generate sampled tree from SAtree
SAt_sampled = FossilSim::sampled.tree.from.combined(tree = SAt$tree,
                                                    sampled_tips = sampled_extant_tips)
plot(SAt_sampled,
     main = "Sampled ancestor tree, sampled")
# simple Mk2 model
n_states = 2
br_rate = 0.5
trait_num = 10
trait_to_display = 2
morpho_mk = MorphSim::sim.morpho(k = n_states,
                                 time.tree = tree,
                                 br.rates = br_rate,
                                 trait.num = trait_num, 
                                 variable = FALSE,
                                 fossil = fossils2)
MorphSim::plotMorphoGrid(morpho_mk)
plot(morpho_mk,
     trait = trait_to_display)
title(main = paste0("Trait no. ", trait_to_display))
# export character matrix
MorphSim::write.morpho(data = morpho_mk,
                       file = "char_mat_temp.nex",
                       type = "matrix",
                       all = TRUE,
                       reconstructed = TRUE)
# export taxon ages
MorphSim::write.morpho(data = morpho_mk,
                       file = "ages_temp.tsv",
                       type = "ages",
                       all = TRUE,
                       reconstructed = TRUE)

# export the true tree 
MorphSim::write.morpho(data = morpho_mk,
                       file = "true_tree.tre",
                       type = "tree",
                       all = TRUE,
                       reconstructed = TRUE)
