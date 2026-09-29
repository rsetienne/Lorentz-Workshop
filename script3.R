library(ape)
library(MorphSim)
set.seed(12)
tree_list <- sim.bd.taxa(
  n        = 10,   # number of extant tips
  numbsim  = 1,    # number of trees to simulate
  lambda   = 1,    # speciation rate
  mu       = 0.5,  # extinction rate
  complete = TRUE # if set to false drop all extinct lineages and keep only extant tips
)
tree <- tree_list[[1]]
# total tip count, including BOTH extant and extinct lineages
# (this will be > 10, since n only controls extant tips)
Ntip(tree)
# this is the check that actually confirms n = 10 was enforced:
# count tips whose distance from the root equals the tree's max depth
# (i.e. tips sampled at the present)
node_depths <- node.depth.edgelength(tree)
tip_depths  <- node_depths[1:Ntip(tree)]
sum(abs(tip_depths - max(tip_depths)) < 1e-8)  
# with complete = TRUE, the tree should generally NOT be ultrametric,
# since extinct tips end before the present
is.ultrametric(tree)
# TreeSim always produces a strictly bifurcating tree so this
# should be TRUE regardless of complete = TRUE/FALSE
is.binary(tree)
plot(tree, cex = 0.8, no.margin = TRUE)
morpho_mk <- sim.morpho(k = 3, time.tree = tree, br.rates = 0.1, trait.num = 15)
plotMorphoGrid(morpho_mk)
morpho_mkv <- sim.morpho(k = 3, time.tree = tree, br.rates = 0.1,
                         trait.num = 15, variable = TRUE)
plotMorphoGrid(morpho_mkv)
morpho_full <- sim.morpho(k = 3, time.tree = tree, br.rates = 0.1,
                          trait.num = 15, variable = TRUE, full.states = TRUE)
plotMorphoGrid(morpho_full)
# tolerate autapomorphies
morpho_pars_standard <- sim.morpho(k = 3, time.tree = tree, br.rates = 0.1,
                                   trait.num = 15, variable = TRUE,
                                   parsimony = "standard")
plotMorphoGrid(morpho_pars_standard)
# every state must occur in 2+ taxa = stricter
morpho_pars_strict <- sim.morpho(k = 3, time.tree = tree, br.rates = 0.1,
                                 trait.num = 15, variable = TRUE,
                                 parsimony = "strict")
plotMorphoGrid(morpho_pars_strict)
morpho_acrv <- sim.morpho(k = 3, time.tree = tree, br.rates = 0.1, trait.num = 15,
                          variable = TRUE, ACRV = "gamma", alpha.gamma = 1, ACRV.ncats = 4)

plotMorphoGrid(morpho_acrv)
ord_Q <- matrix(c(-0.5, 0.5, 0.0,
                  0.3333333, -0.6666667, 0.3333333,
                  0.0, 0.5, -0.5),
                nrow = 3, byrow = TRUE)
morpho_ordered <- sim.morpho(k = 3, time.tree = tree, br.rates = 0.1, trait.num = 15,
                             variable = TRUE, ACRV = "gamma", alpha.gamma = 1,
                             ACRV.ncats = 4, define.Q = ord_Q)

plot(morpho_ordered, trait = 1, box.cex = 2)
part_data <- sim.morpho(
  time.tree  = tree,
  k          = c(2, 3, 4),
  trait.num  = 20,
  partition  = c(10, 5, 5),
  br.rates   = 0.1
)

# model specification per partition
part_data$model$Specified
# partition 1: ordered characters with a custom Q matrix
partition_1 <- sim.morpho(time.tree = tree, k = 3, trait.num = 10,
                          br.rates = 0.1, define.Q = ord_Q)

# partition 2: unordered binary characters with MkV + gamma ACRV
partition_2 <- sim.morpho(time.tree = tree, k = 2, trait.num = 10,
                          br.rates = 0.1, variable = TRUE,
                          ACRV = "gamma", alpha.gamma = 1, ACRV.ncats = 4)

combined <- combine.morpho(partition_1, partition_2)

# 20 traits total, model info from both partitions preserved
length(combined$sequences$tips[[1]])
combined$model$Specified
# hard tissue characters (partition 1): 10% missing
# soft tissue characters (partition 2): 70% missing
missing_part <- sim.missing.data(
  data = combined, method = "partition",
  seq = "tips", probability = c(0.1, 0.7)
)
missing_extinct <- sim.missing.data(
  data = combined , method = "extinct",
  seq = "tips", probability = 0.5
)
# full tree
write.morpho(missing_part, file = "tree.tre", type = "tree",
             reconstructed = FALSE)
# character matrix
write.morpho(missing_part, file = "matrix.nex", type = "matrix", reconstructed = FALSE)
# fossil occurrence ages (compatible with RevBayes)
write.morpho(missing_part, file = "ages.tsv", type = "ages")
# ages with associated uncertainty
write.morpho(missing_part, file = "ages.tsv", type = "ages",
             uncertainty = 2)