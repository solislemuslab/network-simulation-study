library(phylolm)

######
# 0.05
######

setwd("0.05")

quartetCF = read.csv("CFs_cropped_tree0.05.newick_gene_trees_coal_unit.csv")
dat = quartetCF[, c(1:4, 5, 8, 11)]
for (i in 1:4){ dat[,i] = factor(dat[,i])}

# this comes from the best net in snaq_output_h0_CFs_cropped_tree0.05.newick_gene_trees_coal_unit.csv.out
h0string = "(t10_1,((t18_1,(((t12_1,t1_1):0.23424382680566322,t23_1):0.46713929634800333,((t2_1,(t5_1,t6_1):0.16132571918236527):0.009386437453449109,(t8_1,t9_1):0.15407167604261643):0.4201490993715868):0.12205618160098526):0.15868257483839623,t27_1):0.6585743839145971,t20_1);"
h0tree = read.tree(text = h0string)

h0treeprelim = test.tree.preparation(dat, h0tree)
Ntaxa = length(h0tree$tip.label) # 15 of course
h0tree.internal.edges = which(h0tree$edge[,2] > Ntaxa)
h0tree.internal.edges # indices of internal edges: those that need a length
res <- test.one.species.tree(dat,h0tree,h0treeprelim,edge.keep=h0tree.internal.edges)

######
# 0.10
######

setwd("../0.1")

quartetCF = read.csv("CFs_cropped_tree0.1.newick_gene_trees_coal_unit.csv")

######
# 0.50
######

setwd("../0.5")

quartetCF = read.csv("CFs_cropped_tree0.5.newick_gene_trees_coal_unit.csv")

######
# 1.00
######

setwd("../1.0")

quartetCF = read.csv("CFs_cropped_tree1.0.newick_gene_trees_coal_unit.csv")

######
# 2.00
######

setwd("../2.0")

quartetCF = read.csv("CFs_cropped_tree2.0.newick_gene_trees_coal_unit.csv")
