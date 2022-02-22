numbsim1 <- as.numeric(commandArgs(trailingOnly=TRUE)[1])
ntips <- as.numeric(commandArgs(trailingOnly=TRUE)[2])
lambda <- as.numeric(commandArgs(trailingOnly=TRUE)[3])
mu <- as.numeric(commandArgs(trailingOnly=TRUE)[4])
nu <- as.numeric(commandArgs(trailingOnly=TRUE)[5])
hybpropsr <- commandArgs(trailingOnly=TRUE)[6]
num_gene_tree <- as.numeric(commandArgs(trailingOnly=TRUE)[7])
genseed <- as.numeric(commandArgs(trailingOnly=TRUE)[8])


print(numbsim1)
print(ntips)
print(lambda)
print(mu)
print(nu)
print(hybpropsr)
print(num_gene_tree)
print(genseed)
