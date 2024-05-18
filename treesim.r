library("TreeSim")
library(argparser, quietly=TRUE)

p = arg_parser("Simulate a species tree with a given height and leaf number. Other parameters are taken from Molloy and Warnow 2018")
p = add_argument(p, "--l", help="number of leaves", type="integer")
p = add_argument(p, "--h", help="tree height", type="float")
argv = parse_args(p)

br = 10**(-7)
dr = 0
tree = sim.bd.taxa.age(argv$l, 1, br, dr, 1, argv$h, T)
newick = write.tree(tree[[1]])
print(newick)
