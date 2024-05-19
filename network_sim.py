import argparse
import os
import random
from netop import Network
from treeop import Node, Tree, str2tree, getlabs
import subprocess


class TreeDepth(Tree):
    def __init__(self, tup):
        Tree.__init__(self, tup)
        self.root = NodeDepth(tup, None)
        self.nodes = self.root.nodes()
        for i, n in enumerate(self.nodes):
            n.num = i
            n.artificial = False
        if self.srclist:
            l = self.root
            for i in range(len(self.srclist) - 1):
                l.artificial = True
                l = l.l


class NodeDepth(Node):
    def __init__(self, tup, par):
        Node.__init__(self, tup, par)
        if self.parent:
            self.depth = self.parent.depth + self.branchlength
        else:
            self.depth = 0
        if self.src:
            self.c = [NodeDepth(t, self) for t in self.src]

    def __str__(self):
        if self.label == ":":
            dlabel = self.label + str(self.branchlength)
        elif self.label == "":
            dlabel = ""
        else:
            dlabel = self.label + ":" + str(self.branchlength)
        if self.leaf():
            return dlabel
        return "(" + ",".join(str(c) for c in self.c) + ")" + dlabel


def draw_branches(branches, checked):
    branch1, branch2 = None, None
    while branch1 == branch2 or ((branch1, branch2) in checked and (branch2, branch1) in checked):
        branch1 = branches[random.randint(1, len(branches) - 1)]  # from
        branch2 = branches[random.randint(1, len(branches) - 1)]  # to
    if (branch1, branch2) in checked:
        branch1, branch2 = branch2, branch1
        time1 = random.uniform(branch1.parent.depth, branch2.depth)
        time2 = random.uniform(time1, branch2.depth)
    elif (branch2, branch1) in checked:
        time1 = random.uniform(branch1.parent.depth, branch2.depth)
        time2 = random.uniform(time1, branch2.depth)
    else:
        b1_rg = (branch1.parent.depth, branch1.depth)
        b2_rg = (branch2.parent.depth, branch2.depth)
        time1 = random.uniform(b1_rg[0], b1_rg[1])
        time2 = random.uniform(b2_rg[0], b2_rg[1])
        if time1 > time2:
            branch2, branch1 = branch1, branch2
            time1, time2 = time2, time1
    return branch1, branch2, time1, time2


def insert_leaf(n, nodes, rlabel, time1, time2):
    np = n.parent
    ap = NodeDepth(([], [":", str(time1 - np.depth)]), np)
    np.c.remove(n)
    np.c.append(ap)
    nodes.append(ap)
    a = NodeDepth(([], [rlabel, ":", str(time2 - time1)]), ap)
    oldn = n.branchlength
    n.branchlength = n.depth - time1
    n.parent = ap
    ap.c = [a, n]
    nodes.append(a)
    return oldn, nodes, n, a, ap, np


def insert_internal(m, v, r, time2):
    mp = m.parent
    b = NodeDepth(([], [r, ":", str(time2 - m.parent.depth)]), mp)
    mp.c.remove(m)
    mp.c.append(b)
    oldm = m.branchlength
    m.branchlength = m.depth - time2
    b.c = [m]
    b.parent = mp
    m.parent = b
    v.append(b)
    return oldm, v, m, b, mp


def delete_reticulation(n, m, oldn, oldm, v, np, mp, a, b, ap):
    n.branchlength = oldn
    m.branchlength = oldm
    v.remove(b)
    m.parent = mp
    mp.c.remove(b)
    mp.c.append(m)
    v.remove(a)
    v.remove(ap)
    np.c.append(n)
    np.c.remove(ap)
    n.parent = np
    return n, m, np, mp, v


def add_reticulations(tree, rnum, repeat=False):
    adj = str2tree(tree)
    t = TreeDepth(adj)
    nodes = t.nodes.copy()
    reticulations = ['#' + i for i in getlabs(ord('A'), ord('Z'), rnum)]
    inserted = []  # reticulation ids
    checked = []  # pairs of edges
    npairs = (len(nodes) - 1) * (len(nodes) - 2)  # number of all pairs (from,to)
    while reticulations:
        r = reticulations.pop()
        inserted.append(r)
        while len(checked) < npairs:
            n, m, time1, time2 = draw_branches(nodes, checked)
            if n not in inserted:
                oldn, nodes, n, a, ap, np = insert_leaf(n, nodes, r, time1, time2)
                oldm, nodes, m, b, mp = insert_internal(m, nodes, r, time2)
                s = str(t)
                net = Network(str2tree(str(t)))
                if net.isdag() and net.treechild():
                    checked.append((n, m))
                    checked.append((m, n))
                    break
                n, m, np, mp, v = delete_reticulation(n, m, oldn, oldm, nodes, np, mp, a, b, ap)
                if (n, m) not in checked:
                    checked.append((n, m))
                if n.depth < m.parent.depth:
                    if (m, n) not in checked:
                        checked.append((m, n))

        else:
            if repeat:
                return add_reticulations(tree, rnum, repeat)
            else:
                raise ValueError("Can't add any more reticulations")
    return s


def simulate(leaves, ret, height):
    process = subprocess.Popen(f"Rscript treesim.r --l {leaves} --h {height}", shell=True, stdout=subprocess.PIPE)
    (output, error) = process.communicate()
    output = output.decode().split()
    tree_index = output.index("[1]") + 1
    newick = output[tree_index][1:-2]
    if newick[-1] == "0":
        newick = newick[:-2]
    if ret < leaves:
        network = add_reticulations(newick, ret, True)
    else:
        network = add_reticulations(newick, ret)
    return network


def save_nexus(newlist, path):
    with open(path, "w") as f:
        f.write("#NEXUS\nbegin trees;\n")
        for i in range(len(newlist)):
            f.write("\ttree " + str(i) + "=" + newlist[i] + ";\n")
        f.write("end;")


def parse_input():
    parser = argparse.ArgumentParser(
        description='Generate species networks with a given number of leaves and reticulations, and generate their displayed trees')
    parser.add_argument('-o', required=True, help='name of the output folder')
    parser.add_argument('-r', type=int, required=True, help='number of reticulations')
    parser.add_argument('-l', type=int, required=True, help='number of leaves in species network')
    parser.add_argument('-n', type=int, required=True, help='number of species networks to generate')
    parser.add_argument('-d', type=int, required=True, help='number of displayed trees to generate per each network')
    parser.add_argument('-ht', type=int, required=True, help='height of the network in generations')
    parsed = parser.parse_args()
    return parsed


def main():
    args = parse_input()
    os.mkdir(args.o)
    for n in range(1, args.n + 1):
        reppath = args.o + "/" + str(n)
        network_str = simulate(args.l, args.r, args.ht)
        os.mkdir(reppath)

        nwpath = reppath + "/" + "network"
        nwf = open(nwpath, "w")
        nwf.write(network_str)
        nwf.close()

        displayed_trees = []
        network = Network(str2tree(network_str))
        for i in range(args.d):
            dtree_nr = random.randint(0, 2 ** args.r - 1)
            dtree = network.displayedtreebyid(dtree_nr)
            displayed_trees.append(dtree)
        dpath = reppath + "/" + "displayed_trees"
        save_nexus(displayed_trees, dpath)


if __name__ == "__main__":
    main()
