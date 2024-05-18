import sys
from itertools import product
from collections import Counter
from random import randint, shuffle
from timeit import default_timer as timer
from treeop import Tree, Node, str2tree, getlabs, randtreestr, compcostsmp

RET_LEFTUSED = 1
RET_RIGHTUSED = 2
RET_BOTHUSED = RET_LEFTUSED | RET_RIGHTUSED

RET_C1 = 1
RET_C2 = 2
RET_C3 = 4
RET_C4 = 8
RET_C5 = 16
RET_CMASK = 31


class Network(Tree):
    def __init__(self, tup):
        Tree.__init__(self, tup)

        # recognize reticulations
        dlf = {}
        din = {}
        err = 0
        for i, n in enumerate(self.nodes):
            if n.label and n.label[0] == "#":
                retid = n.label[1:]
                n.reticulation = 1
                if n.leaf():
                    if retid in dlf:
                        print("Reticulation id <%s> already defined" % retid)
                        err = 1
                    dlf[retid] = n
                    n.reticulationleaf = 1
                else:
                    if retid in din:
                        print("Reticulation id <%s> already defined" % retid)
                        err = 1
                    din[retid] = n
                n.retid = retid
            else:
                n.reticulation = 0
                n.retid = 0

        if set(dlf.keys()) != set(din.keys()) or err:
            for k in din:
                if k not in dlf:
                    print("Missing reticulation label <%s> in leaves" % k)
            for k in dlf:
                if k not in din:
                    print("Missing reticulation label <%s> in internal nodes" % k)

            sys.exit(-1)

        self.reticulations = []
        for retnum, retid in enumerate(dlf, 1):
            l = dlf[retid]
            i = din[retid]
            lpar = l.parent
            ipar = i.parent
            lpar.c.remove(l)
            ipar.c.remove(i)

            lpar.c.insert(0, i)  # inserted at 0
            ipar.c.insert(0, i)  # inserted at 0

            i.lftparent = ipar
            i.rghparent = lpar
            i.retnum = retnum

            if i.branchlengthset and l.branchlengthset:  # set only if both are defined
                i.branchlengthset = True
                i.branchlength = (i.branchlength, l.branchlength)
            else:
                i.branchlengthset = False

            self.reticulations.append(i)
            self.nodes.remove(l)

        for n in self.nodes:
            n._setcluster()  # reconstruct clusters

    def isdag(self):
        return sorttop([(n.num, c.num) for n in self.nodes for c in n.c])

    def istimeconsistent(self):
        lftparent_nums = []
        edges = []

        # glue reticulation parents
        for r in self.reticulations:
            lftparent_nums.append(r.lftparent.num)
            r.lftparent.num = r.rghparent.num

        # create edge representation of graph
        for n in self.nodes:
            if n == self.root:
                pass
            elif n.reticulation:
                edges.append((n.lftparent.num, n.num))
            else:
                edges.append((n.parent.num, n.num))

        # unglue reticulation parents
        for i, n in enumerate(self.reticulations):
            n.lftparent.num = lftparent_nums[i]

        return bool(sorttop(edges))

    def ptab(self):
        print("=" * 80)
        for n in self.nodes:
            print(n.num, n.reticulation, n.retid, end='')
            if n.reticulation:
                if hasattr(n, "lftparent"):
                    print("rtp:%d,%d" % (n.lftparent.num, n.rghparent.num), end='')
            elif n.parent:
                print("par:%d" % n.parent.num, end='')
            print('==<<', n.netrepr(), ">>", end='')
            if n.leaf():
                print("LEAF", end='')
            else:
                print("Children=", end='')
                for c in n.c:
                    print(c.num, end=' ')
            print()

    def todotfile(self, f, nodeprefix=""):
        for n in self.nodes:
            comments = "\n".join(n.comments)
            if comments:
                comments = "\n" + comments
            f.write(nodeprefix)
            if n.leaf():
                f.write("%d [label=\"%s %d%s\"];\n" % (n.num, n, n.num, comments))
            else:
                if n.reticulation:
                    f.write("%d [shape=box,color=red, label=\"%d #%s%s\"];\n"
                            % (n.num, n.num, n.retid, comments))
                else:
                    f.write("%d [label=\"%d%s\"];\n" % (n.num, n.num, comments))
        for n in self.nodes:
            if n.reticulation:
                if hasattr(n, "lftparent"):
                    f.write("%s%d -> %s%d [color=green,label=\"%sl\"];\n" % (
                        nodeprefix, n.lftparent.num, nodeprefix, n.num, n.retid))
                    f.write("%s%d -> %s%d [color=blue,label=\"%sr\"];\n" % (
                        nodeprefix, n.rghparent.num, nodeprefix, n.num, n.retid))
            elif n.parent:
                f.write("%s%d -> %s%d;\n" % (nodeprefix, n.parent.num, nodeprefix, n.num))

    def get_leaves(self):
        """Get nodes of 1-indegree and 0-outdegree."""
        return [node for node in self.nodes if node.leaf()]

    def get_labels(self):
        """Get labels of all leaves."""
        return [node.label for node in self.get_leaves()]

    def get_inner_nodes(self):
        """Get all nodes (tree and reticulation) except of leaves."""
        return [node for node in self.nodes if not node.leaf()]

    def get_inner_tree_nodes(self):
        """Get nodes of 1-indegree and 2-outdegree plus root."""
        return [node for node in self.get_inner_nodes() if not node.reticulation]

    def valid_binary(self):
        """
        Sanity check whether network is DAG, bijective, binary nodes and
        reticulations, no parallel edges, unique root. Some of these properties
        are partially assured by the class constructor.

        Returns:
            bool for satisfying conditions
        """

        # check if network directed acyclic
        if not self.isdag():
            return False

        # check if leaves are bijective
        labels = self.get_labels()
        if len(labels) != len(set(labels)):
            return False

        # check if binary
        reticulations = self.reticulations
        for node in reticulations:
            if len(node.c) != 1:
                return False

        inner_tree_nodes = self.get_inner_tree_nodes()
        for node in inner_tree_nodes:
            if len(node.c) != 2:
                return False

        # check for parallel edges
        for node in reticulations:
            if node.lftparent == node.rghparent:
                return False

        for node in inner_tree_nodes:
            if node.c[0] == node.c[1]:
                return False

        return True

    def treechild(self):
        """
        Check if network is a valid binary network and belongs to the
        Tree-Child class.

        Returns:
            bool for satisfying conditions
        """

        # check if valid network

        if not self.valid_binary():
            return False

        # check if each inner node has a tree or leaf child node

        for node in self.nodes:
            if not node.leaf() and all(child in self.reticulations for child in node.c):
                return False

        return True

    def type1net(self):
        """
        Check if network is a valid binary network such that
        no node has >= two reticulation parents

        Returns:
            bool for satisfying conditions
        """

        # check if valid network

        if not self.valid_binary():
            return False

        # check if each inner node has a tree or leaf child node

        for node in self.nodes:
            if not node.leaf() and len(node.c) > 1 and all(child in self.reticulations for child in node.c):
                return False

        return True

    def __repr__(self):
        return self.root.netrepr()

    def __str__(self):
        return self.root.netrepr()

    def robinsonfouldsnet(self, net):
        """
        Compute network version of RF distance between two networks.
        Normalizing divider: 2(n-2)*2^k
        """

        disp_trees1 = [Tree(str2tree(t)) for t in self.displayedtrees()]
        disp_trees2 = [Tree(str2tree(t)) for t in net.displayedtrees()]

        counter1 = Counter([n.cluster for tree in disp_trees1 for n in tree.nodes])
        counter2 = Counter([n.cluster for tree in disp_trees2 for n in tree.nodes])

        dif1, dif2 = counter1 - counter2, counter2 - counter1

        return sum(dif1.values()) + sum(dif2.values())

    def displayedtreebyid(self, displayedtreeid):
        """
        Return display tree via id
        """
        ign = 0
        fmt = "{0:0%db}" % len(self.reticulations)  # format with leading zeros

        try:
            return self._bltree(
                *self._displtree(dict(zip(self.reticulations, map(int, fmt.format(displayedtreeid)))), self.root, None))
        except Exception:
            ign += 1

        if ign:
            print(f"{ign} tree(s) ignored. Is our network tree-child?", file=sys.stderr)

    # insert bl in present
    def _bltree(self, t, blset, branchlength):
        return t + (f":{branchlength}" if blset else "")

    # generic display tree generator via dictionary of ret. usages
    # given vector of reticulation usages b, a node n and parent node par from which n is reached
    # returns a tuple: (display tree encoded in string, branchlengthset_flag, branchlength )
    def _displtree(self, b, n, par):

        if n.branchlengthset:
            bl = n.branchlength
        else:
            bl = 0

        if n.leaf():
            return n.label, n.branchlengthset, bl

        if n in b:  # reticulation;  aggregate branch lengths
            t = None
            if b[n]:
                if n.lftparent == par:
                    r = self._displtree(b, n.c[0], n)
                    if not r:
                        return None
                    t, blsetc, blc = r
                    if n.branchlengthset:
                        bl = bl[0] + blc
            else:
                if n.rghparent == par:
                    r = self._displtree(b, n.c[0], n)
                    if not r:
                        return None
                    t, blsetc, blc = r
                    if n.branchlengthset:
                        bl = bl[1] + blc

            if t:
                if n.branchlengthset and blsetc:
                    return t, True, bl
                return t, False, 0

            return None
        else:
            l = [self._displtree(b, c, n) for c in n.c]
            l = [s for s in l if s]

            if not l:
                return ''

            if len(l) == 1:
                t, blsetc, blc = l[0]
                if n.branchlengthset and blsetc:
                    return t, True, bl + blc
                return t, False, 0
            s = ",".join(self._bltree(*t) for t in l)

            return "(" + s + ")", n.branchlengthset, bl

    def displayedtrees(self):

        ign = 0

        for b in product([0, 1], repeat=len(self.reticulations)):
            try:
                yield self._bltree(*self._displtree(dict(zip(self.reticulations, b)), self.root, None))
            except Exception:
                ign += 1
        if ign:
            print(f"{ign} tree(s) ignored. Is your network tree-child?", file=sys.stderr)

    def retmindup(self, gt, **opt):
        # custom multiplier for custom hashing
        hash_val = max(len(gt.nodes), len(self.nodes) + len(self.reticulations)) + 1

        # with the way we do hashes max value we store is smaller than 2 * (hash_val ** 2)
        max_size = 2 * hash_val ** 2

        # Introduce custom hashing for tuples
        # Our node numbers will be small, so we can introduce no-collision tuple hashing
        # This speeds up dict operations immensely
        # Assumption - hash_val > n1 and hash_val > n2
        def hash_tuple(n1, n2, n3=-1):
            if n3 == -1:
                return n1 * hash_val + n2
            # since n3 is either 0 or 1, changing hashing order to n3, n1, n2 reduces hash space
            return n3 * hash_val ** 2 + n1 * hash_val + n2

        # initialize lists with None, since our costs may be equal to zero
        deltav = [None] * max_size
        deltaretusage = [None] * max_size
        deltaupv = [None] * max_size
        deltaupretusage = [None] * max_size

        INFTY = 1e100

        # find mapping of leaves from gt to self
        labels = {}
        for leaf in self.leaves():
            labels[leaf.label] = leaf

        for n in gt.leaves():
            if n.label in labels:
                n.map = labels[n.label]
            else:
                raise Exception("Gene tree leaf", n, "cannot be mapped to species tree")

        def delta(tree_node, net_node):
            hsh = hash_tuple(tree_node.num, net_node.num)
            if deltav[hsh]:
                return deltav[hsh], deltaretusage[hsh]

            if tree_node.leaf():
                retusage = 0
                if tree_node.map == net_node:
                    res = 0
                else:
                    if net_node.reticulation:
                        res, retusage = delta(tree_node, net_node.c[0])
                    else:
                        res = INFTY
            else:
                tree_c0, tree_c1 = tree_node.c[0], tree_node.c[1]
                res0, retusage0 = delta(tree_c0, net_node)
                res1, retusage1 = deltaup(tree_c1, net_node)
                res = res0 + res1 + 1
                retusage = retusage0 | retusage1

                res0, retusage0 = delta(tree_c1, net_node)
                res1, retusage1 = deltaup(tree_c0, net_node)
                if res > res0 + res1 + 1:
                    res = res0 + res1 + 1
                    retusage = retusage0 | retusage1

                if not net_node.leaf():
                    if not net_node.reticulation:
                        applied = False
                        res0, retusage0 = deltaup(tree_c0, net_node.c[0])
                        res1, retusage1 = deltaup(tree_c1, net_node.c[1])
                        if res0 + res1 < res:
                            res = res0 + res1
                            retusage = retusage0 | retusage1
                            applied = True
                        res0, retusage0 = deltaup(tree_c1, net_node.c[0])
                        res1, retusage1 = deltaup(tree_c0, net_node.c[1])
                        if res0 + res1 < res:
                            res = res0 + res1
                            retusage = retusage0 | retusage1
                            applied = True
                        if applied:
                            for ch in net_node.c:
                                if ch.reticulation:
                                    if ch.lftparent == net_node:
                                        retusage |= ch.left_reticulation_used()
                                    else:
                                        retusage |= ch.right_reticulation_used()
                    else:
                        res0, retusage0 = deltaup(tree_c0, net_node.c[0])
                        res1, retusage1 = deltaup(tree_c1, net_node.c[0])
                        if res > res0 + res1 + 1:
                            res = res0 + res1 + 1
                            retusage = retusage0 | retusage1
                else:
                    res = INFTY
                    retusage = 0
            deltav[hsh] = res
            deltaretusage[hsh] = retusage
            return res, retusage

        def deltaup(tree_node, net_node):
            hsh = hash_tuple(tree_node.num, net_node.num)
            if deltaupv[hsh]:
                return deltaupv[hsh], deltaupretusage[hsh]

            res, retusage = delta(tree_node, net_node)

            if net_node.reticulation:
                res0, retusage0 = deltaup(tree_node, net_node.c[0])
                if res > res0:
                    res = res0
                    retusage = retusage0
            elif not net_node.leaf():
                net_c0, net_c1 = net_node.c
                res0, retusage0 = deltaup(tree_node, net_c0)
                res1, retusage1 = deltaup(tree_node, net_c1)
                res01 = min(res0, res1)

                if res > res01:
                    # min achieved for a child s0 or s1
                    used = notused = None
                    res = res01
                    if res == res0:
                        if net_c0.reticulation:
                            used = net_c0
                        if net_c1.reticulation:
                            notused = net_c1
                        retusage = retusage0
                    elif res == res1:
                        if net_c1.reticulation:
                            used = net_c1
                        if net_c0.reticulation:
                            notused = net_c0
                        retusage = retusage1
                    if used:
                        if net_node == used.lftparent:
                            retusage |= used.left_reticulation_used()
                        else:
                            retusage |= used.right_reticulation_used()

                    if notused:
                        if net_node == notused.lftparent:
                            retusage |= notused.right_reticulation_used()  # right must be used
                        else:
                            retusage |= notused.left_reticulation_used()  # left must be used

            deltaupv[hsh] = res
            deltaupretusage[hsh] = retusage
            return res, retusage

        res = min(delta(gt.root, s)[0] for s in self.nodes)

        opt_root = None
        opt_map = None

        for s in self.nodes:
            hsh = hash_tuple(gt.root.num, s.num)
            if deltav[hsh] == res:
                opt_root = s
                opt_map = {}
                for retic in self.reticulations:
                    opt_map[retic] = retic.get_node_retusage(deltaretusage[hsh])
                break

        if 'additional_return' in opt:
            return res, opt_root, opt_map

        return res

    # DP main for DC-non-classic
    def retmindc(self, gt, **opt):

        # custom multiplier for custom hashing
        hash_val = max(len(gt.nodes), len(self.nodes) + len(self.reticulations)) + 1

        # with the way we do hashes max value we store is smaller than 2 * (hash_val ** 2)
        max_size = 2 * hash_val ** 2

        # initialize lists with None, since our costs may be equal to zero
        deltav = [None] * max_size
        deltaretusage = [None] * max_size
        deltaupv = [None] * max_size
        deltaupoptval = [None] * max_size
        deltaupretusage = [None] * max_size

        INFTY = 1e100
        # find mapping of leaves from gt to self
        labels = {}
        for leaf in self.leaves():
            labels[leaf.label] = leaf

        for n in gt.leaves():
            if n.label in labels:
                n.map = labels[n.label]
            else:
                raise Exception("Gene tree leaf", n, "cannot be mapped to species tree")

        # Introduce custom hashing for tuples
        # Our node numbers will be small, so we can introduce no-collision tuple hashing
        # This speeds up dict operations immensely
        # Assumption - hash_val > n1 and hash_val > n2
        def hash_tuple(n1, n2, n3=-1):
            if n3 == -1:
                return n1 * hash_val + n2
            # since n3 is either 0 or 1, changing hashing order to n3, n1, n2 reduces hash space
            return n3 * hash_val ** 2 + n1 * hash_val + n2

        def delta(g, s):
            hsh = hash_tuple(g.num, s.num)
            if deltav[hsh]:
                return deltav[hsh], deltaretusage[hsh]

            if g.leaf():
                retusage = 0
                if g.map == s:
                    res = 0
                else:
                    res = INFTY
            else:
                res0, retusage0 = deltaup(g.c[0], s, 1)
                res1, retusage1 = deltaup(g.c[1], s, 1)
                res = res0 + res1
                retusage = retusage0 | retusage1

            deltav[hsh] = res
            deltaretusage[hsh] = retusage
            return res, retusage

        def deltaup(g, s, first):
            hsh = hash_tuple(g.num, s.num, first)
            if deltaupv[hsh]:
                return deltaupv[hsh], deltaupretusage[hsh]

            if len(s.c) == 1:  # s is a reticulation
                sc = s.c[0]
                res, retusage = deltaup(g, sc, 0)
                if sc.reticulation:  # special new case nonTC1
                    if s == sc.lftparent:
                        retusage |= sc.left_reticulation_used()
                    else:
                        retusage |= sc.right_reticulation_used()
                    optval = ("D^", sc)
                    # end of special case nonTC1
                else:
                    res += 1
                    optval = ("D^", sc)

            elif s.leaf():

                res, retusage = delta(g, s)
                optval = ("D", s)

            else:
                res, retusage = delta(g, s)
                optval = ("D", s)
                # s has 2 children s0 and s1
                s0, s1 = s.c

                # two variants
                du0, retusage0 = deltaup(g, s0, 0)
                du1, retusage1 = deltaup(g, s1, 0)

                if first:
                    # ignore edge <s,s0/s1> if s0/s1 is a reticulation
                    res0 = 1 - s0.reticulation + du0
                    res1 = 1 - s1.reticulation + du1
                    res = min(res, res0, res1)
                    sc = None
                    if res == res0:
                        sc, retusagec, optval = s0, retusage0, ("D^", s0)
                    elif res == res1:
                        sc, retusagec, optval = s1, retusage1, ("D^", s1)

                    # res is min
                    # problem: multiple min's
                    if sc:
                        # min is achieved by a kid
                        if sc.reticulation:
                            # kid is a reticulation node
                            retusage = retusagec
                            if s == sc.lftparent:
                                retusage |= sc.left_reticulation_used()
                            else:
                                retusage |= sc.right_reticulation_used()
                        else:
                            retusage = retusagec
                else:
                    res0, res1 = du0, du1
                    if s0.reticulation == s1.reticulation == 0:
                        # adjust if a child is not reticulation
                        res0 += 1
                        res1 += 1
                    elif s0.reticulation and s1.reticulation:
                        # not allowed in tree-child networks
                        raise Exception("Tree child network expected")

                    res01 = min(res0, res1)
                    if res > res01:
                        # min achieved for a child s0 or s1
                        used = notused = None
                        res = res01
                        if res == res0:
                            optval = ("D^", s0)
                            if s0.reticulation:
                                used = s0
                            if s1.reticulation:
                                notused = s1
                            retusage = retusage0
                        elif res == res1:
                            optval = ("D^", s1)
                            if s1.reticulation:
                                used = s1
                            if s0.reticulation:
                                notused = s0
                            retusage = retusage1
                        if used:
                            if s == used.lftparent:
                                retusage |= used.left_reticulation_used()
                            else:
                                retusage |= used.right_reticulation_used()

                        if notused:
                            if s == notused.lftparent:
                                retusage |= notused.right_reticulation_used()  # right must be used
                            else:
                                retusage |= notused.left_reticulation_used()  # left must be used

            deltaupv[hsh] = res
            deltaupretusage[hsh] = retusage
            deltaupoptval[hsh] = optval

            return res, retusage

        def ncet(v):
            c = ''
            if v & RET_LEFTUSED:
                c += "l"
            if v & RET_RIGHTUSED:
                c += "r"
            if c:
                return c
            return '.'

        def nc(dr):
            if dr:
                return "[" + ", ".join("%s%s" % (r.retid, ncet(dr[r])) for r in dr) + "]"
            return " []"

        def pptree(t, ar):
            add = " num=%d" % t.num
            if ar:
                tb = ''
                for g, s in sorted(deltav):
                    if g == t and deltav[g, s] < INFTY:
                        tb += " s%d=%d%s" % (s.num, deltav[g, s], nc(deltaretusage[g, s]))
                add += " tabv='%s'" % tb.strip()
                if not t.leaf():
                    tb = "|"
                else:
                    tb = ''
                sep = "|" if t.leaf() else " "
                tbs = ''
                for g, s, i in sorted(deltaupv):
                    if g == t and deltaupv[g, s, i] < INFTY:
                        if not t.leaf() and len(tbs) > 70:
                            tb += tbs + "|"
                            tbs = '  '
                        tbs += "%ss%s:%d.%d%s" % (sep, s.num, i, deltaupv[g, s, i], nc(deltaupretusage[g, s, i]))
                tb += tbs
                add += " tabup='%s'" % tb.strip()

            if t.leaf():
                return t.label + add
            return "(" + ",".join(pptree(c, ar) for c in t.c) + ")" + add

        def ppstree(t):
            add = " num=%d" % t.num
            if t.parent:
                add += " :%d" % (t.dagdepth() - t.parent.dagdepth())

            if t.leaf():
                return t.smplabel() + add
            return ""

        # optimal cost
        res = min(delta(gt.root, s)[0] for s in self.nodes)

        resup = res - len(gt.nodes) + 1  # DC classic

        if "dpdebug" in opt:
            print("&g " + pptree(gt.root, 0) + " dcdp=%d" % res + " dcclassic=%d" % resup)

            for t in gt.nodes:
                print(t.num, t)
                for g, s in sorted(deltav):
                    if g == t and deltav[g, s] < INFTY:
                        print("   D : %s s%d=%d %s" % (s, s.num, deltav[g, s], nc(deltaretusage[g, s])))

                for g, s, i in sorted(deltaupv):
                    if g == t and deltaupv[g, s, i] < INFTY:
                        print("   D^: %s s%d.%d=%d %s" % (s, s.num, i, deltaupv[g, s, i], nc(deltaupretusage[g, s, i])),
                              end='')
                        if not g.leaf() and deltaupoptval[g, s, i]:
                            print(" from: %s s%s" % (deltaupoptval[g, s, i][0], deltaupoptval[g, s, i][1].num), end='')
                        print()

        opt_root = None
        opt_map = None

        # print RET usage
        for s in self.nodes:
            hsh = hash_tuple(gt.root.num, s.num)
            if deltav[hsh] == res:
                opt_root = s
                opt_map = {}
                for retic in self.reticulations:
                    opt_map[retic] = retic.get_node_retusage(deltaretusage[hsh])
                if "optimalmap" in opt:
                    print("Optimal map: s%d dc=%d retusage=%s" % (s.num, res, nc(deltaretusage[gt.root, s])))

        # single print of dp arrays
        def pdparrays():
            for g, s in sorted(deltav):
                if deltav[g, s] < INFTY:
                    print(" D[", g, "-->", s.num, "]=", deltav[g, s], deltaretusage[g, s], file=sys.stderr)

            for g, s, i in sorted(deltaupv):
                if deltaupv[g, s, i] < INFTY:
                    print("D^[", g, "-->", s.num, i, "]=", deltaupv[g, s, i], deltaupretusage[g, s, i], file=sys.stderr)

            for s in self.nodes:
                if res == delta(gt.root, s):
                    print("#retDC", gt.root, "-->", s.num, " RET:", deltaretusage[gt.root, s])

        # add info to comments
        if 'dotspec' in opt:
            for gid in opt['dotspec']:
                for g in gt.nodes:
                    if g.num == gid:
                        for s in self.nodes:
                            if (g, s) in deltav:
                                if deltav[g, s] < INFTY:
                                    s.comments.append("D[%d]=%d %s" % (gid, deltav[g, s], nc(deltaretusage[g, s])))
                            for i in (0, 1):
                                if (g, s, i) in deltaupv and deltaupv[g, s, i] < INFTY:
                                    s.comments.append(
                                        "D^[%d,%d]=%d %s" % (gid, i, deltaupv[g, s, i], nc(deltaupretusage[g, s, i])))

        if 'additional_return' in opt:
            return res, opt_root, opt_map

        return res

    def branch_and_bound(self, g, use_already_set_nodes=False, branch_limit=1000000000000,
                         return_stats=False, cost_func='C'):
        """
        Find the best, complete solution for minimal displayed tree problem
        Dynamic programming solution (mindc method) allows some nodes to use both parents in tree embedding
        We use branch-and-bound method to set one parent for a non-resolved node and recurse to get a full solution
        Returns the minimal cost
        """

        def set_parent(node, parent, opt_root):
            """
            Set parent of reticulation node to parent
            Contract the network in place, so that it doesn't have unnecessary nodes and edges
            """

            def contract(node, parent, opt_root):
                """
                After choosing a parent for a node we might end up with internal nodes with in- and out- degree equal 1
                Those nodes should be contracted, which might lead to another node needing contraction and so on
                Invariant: node has a single child (reticulation node or we've already removed one of nodes' children)
                """
                if parent is None or node not in opt_root.nodes():
                    return
                i = parent.c.index(node)
                parent.c.remove(node)

                if not hasattr(node, 'removed') and node.parent == parent:
                    # If a parent is a "proper" one
                    # And we haven't marked the node as removed (only done in else below)
                    # We can simply contract the node and reroute its only child
                    if node.c[0].reticulation:
                        if node.c[0].lftparent == node:
                            node.c[0].lftparent = node.parent
                        else:
                            node.c[0].rghparent = node.parent
                    else:
                        node.c[0].parent = node.parent

                    node.parent.c.insert(i, node.c[0])
                else:
                    if parent.reticulation:
                        # If the parent is a reticulation we need to remove it completely
                        # Thus, setting the removed attribute and recursing for both  parents
                        parent.removed = True
                        contract(parent, parent.lftparent, opt_root)
                        contract(parent, parent.rghparent, opt_root)
                        del parent.removed

                    else:
                        # Parent node now has a single child, so we need to contract further
                        contract(parent, parent.parent, opt_root)

            node.parent = parent
            contract(node, node.rghparent, opt_root)
            contract(node, node.lftparent, opt_root)

        def copy(n):
            d = Network(str2tree(n.netrepr()))
            return d

        start_time = timer()
        assert branch_limit >= 0
        # Pick starting cost by comparing with one of displayed trees
        # Variable holds a cost for a complete tree
        best_cost = None
        for t in self.displayedtrees():
            best_cost = compcostsmp(g, Tree(str2tree(t)), cost_func)
            break
        if best_cost is None:
            raise Exception('No displayed trees found in network.')
        branching_count = 0
        max_depth = 0
        best_depth = 0
        dp_called = 0
        # BFS-like traversal, browsing reticulations from top to bottom
        queue = [(copy(self.root), branch_limit)]
        while queue:
            u, branch_left = queue.pop(0)
            if cost_func == 'C':
                cost, _, opt_map = u.retmindc(g, additional_return=True)
            elif cost_func == 'D':
                cost, _, opt_map = u.retmindup(g, additional_return=True)
            else:
                raise Exception("Unknown cost, expected C or D")
            opt_root = u.root

            depth = branch_limit - branch_left
            dp_called += 1
            max_depth = max(max_depth, depth)

            # If a cost for a potentially incomplete tree is higher than the one we achieved,
            # we can stop at that point
            if cost >= best_cost:
                continue

            # Solution has all reticulation parents set, we can use the cost
            if not any([value == RET_BOTHUSED for value in opt_map.values()]) or branch_left == 0:
                best_cost = cost
                best_depth = depth
                continue

            # Heuristic - we try to use as much as possible from our incomplete solution
            # to reduce computation time
            if use_already_set_nodes:
                # Set nodes that already have their parents known
                for retic, val in opt_map.items():
                    if val == RET_RIGHTUSED:
                        set_parent(retic, retic.rghparent, opt_root)
                    elif val == RET_LEFTUSED:
                        set_parent(retic, retic.lftparent, opt_root)

            # Since our opt_root might not be the initial root of the tree,
            # we need to set parents of reticulations
            # that only have one parent in a new subnetwork
            for retic, val in opt_map.items():
                if val == RET_BOTHUSED:
                    if retic.rghparent not in opt_root.nodes():
                        set_parent(retic, retic.lftparent, opt_root)
                        opt_map[retic] = RET_LEFTUSED
                    elif retic.lftparent not in opt_root.nodes():
                        set_parent(retic, retic.rghparent, opt_root)
                        opt_map[retic] = RET_RIGHTUSED

            branching_count += 1
            cp1 = copy(opt_root)
            cp2 = copy(opt_root)

            # If parent is not known, set both values in different copies and branch
            for retic, val in opt_map.items():
                if val == RET_BOTHUSED:
                    # find retic in cp1
                    x = [str(i) for i in cp1.nodes]
                    if str(retic) in x:
                        node = cp1.nodes[x.index(str(retic))]
                        set_parent(node, node.rghparent, cp1.root)
                        break

            for retic, val in opt_map.items():
                if val == RET_BOTHUSED:
                    x = [str(i) for i in cp2.nodes]
                    if str(retic) in x:
                        node = cp2.nodes[x.index(str(retic))]
                        set_parent(node, node.lftparent, cp2.root)
                        break
            # Workaround for cases when root has only one child
            if len(cp1.root.c) == 1:
                cp1.root = cp1.root.c[0]
            if len(cp2.root.c) == 1:
                cp2.root = cp2.root.c[0]

            queue.append((copy(cp1.root), branch_left - 1))
            queue.append((copy(cp2.root), branch_left - 1))
        end_time = timer()
        if return_stats:
            stats = {
                'dp_called': dp_called,
                'branching_count': branching_count,
                'max_depth': max_depth,
                'best_depth': best_depth,
                'time_in_seconds': end_time - start_time
            }
            return best_cost, stats
        return best_cost


# input: list of directed edges, must be non-empty
# return: None if cycle or empty graph
#         top. sort otherwise
def sorttop(e):
    l = []
    d1 = {}
    d2 = {}
    for x, y in e:
        d1.setdefault(x, []).append(y)
        d2.setdefault(y, []).append(x)
    roots = list(set(d1.keys()).difference(d2.keys()))
    while roots:
        n = roots.pop()
        l.append(n)
        if n not in d1:
            continue
        ms = d1[n][:]
        for m in ms:
            d1[n].remove(m)
            d2[m].remove(n)
            if not d2[m]:
                roots.append(m)
                d2.pop(m)
            if not d1[n]:
                d1.pop(n)
    if d1 or d2:
        return None
    return l


def addretstr(s, reticulations, skip=0, networktype=0, time_consistent=False):
    """
    Exhaustively add given number of reticulations to a tree or tree
    representation of network.

    networktype=0 -> treechild
    networktype=1 -> nontreechild type 1
    networktype=2 -> general

    skip is how many reticulations label to skip, e.g. skip=2 omits 'A' and 'B'
    time_consistent=True for networks suitable for HGT model

    Returns None if the network cannot be constructed.
    """

    treechild = networktype == 0

    t = Tree(str2tree(s))
    v = t.nodes.copy()

    reticulations = ['#' + i for i in getlabs(ord('A'), ord('Z'), reticulations + skip)]
    reticulations = reticulations[skip:]
    inserted = []

    # quite ugly code
    while reticulations:

        r = reticulations.pop()
        inserted.append(r)

        # insert a leaf labelled <r>

        while True:
            n = v[randint(0, len(v) - 1)]  # root is allowed
            if not treechild or n.label not in inserted:
                break

        np = n.parent
        ap = Node(([], []), np)
        a = Node(([], [r]), ap)

        if np:
            np.c.remove(n)
            np.c.append(ap)
            v.append(ap)
        else:
            t.root = ap
            v.insert(0, ap)  # new root

        n.parent = ap
        ap.c = [a, n]
        v.append(a)

        vshuffled = v[1:]  # skip the root
        shuffle(vshuffled)

        while vshuffled:

            # insert internal node with the label <r>; avoid cycles
            m = vshuffled.pop()

            mp = m.parent
            b = Node(([], [r]), mp)
            mp.c.remove(m)
            mp.c.append(b)
            b.c = [m]
            b.parent = mp
            m.parent = b
            v.append(b)  # not the last

            # check if dag
            s = str(t)
            net = Network(str2tree(s))  # can be done without using network (todo)
            if net.isdag():
                # dag is found
                # check if treechild if needed
                if networktype == 2 or treechild and net.treechild() or networktype == 1 \
                        and net.type1net():
                    # check if time_consistent if needed
                    if not time_consistent or time_consistent and net.istimeconsistent():
                        break  # OK; next reticulation is OK, end while

            # clean and try again
            v.remove(b)
            m.parent = mp
            mp.c.remove(b)
            mp.c.append(m)

        else:

            # clean; this position of leaf labelled <r> wasn't sucessful
            v.remove(a)
            v.remove(ap)
            if np:
                np.c.append(n)
                np.c.remove(ap)
            else:
                t.root = v[0]
            n.parent = np

            # try again with different position of leaf reticulation
            reticulations.append(r)

    return s  # return tree representation of the dag


def randdagstr(leaves, reticulations, networktype=0):
    """
    Return a network with given number of leaves and reticulations.
    If treechild==0 the network has tree-child property.
    networktype=1 -> nontreechild type 1
    networktype=2 -> general

    Returns None if the network cannot be constructed
    """

    if networktype == 0 and reticulations >= leaves:
        return None

    s = randtreestr(leaves)
    if not reticulations:
        return s
    return addretstr(s, reticulations, networktype=networktype)


import sys
from itertools import product
from collections import Counter
from random import randint, shuffle
from timeit import default_timer as timer
from treeop import Tree, Node, str2tree, getlabs, randtreestr, compcostsmp

RET_LEFTUSED = 1
RET_RIGHTUSED = 2
RET_BOTHUSED = RET_LEFTUSED | RET_RIGHTUSED

RET_C1 = 1
RET_C2 = 2
RET_C3 = 4
RET_C4 = 8
RET_C5 = 16
RET_CMASK = 31


class Network(Tree):
    def __init__(self, tup):
        Tree.__init__(self, tup)

        # recognize reticulations
        dlf = {}
        din = {}
        err = 0
        for i, n in enumerate(self.nodes):
            if n.label and n.label[0] == "#":
                retid = n.label[1:]
                n.reticulation = 1
                if n.leaf():
                    if retid in dlf:
                        print("Reticulation id <%s> already defined" % retid)
                        err = 1
                    dlf[retid] = n
                    n.reticulationleaf = 1
                else:
                    if retid in din:
                        print("Reticulation id <%s> already defined" % retid)
                        err = 1
                    din[retid] = n
                n.retid = retid
            else:
                n.reticulation = 0
                n.retid = 0

        if set(dlf.keys()) != set(din.keys()) or err:
            for k in din:
                if k not in dlf:
                    print("Missing reticulation label <%s> in leaves" % k)
            for k in dlf:
                if k not in din:
                    print("Missing reticulation label <%s> in internal nodes" % k)

            sys.exit(-1)

        self.reticulations = []
        for retnum, retid in enumerate(dlf, 1):
            l = dlf[retid]
            i = din[retid]
            lpar = l.parent
            ipar = i.parent
            lpar.c.remove(l)
            ipar.c.remove(i)

            lpar.c.insert(0, i)  # inserted at 0
            ipar.c.insert(0, i)  # inserted at 0

            i.lftparent = ipar
            i.rghparent = lpar
            i.retnum = retnum

            if i.branchlengthset and l.branchlengthset:  # set only if both are defined
                i.branchlengthset = True
                i.branchlength = (i.branchlength, l.branchlength)
            else:
                i.branchlengthset = False

            self.reticulations.append(i)
            self.nodes.remove(l)

        for n in self.nodes:
            n._setcluster()  # reconstruct clusters

    def isdag(self):
        return sorttop([(n.num, c.num) for n in self.nodes for c in n.c])

    def istimeconsistent(self):
        lftparent_nums = []
        edges = []

        # glue reticulation parents
        for r in self.reticulations:
            lftparent_nums.append(r.lftparent.num)
            r.lftparent.num = r.rghparent.num

        # create edge representation of graph
        for n in self.nodes:
            if n == self.root:
                pass
            elif n.reticulation:
                edges.append((n.lftparent.num, n.num))
            else:
                edges.append((n.parent.num, n.num))

        # unglue reticulation parents
        for i, n in enumerate(self.reticulations):
            n.lftparent.num = lftparent_nums[i]

        return bool(sorttop(edges))

    def ptab(self):
        print("=" * 80)
        for n in self.nodes:
            print(n.num, n.reticulation, n.retid, end='')
            if n.reticulation:
                if hasattr(n, "lftparent"):
                    print("rtp:%d,%d" % (n.lftparent.num, n.rghparent.num), end='')
            elif n.parent:
                print("par:%d" % n.parent.num, end='')
            print('==<<', n.netrepr(), ">>", end='')
            if n.leaf():
                print("LEAF", end='')
            else:
                print("Children=", end='')
                for c in n.c:
                    print(c.num, end=' ')
            print()

    def todotfile(self, f, nodeprefix=""):
        for n in self.nodes:
            comments = "\n".join(n.comments)
            if comments:
                comments = "\n" + comments
            f.write(nodeprefix)
            if n.leaf():
                f.write("%d [label=\"%s %d%s\"];\n" % (n.num, n, n.num, comments))
            else:
                if n.reticulation:
                    f.write("%d [shape=box,color=red, label=\"%d #%s%s\"];\n"
                            % (n.num, n.num, n.retid, comments))
                else:
                    f.write("%d [label=\"%d%s\"];\n" % (n.num, n.num, comments))
        for n in self.nodes:
            if n.reticulation:
                if hasattr(n, "lftparent"):
                    f.write("%s%d -> %s%d [color=green,label=\"%sl\"];\n" % (
                        nodeprefix, n.lftparent.num, nodeprefix, n.num, n.retid))
                    f.write("%s%d -> %s%d [color=blue,label=\"%sr\"];\n" % (
                        nodeprefix, n.rghparent.num, nodeprefix, n.num, n.retid))
            elif n.parent:
                f.write("%s%d -> %s%d;\n" % (nodeprefix, n.parent.num, nodeprefix, n.num))

    def get_leaves(self):
        """Get nodes of 1-indegree and 0-outdegree."""
        return [node for node in self.nodes if node.leaf()]

    def get_labels(self):
        """Get labels of all leaves."""
        return [node.label for node in self.get_leaves()]

    def get_inner_nodes(self):
        """Get all nodes (tree and reticulation) except of leaves."""
        return [node for node in self.nodes if not node.leaf()]

    def get_inner_tree_nodes(self):
        """Get nodes of 1-indegree and 2-outdegree plus root."""
        return [node for node in self.get_inner_nodes() if not node.reticulation]

    def valid_binary(self):
        """
        Sanity check whether network is DAG, bijective, binary nodes and
        reticulations, no parallel edges, unique root. Some of these properties
        are partially assured by the class constructor.

        Returns:
            bool for satisfying conditions
        """

        # check if network directed acyclic
        if not self.isdag():
            return False

        # check if leaves are bijective
        labels = self.get_labels()
        if len(labels) != len(set(labels)):
            return False

        # check if binary
        reticulations = self.reticulations
        for node in reticulations:
            if len(node.c) != 1:
                return False

        inner_tree_nodes = self.get_inner_tree_nodes()
        for node in inner_tree_nodes:
            if len(node.c) != 2:
                return False

        # check for parallel edges
        for node in reticulations:
            if node.lftparent == node.rghparent:
                return False

        for node in inner_tree_nodes:
            if node.c[0] == node.c[1]:
                return False

        return True

    def treechild(self):
        """
        Check if network is a valid binary network and belongs to the
        Tree-Child class.

        Returns:
            bool for satisfying conditions
        """

        # check if valid network

        if not self.valid_binary():
            return False

        # check if each inner node has a tree or leaf child node

        for node in self.nodes:
            if not node.leaf() and all(child in self.reticulations for child in node.c):
                return False

        return True

    def type1net(self):
        """
        Check if network is a valid binary network such that
        no node has >= two reticulation parents

        Returns:
            bool for satisfying conditions
        """

        # check if valid network

        if not self.valid_binary():
            return False

        # check if each inner node has a tree or leaf child node

        for node in self.nodes:
            if not node.leaf() and len(node.c) > 1 and all(child in self.reticulations for child in node.c):
                return False

        return True

    def __repr__(self):
        return self.root.netrepr()

    def __str__(self):
        return self.root.netrepr()

    def robinsonfouldsnet(self, net):
        """
        Compute network version of RF distance between two networks.
        Normalizing divider: 2(n-2)*2^k
        """

        disp_trees1 = [Tree(str2tree(t)) for t in self.displayedtrees()]
        disp_trees2 = [Tree(str2tree(t)) for t in net.displayedtrees()]

        counter1 = Counter([n.cluster for tree in disp_trees1 for n in tree.nodes])
        counter2 = Counter([n.cluster for tree in disp_trees2 for n in tree.nodes])

        dif1, dif2 = counter1 - counter2, counter2 - counter1

        return sum(dif1.values()) + sum(dif2.values())

    def displayedtreebyid(self, displayedtreeid):
        """
        Return display tree via id
        """
        ign = 0
        fmt = "{0:0%db}" % len(self.reticulations)  # format with leading zeros

        try:
            return self._bltree(
                *self._displtree(dict(zip(self.reticulations, map(int, fmt.format(displayedtreeid)))), self.root, None))
        except Exception:
            ign += 1

        if ign:
            print(f"{ign} tree(s) ignored. Is our network tree-child?", file=sys.stderr)

    # insert bl in present
    def _bltree(self, t, blset, branchlength):
        return t + (f":{branchlength}" if blset else "")

    # generic display tree generator via dictionary of ret. usages
    # given vector of reticulation usages b, a node n and parent node par from which n is reached
    # returns a tuple: (display tree encoded in string, branchlengthset_flag, branchlength )
    def _displtree(self, b, n, par):

        if n.branchlengthset:
            bl = n.branchlength
        else:
            bl = 0

        if n.leaf():
            return n.label, n.branchlengthset, bl

        if n in b:  # reticulation;  aggregate branch lengths
            t = None
            if b[n]:
                if n.lftparent == par:
                    r = self._displtree(b, n.c[0], n)
                    if not r:
                        return None
                    t, blsetc, blc = r
                    if n.branchlengthset:
                        bl = bl[0] + blc
            else:
                if n.rghparent == par:
                    r = self._displtree(b, n.c[0], n)
                    if not r:
                        return None
                    t, blsetc, blc = r
                    if n.branchlengthset:
                        bl = bl[1] + blc

            if t:
                if n.branchlengthset and blsetc:
                    return t, True, bl
                return t, False, 0

            return None
        else:
            l = [self._displtree(b, c, n) for c in n.c]
            l = [s for s in l if s]

            if not l:
                return ''

            if len(l) == 1:
                t, blsetc, blc = l[0]
                if n.branchlengthset and blsetc:
                    return t, True, bl + blc
                return t, False, 0
            s = ",".join(self._bltree(*t) for t in l)

            return "(" + s + ")", n.branchlengthset, bl

    def displayedtrees(self):

        ign = 0

        for b in product([0, 1], repeat=len(self.reticulations)):
            try:
                yield self._bltree(*self._displtree(dict(zip(self.reticulations, b)), self.root, None))
            except Exception:
                ign += 1
        if ign:
            print(f"{ign} tree(s) ignored. Is your network tree-child?", file=sys.stderr)

    def retmindup(self, gt, **opt):
        # custom multiplier for custom hashing
        hash_val = max(len(gt.nodes), len(self.nodes) + len(self.reticulations)) + 1

        # with the way we do hashes max value we store is smaller than 2 * (hash_val ** 2)
        max_size = 2 * hash_val ** 2

        # Introduce custom hashing for tuples
        # Our node numbers will be small, so we can introduce no-collision tuple hashing
        # This speeds up dict operations immensely
        # Assumption - hash_val > n1 and hash_val > n2
        def hash_tuple(n1, n2, n3=-1):
            if n3 == -1:
                return n1 * hash_val + n2
            # since n3 is either 0 or 1, changing hashing order to n3, n1, n2 reduces hash space
            return n3 * hash_val ** 2 + n1 * hash_val + n2

        # initialize lists with None, since our costs may be equal to zero
        deltav = [None] * max_size
        deltaretusage = [None] * max_size
        deltaupv = [None] * max_size
        deltaupretusage = [None] * max_size

        INFTY = 1e100

        # find mapping of leaves from gt to self
        labels = {}
        for leaf in self.leaves():
            labels[leaf.label] = leaf

        for n in gt.leaves():
            if n.label in labels:
                n.map = labels[n.label]
            else:
                raise Exception("Gene tree leaf", n, "cannot be mapped to species tree")

        def delta(tree_node, net_node):
            hsh = hash_tuple(tree_node.num, net_node.num)
            if deltav[hsh]:
                return deltav[hsh], deltaretusage[hsh]

            if tree_node.leaf():
                retusage = 0
                if tree_node.map == net_node:
                    res = 0
                else:
                    if net_node.reticulation:
                        res, retusage = delta(tree_node, net_node.c[0])
                    else:
                        res = INFTY
            else:
                tree_c0, tree_c1 = tree_node.c[0], tree_node.c[1]
                res0, retusage0 = delta(tree_c0, net_node)
                res1, retusage1 = deltaup(tree_c1, net_node)
                res = res0 + res1 + 1
                retusage = retusage0 | retusage1

                res0, retusage0 = delta(tree_c1, net_node)
                res1, retusage1 = deltaup(tree_c0, net_node)
                if res > res0 + res1 + 1:
                    res = res0 + res1 + 1
                    retusage = retusage0 | retusage1

                if not net_node.leaf():
                    if not net_node.reticulation:
                        applied = False
                        res0, retusage0 = deltaup(tree_c0, net_node.c[0])
                        res1, retusage1 = deltaup(tree_c1, net_node.c[1])
                        if res0 + res1 < res:
                            res = res0 + res1
                            retusage = retusage0 | retusage1
                            applied = True
                        res0, retusage0 = deltaup(tree_c1, net_node.c[0])
                        res1, retusage1 = deltaup(tree_c0, net_node.c[1])
                        if res0 + res1 < res:
                            res = res0 + res1
                            retusage = retusage0 | retusage1
                            applied = True
                        if applied:
                            for ch in net_node.c:
                                if ch.reticulation:
                                    if ch.lftparent == net_node:
                                        retusage |= ch.left_reticulation_used()
                                    else:
                                        retusage |= ch.right_reticulation_used()
                    else:
                        res0, retusage0 = deltaup(tree_c0, net_node.c[0])
                        res1, retusage1 = deltaup(tree_c1, net_node.c[0])
                        if res > res0 + res1 + 1:
                            res = res0 + res1 + 1
                            retusage = retusage0 | retusage1
                else:
                    res = INFTY
                    retusage = 0
            deltav[hsh] = res
            deltaretusage[hsh] = retusage
            return res, retusage

        def deltaup(tree_node, net_node):
            hsh = hash_tuple(tree_node.num, net_node.num)
            if deltaupv[hsh]:
                return deltaupv[hsh], deltaupretusage[hsh]

            res, retusage = delta(tree_node, net_node)

            if net_node.reticulation:
                res0, retusage0 = deltaup(tree_node, net_node.c[0])
                if res > res0:
                    res = res0
                    retusage = retusage0
            elif not net_node.leaf():
                net_c0, net_c1 = net_node.c
                res0, retusage0 = deltaup(tree_node, net_c0)
                res1, retusage1 = deltaup(tree_node, net_c1)
                res01 = min(res0, res1)

                if res > res01:
                    # min achieved for a child s0 or s1
                    used = notused = None
                    res = res01
                    if res == res0:
                        if net_c0.reticulation:
                            used = net_c0
                        if net_c1.reticulation:
                            notused = net_c1
                        retusage = retusage0
                    elif res == res1:
                        if net_c1.reticulation:
                            used = net_c1
                        if net_c0.reticulation:
                            notused = net_c0
                        retusage = retusage1
                    if used:
                        if net_node == used.lftparent:
                            retusage |= used.left_reticulation_used()
                        else:
                            retusage |= used.right_reticulation_used()

                    if notused:
                        if net_node == notused.lftparent:
                            retusage |= notused.right_reticulation_used()  # right must be used
                        else:
                            retusage |= notused.left_reticulation_used()  # left must be used

            deltaupv[hsh] = res
            deltaupretusage[hsh] = retusage
            return res, retusage

        res = min(delta(gt.root, s)[0] for s in self.nodes)

        opt_root = None
        opt_map = None

        for s in self.nodes:
            hsh = hash_tuple(gt.root.num, s.num)
            if deltav[hsh] == res:
                opt_root = s
                opt_map = {}
                for retic in self.reticulations:
                    opt_map[retic] = retic.get_node_retusage(deltaretusage[hsh])
                break

        if 'additional_return' in opt:
            return res, opt_root, opt_map

        return res

    # DP main for DC-non-classic
    def retmindc(self, gt, **opt):

        # custom multiplier for custom hashing
        hash_val = max(len(gt.nodes), len(self.nodes) + len(self.reticulations)) + 1

        # with the way we do hashes max value we store is smaller than 2 * (hash_val ** 2)
        max_size = 2 * hash_val ** 2

        # initialize lists with None, since our costs may be equal to zero
        deltav = [None] * max_size
        deltaretusage = [None] * max_size
        deltaupv = [None] * max_size
        deltaupoptval = [None] * max_size
        deltaupretusage = [None] * max_size

        INFTY = 1e100
        # find mapping of leaves from gt to self
        labels = {}
        for leaf in self.leaves():
            labels[leaf.label] = leaf

        for n in gt.leaves():
            if n.label in labels:
                n.map = labels[n.label]
            else:
                raise Exception("Gene tree leaf", n, "cannot be mapped to species tree")

        # Introduce custom hashing for tuples
        # Our node numbers will be small, so we can introduce no-collision tuple hashing
        # This speeds up dict operations immensely
        # Assumption - hash_val > n1 and hash_val > n2
        def hash_tuple(n1, n2, n3=-1):
            if n3 == -1:
                return n1 * hash_val + n2
            # since n3 is either 0 or 1, changing hashing order to n3, n1, n2 reduces hash space
            return n3 * hash_val ** 2 + n1 * hash_val + n2

        def delta(g, s):
            hsh = hash_tuple(g.num, s.num)
            if deltav[hsh]:
                return deltav[hsh], deltaretusage[hsh]

            if g.leaf():
                retusage = 0
                if g.map == s:
                    res = 0
                else:
                    res = INFTY
            else:
                res0, retusage0 = deltaup(g.c[0], s, 1)
                res1, retusage1 = deltaup(g.c[1], s, 1)
                res = res0 + res1
                retusage = retusage0 | retusage1

            deltav[hsh] = res
            deltaretusage[hsh] = retusage
            return res, retusage

        def deltaup(g, s, first):
            hsh = hash_tuple(g.num, s.num, first)
            if deltaupv[hsh]:
                return deltaupv[hsh], deltaupretusage[hsh]

            if len(s.c) == 1:  # s is a reticulation
                sc = s.c[0]
                res, retusage = deltaup(g, sc, 0)
                if sc.reticulation:  # special new case nonTC1
                    if s == sc.lftparent:
                        retusage |= sc.left_reticulation_used()
                    else:
                        retusage |= sc.right_reticulation_used()
                    optval = ("D^", sc)
                    # end of special case nonTC1
                else:
                    res += 1
                    optval = ("D^", sc)

            elif s.leaf():

                res, retusage = delta(g, s)
                optval = ("D", s)

            else:
                res, retusage = delta(g, s)
                optval = ("D", s)
                # s has 2 children s0 and s1
                s0, s1 = s.c

                # two variants
                du0, retusage0 = deltaup(g, s0, 0)
                du1, retusage1 = deltaup(g, s1, 0)

                if first:
                    # ignore edge <s,s0/s1> if s0/s1 is a reticulation
                    res0 = 1 - s0.reticulation + du0
                    res1 = 1 - s1.reticulation + du1
                    res = min(res, res0, res1)
                    sc = None
                    if res == res0:
                        sc, retusagec, optval = s0, retusage0, ("D^", s0)
                    elif res == res1:
                        sc, retusagec, optval = s1, retusage1, ("D^", s1)

                    # res is min
                    # problem: multiple min's
                    if sc:
                        # min is achieved by a kid
                        if sc.reticulation:
                            # kid is a reticulation node
                            retusage = retusagec
                            if s == sc.lftparent:
                                retusage |= sc.left_reticulation_used()
                            else:
                                retusage |= sc.right_reticulation_used()
                        else:
                            retusage = retusagec
                else:
                    res0, res1 = du0, du1
                    if s0.reticulation == s1.reticulation == 0:
                        # adjust if a child is not reticulation
                        res0 += 1
                        res1 += 1
                    elif s0.reticulation and s1.reticulation:
                        # not allowed in tree-child networks
                        raise Exception("Tree child network expected")

                    res01 = min(res0, res1)
                    if res > res01:
                        # min achieved for a child s0 or s1
                        used = notused = None
                        res = res01
                        if res == res0:
                            optval = ("D^", s0)
                            if s0.reticulation:
                                used = s0
                            if s1.reticulation:
                                notused = s1
                            retusage = retusage0
                        elif res == res1:
                            optval = ("D^", s1)
                            if s1.reticulation:
                                used = s1
                            if s0.reticulation:
                                notused = s0
                            retusage = retusage1
                        if used:
                            if s == used.lftparent:
                                retusage |= used.left_reticulation_used()
                            else:
                                retusage |= used.right_reticulation_used()

                        if notused:
                            if s == notused.lftparent:
                                retusage |= notused.right_reticulation_used()  # right must be used
                            else:
                                retusage |= notused.left_reticulation_used()  # left must be used

            deltaupv[hsh] = res
            deltaupretusage[hsh] = retusage
            deltaupoptval[hsh] = optval

            return res, retusage

        def ncet(v):
            c = ''
            if v & RET_LEFTUSED:
                c += "l"
            if v & RET_RIGHTUSED:
                c += "r"
            if c:
                return c
            return '.'

        def nc(dr):
            if dr:
                return "[" + ", ".join("%s%s" % (r.retid, ncet(dr[r])) for r in dr) + "]"
            return " []"

        def pptree(t, ar):
            add = " num=%d" % t.num
            if ar:
                tb = ''
                for g, s in sorted(deltav):
                    if g == t and deltav[g, s] < INFTY:
                        tb += " s%d=%d%s" % (s.num, deltav[g, s], nc(deltaretusage[g, s]))
                add += " tabv='%s'" % tb.strip()
                if not t.leaf():
                    tb = "|"
                else:
                    tb = ''
                sep = "|" if t.leaf() else " "
                tbs = ''
                for g, s, i in sorted(deltaupv):
                    if g == t and deltaupv[g, s, i] < INFTY:
                        if not t.leaf() and len(tbs) > 70:
                            tb += tbs + "|"
                            tbs = '  '
                        tbs += "%ss%s:%d.%d%s" % (sep, s.num, i, deltaupv[g, s, i], nc(deltaupretusage[g, s, i]))
                tb += tbs
                add += " tabup='%s'" % tb.strip()

            if t.leaf():
                return t.label + add
            return "(" + ",".join(pptree(c, ar) for c in t.c) + ")" + add

        def ppstree(t):
            add = " num=%d" % t.num
            if t.parent:
                add += " :%d" % (t.dagdepth() - t.parent.dagdepth())

            if t.leaf():
                return t.smplabel() + add
            return ""

        # optimal cost
        res = min(delta(gt.root, s)[0] for s in self.nodes)

        resup = res - len(gt.nodes) + 1  # DC classic

        if "dpdebug" in opt:
            print("&g " + pptree(gt.root, 0) + " dcdp=%d" % res + " dcclassic=%d" % resup)

            for t in gt.nodes:
                print(t.num, t)
                for g, s in sorted(deltav):
                    if g == t and deltav[g, s] < INFTY:
                        print("   D : %s s%d=%d %s" % (s, s.num, deltav[g, s], nc(deltaretusage[g, s])))

                for g, s, i in sorted(deltaupv):
                    if g == t and deltaupv[g, s, i] < INFTY:
                        print("   D^: %s s%d.%d=%d %s" % (s, s.num, i, deltaupv[g, s, i], nc(deltaupretusage[g, s, i])),
                              end='')
                        if not g.leaf() and deltaupoptval[g, s, i]:
                            print(" from: %s s%s" % (deltaupoptval[g, s, i][0], deltaupoptval[g, s, i][1].num), end='')
                        print()

        opt_root = None
        opt_map = None

        # print RET usage
        for s in self.nodes:
            hsh = hash_tuple(gt.root.num, s.num)
            if deltav[hsh] == res:
                opt_root = s
                opt_map = {}
                for retic in self.reticulations:
                    opt_map[retic] = retic.get_node_retusage(deltaretusage[hsh])
                if "optimalmap" in opt:
                    print("Optimal map: s%d dc=%d retusage=%s" % (s.num, res, nc(deltaretusage[gt.root, s])))

        # single print of dp arrays
        def pdparrays():
            for g, s in sorted(deltav):
                if deltav[g, s] < INFTY:
                    print(" D[", g, "-->", s.num, "]=", deltav[g, s], deltaretusage[g, s], file=sys.stderr)

            for g, s, i in sorted(deltaupv):
                if deltaupv[g, s, i] < INFTY:
                    print("D^[", g, "-->", s.num, i, "]=", deltaupv[g, s, i], deltaupretusage[g, s, i], file=sys.stderr)

            for s in self.nodes:
                if res == delta(gt.root, s):
                    print("#retDC", gt.root, "-->", s.num, " RET:", deltaretusage[gt.root, s])

        # add info to comments
        if 'dotspec' in opt:
            for gid in opt['dotspec']:
                for g in gt.nodes:
                    if g.num == gid:
                        for s in self.nodes:
                            if (g, s) in deltav:
                                if deltav[g, s] < INFTY:
                                    s.comments.append("D[%d]=%d %s" % (gid, deltav[g, s], nc(deltaretusage[g, s])))
                            for i in (0, 1):
                                if (g, s, i) in deltaupv and deltaupv[g, s, i] < INFTY:
                                    s.comments.append(
                                        "D^[%d,%d]=%d %s" % (gid, i, deltaupv[g, s, i], nc(deltaupretusage[g, s, i])))

        if 'additional_return' in opt:
            return res, opt_root, opt_map

        return res

    def branch_and_bound(self, g, use_already_set_nodes=False, branch_limit=1000000000000,
                         return_stats=False, cost_func='C'):
        """
        Find the best, complete solution for minimal displayed tree problem
        Dynamic programming solution (mindc method) allows some nodes to use both parents in tree embedding
        We use branch-and-bound method to set one parent for a non-resolved node and recurse to get a full solution
        Returns the minimal cost
        """

        def set_parent(node, parent, opt_root):
            """
            Set parent of reticulation node to parent
            Contract the network in place, so that it doesn't have unnecessary nodes and edges
            """

            def contract(node, parent, opt_root):
                """
                After choosing a parent for a node we might end up with internal nodes with in- and out- degree equal 1
                Those nodes should be contracted, which might lead to another node needing contraction and so on
                Invariant: node has a single child (reticulation node or we've already removed one of nodes' children)
                """
                if parent is None or node not in opt_root.nodes():
                    return
                i = parent.c.index(node)
                parent.c.remove(node)

                if not hasattr(node, 'removed') and node.parent == parent:
                    # If a parent is a "proper" one
                    # And we haven't marked the node as removed (only done in else below)
                    # We can simply contract the node and reroute its only child
                    if node.c[0].reticulation:
                        if node.c[0].lftparent == node:
                            node.c[0].lftparent = node.parent
                        else:
                            node.c[0].rghparent = node.parent
                    else:
                        node.c[0].parent = node.parent

                    node.parent.c.insert(i, node.c[0])
                else:
                    if parent.reticulation:
                        # If the parent is a reticulation we need to remove it completely
                        # Thus, setting the removed attribute and recursing for both  parents
                        parent.removed = True
                        contract(parent, parent.lftparent, opt_root)
                        contract(parent, parent.rghparent, opt_root)
                        del parent.removed

                    else:
                        # Parent node now has a single child, so we need to contract further
                        contract(parent, parent.parent, opt_root)

            node.parent = parent
            contract(node, node.rghparent, opt_root)
            contract(node, node.lftparent, opt_root)

        def copy(n):
            d = Network(str2tree(n.netrepr()))
            return d

        start_time = timer()
        assert branch_limit >= 0
        # Pick starting cost by comparing with one of displayed trees
        # Variable holds a cost for a complete tree
        best_cost = None
        for t in self.displayedtrees():
            best_cost = compcostsmp(g, Tree(str2tree(t)), cost_func)
            break
        if best_cost is None:
            raise Exception('No displayed trees found in network.')
        branching_count = 0
        max_depth = 0
        best_depth = 0
        dp_called = 0
        # BFS-like traversal, browsing reticulations from top to bottom
        queue = [(copy(self.root), branch_limit)]
        while queue:
            u, branch_left = queue.pop(0)
            if cost_func == 'C':
                cost, _, opt_map = u.retmindc(g, additional_return=True)
            elif cost_func == 'D':
                cost, _, opt_map = u.retmindup(g, additional_return=True)
            else:
                raise Exception("Unknown cost, expected C or D")
            opt_root = u.root

            depth = branch_limit - branch_left
            dp_called += 1
            max_depth = max(max_depth, depth)

            # If a cost for a potentially incomplete tree is higher than the one we achieved,
            # we can stop at that point
            if cost >= best_cost:
                continue

            # Solution has all reticulation parents set, we can use the cost
            if not any([value == RET_BOTHUSED for value in opt_map.values()]) or branch_left == 0:
                best_cost = cost
                best_depth = depth
                continue

            # Heuristic - we try to use as much as possible from our incomplete solution
            # to reduce computation time
            if use_already_set_nodes:
                # Set nodes that already have their parents known
                for retic, val in opt_map.items():
                    if val == RET_RIGHTUSED:
                        set_parent(retic, retic.rghparent, opt_root)
                    elif val == RET_LEFTUSED:
                        set_parent(retic, retic.lftparent, opt_root)

            # Since our opt_root might not be the initial root of the tree,
            # we need to set parents of reticulations
            # that only have one parent in a new subnetwork
            for retic, val in opt_map.items():
                if val == RET_BOTHUSED:
                    if retic.rghparent not in opt_root.nodes():
                        set_parent(retic, retic.lftparent, opt_root)
                        opt_map[retic] = RET_LEFTUSED
                    elif retic.lftparent not in opt_root.nodes():
                        set_parent(retic, retic.rghparent, opt_root)
                        opt_map[retic] = RET_RIGHTUSED

            branching_count += 1
            cp1 = copy(opt_root)
            cp2 = copy(opt_root)

            # If parent is not known, set both values in different copies and branch
            for retic, val in opt_map.items():
                if val == RET_BOTHUSED:
                    # find retic in cp1
                    x = [str(i) for i in cp1.nodes]
                    if str(retic) in x:
                        node = cp1.nodes[x.index(str(retic))]
                        set_parent(node, node.rghparent, cp1.root)
                        break

            for retic, val in opt_map.items():
                if val == RET_BOTHUSED:
                    x = [str(i) for i in cp2.nodes]
                    if str(retic) in x:
                        node = cp2.nodes[x.index(str(retic))]
                        set_parent(node, node.lftparent, cp2.root)
                        break
            # Workaround for cases when root has only one child
            if len(cp1.root.c) == 1:
                cp1.root = cp1.root.c[0]
            if len(cp2.root.c) == 1:
                cp2.root = cp2.root.c[0]

            queue.append((copy(cp1.root), branch_left - 1))
            queue.append((copy(cp2.root), branch_left - 1))
        end_time = timer()
        if return_stats:
            stats = {
                'dp_called': dp_called,
                'branching_count': branching_count,
                'max_depth': max_depth,
                'best_depth': best_depth,
                'time_in_seconds': end_time - start_time
            }
            return best_cost, stats
        return best_cost


# input: list of directed edges, must be non-empty
# return: None if cycle or empty graph
#         top. sort otherwise
def sorttop(e):
    l = []
    d1 = {}
    d2 = {}
    for x, y in e:
        d1.setdefault(x, []).append(y)
        d2.setdefault(y, []).append(x)
    roots = list(set(d1.keys()).difference(d2.keys()))
    while roots:
        n = roots.pop()
        l.append(n)
        if n not in d1:
            continue
        ms = d1[n][:]
        for m in ms:
            d1[n].remove(m)
            d2[m].remove(n)
            if not d2[m]:
                roots.append(m)
                d2.pop(m)
            if not d1[n]:
                d1.pop(n)
    if d1 or d2:
        return None
    return l


def addretstr(s, reticulations, skip=0, networktype=0, time_consistent=False):
    """
    Exhaustively add given number of reticulations to a tree or tree
    representation of network.

    networktype=0 -> treechild
    networktype=1 -> nontreechild type 1
    networktype=2 -> general

    skip is how many reticulations label to skip, e.g. skip=2 omits 'A' and 'B'
    time_consistent=True for networks suitable for HGT model

    Returns None if the network cannot be constructed.
    """

    treechild = networktype == 0

    t = Tree(str2tree(s))
    v = t.nodes.copy()

    reticulations = ['#' + i for i in getlabs(ord('A'), ord('Z'), reticulations + skip)]
    reticulations = reticulations[skip:]
    inserted = []

    # quite ugly code
    while reticulations:

        r = reticulations.pop()
        inserted.append(r)

        # insert a leaf labelled <r>

        while True:
            n = v[randint(0, len(v) - 1)]  # root is allowed
            if not treechild or n.label not in inserted:
                break

        np = n.parent
        ap = Node(([], []), np)
        a = Node(([], [r]), ap)

        if np:
            np.c.remove(n)
            np.c.append(ap)
            v.append(ap)
        else:
            t.root = ap
            v.insert(0, ap)  # new root

        n.parent = ap
        ap.c = [a, n]
        v.append(a)

        vshuffled = v[1:]  # skip the root
        shuffle(vshuffled)

        while vshuffled:

            # insert internal node with the label <r>; avoid cycles
            m = vshuffled.pop()

            mp = m.parent
            b = Node(([], [r]), mp)
            mp.c.remove(m)
            mp.c.append(b)
            b.c = [m]
            b.parent = mp
            m.parent = b
            v.append(b)  # not the last

            # check if dag
            s = str(t)
            net = Network(str2tree(s))  # can be done without using network (todo)
            if net.isdag():
                # dag is found
                # check if treechild if needed
                if networktype == 2 or treechild and net.treechild() or networktype == 1 \
                        and net.type1net():
                    # check if time_consistent if needed
                    if not time_consistent or time_consistent and net.istimeconsistent():
                        break  # OK; next reticulation is OK, end while

            # clean and try again
            v.remove(b)
            m.parent = mp
            mp.c.remove(b)
            mp.c.append(m)

        else:

            # clean; this position of leaf labelled <r> wasn't sucessful
            v.remove(a)
            v.remove(ap)
            if np:
                np.c.append(n)
                np.c.remove(ap)
            else:
                t.root = v[0]
            n.parent = np

            # try again with different position of leaf reticulation
            reticulations.append(r)

    return s  # return tree representation of the dag


def randdagstr(leaves, reticulations, networktype=0):
    """
    Return a network with given number of leaves and reticulations.
    If treechild==0 the network has tree-child property.
    networktype=1 -> nontreechild type 1
    networktype=2 -> general

    Returns None if the network cannot be constructed
    """

    if networktype == 0 and reticulations >= leaves:
        return None

    s = randtreestr(leaves)
    if not reticulations:
        return s
    return addretstr(s, reticulations, networktype=networktype)
