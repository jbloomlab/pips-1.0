"""Module establishing classes for representing phylogenetic trees.

Written by Jesse Bloom, 2007."""



import sys
import re


def TallySiteChangesAlongTree(t):
    """Counts the number of times a site has changed along a reconstructed tree.

    Takes as input a single argument 't' which shuld represent a 'tree.Tree'
        object.  This true must be rooted, and must have sequence assigned to
        both all tip nodes and all interior nodes.  That is, the interior
        sequences must already be reconstructed by some method.  These sequences
        should all be the same length (i.e., they should be aligned).
    This method traverses along the tree starting at the root.  For each site,
        it counts the number of times that the identity of that site has
        changed, and what these changes have been.  It returns a dictionary 'tally'
        keyed by sequence numbers (ranging from 1 to the length of the sequence).
        Entry 'tally[i]' is the 2-tuple '(n, changes)'.  'n' is the total number
        of times site 'i' is inferred to have changed identities.  'changes' is
        itself a dictinary.  For each transition from x to y at site i, 
        'changes[(x, y)]' counts the number of times that transition occurred."""
    assert isinstance(t, Tree)
    root = t.GetRoot()
    assert root.root
    if root.tip:
        raise ValueError("Root is also tip, indicating a null tree.")
    if not root.seq:
        raise ValueError("No sequence defined for root node.")
    seqlength = len(root.seq)
    tally = {}
    for i in range(1, seqlength + 1):
        tally[i] = (0, {})
    
    def _RecursivelyTally(tally, ancestorseq, node, seqlength):
        """Recursive internal subroutine for tallying changes."""
        if not node.seq:
            raise ValueError("Node has no sequence.")
        if len(node.seq) != seqlength:
            raise ValueError("Node sequence has different length.")
        for i in range(seqlength):
            if ancestorseq[i] != node.seq[i]:
                (n, itally) = tally[i + 1]
                n += 1
                tup = (ancestorseq[i], node.seq[i])
                if tup in itally:
                    itally[tup] += 1
                else:
                    itally[tup] = 1
                tally[i + 1] = (n, itally)
        if node.tip:
            return
        else:
            _RecursivelyTally(tally, node.seq, node.rightdescendent, seqlength)
            _RecursivelyTally(tally, node.seq, node.leftdescendent, seqlength)

    # Start recursively tallying down the tree
    _RecursivelyTally(tally, root.seq, root.rightdescendent, seqlength)
    _RecursivelyTally(tally, root.seq, root.leftdescendent, seqlength)
    return tally



def SplitNewick(newick_tree):
    """Splits a bifurcating phylogenetic tree into its next two branches.

    Takes as input a string specifying a sub-tree of a Newick tree.
        If the string specifies the full tree, the trailing semicolon
        and the parentheses surrounding the root must already be
        removed.  In other words, the tree
        '(((One:0.2,Two:0.3):0.3,(Three:0.5,Four:0.3):0.2):0.3,Five:0.7):0.0;'
        is not acceptable, but the tree
        '((One:0.2,Two:0.3):0.3,(Three:0.5,Four:0.3):0.2):0.3,Five:0.7'
        is acceptable.  The same limitations to the Newick tree apply
        as are listed in the documentation string for the '__init__' 
        method of 'Tree'.
    Returns the 2-tuple '(left_branch, right_branch)'.  The elements of this
        tuple specify the left and right branches of the tree, respectively.
        Each branch is itself a 3-tuple of the form '(branch, length, tip)'.
        'branch' is a string giving the subtree if this branch is a subtree,
        or the name of the node if it is a tip.  'length' is a number giving
        the length of the branch leading to this tip or subtree.  'tip' is
        'True' if the branch is a tip, or 'False' otherwise.  For example,
        calling this method on the example subtree given above would return
        'right_branch' as '("Five", 0.7, True)' and 'left_branch' as
        '("(One:0.2,Two:0.3):0.3,(Three:0.5,Four:0.2):0.3", 0.3, False)'."""
    assert isinstance(newick_tree, str)
    # We find the comma that is not enclosed in any parentheses.  There should
    # be exactly one such comma.
    nopen = 0 # 'nopen' is the number of open parentheses
    comma_position = []
    for i in range(len(newick_tree)):
        if newick_tree[i] == '(':
            nopen += 1
        elif newick_tree[i] == ')':
            nopen -= 1
            if nopen < 0:
                raise ValueError("Negative parentheses count at position %d of %d in tree %s" % (i, len(newick_tree), newick_tree))
        elif newick_tree[i] == ',':
            if nopen == 0:
                comma_position.append(i)
    if len(comma_position) == 0:
        raise ValueError, "There are no commas not enclosed in parentheses in %s." % newick_tree
    elif len(comma_position) > 1:
        raise ValueError, "There are multiple commas not enclosed in parentheses in %s." % newick_tree
    icomma = comma_position[0]
    (branch_l, branch_r) = (newick_tree[ 0 : icomma], newick_tree[icomma + 1 : ])
    name_length_match = re.compile('^(?P<name>.*):(?P<length>[\d\.]+)$')
    left_branch = []
    right_branch = []
    for (in_branch, out_branch) in [(branch_l, left_branch), (branch_r, right_branch)]:
        m = name_length_match.search(in_branch)
        if not m:
            raise ValueError, "Could not parse %s." % in_branch
        length = float(m.group('length'))
        branch = m.group('name').strip()
        if '(' in branch:
            assert branch[0] == '(' and branch[-1] == ')'
            branch = branch[1 : -1]
            tip = False
        else:
            tip = True
        out_branch += [branch, length, tip]
    return (tuple(left_branch), tuple(right_branch))



def RecursivelyAssignDescendents(node, newick_tree, scalebranchlength, tipnames_sequences_dict):
    """Recursively assigns the descendents to a node in a phylogenetic tree.

    On input, 'node' is a Node that is not a tip and has no descendents
        assigned.  'newick_tree' is a string giving the Newick tree
        that has this node as its root.  The same restrictions apply
        to 'newick_tree' as those mentioned in the '__init__' method
        of 'Tree'.  
    This method assigns the descendent nodes for 'node', and the same
        for these descendent nodes, until the tips of the tree have
        been reached.  The tips of the trees have names assigned.  Branch
        lengths are assigned to all nodes.
    'scalebranchlength' is an optional argument that allows us to rescale
        all branch lengths.  If it has a value of 'None', then
        rescaling is done.  If it has another value, then it should be a number
        greater than zero.  In this case, all branch lengths specified by 
        'newick_tree' are multiplied by 'scalebranchlength'.
    'tipnames_sequences_dict' has a similar meaning as described in the '__init__'
        method of 'Tree'.  It is used to assign sequences to the tip nodes.  It
        can still be 'None', meaning that no sequences are assigned.  Otherwise
        it is a dictionary, with the keys being the sequence names and the
        values being the sequences.  There must be a key for every tip node.
    """
    if scalebranchlength == None:
        scalebranchlength = 1.0 # No scaling is the same as scaling by one
    else:
        assert 0 < scalebranchlength and isinstance(scalebranchlength, (int, float))
        pass # we scale by whatever number is set
    assert isinstance(node, Node) and not node.tip
    (left_branch, right_branch) = SplitNewick(newick_tree)
    # set left branch 
    if left_branch[2]: # left branch is a tip
        name = left_branch[0]
        try:
            seq = tipnames_sequences_dict[name]
        except KeyError:
            raise ValueError, "No sequence assigned for %s in %s." % (name, str(tipnames_sequences_dict))
        except TypeError:
            if tipnames_sequences_dict == None:
                seq = None
            else:
                raise
        left_node = Node(tip = True, ancestor = node, name = name, seq = seq, ancestorbranch = left_branch[1] * scalebranchlength)
    else: # left branch is not a tip, proceed recursively
        left_node = Node(ancestor = node, ancestorbranch = left_branch[1] * scalebranchlength)
        RecursivelyAssignDescendents(left_node, left_branch[0], scalebranchlength, tipnames_sequences_dict)
    node.leftdescendent = left_node
    node.leftbranch = left_branch[1] * scalebranchlength
    # set right branch 
    if right_branch[2]: # left branch is a tip
        name = right_branch[0]
        try:
            seq = tipnames_sequences_dict[name]
        except KeyError:
            raise ValueError, "No sequence assigned for %s." % name
        except TypeError:
            if tipnames_sequences_dict == None:
                seq = None
            else:
                raise
        right_node = Node(tip = True, ancestor = node, name = name, seq = seq, ancestorbranch = right_branch[1] * scalebranchlength)
    else: # right branch is not a tip, proceed recursively
        right_node = Node(ancestor = node, ancestorbranch = right_branch[1] * scalebranchlength)
        RecursivelyAssignDescendents(right_node, right_branch[0], scalebranchlength, tipnames_sequences_dict)
    node.rightdescendent = right_node
    node.rightbranch = right_branch[1] * scalebranchlength



def RecursivelySetNumbers(node, number):
    """Function for setting the numbers for nodes of a tree.
   
    Calling this function on the root node of a tree will set unique 
        numbers for all nodes of the tree.
    Calling this function with the root node and 'number' set to zero will
        return an integer giving the number of nodes in a tree (tip and
        internal nodes combined).  Typically, that is how you will use
        the function: 'RecursivelySetNumbers(tree.GetRoot(), 0)'
    On call, 'node' should be a 'Node'.
    'number' should be an integer that is the number that is set for 'node'.
    Sets a number for 'node', and also recursively sets numbers for any of
        its descendents.  The returned value is the next number (number + 1)
        that is not set for any node.
    """
    node.number = number
    number += 1 # number of next node
    if not node.tip:
        number = RecursivelySetNumbers(node.rightdescendent, number)
        number = RecursivelySetNumbers(node.leftdescendent, number)
    return number



def RecursivelySetTipAndInternalNumbers(node, tip_number, internal_number):
    """Function for setting numbers for tip and internal nodes.

    The node numbers are the values accessed by 'node.number'.  Any
        existing values for these numbers are cleared.
    Typically this function will initially be called with 'node' as the
        root node of a tree, and 'tip_number' and 'internal_number' both
        equal to zero.  In this case, all of the tip nodes in the tree
        will be assigned integer numbers going from 0, 1, ...
        All of the inernal nodes in the tree will also be assigned integer
        numbers going from 0, 1, ...  Crucially, the numbers for the internal
        nodes are assigned so that if 'nodej' is a descendent of 'nodei', then
        'nodej.number < nodei.number' for all nodes.  This leads to the implication
        that the root node will have the highest internal node number.  Any existing 
        numbers for the nodes are cleared.  Note that it is NOT the case that all nodes
        in the tree have unique numbers; all tip nodes have unique numbers,
        and all internal nodes have unique numbers.  The returned value is
        the 2-tuple (n_tips, n_internal) where 'n_tips' is the number of
        tip nodes and 'n_internal' is the number of internal nodes.
    The above paragraph describes how the function will typically function.  In its
        implementation, this function works by recursively calling itself.  In
        these recursive calls, the values of 'tip_number' and 'internal_number' are
        incremented as nodes are assigned numbers -- the value of each of these
        corresponds to the next number that will be assigned to the next encountered
        tip or internal node.
    """
    if node.tip:
        node.number = tip_number
        return (tip_number + 1, internal_number)
    else:
        (tip_number, internal_number) = RecursivelySetTipAndInternalNumbers(node.rightdescendent, tip_number, internal_number)
        (tip_number, internal_number) = RecursivelySetTipAndInternalNumbers(node.leftdescendent, tip_number, internal_number)
        node.number = internal_number
        return (tip_number, internal_number + 1)



def ListsOfNodes(node, tip_list, internal_list):
    """Function that returns the list of all nodes in a tree.

    Typically you will call this function with 'node' equal to the root node
        of a tree, and 'tip_list' and 'internal_list' as variables set to empty lists.
        On return, 'tip_list' is a list of all tip nodes in the tree with 'node' as its 
        root node.  'internal_list' is a list of all internal nodes in the tree with
        'node' as its root.
    If the method is called with 'node' equal to something other than the root node,
        then the lists are of all nodes in the subtree with root 'node'.
    If the method is called with 'tip_list' and 'internal_list' as non-empty lists,
        then the tip and internal nodes are merely appended to the existing entries
        of these lists.
    This function works by recursively calling itself.
    """
    if node.tip:
        tip_list.append(node)
    else:
        internal_list.append(node)
        ListsOfNodes(node.rightdescendent, tip_list, internal_list)
        ListsOfNodes(node.leftdescendent, tip_list, internal_list)



def QueryTree(tree, nodelocation):
    """Gets information about a particular node of a tree.

    'tree' is a 'Tree' object specifying a phylogenetic tree.
    'nodelocation' is a string giving the location of the node that
        we want to query from the tree.  It is specified in terms of
        tracing the tree forward from its root node.  A node location
        of '' (empty string) refers to the root node.  A node location
        of 'right' refers to the right descendent of the root node.
        A node location of 'right, left, left' refers to the node reached
        by beginning at the root node, taking the right descendent node,
        then this node's left descendent, then that node's left descendent.
    If 'nodelocation' fails to specify an existing node of the tree, then
        an exception is raised.
    Otherwise, the function returns the 3-tuple:
        '(tip, name, sequence)'
        'tip' is 'True' if the node is a tip node, and 'False' otherwise.
        'name' is the string giving the name of the node, or 'None' if
            no name is assigned to the node.
        'sequence' is a string giving the sequence assigned to the node,
            or 'None' if no sequence is assigned to this node.
    """
    assert isinstance(tree, Tree)
    assert isinstance(nodelocation, string)
    nodelocation = [entry.strip() for entry in nodelocation.split(',')]
    node = tree.GetRoot()
    for direction in nodelocation:
        if node.tip:
            raise ValueError("'nodelocation' is specifying a descendent to a tip node with a %s" % direction)
        if direction == 'right':
            node = node.rightdescendent
        elif direction == 'left':
            node = node.leftdescendent
        else:
            raise ValueError("Invalid descendent direction of %s.  Must be 'right' or 'left'." % direction)
    return (node.tip, node.name, node.sequence)



def PrintAllTreeInformation(tree, f=sys.stdout):
    """Prints information about all nodes, names, and sequences in a tree.

    'tree' is a Tree object specifying a phylogenetic tree.
    'f' is a writeable filelike object.  By default, it is set to 'sys.stdout'.
        You can also set it to some other writeable file.
    This function writes a FASTA file containing an entry for every node
        in the tree.  The headers begin with the node's location in the tree,
        then whether this is an interior or tip node, and finally the node's
        assigned name.  The sequence is the node's assigned sequence, or the
        string "None" if no sequence is assigned.  The location is specified
        as "ROOT" if the node is a root node, and otherwise as a string of
        comma delimited "LEFT" and "RIGHT" characters indicating how the node
        is reached from the root node.  If no name is assigned for the node,
        the name is printed as "None". For example, the function might write:
        >ROOT INTERIOR None 
        None
        >LEFT INTERIOR None
        None
        >LEFT,LEFT TIP A/WSN/33 (H1N1)
        MTAGCKL
        >LEFT,RIGHT TIP A/PR/8/34 (H1N1)
        MTAGNL
        >RIGHT TIP A/Aichi/1968 (H3N2)
        MSGGNL
    """
    root = tree.GetRoot()
    assert not root.tip
    f.write(">ROOT INTERIOR %s\n%s\n" % (root.name, root.seq))
    # create function for recursively writing out node information
    def _RecursivelyWriteNode(node, location):
        """Recursive function for printing node descriptions."""
        if node.tip:
            f.write(">%s TIP %s\n%s\n" % (location, node.name, node.seq))
        else:
            f.write(">%s INTERIOR %s\n%s\n" % (location, node.name, node.seq))
            _RecursivelyWriteNode(node.leftdescendent, "%s,LEFT" % location)
            _RecursivelyWriteNode(node.rightdescendent, "%s,RIGHT" % location)
    # write out all of the descendent nodes
    _RecursivelyWriteNode(root.leftdescendent, 'LEFT')
    _RecursivelyWriteNode(root.rightdescendent, 'RIGHT')



class Tree(object):
    """Class for representing a phylogenetic tree.

    Currently only competent for representing a rooted bifurcating tree."""

    def __init__(self, newick_tree, scalebranchlength=None, tipnames_sequences=None):
        """Creates a phylogenetic tree.

        This class can currently represent only rotted bifurcating trees.
        On call, 'newick_tree' should be a string specifying a tree in Newick
            format with branch lengths.  The tree tips must have non-empty names,
            and the names cannot contain spaces.  Names cannot be duplicated.
            Currently, only certain aspects of the full Newick format
            (http://evolution.gs.washington.edu/phylip/newick_doc.html) are
            handled by this method.  Specifically, there is no mechanism for:
            * handling names for internal nodes
            * handling quoted labels
            * converting underscores to blanks
            * handling newlines
            * handling comments
            * handling the absence of branch lengths.
            It is optional whether the branch length is specified for the root
                node.  If it is specified, it is ignored.
            Any spaces and line breaks are parsed out of 'newick_tree'.
        'scalebranchlength' is an optional argument that allows us to rescale
            all branch lengths.  If it has its default value of 'None', then no
            rescaling is done.  If it has another value, then it should be a number
            greater than zero.  In this case, all branch lengths specified by 
            'newick_tree' are multiplied by 'scalebranchlength'.
        'tipnames_sequences' is an optional argument that allows us to specify
            sequences for the tips of the tree.  By default, it is 'None' which
            means that no sequences are added.  If it has another value, it
            must be a list of 2-tuples of the format (name, seq).  'name' gives
            the name of the tip nodes, and there must be a unique value of 'name'
            for each 'name' specified for the tree tips.  The sequence for the
            tip node with name 'name' is set to 'seq'.
        Each node of the tree is assigned a unique integer identifier n, with
            0 <= n < number of nodes in tree.
        """
        assert isinstance(newick_tree, str)
        newick_tree = ''.join(newick_tree.split()) # get rid of any spaces
        assert newick_tree
        assert scalebranchlength == None or (isinstance(scalebranchlength, (int, float)) and scalebranchlength > 0)
        if tipnames_sequences == None: 
            tipnames_sequences_dict = None
        else:
            assert isinstance(tipnames_sequences, list)
            tipnames_sequences_dict = {}
            for (name, seq) in tipnames_sequences:
                if name in tipnames_sequences_dict:
                    raise ValueError, "Duplicate name of %s." % name
                tipnames_sequences_dict[name] = seq
        self._newick = newick_tree # string giving Newick format
        # strip away the space and semicolon do some error checking on the Newick tree 
        newick_tree = newick_tree.strip()
        if newick_tree[-1] != ';':
            raise ValueError, "Newick tree does not end in semicolon."
        else: # strip away semicolon and the branch length for the root
            ifinal = len(newick_tree) - 1
            while ifinal > 0 and newick_tree[ifinal] != ')':
                ifinal -= 1
            newick_tree = newick_tree[ : ifinal + 1] 
        if not (newick_tree.count('(') == newick_tree.count(')') > 0):
            raise ValueError, "Invalid number of begin and end parentheses in Newick tree."
        # initialize the tree root
        assert newick_tree[0] == '(' and newick_tree[-1] == ')'
        newick_tree = newick_tree[1 : -1]
        self._root = Node(root = True)
        RecursivelyAssignDescendents(self._root, newick_tree, scalebranchlength=scalebranchlength, tipnames_sequences_dict=tipnames_sequences_dict)
        RecursivelySetNumbers(self._root, 0)
    
    def GetRoot(self):
        """Returns the Node object that is the root of the tree."""
        return self._root
    
    def GetNewickTree(self):
        """Returns the Newick string representing the tree."""
        return self._newick



class Node(object):
    """Class for representing nodes of a phylogenetic tree.
    
    Currently only competent for representing nodes on a bifurcating tree.
    A node object 'n' has the following public attributes:
        n.tip - True if the node is a tip, False otherwise.
        n.root - True if the node is a root node, False otherwise
        n.number - an integer that can be assigned to the node as an identifier.
            It is 'None' of no integer has been set.  Typically, this value
            is used to uniquely identify nodes in a tree, and is set during
            construction of the tree.
        n.name - a string giving the name of a node, or 'None' if no name exists.
        n.seq - a string giving the (protein or DNA) sequence corresponding 
            to a node, or 'None' if there is no sequence.  Should only have a 
            value other than 'None' if n.tip is True.
        n.ancestor - the node that is the ancestor of this node, or 'None' if
            the node has no ancestor.  Typically there will be no ancestor
            if n.root is True.
        n.rightdescendent - the node that is the right descendent of this node,
            or 'None' if there is no right descendent.  Typically there will be
            no right descendent if n.tip is True.
        n.leftdescendent - like n.rightdescendent, but for the left descendent.
        n.ancestorbranch is the length of the branch to the ancestor, or 'None'
            if no such length is set.
        n.rightbranch is the length to the branch to the right descendent, or 'None'
            if no such length is set.
        n.leftbranch is like n.rightbranch, but for the left ancestor.
    """
    
    def __init__(self, root = False, tip = False, name = None, seq = None, ancestor = None, leftdescendent = None, rightdescendent = None, ancestorbranch = None, leftbranch = None, rightbranch = None, number = None):
        """Creates an instance of 'Node'.

        There are no required arguments.  There are the following optional arguments:
        'root' should be 'True' iff this is the root of the tree; 'False' by default.
        'tip' should be 'True' iff this is a tip of the tree; 'False' by default.
        'name' can be a string giving the name of the node, if such a name exists.
            By default is 'None'.
        'seq' can be a string giving the gene/protein sequence of the node, but 
            only if 'tip' is 'True'.  By default is 'None'.
        'ancestor' can be another 'Node' object which is the ancestor of this node.
            By default is 'None', and must be 'None' if 'root' is True.
        'leftdescendent' can be another 'Node' object which is the left descendent of 
            this node.  By default is 'None', and must be 'None' if 'tip' is True.
        'rightdescendent' is like 'leftdescendent', but is the right descendent.
        'ancestorbranch' can be a number giving the length of the branch back to the
            ancestor of this node.  By default is 'None', and must be 'None' if 
            'ancestor' is none.
        'leftbranch' can be a number giving the length of the left branch of the
            tree.  By default is 'None', and must be 'None' if 'leftdescendent' is 'None'.
        'rightbranch' is like 'leftbranch', but for the right branch.
        'number' is an integer that can be assigned to the node as an identifier.  It
            is 'None' by default. 
        """
        assert isinstance(root, bool)
        self.root = root
        assert isinstance(tip, bool)
        self.tip = tip
        assert name == None or isinstance(name, str)
        self.name = name
        assert seq == None or (self.tip and isinstance(seq, str))
        self.seq = seq
        assert ancestor == None or (not self.root and isinstance(ancestor, Node))
        self.ancestor = ancestor
        assert leftdescendent == None or (not self.tip and isinstance(leftdescendent, Node))
        assert rightdescendent == None or (not self.tip and isinstance(rightdescendent, Node))
        self.leftdescendent = leftdescendent
        self.rightdescendent = rightdescendent
        assert ancestorbranch == None or isinstance(ancestorbranch, (int, float))
        assert rightbranch == None or isinstance(rightbranch, (int, float))
        assert leftbranch == None or isinstance(leftbranch, (int, float))
        assert not (self.root and ancestorbranch != None)
        assert not (self.tip and rightbranch != None)
        assert not (self.tip and leftbranch != None)
        self.ancestorbranch = ancestorbranch
        self.rightbranch = rightbranch
        self.leftbranch = leftbranch
        assert number == None or isinstance(number, int)
        self.number = number
