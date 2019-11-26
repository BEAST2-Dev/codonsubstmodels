package beast.evolution.tree;

import beast.core.Input;
import beast.core.StateNode;
import beast.core.util.Log;
import beast.evolution.alignment.CodonAlignment;

import java.io.PrintStream;
import java.util.Arrays;

/**
 * The wrapper class to store all node states including tips.
 * {@link NodeStates} is referred to a node, and its array are sites.
 * and the state starts from 0. <br>
 *
 * <code>nodeNr</code> is the node index number starting from 0,
 * and for internal nodes, they are ranged from the number of leaf nodes
 * to the total nodes - 1.<br>
 */
public class NodesStatesAndTree extends NodesStates {

    final public Input<TreeInterface> treeInput = new Input<>("tree",
            "phylogenetic beast.tree to map sequence data correctly to the tips", Input.Validate.REQUIRED);

    final public Input<String> initINSInput = new Input<>("initINS",
            "whether and how to initialise the states of internal nodes. Choose 'random' or 'parsimony'",
            Input.Validate.OPTIONAL);

    final public Input<Long> initINSRandomSeedInput = new Input<>("seed",
            "seed to initialise the states of internal nodes if choosing 'random' method",
            new Long(777));


    protected TreeInterface tree;

    /**
     * init class and states in tips given {@link CodonAlignment}.
     * internal nodes states need to init after this.
     * @param codonAlignment
     */
    public NodesStatesAndTree(CodonAlignment codonAlignment, TreeInterface tree) {
        super(codonAlignment);
        this.tree = tree;

        // sanity check: alignment should have same #taxa as tree
        if (codonAlignment.getTaxonCount() != tree.getLeafNodeCount())
            throw new IllegalArgumentException("The number of nodes in the tree does not match the number of sequences");
        // nodeCount == tree.getNodeCount()
        if (getNodeCount() != tree.getNodeCount())
            throw new IllegalArgumentException("The dimension of nodes states should equal to " +
                    "the number of nodes in the tree !\n" + getNodeCount() + " != " + tree.getNodeCount());
        // get***Count uses nodeCount
        assert getTipsCount() == tree.getLeafNodeCount();
        assert getInternalNodeCount() == tree.getInternalNodeCount();

        Log.info.println("  " + codonAlignment.toString(true));

        // init from constructor
        nodesStates = new NodeStates[nodeCount];

        // set tips states
        initTipStates(codonAlignment, tree);
        // call initINS
    }

    @Override
    public void initAndValidate() {
        // need data type, site count, taxa count
        CodonAlignment codonAlignment = codonAlignmentInput.get();
        TreeInterface tree = treeInput.get();
        // init params
        NodesStatesAndTree ns = new NodesStatesAndTree(codonAlignment, tree);

        // internal nodes
        initINS(initINSInput.get(), initINSRandomSeedInput.get());
    }

    // set tips states to nodesStates[] by nodeNr indexing
    protected void initTipStates(CodonAlignment codonAlignment, TreeInterface tree) {
        // tips
        for (int i=0; i < getTipsCount() ; i++) {
            Node tip = tree.getNode(i);
            // use nodeNr to map the index of nodesStates[]
            nodesStates[i] = new NodeStates(tip, codonAlignment);
        }
    }

    // "random", "parsimony"
    public void initINS(String initMethod, long seed) {
        // internal nodes

        if ("random".equalsIgnoreCase(initMethod)) {
            int[][] inStates = initINStatesRandom(getInternalNodeCount(),
                    getStateCount(), getSiteCount(), seed);

            for (int i=0; i < getInternalNodeCount(); i++) {
                int nR = i + getTipsCount();
                nodesStates[nR] = new NodeStates(nR, inStates[i], getStateCount());
            }
        } else if ("parsimony".equalsIgnoreCase(initMethod)) {
            throw new UnsupportedOperationException("in dev");

        } else {
            throw new IllegalArgumentException("No method selected to initialise the states at internal nodes !");
        }

    }

    /**
     * Parsimony to init states Snell & Childress 1987.
     * Equally to choose a state from the ambiguous set.
     */
    public void initINStatesParsimony() {
        throw new UnsupportedOperationException();

    }

    /**
     * Assuming rootIndex = tree.getNodeCount() - 1,
     * and then validate rootIndex == tree.getRoot().getNr()
     * @return {@link TreeInterface}
     */
    public TreeInterface getValidTree() {
        // exclude root node, branches = nodes - 1
        final int rootIndex = tree.getNodeCount() - 1;
        if (rootIndex != tree.getRoot().getNr())
            throw new RuntimeException("Invalid root index " + tree.getRoot().getNr() + " != " + rootIndex);

        return tree;
    }

    public int getRootIndex() {
        int rootIndex = tree.getNodeCount() - 1;
        assert rootIndex == tree.getRoot().getNr();
        return rootIndex;
    }

//    public TreeInterface getTree() {
//        return tree;
//    }


    @Override
    public void log(long sampleNr, PrintStream out) {
        super.log(sampleNr, out);
    }


    //TODO full copy or just currentMatrixIndex?
    @Override
    public StateNode copy() {
        try {
            @SuppressWarnings("unchecked")
            final NodesStatesAndTree copy = (NodesStatesAndTree) this.clone();
            for (int i = 0; i < nodesStates.length; i++)
                copy.nodesStates[i] = (NodeStates) nodesStates[i].copy();

            //storedStates = copy.storedStates.clone();
            copy.nodeIsDirty = new boolean[getNodeCount()];
            return copy;
        } catch (Exception e) {
            e.printStackTrace();
        }
        return null;

    }

    @Override
    public void assignTo(StateNode other) {
        @SuppressWarnings("unchecked")
        final NodesStatesAndTree copy = (NodesStatesAndTree) other;
        copy.setID(getID());
        copy.index = index;
        for (int i = 0; i < nodesStates.length; i++)
            nodesStates[i].assignTo(copy.nodesStates[i]);
        //storedStates = copy.storedStates.clone();
        copy.lower = lower;
        copy.upper = upper;
        copy.nodeIsDirty = new boolean[getNodeCount()];
    }

    @Override
    public void assignFrom(StateNode other) {
        @SuppressWarnings("unchecked")
        final NodesStatesAndTree source = (NodesStatesAndTree) other;
        setID(source.getID());
        for (int i = 0; i < nodesStates.length; i++)
            nodesStates[i].assignFrom(source.nodesStates[i]);
        lower = source.lower;
        upper = source.upper;
        nodeIsDirty = new boolean[source.getNodeCount()];
    }

    @Override
    public void assignFromFragile(StateNode other) {
        @SuppressWarnings("unchecked")
        final NodesStatesAndTree source = (NodesStatesAndTree) other;
        for (int i = 0; i < nodesStates.length; i++)
            nodesStates[i].assignFromFragile(source.nodesStates[i]);
        Arrays.fill(nodeIsDirty, false);
    }

    @Override
    public void fromXML(org.w3c.dom.Node node) {
        throw new UnsupportedOperationException("in dev");
    }

    @Override
    public int scale(double scale) {
        throw new UnsupportedOperationException("in dev");
    }


}
