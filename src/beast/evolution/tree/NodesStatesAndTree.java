package beast.evolution.tree;

import beast.core.Input;
import beast.core.StateNode;
import beast.core.util.Log;
import beast.evolution.alignment.CodonAlignment;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.Random;

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


    // all nodes, use nodeNr to map the index of nodesStates[]
    protected  NodeStates[] nodesStates;

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
     * random states given genetic code
     */
    public int[][] initINStatesRandom(final int internalNodeCount, final int stateCount, final int siteCount, long seed) {

        Log.info("Random generate codon states using " + getGeneticCode().getDescription() +
                " for " + internalNodeCount + " internal nodes, " + getSiteCount() + " codon, seed = " + seed);

        Random generator;
        if (seed > 0)
            generator = new Random(seed);
        else
            generator = new Random();

        // states[internal nodes][sites], where i from 0 to internalNodeCount-1
        int[][] inStates = new int[internalNodeCount][siteCount];
        for (int i=0; i < inStates.length; i++) {
            for (int j = 0; j < inStates[0].length; j++) {
                // 0 - 63
                inStates[i][j] = (int) (generator.nextDouble() * stateCount);
                // skip stop codon states, such as vertebrateMitochondrial: 8  10  48  50
                while (getGeneticCode().isStopCodon(inStates[i][j]))
                    inStates[i][j] = (int) (generator.nextDouble() * stateCount);
            }
        }
        return inStates;
    }

    /**
     * Parsimony to init states. Equally to choose a state from the ambiguous set.
     */
    public void initINStatesParsimony() {
        throw new UnsupportedOperationException();
        //traverse 1
        //traverse 2
    }

    /**
     * cleans up and deallocates arrays.
     */
    @Override
    public void finalize() throws Throwable {
        for (NodeStates nodeStates : nodesStates)
            nodeStates.finalize();
        nodeIsDirty = null;
    }


    //TODO store restore a site
    @Override
    protected void store() {
        for (NodeStates nodeStates : nodesStates)
            nodeStates.store();
        //TODO need nodeIsDirty[]?
    }

    @Override
    public void restore() {
        for (NodeStates nodeStates : nodesStates)
            nodeStates.restore();
    }

//    @Override
//    public void unstore() {
//        System.arraycopy(storedMatrixIndex, 0, currentMatrixIndex, 0, getNodeCount());
//    }


    /**
     * indicate that the states matrix for node nodeIndex is to be changed *
     */
    public void setStatesForUpdate(int nodeIndex) {
        nodesStates[nodeIndex].setStatesForUpdate(); // 0 or 1
    }

    /**
     * Get an internal node states. The state starts from 0.
     *
     * @param nodeIndex the node index :<br>
     *      *                Leaf nodes are number 0 to <code>(nodeIndex+1)/2</code>;
     *      *                The root node is always numbered <code>nodeIndex-1</code>.
     * @return the codon states of the internal node sequence.
     */
    public int[] getNodeStates(final int nodeIndex) {
        // internal node index starts from getTipsCount();
        return nodesStates[nodeIndex].getNodeStates();
    }

    /**
     * Set the states to an internal node given its node index.
     * If {@link NodeStates} not init, then create with the given states.
     *
     * @param nodeIndex the node index :<br>
     *      *                Leaf nodes are number 0 to <code>(nodeIndex+1)/2</code>;
     *      *                The root node is always numbered <code>nodeIndex-1</code>.
     * @param states int[] states
     */
    public void setNodeStates(final int nodeIndex, final int[] states) {
        if (nodesStates[nodeIndex] == null) {
            // internal node not init
            nodesStates[nodeIndex] = new NodeStates(nodeIndex, states, getStateCount());
        }

        // internal node index starts from getTipsCount();
        setValue(nodeIndex, states);
    }

    /**
     * modify setValue for 2d matrix to take int[].
     *
     * @param nodeIndex the start index of parameter to set to the flattened matrix,
     *                 the values of the parameter is given by an int[].
     * @param vals     int[] values of the parameter
     */
    protected void setValue(final int nodeIndex, final int[] vals) {
        startEditing(null);

        nodesStates[nodeIndex].setValue(vals);
        nodeIsDirty[nodeIndex] = true;
//        m_nLastDirty = nodeIndex;
    }

    /**
     * get the sites at a site of all internal nodes.
     *
     * @param codonNr the codon site index.
     * @return
     */
    public int[] getSites(final int codonNr) {
        int[] col = new int[getNodeCount()];
        for (int i = 0; i < getNodeCount(); i++)
            col[i] = getASite(i, codonNr);
        return col;
    }

    /**
     * Get a codon state from the site at the node.
     * The state starts from 0.
     * <code>matrix[i,j] = values[i * minorDimension + j]</code>
     *
     * @param nodeIndex  the node index :<br>
     *                Leaf nodes are number 0 to <code>(nodeIndex+1)/2</code>;
     *                The root node is always numbered <code>nodeIndex-1</code>.
     * @param codonNr the site index.
     * @return
     */
    public int getASite(final int nodeIndex, final int codonNr) {
        return nodesStates[nodeIndex].getASite(codonNr);
    }


//    /**
//     * Set a codon state to the site at the internal node.
//     * The state starts from 0.
//     * <code>matrix[i,j] = values[i * minorDimension + j]</code>
//     *
//     * @param nodeNr the node index <code>i = nodeNr - internalNodeCount - 1</code><br>
//     *               Leaf nodes are number 0 to <code>leafnodes-1</code>;
//     *               Internal nodes are numbered  <code>leafnodes</code> up to <code>nodes-1</code>;
//     *               The root node is always numbered <code>nodes-1</code>.
//     * @param codonNr the site index.
//     */
//    public void setASite(final int nodeNr, final int codonNr, final int value) {
//        // internal node index nodeNr starts from getTaxonCount();
//        assert nodeNr > getInternalNodeCount() && nodeNr < getNodeCount();
//        // convert nodeNr into 2d matrix row index
//        int rowIndex = nodeNr - getInternalNodeCount() - 1;
//
//        setValue(rowIndex, codonNr, value);
//    }

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

    public int getTipsCount() {
        return (getNodeCount() + 1) / 2;
    }

    public int getInternalNodeCount() {
        return getTipsCount() - 1;
    }


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
