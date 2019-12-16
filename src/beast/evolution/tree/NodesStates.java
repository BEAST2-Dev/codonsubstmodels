package beast.evolution.tree;

import beast.core.Input;
import beast.core.StateNode;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.CodonAlignment;
import beast.evolution.datatype.Codon;
import beast.evolution.datatype.GeneticCode;
import beast.util.Randomizer;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;

/**
 * Store states for all nodes including tips in the array of {@link NodeStates},
 * which contains an array of states at all sites.
 * The state is an integer starting from 0. <br>
 * The array index is <code>nodeNr</code> which is the node index
 * used in {@link TreeInterface}.
 * The tip index starts from 0, and the index of an internal node is ranged
 * from the number of leaf nodes to the total nodes - 1.
 * Root index is always the last number.<br>
 * @author Walter Xie
 */
public class NodesStates extends StateNode {

    final public Input<CodonAlignment> codonAlignmentInput = new Input<>("data",
            "codon alignment to initialise the 2-d matrix of nodes (including tips) states",
            Input.Validate.REQUIRED);

    final public Input<TreeInterface> treeInput = new Input<>("tree",
            "phylogenetic beast.tree to map states data correctly to the nodes", Input.Validate.REQUIRED);

    final public Input<String> initINSInput = new Input<>("initINS",
            "whether and how to initialise the states of internal nodes. Choose 'random' or 'parsimony'",
            "parsimony", Input.Validate.OPTIONAL);

    /**
     * upper & lower bound These are located before the inputs (instead of
     * after the inputs, as usual) so that valuesInput can determines the
     * class
     */
    protected int upper;
    protected int lower;


    Codon codonDataType; // getStateCount
    int siteCount;
    int nodeCount; // == nodesStates.length

    // all nodes, use nodeNr to map the index of nodesStates[]
    protected  NodeStates[] nodesStates;
//    protected  NodeStates[] storedNodesStates;

    /**
     * isDirty array flags for each node
     */
    protected boolean[] nodeIsDirty;

    //for XML parser
    public NodesStates() { }

    public NodesStates(CodonAlignment codonAlignment, TreeInterface tree) {
        initByAlignment(codonAlignment);
        validateTree(tree);
        setTipsStates(codonAlignment, tree);
    }

    /**
     * Init class and states in tips given {@link CodonAlignment},
     * and generate states for internal nodes.
     * @param codonAlignment   data
     * @param tree             tree
     * @param initMethod       method to generate states for internal nodes.
     *                         see {@link #initInternalNodesStates(String, TreeInterface)}
     */
    public NodesStates(CodonAlignment codonAlignment, TreeInterface tree, String initMethod) {
        this(codonAlignment, tree);
        initInternalNodesStates(initMethod, tree);
        validateNodeStates();
    }

    @Override
    public void initAndValidate() {
        // need data type, site count, taxa count
        CodonAlignment codonAlignment = codonAlignmentInput.get();
        // cannot init params by constructor
        initByAlignment(codonAlignment);

        TreeInterface tree = treeInput.get();
        validateTree(tree);
        setTipsStates(codonAlignment, tree);

        // get BEAST seed long seed = Randomizer.getSeed();
        String initMethod = initINSInput.get();
        initInternalNodesStates(initMethod, tree);

        // validation
        validateNodeStates();
    }


    public void validateTree(TreeInterface tree) {
        // sanity check: alignment should have same #taxa as tree
        if (getTipsCount() != tree.getLeafNodeCount())
            throw new IllegalArgumentException("The number of tips in the tree does not match the number of sequences");
        // nodeCount == tree.getNodeCount()
        if (getNodeCount() != tree.getNodeCount())
            throw new IllegalArgumentException("The dimension of nodes states should equal to " +
                    "the number of nodes in the tree !\n" + getNodeCount() + " != " + tree.getNodeCount());
        // get***Count uses nodeCount
//        assert getTipsCount() == tree.getLeafNodeCount();
//        assert getInternalNodeCount() == tree.getInternalNodeCount();
    }

    public void validateNodeStates() {
        for (int i=0; i < getNodeCount(); i++) {
            NodeStates node = getNodeStates(i);
            assert node.getNodeNr() == i;

            if (i < getTipsCount()) { // tips ID
                // use it to recognise tips in setter
                if (node.getTipID() == null)
                    throw new IllegalArgumentException("Tip (" + i + ") requires ID taxon name in NodeStates object !");
                if (node.getNodeNr() >= getTipsCount())
                    throw new IllegalArgumentException("Invalid node index at tip : " + i + " cannot >= " + getTipsCount() + " !");
            } else { // internal nodes
                // use it to recognise tips in setter
                if (node.getTipID() != null)
                    throw new IllegalArgumentException("TipID should be null in internal node (" + i + ") NodeStates object !");
                if (node.getNodeNr() < getTipsCount() || node.getNodeNr() >= getNodeCount())
                    throw new IllegalArgumentException("Invalid node index at internal node : " + i +
                            " should range [" + getTipsCount() + ", " + (getNodeCount()-1) + "] !");
            }
        } // end i loop
    }

    // nodeCount = 2 * codonAlignment.getTaxonCount() - 1
    protected void initByAlignment(CodonAlignment codonAlignment) {
        this.codonDataType = codonAlignment.getDataType();
        final int stateCount = getStateCount();
//        assert stateCount == 64;
        // 0 - 63, ignore lowerValueInput upperValueInput
        lower = 0;
        upper = stateCount - 1; // not deal unknown codon

        // nodeCount defines rows of 2d matrix
        nodeCount = 2 * codonAlignment.getTaxonCount() - 1;
        assert nodeCount > 2;

        // siteCount defines cols of 2d matrix
        // siteCount = num of codon = nucleotides / 3, transformation in CodonAlignment
        siteCount = codonAlignment.getSiteCount();

        // init from constructor
        nodesStates = new NodeStates[nodeCount];
//        storedNodesStates = new NodeStates[nodeCount];
        nodeIsDirty = new boolean[nodeCount];

        // fix ID
        if (getID() == null || getID().length() < 1)
            setID(codonAlignment.getID());
    }

    // set tips states to nodesStates[] by nodeNr indexing
    protected void setTipsStates(CodonAlignment codonAlignment, TreeInterface tree) {
        // sanity check: alignment should have same #taxa as tree
        if (codonAlignment.getTaxonCount() != tree.getLeafNodeCount())
            throw new IllegalArgumentException("The number of nodes in the tree does not match the number of sequences");

        Log.info.println("  " + codonAlignment.toString(true));

        // set tips states
        for (Node tip : tree.getExternalNodes()) {
            int nr = tip.getNr();
            // use nodeNr to map the index of nodesStates[]
            nodesStates[nr] = new NodeStates(nr, tip.getID(), codonAlignment);
        }
        // call initINS
    }

    /**
     * Initialise states at internal nodes by the "random" or "parsimony" method.
     * The seed can be fixed by {@link Randomizer#setSeed(long)}.
     * @param initMethod   "random" or "parsimony"
     * @param tree         used for parsimony Fitch (1971) algorithm.
     */
    public void initInternalNodesStates(String initMethod, TreeInterface tree) {
        // internal nodes
        int[][] inStates;
        if ("random".equalsIgnoreCase(initMethod)) {
            inStates = initINStatesRandom(getInternalNodeCount(), getSiteCount(), getStateCount());

        } else if ("parsimony".equalsIgnoreCase(initMethod)) {
            inStates = initINStatesParsimony(tree);

        } else {
            throw new IllegalArgumentException("No method selected to initialise the states at internal nodes !");
        }
        assert inStates.length == getInternalNodeCount() && inStates[0].length == getSiteCount();

        for (int i=0; i < getInternalNodeCount(); i++) {
            int nR = i + getTipsCount();
            nodesStates[nR] = new NodeStates(nR, inStates[i], getStateCount());
        }

    }

    /**
     * Init states using parsimony Fitch (1971) algorithm.
     * Equally to choose a state from the ambiguous set.
     * The seed can be fixed by {@link Randomizer#setSeed(long)}.
     * @param tree         used for parsimony Fitch (1971) algorithm.
     * @return   states[internal nodes][sites], where i from 0 to internalNodeCount-1
     * @see RASParsimony1Site
     */
    public int[][] initINStatesParsimony(TreeInterface tree) {
        // states[internal nodes][sites], where i from 0 to internalNodeCount-1
        int[][] inStates = new int[getInternalNodeCount()][getSiteCount()];

        int[] tipsStates = new int[getTipsCount()];
        for (int k=0; k < getSiteCount(); k++) {
            getTipsStatesAtSite(k, tipsStates);
            RASParsimony1Site ras1Site = new RASParsimony1Site(tipsStates, tree);

            int[] as = ras1Site.reconstructAncestralStates();
            // internal nodes i from 0 to internalNodeCount-1
            for (int i=0; i < getInternalNodeCount(); i++)
                inStates[i][k] = as[i];
        }
        return inStates;
    }


    /**
     * Assuming rootIndex = tree.getNodeCount() - 1,
     * and then validate rootIndex == tree.getRoot().getNr()
     * @return {@link TreeInterface}
     */
//    public TreeInterface getTree() {
//        // exclude root node, branches = nodes - 1
//        final int rootIndex = tree.getNodeCount() - 1;
//        if (rootIndex != tree.getRoot().getNr())
//            throw new RuntimeException("Invalid root index " + tree.getRoot().getNr() + " != " + rootIndex);
//
//        return tree;
//    }

//    public int getRootIndex() {
//        int rootIndex = tree.getNodeCount() - 1;
//        assert rootIndex == tree.getRoot().getNr();
//        return rootIndex;
//    }

//    public TreeInterface getTree() {
//        return tree;
//    }


    /**
     * Randomly generate states at internal nodes given genetic code.
     * The seed can be fixed by {@link Randomizer#setSeed(long)}.
     * @param internalNodeCount  the number of internal nodes at the rows of 2d array.
     * @param siteCount          the number of codon sites at the columns of 2d array.
     * @param stateCount         the maximum state to generate.
     * @return    states[internal nodes][sites], where i from 0 to internalNodeCount-1
     */
    public int[][] initINStatesRandom(final int internalNodeCount, final int siteCount, final int stateCount) {

        Log.info("Randomly generate codon states using " + getGeneticCode().getDescription() + " for " +
                internalNodeCount + " internal nodes, " + getSiteCount() + " codon, seed = " + Randomizer.getSeed());

        // states[internal nodes][sites], where i from 0 to internalNodeCount-1
        int[][] inStates = new int[internalNodeCount][siteCount];
        for (int i=0; i < inStates.length; i++) {
            for (int j = 0; j < inStates[0].length; j++) {
                // 0 - 60/61, no stop codon
                inStates[i][j] = (int) (Randomizer.nextDouble() * stateCount);
            }
        }
        return inStates;
    }


    /**
     * For CodonAlignment, convert triplets string into list of codon states.
     * @param triplets coded nucleotides in string
     * @return the list of codon states. The size should be 1/3 of string length
     */
    public List<Integer> stringToEncoding(String triplets) {
        // remove spaces
        triplets = triplets.replaceAll("\\s", "");
        Alignment alignment = codonAlignmentInput.get();
        List<Integer> sequence = alignment.getDataType().stringToEncoding(triplets);
        if (sequence.size() * 3 != triplets.length())
            throw new IllegalArgumentException("The string of triplets has invalid number ! " +
                    triplets.length() + " != " + sequence.size() + " * 3");
        return sequence;
    }

    /**
     * @param triplets coded nucleotides in string
     * @return the int[] of codon states. The size should be 1/3 of string length
     */
    public int[] stringToStates(String triplets) {
        List<Integer> sequence = stringToEncoding(triplets);

        return sequence.stream().mapToInt(i -> i).toArray();
    }

    /**
     * indicate that the states matrix for node nodeIndex is to be changed *
     */
    public void setStatesForUpdate(int nodeIndex) {
        nodesStates[nodeIndex].setStatesForUpdate(); // 0 or 1
    }

    /**
     * Get an internal node states. The state starts from 0.
     *
     * @param nodeIndex The node index :<br>
     *                  Leaf nodes are indexed from 0 to <code>tipsCount - 1</code>.
     *                  Internal nodes are indexed from <code>tipsCount</code> to <code>nodeCount-1</code>.
     *                  The root node is always indexed <code>nodeCount-1</code>.
     * @return the codon states of the internal node sequence.
     */
    public int[] getStates(final int nodeIndex) {
        // internal node index starts from getTipsCount();
        return nodesStates[nodeIndex].getStates();
    }

    /**
     * Set the states to an internal node given its node index.
     * If {@link NodeStates} not init, then create with the given states.
     *
     * @param states int[] states
     * @param nodeIndex The node index :<br>
     *                  Leaf nodes are indexed from 0 to <code>tipsCount - 1</code>.
     *                  Internal nodes are indexed from <code>tipsCount</code> to <code>nodeCount-1</code>.
     *                  The root node is always indexed <code>nodeCount-1</code>.
     */
    public void setStates(final int[] states, final int nodeIndex) {
        if (nodesStates[nodeIndex] == null) {
            // internal node not init
            nodesStates[nodeIndex] = new NodeStates(nodeIndex, states, getStateCount());
        } else {
            // internal node index starts from getTipsCount();
            nodesStates[nodeIndex].setStates(states);
        }
        nodeIsDirty[nodeIndex] = true;
    }


    /**
     * Get tips states at given 1 site (codon).
     * @param nrOfSite      The codon site index in the alignment
     * @param tipsStates    States array to fill in. Tips nr = [0, TipsCount-1].
     */
    public void getTipsStatesAtSite(int nrOfSite, int[] tipsStates) {
        assert tipsStates.length == getTipsCount();
        int state;
        // BEAST tips nr = [0, TipsCount-1]
        for (int i=0; i < getTipsCount() ; i++) {
            state = getState(i, nrOfSite);
            tipsStates[i] = state;
        }
    }

//    /**
//     * get the sites at a site of all nodes.
//     *
//     * @param codonNr the codon site index.
//     * @return
//     */
//    public int[] getStatesAtSite(final int codonNr) {
//        int[] col = new int[getNodeCount()];
//        for (int i = 0; i < getNodeCount(); i++)
//            col[i] = getState(i, codonNr);
//        return col;
//    }

    /**
     * Get a codon state from the site at the node.
     * The state starts from 0.
     * <code>matrix[i,j] = values[i * minorDimension + j]</code>
     *
     * @param nodeIndex The node index :<br>
     *                  Leaf nodes are indexed from 0 to <code>tipsCount - 1</code>.
     *                  Internal nodes are indexed from <code>tipsCount</code> to <code>nodeCount-1</code>.
     *                  The root node is always indexed <code>nodeCount-1</code>.
     * @param codonNr the site index.
     * @return
     */
    public int getState(final int nodeIndex, final int codonNr) {
        return nodesStates[nodeIndex].getState(codonNr);
    }


    /**
     * Set a codon state to the site at the internal node.
     * The state starts from 0.
     * <code>matrix[i,j] = values[i * minorDimension + j]</code>
     *
     * @param state   new state to set
     * @param nodeIndex The node index :<br>
     *                  Leaf nodes are indexed from 0 to <code>tipsCount - 1</code>.
     *                  Internal nodes are indexed from <code>tipsCount</code> to <code>nodeCount-1</code>.
     *                  The root node is always indexed <code>nodeCount-1</code>.
     * @param codonNr the site index.
     */
    public void setState(final int state, final int nodeIndex, final int codonNr) {
        if (nodesStates[nodeIndex] == null)
            throw new IllegalArgumentException("Node (" + nodeIndex + ") states are not initiated !");

        // internal node index starts from getTipsCount();
        nodesStates[nodeIndex].setState(state, codonNr);
        nodeIsDirty[nodeIndex] = true;
    }

    /**
     * @param nodeIndex The node index :<br>
     *                  Leaf nodes are indexed from 0 to <code>tipsCount - 1</code>.
     *                  Internal nodes are indexed from <code>tipsCount</code> to <code>nodeCount-1</code>.
     *                  The root node is always indexed <code>nodeCount-1</code>.
     * @return if the states have been changed,
     *         but always false if node is a tip,
     *         where <code>nodeIndex < {@link #getTipsCount()}</code>.
     */
    public boolean isNodeStatesDirty(final int nodeIndex) {
        if (nodeIndex < getTipsCount())
            return false;
        return nodeIsDirty[nodeIndex];
    }


//    /**
//     * Internal node states changed
//     * @return
//     */
//    @Override
//    protected boolean requiresRecalculation() {
//        // internal nodes only
//        for (int i = getTipsCount(); i < getNodeCount(); i++) {
//            if (isNodeStatesDirty(i)) return true;
//        }
//        return false;
//        //return super.requiresRecalculation();
//    }


    public NodeStates getNodeStates(final int nodeIndex) {
        return nodesStates[nodeIndex];
    }

    public Codon getCodonDataType() {
        return codonDataType;
    }

    public int getStateCount() {
        return codonDataType.getStateCount();
    }

    public GeneticCode getGeneticCode() {
        return codonDataType.getGeneticCode();
    }

    public int getSiteCount() {
        return siteCount;
    }

    /**
     * equal to the length of {@link NodeStates}[],
     * which is initialised by 2 * {@link CodonAlignment#getTaxonCount()} - 1.
     * @return    number of nodes that have states
     */
    public int getNodeCount() {
        return nodeCount;
    }

    /**
     * ({@link #getNodeCount()} + 1) / 2
     */
    public int getTipsCount() {
        return (getNodeCount() + 1) / 2;
    }

    /**
     * {@link #getTipsCount()} - 1
     */
    public int getInternalNodeCount() {
        return getTipsCount() - 1;
    }


    //******* MCMC StateNode *******
    @Override
    public void log(long sampleNr, PrintStream out) {
        super.log(sampleNr, out);
    }

    /**
     * all tree nodes at the moment
     */
    @Override
    public String toString() {
        final StringBuilder buf = new StringBuilder();
        buf.append(getID()).append("[").
                append(" ").append(getNodeCount());
        buf.append(" ").append(getSiteCount());
        buf.append("] ");
        buf.append("(").append(lower).append(",").append(upper).append("): ");
//        for (final int value : states) {
//            buf.append(value).append(" ");
//        }
        return buf.toString();
    }

    @Override
    public void init(PrintStream out) {
        final int nodeCount = getNodeCount();
        if (nodeCount == 1) {
            out.print(getID() + "\t");
        } else {
            for (int i = 0; i < nodeCount; i++) {
                out.print(getID() + (i + 1) + "\t");
            }
        }
    }

    @Override
    public void close(PrintStream out) {
// nothing to do
    }

    @Override
    public int getDimension() { // use getInternalNodeCount() & getSiteCount()
        throw new UnsupportedOperationException("Unsupported");
    }

    @Override
    public double getArrayValue(int dim) { // use getASite(int, int)
        throw new UnsupportedOperationException("Unsupported");
    }

    /**
     * this will only set internal nodes dirty in array {@link #nodeIsDirty}
     * @param isDirty
     */
    @Override
    public void setEverythingDirty(boolean isDirty) {
        setSomethingIsDirty(isDirty);
//        Arrays.fill(nodeIsDirty, isDirty);
        for (int i=getTipsCount(); i < getNodeCount(); i++)
            nodeIsDirty[i] = true;
    }

    //cleans up and deallocates arrays
    @Override
    public void finalize() throws Throwable {
        for (NodeStates nodeStates : nodesStates)
            nodeStates.finalize();
        nodeIsDirty = null;
    }


    /**
     * only call internal nodes {@link NodeStates#store()}
     */
    //TODO store/restore per site
    @Override
    protected void store() {
        // internal nodes only
        for (int i = getTipsCount(); i < getNodeCount(); i++){
            //TODO need nodeIsDirty[]?
//            storedNodesStates[i] = nodesStates[i].copy();
            nodesStates[i].store();
        }



    }

    /**
     * only call internal nodes {@link NodeStates#restore()}
     */
    @Override
    public void restore() {
        // internal nodes only
        for (int i = getTipsCount(); i < getNodeCount(); i++){
            //TODO need nodeIsDirty[]?
//            nodesStates[i] = storedNodesStates[i].copy();
            nodesStates[i].restore();
        }

        hasStartedEditing = false;
    }

//    @Override
//    public void unstore() {
//        System.arraycopy(storedMatrixIndex, 0, currentMatrixIndex, 0, getNodeCount());
//    }




    //******* use those in NodesStates not these below *******

    //TODO full copy or just currentMatrixIndex?
    @Override
    public StateNode copy() {
        try {
            @SuppressWarnings("unchecked")
            final NodesStates copy = (NodesStates) this.clone();
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
        final NodesStates copy = (NodesStates) other;
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
        final NodesStates source = (NodesStates) other;
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
        final NodesStates source = (NodesStates) other;
        for (int i = 0; i < nodesStates.length; i++)
            nodesStates[i].assignFromFragile(source.nodesStates[i]);
        Arrays.fill(nodeIsDirty, false);
    }

    @Override
    public void fromXML(org.w3c.dom.Node node) {
        throw new UnsupportedOperationException("Unsupported");
    }

    @Override
    public int scale(double scale) {
        throw new UnsupportedOperationException("Unsupported");
    }

}
