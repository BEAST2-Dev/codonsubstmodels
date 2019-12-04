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
 * Code to deal with node states without requiring tree.
 *
 * @author Walter Xie
 */
public class NodesStates extends StateNode {

    final public Input<CodonAlignment> codonAlignmentInput = new Input<>("codon",
            "codon alignment to initialise the 2-d matrix of nodes (including tips) states",
            Input.Validate.REQUIRED);

    /**
     * upper & lower bound These are located before the inputs (instead of
     * after the inputs, as usual) so that valuesInput can determines the
     * class
     */
    protected int upper;
    protected int lower;


    final Codon codonDataType; // getStateCount
    final int siteCount;
    final int nodeCount; // == nodesStates.length

    // all nodes, use nodeNr to map the index of nodesStates[]
    protected  NodeStates[] nodesStates;

    // isDirty flags for each node
    protected boolean[] nodeIsDirty;


    public NodesStates(CodonAlignment codonAlignment) {
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

        nodeIsDirty = new boolean[nodeCount];

        // fix ID
        if (getID() == null || getID().length() < 1)
            setID(codonAlignment.getID());
    }


    @Override
    public void initAndValidate() {
        // init in NodesStatesAndTree using its constructor
        throw new IllegalArgumentException("Not used");
//        // need data type, site count, taxa count
//        CodonAlignment codonAlignment = codonAlignmentInput.get();
//        // initParam
//        NodesStates nodesStates = new NodesStates(codonAlignment);
    }

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

    @Override
    public void setEverythingDirty(boolean isDirty) {
        setSomethingIsDirty(isDirty);
        Arrays.fill(nodeIsDirty, isDirty);
    }

    //cleans up and deallocates arrays
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

    //******* use those in NodesStatesAndTree not these below *******

    @Override
    public StateNode copy() {
        throw new UnsupportedOperationException("Unsupported");
    }

    @Override
    public void assignTo(StateNode other) {
        throw new UnsupportedOperationException("Unsupported");
    }

    @Override
    public void assignFrom(StateNode other) {
        throw new UnsupportedOperationException("Unsupported");
    }

    @Override
    public void assignFromFragile(StateNode other) {
        throw new UnsupportedOperationException("Unsupported");
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
