package beast.evolution.tree;

import beast.core.Input;
import beast.core.StateNode;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.CodonAlignment;
import beast.evolution.datatype.GeneticCode;
import org.w3c.dom.Node;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * The large 2-d matrix to store all node states including tips.
 * Rows are nodes, cols are sites.
 * The values are flattened into 1d <code>int[]</code>,
 * and the state starts from 0. <br>
 *
 * <code>minorDimension</code> is the number of codons.<br>
 * <code>internalNodeCount</code> is the number of internal nodes.<br>
 * <code>nodeNr</code> is the node index number starting from 0,
 * and for internal nodes, they are ranged from the number of leaf nodes
 * to the total nodes - 1.<br>
 * So use the following formula to convert index:<br>
 * <code>rowIndex = nodeNr - internalNodeCount - 1</code><br>
 */
public class NodesStates extends StateNode {

    final public Input<CodonAlignment> codonAlignmentInput = new Input<>("codon",
            "codon alignment to initialise the 2-d matrix of nodes (including tips) states",
            Input.Validate.REQUIRED);

    final public Input<String> initINInput = new Input<>("initIN",
            "whether and how to initialise the states of internal nodes. Choose 'random' or 'parsimony'",
            Input.Validate.OPTIONAL);

    final public Input<Long> initINSeedInput = new Input<>("seed",
            "seed to initialise the states of internal nodes if choosing 'random' method",
            Input.Validate.OPTIONAL);

    final GeneticCode geneticCode;

    // 1st dimension is matrix index (current, stored),
    // 2nd is node index,
    // 3rd is nrOfSites
    // states[matrix index][node index][sites]
    protected int[][][] states;
    protected int[][][] storedStates;
    // store the matrix index, instead of different matrices
    protected int[] currentMatrixIndex; // node count
    protected int[] storedMatrixIndex;

    /**
     * upper & lower bound These are located before the inputs (instead of
     * after the inputs, as usual) so that valuesInput can determines the
     * class
     */
    protected int upper;
    protected int lower;
    /**
     * isDirty flags for each node
     */
    protected boolean[] nodeIsDirty;


    /**
     * init class and states in tips given {@link CodonAlignment}.
     * @param codonAlignment
     */
    public NodesStates(CodonAlignment codonAlignment) {
        final int stateCount = codonAlignment.getDataType().getStateCount();
        assert stateCount == 64;

        this.geneticCode = codonAlignment.getDataType().getGeneticCode();

        // nodeCount defines rows of 2d matrix
        final int nodeCount = 2 * codonAlignment.getTaxonCount() - 1;
        assert nodeCount > 2;

        // siteCount defines cols of 2d matrix
        // siteCount = num of codon = nucleotides / 3, transformation in CodonAlignment
        final int siteCount = codonAlignment.getSiteCount();

        // init from constructor
        initParam(stateCount, nodeCount, siteCount);
        // set tips states after initParam
        setTipStates(codonAlignment);

    }

    /**
     * The large 2-d matrix to store node states.
     *
     * @param stateCount        to define upper
     * @param nodeCount         to define rows of 2d matrix
     * @param siteCount         to define cols of 2d matrix
     */
    public NodesStates(int stateCount, int nodeCount, int siteCount, GeneticCode geneticCode) {
        this.geneticCode = geneticCode;
        initParam(stateCount, nodeCount, siteCount);
    }

    @Override
    public void initAndValidate() {
        // need data type, site count, taxa count
        CodonAlignment codonAlignment = codonAlignmentInput.get();
        final int stateCount = codonAlignment.getDataType().getStateCount();
        assert stateCount == 64;

        GeneticCode geneticCode = codonAlignment.getDataType().getGeneticCode();

        // nodeCount defines rows of 2d matrix
        final int nodeCount = 2 * codonAlignment.getTaxonCount() - 1;
        assert nodeCount > 2;

        // siteCount defines cols of 2d matrix
        // siteCount = num of codon = nucleotides / 3, transformation in CodonAlignment
        final int siteCount = codonAlignment.getSiteCount();

        // init from constructor
        NodesStates nodesStates = new NodesStates(stateCount, nodeCount, siteCount, geneticCode);
        // set tips states after initParam
        nodesStates.setTipStates(codonAlignment);

        //TODO init internal nodes states
        if ("random".equalsIgnoreCase(initINInput.get())) {
            long seed = initINSeedInput.get();
            nodesStates.initINSates(seed);
        } else if ("parsimony".equalsIgnoreCase(initINInput.get())) {
            nodesStates.initINSates(-1);
        } else {
            throw new IllegalArgumentException("The internal nodes states have to be initialised !");
        }
    }

    // stateCount -> upper, nodeCount * siteCount -> 2d matrix
    private void initParam(int stateCount, int nodeCount, int siteCount) {
        // siteCount = num of codons, overwrite in CodonAlignment /= 3
        states = new int[2][nodeCount][siteCount];
        storedStates = new int[2][nodeCount][siteCount];

        currentMatrixIndex = new int[nodeCount];
        storedMatrixIndex = new int[nodeCount];

        // 0 - 63, ignore lowerValueInput upperValueInput
        lower = 0;
        upper = stateCount - 1; // not deal unknown codon

        nodeIsDirty = new boolean[nodeCount];

        Log.info.println("Create nodes states matrix : " + nodeCount + " * " + siteCount);


    }

    // set tips states after initParam
    public void setTipStates(CodonAlignment codonAlignment) {
        final int tipsCount = codonAlignment.getTaxonCount();

        for (int taxonIndex = 0; taxonIndex < tipsCount; taxonIndex++) {
            // no patterns
            List<Integer> statesList = codonAlignment.getCounts().get(taxonIndex);

            assert statesList.size() == codonAlignment.getSiteCount();
            // Java 8
            states[currentMatrixIndex[taxonIndex]][taxonIndex] = statesList.stream().mapToInt(i->i).toArray();
        }
    }


    // init internal nodes states
    public void initINSates(long seed) {
        // init internal nodes states
        if (seed < 0) {
            initINStatesParsimony();
        } else {
            initINStatesRandom(getInternalNodeCount(), getSiteCount(), geneticCode, seed);
        }
    }

    /**
     * random states given genetic code
     */
    protected void initINStatesRandom(int internalNodeCount, int siteCount, GeneticCode geneticCode, long seed) {
        final int stateCount = 64;
        Random generator = new Random(seed);
        // vertebrateMitochondrial: 8   10   48   50
//        int[] stopCodons = geneticCode.getStopCodonStates();

        Log.info("Random generate codon states using " + geneticCode.getDescription() +
                " for " + internalNodeCount + " internal nodes, " + siteCount + " codon, seed = " + seed);

        // internal nodes, i from NrTips to NrRoot
        for (int i=getTipsCount(); i < getNodeCount() ; i++) {
            // states[matrix index][node index][sites]
            int[] inStates = states[currentMatrixIndex[i]][i];
            for (int j=0; j < inStates.length; j++) {
                // 0 - 63
                inStates[j] = (int)(generator.nextDouble() * stateCount);
                // skip stop codon states
                while(geneticCode.isStopCodon(inStates[j]))
                    inStates[j] = (int)(generator.nextDouble() * stateCount);
            }
        }// end i loop

    }

    /**
     * Parsimony to init states. Equally to choose a state from the ambiguous set.
     */
    protected void initINStatesParsimony() {

        throw new UnsupportedOperationException();

        //traverse 1


        //traverse 2


    }




    /**
     * indicate that the states matrix for node nodeIndex is to be changed *
     */
    public void setStatesForUpdate(int nodeIndex) {
        currentMatrixIndex[nodeIndex] = 1 - currentMatrixIndex[nodeIndex]; // 0 or 1
    }



    /**
     * cleans up and deallocates arrays.
     */
    @Override
    public void finalize() throws Throwable {
        states = null;
        storedStates = null;
        nodeIsDirty = null;
        currentMatrixIndex = null;
        storedMatrixIndex = null;
    }


    //TODO store restore a site
    @Override
    protected void store() {
        System.arraycopy(currentMatrixIndex, 0, storedMatrixIndex, 0, getNodeCount());
    }

    @Override
    public void restore() {
        // Rather than copying the stored stuff back, just swap the pointers...
        int[] tmp1 = currentMatrixIndex;
        currentMatrixIndex = storedMatrixIndex;
        storedMatrixIndex = tmp1;

        hasStartedEditing = false;
        if (nodeIsDirty.length != getNodeCount()) {
            nodeIsDirty = new boolean[getNodeCount()];
        }

    }

//    @Override
//    public void unstore() {
//        System.arraycopy(storedMatrixIndex, 0, currentMatrixIndex, 0, getNodeCount());
//    }


    public int getTipsCount() {
        return (states[currentMatrixIndex[0]].length + 1) / 2;
    }

    public int getInternalNodeCount() {
        return getTipsCount() - 1;
    }

    public int getNodeCount() {
        return states[currentMatrixIndex[0]].length;
    }

    public int getSiteCount() {
        return states[currentMatrixIndex[0]][0].length;
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
        return states[currentMatrixIndex[nodeIndex]][nodeIndex];
    }

    /**
     * Set the states to an internal node.
     * The node index has to convert to the array index before setValue.
     *
     * @param nodeIndex the node index :<br>
     *      *                Leaf nodes are number 0 to <code>(nodeIndex+1)/2</code>;
     *      *                The root node is always numbered <code>nodeIndex-1</code>.
     * @param states int[]
     */
    public void setNodeStates(final int nodeIndex, final int[] states) {
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

        System.arraycopy(vals, 0, this.states[currentMatrixIndex[nodeIndex]][nodeIndex], 0, vals.length);
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
            col[i] = states[currentMatrixIndex[i]][i][codonNr];
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
        return states[currentMatrixIndex[nodeIndex]][nodeIndex][codonNr];
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
    public void log(long sampleNr, PrintStream out) {
        super.log(sampleNr, out);
    }

    @Override
    public void close(PrintStream out) {
        // nothing to do
    }

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
    public void setEverythingDirty(boolean isDirty) {
        setSomethingIsDirty(isDirty);
        Arrays.fill(nodeIsDirty, isDirty);
    }

    @Override
    public StateNode copy() {
        try {
            @SuppressWarnings("unchecked")
            final NodesStates copy = (NodesStates) this.clone();
//            copy.states = states.clone(); // always [2][][]
            for (int i = 0; i < getNodeCount(); i++) {
                copy.states[0][i] = states[0][i].clone();
                copy.states[1][i] = states[1][i].clone();
            }

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
        for (int i = 0; i < getNodeCount(); i++) {
            copy.states[0][i] = states[0][i].clone();
            copy.states[1][i] = states[1][i].clone();
        }
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
//        states = source.states.clone();
//        storedStates = source.storedStates.clone();
        for (int i = 0; i < getNodeCount(); i++) {
            states[0][i] = source.states[0][i].clone();
            states[1][i] = source.states[1][i].clone();
            storedStates[0][i] = source.storedStates[0][i].clone();
            storedStates[1][i] = source.storedStates[1][i].clone();
        }
        lower = source.lower;
        upper = source.upper;
        nodeIsDirty = new boolean[source.getNodeCount()];
    }

    @Override
    public void assignFromFragile(StateNode other) {
        @SuppressWarnings("unchecked")
        final NodesStates source = (NodesStates) other;
        for (int i = 0; i < getNodeCount(); i++) {
            System.arraycopy(source.states[0][i], 0, states[0][i], 0, source.states[0][i].length);
            System.arraycopy(source.states[1][i], 0, states[1][i], 0, source.states[1][i].length);
        }
        Arrays.fill(nodeIsDirty, false);
    }

    @Override
    public void fromXML(Node node) {
        throw new UnsupportedOperationException("in dev");
    }

    @Override
    public int scale(double scale) {
        throw new UnsupportedOperationException("in dev");
    }

    @Override
    public int getDimension() { // use getInternalNodeCount() & getSiteCount()
        throw new UnsupportedOperationException("Unsupported");
    }

    @Override
    public double getArrayValue(int dim) { // use getASite(int, int)
        throw new UnsupportedOperationException("Unsupported");
    }

}
