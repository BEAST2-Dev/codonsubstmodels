package beast.evolution.tree;

import beast.core.StateNode;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.CodonAlignment;
import beast.evolution.datatype.DataType;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;

/**
 * The large array to store one node states (can be a tip).<br>
 *
 * <code>nodeNr</code> used for all indexing system
 * is the node index number starting from 0 to root index (nodeCount - 1).
 * For internal nodes, they are ranged from the number of leaf nodes
 * to the nodeCount - 1.<br>
 *
 * Not use {@link beast.core.parameter.IntegerParameter}
 * because of both <code>arraycopy</code> and storing Object are slow.
 * Apple the index trick implemented in {@link beast.evolution.likelihood.BeerLikelihoodCore}
 * to avoid <code>arraycopy</code>.<br>
 *
 * No {@link beast.core.Input} here.
 * Init and validate in {@link NodesStates#initAndValidate()}.
 */
public class NodeStates extends StateNode {

    /**
     * The node index :<br>
     * Leaf nodes are indexed from 0 to <code>tipsCount - 1</code>.
     * Internal nodes are indexed from <code>tipsCount</code> to <code>nodeCount-1</code>.
     * The root node is always indexed <code>nodeCount-1</code>.
     */
    private final int nodeNr; // used to map to index of DABranchLikelihoodCore[]
    /**
     * taxon name, null if not tip. use tipID to recognise tips
     */
    private final String tipID;
//    final GeneticCode geneticCode;

    /**
     * states[matrix index][sites]
     * 1st dimension is matrix index (current, stored),
     * 2nd is nrOfSites.
     */
    protected int[][] states;
    /**
     * store the matrix index, instead of different matrices
     */
    protected int currentMatrixIndex;
    protected int storedMatrixIndex;

    /**
     * upper & lower bound These are located before the inputs (instead of
     * after the inputs, as usual) so that valuesInput can determines the
     * class
     */
    protected int upper;
    protected int lower;
    /**
     * isDirty flags for each site
     */
    protected boolean[] siteIsDirty; // TODO setASiteDirty


    /**
     * Set states to internal nodes and initialise parameters.
     * Call {@link #setStates(int[])}.
     * @param nodeNr       nr from {@link Node#getNr()}
     * @param states       states[]
     * @param stateCount   maximum state, e.g. 60
     */
    public NodeStates(int nodeNr, int[] states, int stateCount) {
        this.nodeNr = nodeNr;
        this.tipID = null;

        // init from constructor
        initParam(stateCount, states.length);
        // set internal nodes or tips states after initParam
        setStates(states);
    }

    /**
     * set states to tips given {@link CodonAlignment}.
     * Because nodeNr != taxonIndex, call {@link #setTipStates(String, CodonAlignment)}.
     * @param nodeNr       nr from {@link Node#getNr()}
     * @param tipID              Leaf {@link Node#getID()}, also taxon name
     * @param codonAlignment    {@link CodonAlignment}
     */
    public NodeStates(int nodeNr, String tipID, CodonAlignment codonAlignment) {
        this.nodeNr = nodeNr;
        this.tipID = tipID;

        final int stateCount = codonAlignment.getDataType().getStateCount();

        // siteCount = num of codon = nucleotides / 3, transformation in CodonAlignment
        final int siteCount = codonAlignment.getSiteCount();

        initParam(stateCount, siteCount);
        // Note: nodeNr != taxonIndex
        setTipStates(tipID, codonAlignment);
    }

    @Override
    public void initAndValidate() {
        //init validate in NodesStates.
    }

    // stateCount -> upper, siteCount -> array
    private void initParam(int stateCount, int siteCount) {
        // siteCount = num of codons, overwrite in CodonAlignment /= 3
        states = new int[2][siteCount];

//        currentMatrixIndex = 0;
//        storedMatrixIndex = 0;

        // 0 - 60/61, ignore lowerValueInput upperValueInput
        lower = 0;
        upper = stateCount - 1; // not deal unknown codon

        siteIsDirty = new boolean[siteCount];

//        Log.info.println("Create node states array length = " + siteCount);

    }

    /**
     * set tips states given {@link CodonAlignment}. <br>
     * Note: nodeNr != taxonIndex.<br>
     * The sequence index (taxonIndex) in {@link CodonAlignment} is different to
     * the node index (nodeNr) in {@link Tree}.
     *
     * @param tipID              Leaf {@link Node#getID()}, also taxon name
     * @param codonAlignment    {@link CodonAlignment}
     */
    public void setTipStates(String tipID, CodonAlignment codonAlignment) {
        // make sure to use taxonIndex to getCounts, nodeNr not mapping to the order in List<Sequence>
        int taxonIndex = getTaxonIndex(tipID, codonAlignment);

        // no patterns
        List<Integer> statesList = codonAlignment.getCounts().get(taxonIndex);
        assert statesList.size() == states[currentMatrixIndex].length;
        // Java 8
        states[currentMatrixIndex] = statesList.stream().mapToInt(i->i).toArray();
    }

    /**
     * Copied from
     * @param taxon the taxon name as a string
     * @param data the alignment
     * @return the taxon index of the given taxon name for accessing its sequence data in the given alignment,
     *         throw RuntimeException when -1 if the taxon is not in the alignment.
     */
    public int getTaxonIndex(String taxon, Alignment data) {
        int taxonIndex = data.getTaxonIndex(taxon);
        if (taxonIndex == -1) {
            if (taxon.startsWith("'") || taxon.startsWith("\"")) {
                taxonIndex = data.getTaxonIndex(taxon.substring(1, taxon.length() - 1));
            }
            if (taxonIndex == -1) {
                throw new RuntimeException("Could not find sequence " + taxon + " in the alignment");
            }
        }
        return taxonIndex;
    }

    public int getNodeNr() {
        return nodeNr;
    }

    /**
     * @return   node index in BEAST tree
     */
    @Override
    public String getID() {
        return Integer.toString(getNodeNr());
    }

    /**
     * @return  null if not tip, otherwise return taxon name
     */
    public String getTipID() {
        return tipID;
    }

    /**
     * indicate that the states matrix for node nodeIndex is to be changed *
     */
    public void setStatesForUpdate() {
        if (tipID != null)
            throw new UnsupportedOperationException("Can not sample tip states yet !");
        currentMatrixIndex = 1 - currentMatrixIndex; // 0 or 1
    }



    /**
     * cleans up and deallocates arrays.
     */
    @Override
    public void finalize() throws Throwable {
        states = null;
        siteIsDirty = null;
        currentMatrixIndex = 0;
        storedMatrixIndex = 0;
    }


    //TODO store/restore per site
    @Override
    protected void store() {
        if (tipID != null)
            throw new UnsupportedOperationException("Can not sample tip states yet !");

        storedMatrixIndex = currentMatrixIndex;
    }

    @Override
    public void restore() {
        if (tipID != null)
            throw new UnsupportedOperationException("Can not sample tip states yet !");

        // Rather than copying the stored stuff back, just swap the pointers...
        currentMatrixIndex = storedMatrixIndex;

        hasStartedEditing = false;
        if (siteIsDirty.length != getSiteCount()) {
            siteIsDirty = new boolean[getSiteCount()];
        }

    }

//    @Override
//    public void unstore() {
//        currentMatrixIndex = storedMatrixIndex;
//    }

    public int getSiteCount() {
        return states[currentMatrixIndex].length;
    }

    /**
     * Get an internal node states. The state starts from 0.
     *
     * @return the codon states of the internal node sequence.
     */
    public int[] getStates() {
        // internal node index starts from getTipsCount();
        return states[currentMatrixIndex];
    }

    /**
     * Set the states to an internal node.
     * The node index has to convert to the array index before setValue.
     *
     * @param states int[]
     */
    public void setStates(final int[] states) {
        // internal node index starts from getTipsCount();
        startEditing(null);

        System.arraycopy(states, 0, this.states[currentMatrixIndex], 0, states.length);
        Arrays.fill(siteIsDirty, true);
    }

    /**
     * Get a codon state from the site at the node.
     * The state starts from 0.
     *
     * @param codonNr the site index.
     * @return
     */
    public int getState(final int codonNr) {
        return states[currentMatrixIndex][codonNr];
    }


    /**
     * Set a codon state to the site at the internal node.
     * The state starts from 0.
     * <code>matrix[i,j] = values[i * minorDimension + j]</code>
     *
     * @param state   new state to set.
     * @param codonNr the site index.
     */
    public void setState(final int state, final int codonNr) {
        startEditing(null);

        states[currentMatrixIndex][codonNr] = state;
        siteIsDirty[codonNr] = true;
    }


    /**
     * For CodonAlignment, convert triplets string into list of codon states.
     * @param triplets coded nucleotides in string
     * @param dataType to parse triplets
     * @return the list of codon states. The size should be 1/3 of string length
     */
    public List<Integer> stringToEncoding(String triplets, DataType dataType) {
        // remove spaces
        triplets = triplets.replaceAll("\\s", "");
        List<Integer> sequence = dataType.stringToEncoding(triplets);
        if (sequence.size() * 3 != triplets.length())
            throw new IllegalArgumentException("The string of triplets has invalid number ! " +
                    triplets.length() + " != " + sequence.size() + " * 3");
        return sequence;
    }

    /**
     * @param triplets coded nucleotides in string
     * @return the int[] of codon states. The size should be 1/3 of string length
     */
    public int[] stringToStates(String triplets, DataType dataType) {
        List<Integer> sequence = stringToEncoding(triplets, dataType);

        return sequence.stream().mapToInt(i -> i).toArray();
    }


    @Override
    public void init(PrintStream out) {
        out.print(nodeNr + "\t");
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
                append(" ").append(nodeNr);
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
        Arrays.fill(siteIsDirty, isDirty);
    }

    //TODO full copy or just currentMatrixIndex?
    @Override
    public NodeStates copy() {
        try {
            @SuppressWarnings("unchecked")
            final NodeStates copy = (NodeStates) this.clone();
//            copy.states = states.clone(); // always [2][][]
            copy.states[0] = states[0].clone();
            copy.states[1] = states[1].clone();

            //storedStates = copy.storedStates.clone();
            copy.siteIsDirty = new boolean[getSiteCount()];
            return copy;
        } catch (Exception e) {
            e.printStackTrace();
        }
        return null;

    }

    @Override
    public void assignTo(StateNode other) {
        @SuppressWarnings("unchecked")
        final NodeStates copy = (NodeStates) other;
        copy.setID(getID());
        copy.index = index;

        copy.states[0] = states[0].clone();
        copy.states[1] = states[1].clone();

        //storedStates = copy.storedStates.clone();
        copy.lower = lower;
        copy.upper = upper;
        copy.siteIsDirty = new boolean[getSiteCount()];
    }

    @Override
    public void assignFrom(StateNode other) {
        @SuppressWarnings("unchecked")
        final NodeStates source = (NodeStates) other;
        setID(source.getID());
//        states = source.states.clone();
//        storedStates = source.storedStates.clone();

        states[0] = source.states[0].clone();
        states[1] = source.states[1].clone();

        lower = source.lower;
        upper = source.upper;
        siteIsDirty = new boolean[source.getSiteCount()];
    }

    @Override
    public void assignFromFragile(StateNode other) {
        @SuppressWarnings("unchecked")
        final NodeStates source = (NodeStates) other;

        System.arraycopy(source.states[0], 0, states[0], 0, source.states[0].length);
        System.arraycopy(source.states[1], 0, states[1], 0, source.states[1].length);

        Arrays.fill(siteIsDirty, false);
    }

    @Override
    public void fromXML(org.w3c.dom.Node xmlNode) {
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
