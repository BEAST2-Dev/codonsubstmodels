package beast.evolution.tree;

import beast.core.StateNode;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.CodonAlignment;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.GeneticCode;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * The large array to store one node states (can be a tip).
 * Cols are sites,
 * and the state starts from 0. <br>
 *
 * <code>nodeNr</code> is the node index number starting from 0,
 * and for internal nodes, they are ranged from the number of leaf nodes
 * to the total nodes - 1.<br>
 */
public class NodeStates extends StateNode {

    private final int nodeNr;

    final GeneticCode geneticCode;

    // 1st dimension is matrix index (current, stored),
    // 2nd is nrOfSites
    // states[matrix index][sites]
    protected int[][] states;
    protected int[][] storedStates;
    // store the matrix index, instead of different matrices
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
    protected boolean[] siteIsDirty;



    public NodeStates(int nodeNr, int[] states, int stateCount) {
        this.nodeNr = nodeNr;
        this.geneticCode = null;

        // init from constructor
        assert stateCount == 64;
        initParam(stateCount, states.length);
        // set tips states after initParam
        setNodeStates(states);
    }

    // for tips
    public NodeStates(Node tip, CodonAlignment codonAlignment) {
        this.nodeNr = tip.getNr();
        this.geneticCode = codonAlignment.getDataType().getGeneticCode();

        final int stateCount = codonAlignment.getDataType().getStateCount();
        assert stateCount == 64;

        // siteCount = num of codon = nucleotides / 3, transformation in CodonAlignment
        final int siteCount = codonAlignment.getSiteCount();

        initParam(stateCount, siteCount);
        // make sure the taxon maps nodeNr
        int taxonIndex = getTaxonIndex(tip.getID(), codonAlignment);
        setTipStates(codonAlignment, taxonIndex);

    }
    /**
     *
     * @param taxon the taxon name as a string
     * @param data the alignment
     * @return the taxon index of the given taxon name for accessing its sequence data in the given alignment,
     *         throw RuntimeException when -1 if the taxon is not in the alignment.
     */
    private int getTaxonIndex(String taxon, Alignment data) {
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

    // for internal nodes
    public NodeStates(int nodeNr, int stateCount, int siteCount, GeneticCode geneticCode, long seed) {
        this.nodeNr = nodeNr;
        this.geneticCode = geneticCode;

        initParam(stateCount, siteCount);
        initINStatesRandom(stateCount, seed);
    }

    // for internal nodes
    public NodeStates(int nodeNr, int stateCount, int siteCount) {
        this.nodeNr = nodeNr;
        this.geneticCode = null;

        initParam(stateCount, siteCount);
        //TODO
        initINStatesParsimony();
    }

    @Override
    public void initAndValidate() {
        //do nothing
    }

    // stateCount -> upper, siteCount -> array
    private void initParam(int stateCount, int siteCount) {
        // siteCount = num of codons, overwrite in CodonAlignment /= 3
        states = new int[2][siteCount];
        storedStates = new int[2][siteCount];

//        currentMatrixIndex = 0;
//        storedMatrixIndex = 0;

        // 0 - 63, ignore lowerValueInput upperValueInput
        lower = 0;
        upper = stateCount - 1; // not deal unknown codon

        siteIsDirty = new boolean[siteCount];

        Log.info.println("Create node states array length = " + siteCount);


    }

    // set tips states after initParam
    public void setTipStates(CodonAlignment codonAlignment, final int taxonIndex) {
        assert taxonIndex == nodeNr;

        // no patterns
        List<Integer> statesList = codonAlignment.getCounts().get(taxonIndex);
        assert statesList.size() == states[currentMatrixIndex].length;
        // Java 8
        states[currentMatrixIndex] = statesList.stream().mapToInt(i->i).toArray();
    }

    /**
     * random states given genetic code
     */
    public void initINStatesRandom(final int stateCount, final long seed) {

        Log.info("Random generate codon states using " + geneticCode.getDescription() +
                " for internal node " + nodeNr + ", " + getSiteCount() + " codon, seed = " + seed);

        Random generator;
        if (seed > 0)
            generator = new Random(seed);
        else
            generator = new Random();

        // internal nodes, i from NrTips to NrRoot
        // states[matrix index][sites]
        int[] inStates = states[currentMatrixIndex];
        for (int j=0; j < inStates.length; j++) {
            // 0 - 63
            inStates[j] = (int)(generator.nextDouble() * stateCount);
            // skip stop codon states, such as vertebrateMitochondrial: 8  10  48  50
            while(geneticCode.isStopCodon(inStates[j]))
                inStates[j] = (int)(generator.nextDouble() * stateCount);
        }

    }

    /**
     * Parsimony to init states. Equally to choose a state from the ambiguous set.
     */
    public void initINStatesParsimony() {

        throw new UnsupportedOperationException();

        //traverse 1


        //traverse 2


    }

    @Override
    public String getID() {
        return Integer.toString(nodeNr);
    }

    /**
     * indicate that the states matrix for node nodeIndex is to be changed *
     */
    public void setStatesForUpdate() {
        currentMatrixIndex = 1 - currentMatrixIndex; // 0 or 1
    }



    /**
     * cleans up and deallocates arrays.
     */
    @Override
    public void finalize() throws Throwable {
        states = null;
        storedStates = null;
        siteIsDirty = null;
        currentMatrixIndex = 0;
        storedMatrixIndex = 0;
    }


    //TODO store restore a site
    @Override
    protected void store() {
        storedMatrixIndex = currentMatrixIndex;
    }

    @Override
    public void restore() {
        // Rather than copying the stored stuff back, just swap the pointers...
        currentMatrixIndex = storedMatrixIndex;

        hasStartedEditing = false;
        if (siteIsDirty.length != getSiteCount()) {
            siteIsDirty = new boolean[getSiteCount()];
        }

    }

//    @Override
//    public void unstore() {
//        System.arraycopy(storedMatrixIndex, 0, currentMatrixIndex, 0, getNodeCount());
//    }

    public int getSiteCount() {
        return states[currentMatrixIndex].length;
    }

    /**
     * Get an internal node states. The state starts from 0.
     *
     * @return the codon states of the internal node sequence.
     */
    public int[] getNodeStates() {
        // internal node index starts from getTipsCount();
        return states[currentMatrixIndex];
    }

    /**
     * Set the states to an internal node.
     * The node index has to convert to the array index before setValue.
     *
     * @param states int[]
     */
    public void setNodeStates(final int[] states) {
        // internal node index starts from getTipsCount();
        setValue(states);
    }

    /**
     * modify setValue for 2d matrix to take int[].
     *
     * @param vals     int[] values of the parameter
     */
    protected void setValue(final int[] vals) {
        startEditing(null);

        System.arraycopy(vals, 0, this.states[currentMatrixIndex], 0, vals.length);
        Arrays.fill(siteIsDirty, true);
//        m_nLastDirty = nodeIndex;
    }

    /**
     * Get a codon state from the site at the node.
     * The state starts from 0.
     *
     * @param codonNr the site index.
     * @return
     */
    public int getASite(final int codonNr) {
        return states[currentMatrixIndex][codonNr];
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
    public StateNode copy() {
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
        storedStates[0] = source.storedStates[0].clone();
        storedStates[1] = source.storedStates[1].clone();

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
