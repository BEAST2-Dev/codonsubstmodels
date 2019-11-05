package beast.evolution.tree;

import beast.core.Input;
import beast.core.StateNode;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.CodonAlignment;
import org.w3c.dom.Node;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;

/**
 * The large 2-d matrix to store all node states.
 * Rows are internal nodes, cols are sites.
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
@Deprecated
public class InternalNodeStates extends StateNode {

    final public Input<Alignment> dataInput = new Input<>("data",
            "codon alignment to initialise the 2-d matrix of internal node states", Input.Validate.REQUIRED);

    /**
     * the actual values of this parameter
     */
    protected int[][] states;
    protected int[][] storedStates;
    //TODO use storeIndex, 0 to select current, 1 to select stored

    /**
     * upper & lower bound These are located before the inputs (instead of
     * after the inputs, as usual) so that valuesInput can determines the
     * class
     */
    protected int m_fUpper;
    protected int m_fLower;
    /**
     * isDirty flags for each internal node
     */
    protected boolean[] m_bIsDirty;


    // row (x) is internal nodes = leafNodeCount - 1
//    protected int internalNodeCount = -1;

    public InternalNodeStates(int stateCount, final int[][] states) {
        initParam(stateCount, states.length, states[0].length);
        for (int i = 0; i < states.length; i++) {
            System.arraycopy(states[i], 0, this.states[i], 0, states[i].length);
            System.arraycopy(states[i], 0, this.storedStates[i], 0, states[i].length);
        }
    }


    public InternalNodeStates(int internalNodeCount, int siteCount) {
        this(64, internalNodeCount, siteCount);
    }

    /**
     * The large 2-d matrix to store internal node states.
     *
     * @param stateCount        to define upper
     * @param internalNodeCount to define rows of 2d matrix
     * @param siteCount         to define cols of 2d matrix
     */
    public InternalNodeStates(int stateCount, int internalNodeCount, int siteCount) {
        initParam(stateCount, internalNodeCount, siteCount);
    }

    @Override
    public void initAndValidate() {
        // need data type, site count, taxa count
        Alignment alignment = dataInput.get();
        int stateCount = alignment.getDataType().getStateCount();

        if (alignment instanceof CodonAlignment) // 64
            assert stateCount == 64;

        // used to adjust Nr
        int internalNodeCount = alignment.getTaxonCount() - 1;
        assert internalNodeCount > 1;
        // siteCount = num of codons, overwrite in CodonAlignment /= 3
        int siteCount = alignment.getSiteCount();

        // stateCount -> upper, internalNodeCount * siteCount -> 2d matrix
        initParam(stateCount, internalNodeCount, siteCount);
    }

    // stateCount -> upper, internalNodeCount * siteCount -> 2d matrix
    private void initParam(int stateCount, int internalNodeCount, int siteCount) {
        // siteCount = num of codons, overwrite in CodonAlignment /= 3
        states = new int[internalNodeCount][siteCount];
        storedStates = new int[internalNodeCount][siteCount];

        // 0 - 63, ignore lowerValueInput upperValueInput
        m_fLower = 0;
        m_fUpper = stateCount - 1; // not deal unknown codon

        m_bIsDirty = new boolean[internalNodeCount];


//        minorDimension = siteCount;

        // internal nodes = n - 1, length = (n-1)*L
//        this.values = new Integer[internalNodeCount * minorDimension];
//        Arrays.fill(this.values, 0);
//        this.storedValues = values.clone();


//        Log.info.println("Create internal node states matrix : " + internalNodeCount + " * " + minorDimension);


//        T[] valuesString = valuesInput.get().toArray((T[]) Array.newInstance(getMax().getClass(), 0));
//
//        int dimension = Math.max(dimensionInput.get(), valuesString.length);
//        dimensionInput.setValue(dimension, this);
//        values = (T[]) Array.newInstance(getMax().getClass(), dimension);
//        storedValues = (T[]) Array.newInstance(getMax().getClass(), dimension);
//        for (int i = 0; i < values.length; i++) {
//            values[i] = valuesString[i % valuesString.length];
//        }
//
//        m_bIsDirty = new boolean[dimensionInput.get()];
//
//        minorDimension = minorDimensionInput.get();
//        if (minorDimension > 0 && dimensionInput.get() % minorDimension > 0) {
//            throw new IllegalArgumentException("Dimension must be divisible by stride");
//        }
//        this.storedValues = values.clone();

    }

    public int getInternalNodeCount() {
        return states.length; // = getTaxonCount() - 1
    }

    public int getNodeCount() {
        return getInternalNodeCount() * 2 + 1;
    }

    public int getSiteCount() {
        return states[0].length;
    }

    /**
     * Get an internal node states. The state starts from 0.
     *
     * @param nodeNr the node index = <code>nodeNr - internalNodeCount - 1</code><br>
     *               Leaf nodes are number 0 to <code>leafnodes-1</code>;
     *               Internal nodes are numbered  <code>leafnodes</code> up to <code>nodes-1</code>;
     *               The root node is always numbered <code>nodes-1</code>.
     * @return the codon states of the internal node sequence.
     */
    public int[] getNrStates(final int nodeNr) {
        // internal node index nodeNr starts from getTaxonCount();
        assert nodeNr > getInternalNodeCount() && nodeNr < getNodeCount();
        // convert nodeNr into 2d matrix row index
        int rowIndex = nodeNr - getInternalNodeCount() - 1;
        return states[rowIndex];
    }

    /**
     * Set the states to an internal node.
     * The node index has to convert to the array index before setValue.
     *
     * @param nodeNr the node index = <code>nodeNr - internalNodeCount - 1</code><br>
     *               Leaf nodes are number 0 to <code>leafnodes-1</code>;
     *               Internal nodes are numbered  <code>leafnodes</code> up to <code>nodes-1</code>;
     *               The root node is always numbered <code>nodes-1</code>.
     * @param states int[]
     */
    public void setNrStates(final int nodeNr, final int[] states) {
        // internal node index nodeNr starts from getTaxonCount();
        assert nodeNr > getInternalNodeCount() && nodeNr < getNodeCount();
        // convert nodeNr into 2d matrix row index
        int rowIndex = nodeNr - getInternalNodeCount() - 1;

        setValue(rowIndex, states);
    }

    /**
     * modify setValue for 2d matrix to take int[].
     *
     * @param rowIndex the start index of parameter to set to the flattened matrix,
     *                 the values of the parameter is given by an int[].
     * @param vals     int[] values of the parameter
     */
    protected void setValue(final int rowIndex, final int[] vals) {
        startEditing(null);

        System.arraycopy(vals, 0, this.states[rowIndex], 0, vals.length);
        m_bIsDirty[rowIndex] = true;
//        m_nLastDirty = rowIndex;
    }

    /**
     * get the sites at a site of all internal nodes.
     *
     * @param codonNr the codon site index.
     * @return
     */
    public int[] getSites(final int codonNr) {
        int[] col = new int[getInternalNodeCount()];
        for (int i = 0; i < getInternalNodeCount(); i++)
            col[i] = states[i][codonNr];
        return col;
    }

    /**
     * Get a codon state from the site at the internal node.
     * The state starts from 0.
     * <code>matrix[i,j] = values[i * minorDimension + j]</code>
     *
     * @param nodeNr  the node index <code>i = nodeNr - internalNodeCount - 1</code><br>
     *                Leaf nodes are number 0 to <code>leafnodes-1</code>;
     *                Internal nodes are numbered  <code>leafnodes</code> up to <code>nodes-1</code>;
     *                The root node is always numbered <code>nodes-1</code>.
     * @param codonNr the site index.
     * @return
     */
    public int getASite(final int nodeNr, final int codonNr) {
        // internal node index nodeNr starts from getTaxonCount();
        assert nodeNr > getInternalNodeCount() && nodeNr < getNodeCount();
        // convert nodeNr into 2d matrix row index
        int rowIndex = nodeNr - getInternalNodeCount() - 1;
        return states[rowIndex][codonNr];
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
        Alignment alignment = dataInput.get();
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
//TODO store restore a site
    @Override
    protected void store() {
        if (storedStates.length != states.length || storedStates[0].length != states[0].length) {
            storedStates = states.clone();
            for (int i = 0; i < states.length; i++)
                storedStates[i] = states[i].clone();
        } else {
            for (int i = 0; i < states.length; i++)
                System.arraycopy(states[i], 0, storedStates[i], 0, states[i].length);
        }
    }

    @Override
    public void restore() {
        // final T[] tmp = storedValues;
        final int[][] tmp = storedStates.clone();
        for (int i = 0; i < storedStates.length; i++)
            System.arraycopy(storedStates[i], 0, tmp[i], 0, storedStates[i].length);
        // storedValues = values;
        for (int i = 0; i < storedStates.length; i++)
            System.arraycopy(states[i], 0, storedStates[i], 0, states[i].length);
        // values = tmp;
        for (int i = 0; i < states.length; i++)
            System.arraycopy(tmp[i], 0, states[i], 0, tmp[i].length);

        hasStartedEditing = false;
        if (m_bIsDirty.length != states.length) {
            m_bIsDirty = new boolean[states.length];
        }

    }

    @Override
    public void init(PrintStream out) {
        final int innoCount = getInternalNodeCount();
        if (innoCount == 1) {
            out.print(getID() + "\t");
        } else {
            for (int i = 0; i < innoCount; i++) {
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
        buf.append(getID()).append("[").append(states.length);
        buf.append(" ").append(getInternalNodeCount());
        buf.append(" ").append(getSiteCount());
        buf.append("] ");
        buf.append("(").append(m_fLower).append(",").append(m_fUpper).append("): ");
//        for (final int value : states) {
//            buf.append(value).append(" ");
//        }
        return buf.toString();
    }

    @Override
    public void setEverythingDirty(boolean isDirty) {
        setSomethingIsDirty(isDirty);
        Arrays.fill(m_bIsDirty, isDirty);
    }

    @Override
    public StateNode copy() {
        try {
            @SuppressWarnings("unchecked")
            final InternalNodeStates copy = (InternalNodeStates) this.clone();
            copy.states = states.clone();
            for (int i = 0; i < states.length; i++)
                copy.states[i] = states[i].clone();
            //storedStates = copy.storedStates.clone();
            copy.m_bIsDirty = new boolean[states.length];
            return copy;
        } catch (Exception e) {
            e.printStackTrace();
        }
        return null;

    }

    @Override
    public void assignTo(StateNode other) {
        @SuppressWarnings("unchecked")
        final InternalNodeStates copy = (InternalNodeStates) other;
        copy.setID(getID());
        copy.index = index;
        copy.states = states.clone();
        for (int i = 0; i < states.length; i++)
            copy.states[i] = states[i].clone();
        //storedStates = copy.storedStates.clone();
        copy.m_fLower = m_fLower;
        copy.m_fUpper = m_fUpper;
        copy.m_bIsDirty = new boolean[states.length];
    }

    @Override
    public void assignFrom(StateNode other) {
        @SuppressWarnings("unchecked")
        final InternalNodeStates source = (InternalNodeStates) other;
        setID(source.getID());
        states = source.states.clone();
        storedStates = source.storedStates.clone();
        for (int i = 0; i < states.length; i++) {
            states[i] = source.states[i].clone();
            storedStates[i] = source.storedStates[i].clone();
        }
        m_fLower = source.m_fLower;
        m_fUpper = source.m_fUpper;
        m_bIsDirty = new boolean[source.states.length];
    }

    @Override
    public void assignFromFragile(StateNode other) {
        @SuppressWarnings("unchecked")
        final InternalNodeStates source = (InternalNodeStates) other;
        for (int i = 0; i < states.length; i++)
            System.arraycopy(source.states[i], 0, states[i], 0, source.states[i].length);
        Arrays.fill(m_bIsDirty, false);
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
