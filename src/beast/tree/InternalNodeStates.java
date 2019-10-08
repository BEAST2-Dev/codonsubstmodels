package beast.tree;

import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.util.Log;
import beast.evolution.alignment.CodonAlignment;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;

/**
 * The large 2-d matrix to store internal node states.
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
 *
 */
public class InternalNodeStates extends IntegerParameter {

    final public Input<CodonAlignment> dataInput = new Input<>("data",
            "codon alignment to initialise the 2-d matrix of internal node states", Input.Validate.REQUIRED);

    // row (x) is internal nodes = leafNodeCount - 1
    protected int internalNodeCount = -1;

    protected CodonAlignment codonAlignment;

    public InternalNodeStates() {
        dimensionInput.setRule(Input.Validate.FORBIDDEN);
        minorDimensionInput.setRule(Input.Validate.FORBIDDEN);
        lowerValueInput.setRule(Input.Validate.FORBIDDEN);
        upperValueInput.setRule(Input.Validate.FORBIDDEN);
    }


    public InternalNodeStates(int internalNodeCount, int siteCount) {
        this(64, internalNodeCount, siteCount);
    }

    /**
     * The large 2-d matrix to store internal node states.
     * @param stateCount           to define upper
     * @param internalNodeCount    to define rows of 2d matrix
     * @param siteCount            to define cols of 2d matrix
     */
    public InternalNodeStates(int stateCount, int internalNodeCount, int siteCount) {
        this(); // reserve  dimension, minorDimension
        initParam(stateCount, internalNodeCount, siteCount);
    }

    @Override
    public void initAndValidate() {
        // need data type, site count, taxa count
        codonAlignment = dataInput.get();
        int stateCount = codonAlignment.getDataType().getStateCount(); // 64
        assert stateCount == 64;

        // used to adjust Nr
        int internalNodeCount = codonAlignment.getTaxonCount() - 1;
        assert internalNodeCount > 1;
        // L = num of codons, overwrite in CodonAlignment /= 3
        int siteCount = codonAlignment.getSiteCount();

        // stateCount -> upper, internalNodeCount * siteCount -> 2d matrix
        initParam(stateCount, internalNodeCount, siteCount);
    }

    // stateCount -> upper, internalNodeCount * siteCount -> 2d matrix
    private void initParam(int stateCount, int internalNodeCount, int siteCount) {
        this.internalNodeCount = internalNodeCount;

        // 0 - 63, ignore lowerValueInput upperValueInput
        m_fLower = 0;
        m_fUpper = stateCount-1; // not deal unknown codon

        // L = num of codons, overwrite in CodonAlignment /= 3
        minorDimension = siteCount;

        // internal nodes = n - 1, length = (n-1)*L
        this.values = new Integer[internalNodeCount * minorDimension];
        Arrays.fill(this.values, 0);
        this.storedValues = values.clone();

        m_bIsDirty = new boolean[values.length];

        Log.info.println("Create internal node states matrix : " + internalNodeCount + " * " + minorDimension);


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

    /**
     * Get an internal node states. The state starts from 0.
     * @param nodeNr the node index = <code>nodeNr - internalNodeCount - 1</code><br>
     *               Leaf nodes are number 0 to <code>leafnodes-1</code>;
     *               Internal nodes are numbered  <code>leafnodes</code> up to <code>nodes-1</code>;
     *               The root node is always numbered <code>nodes-1</code>.
     * @return the codon states of the internal node sequence.
     */
    public int[] getNrStates(final int nodeNr) {
        assert nodeNr > internalNodeCount; // internalNodeCount = getTaxonCount() - 1;
        int rowIndex = nodeNr - internalNodeCount - 1; // convert nodeNr into 2d matrix row index
        Integer[] row = new Integer[minorDimension];
        // matrix[i,j] = values[i * minorDimension + j]
        System.arraycopy(values, rowIndex * minorDimension, row, 0, row.length);

        return Arrays.stream(row).mapToInt(Integer::intValue).toArray();
    }

    /**
     * Set the states to an internal node.
     * The node index has to convert to the array index before setValue.
     * @param nodeNr the node index = <code>nodeNr - internalNodeCount - 1</code><br>
     *               Leaf nodes are number 0 to <code>leafnodes-1</code>;
     *               Internal nodes are numbered  <code>leafnodes</code> up to <code>nodes-1</code>;
     *               The root node is always numbered <code>nodes-1</code>.
     * @param states  int[]
     */
    public void setNrStates(final int nodeNr, final int[] states) {
        assert nodeNr > internalNodeCount; // internalNodeCount = getTaxonCount() - 1;
        int rowIndex = nodeNr - internalNodeCount - 1; // convert nodeNr into 2d matrix row index
        assert states.length == minorDimension; // minorDimension == cols

        //TODO performance check
        Integer[] vals = Arrays.stream( states ).boxed().toArray( Integer[]::new );

        setValue(rowIndex * minorDimension, vals);
    }

    /**
     * modify setValue for 2d matrix to take int[].
     * @param param  the start index of parameter to set to the flattened matrix,
     *               the values of the parameter is given by an int[].
     * @param vals  int[] values of the parameter
     */
    protected void setValue(int param, Integer[] vals) {
        assert vals.length + param <= this.values.length;

        startEditing(null);

        System.arraycopy(vals, 0, this.values, param, vals.length);
        for (int i = param; i < (param + vals.length); i++)
            m_bIsDirty[i] = true;
//        m_nLastDirty = param;
    }

    /**
     * get the sites at a site of all internal nodes.
     * @param codonNr the codon site index.
     * @return
     */
    public int[] getSites(final int codonNr) {
        assert codonNr < minorDimension;
        int[] col = new int[internalNodeCount];
        for (int i = 0; i < internalNodeCount; i++)
            col[i] = values[i * minorDimension + codonNr];
        return col;
    }

    /**
     * Get a codon state from the site at the internal node.
     * The state starts from 0.
     * <code>matrix[i,j] = values[i * minorDimension + j]</code>
     *
     * @param nodeNr the node index <code>i = nodeNr - internalNodeCount - 1</code><br>
     *               Leaf nodes are number 0 to <code>leafnodes-1</code>;
     *               Internal nodes are numbered  <code>leafnodes</code> up to <code>nodes-1</code>;
     *               The root node is always numbered <code>nodes-1</code>.
     * @param codonNr the site index.
     * @return
     */
    public int getASite(final int nodeNr, final int codonNr) {
        assert nodeNr > internalNodeCount; // internalNodeCount = getTaxonCount() - 1;
        int rowIndex = nodeNr - internalNodeCount - 1; // convert nodeNr into 2d matrix row index
        assert codonNr < minorDimension;
        return values[rowIndex * minorDimension + codonNr];
    }


    /**
     * Set a codon state to the site at the internal node.
     * The state starts from 0.
     * <code>matrix[i,j] = values[i * minorDimension + j]</code>
     *
     * @param nodeNr the node index <code>i = nodeNr - internalNodeCount - 1</code><br>
     *               Leaf nodes are number 0 to <code>leafnodes-1</code>;
     *               Internal nodes are numbered  <code>leafnodes</code> up to <code>nodes-1</code>;
     *               The root node is always numbered <code>nodes-1</code>.
     * @param codonNr the site index.
     */
    public void setASite(final int nodeNr, final int codonNr, final int value) {
        assert nodeNr > internalNodeCount; // internalNodeCount = getTaxonCount() - 1;
        int rowIndex = nodeNr - internalNodeCount - 1; // convert nodeNr into 2d matrix row index
        assert codonNr < minorDimension;

        setValue(rowIndex * minorDimension + codonNr, value);
    }


    @Override
    public void log(long sampleNr, PrintStream out) {
        super.log(sampleNr, out);
    }


    /**
     *
     * @param triplets coded nucleotides in string
     * @return the list of codon states. The size should be 1/3 of string length
     */
    public List<Integer> stringToEncoding(String triplets) {
        // remove spaces
        triplets = triplets.replaceAll("\\s", "");
        List<Integer> sequence = codonAlignment.getDataType().stringToEncoding(triplets);
        if (sequence.size() * 3 != triplets.length())
            throw new IllegalArgumentException("The string of triplets has invalid number ! " +
                    triplets.length() + " != " + sequence.size() + " * 3");
        return sequence;
    }

    /**
     *
     * @param triplets coded nucleotides in string
     * @return the int[] of codon states. The size should be 1/3 of string length
     */
    public int[] stringToStates(String triplets) {
        List<Integer> sequence = stringToEncoding(triplets);

        return sequence.stream().mapToInt(i->i).toArray();
    }

}
