package beast.evolution.tree;

import beast.core.Input;
import beast.core.StateNode;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.CodonAlignment;
import beast.evolution.datatype.Codon;
import beast.evolution.datatype.GeneticCode;

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
    final int nodeCount;

    // isDirty flags for each node
    protected boolean[] nodeIsDirty;


    public NodesStates(CodonAlignment codonAlignment) {
        this.codonDataType = codonAlignment.getDataType();
        final int stateCount = getStateCount();
        assert stateCount == 64;
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

    public int getNodeCount() {
        return nodeCount;
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


    //******* use those in NodesStatesAndTree not these below *******

    @Override
    protected void store() {
        throw new UnsupportedOperationException("Unsupported");
    }

    @Override
    public void restore() {
        throw new UnsupportedOperationException("Unsupported");
    }

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
