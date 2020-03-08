package beast.evolution.tree;

import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.CodonAlignment;
import beast.evolution.datatype.DataType;

import java.io.PrintStream;
import java.util.List;

/**
 * The large array to store one node states (can be a tip).
 * A non-ambiguous state starts from 0 ends to <code>stateCount - 1</code>,
 * such as the universal genetic code has states 0-60. <br>
 *
 * <code>nodeNr</code> used for all indexing system
 * is the node index number starting from 0 to root index (nodeCount - 1).
 * For internal nodes, they are ranged from the number of leaf nodes
 * to the nodeCount - 1.<br>
 *
 * Not use {@link beast.core.parameter.IntegerParameter}, because of performance.
 * https://softwareengineering.stackexchange.com/questions/203970/when-to-use-primitive-vs-class-in-java/203974.
 *
 * No {@link beast.core.Input} here.
 * Init and validate in {@link NodeStatesArray#initAndValidate()}.
 */
public class NodeStates implements Cloneable {// cannot extend StateNode

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
     * states[sites], the int states for this node.
     */
    protected int[] states;
    protected int[] storedStates;

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
//    protected boolean[] siteIsDirty; // TODO setASiteDirty


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

//    @Override
//    public void initAndValidate() {
//        //init validate in NodesStates.
//    }

    // stateCount -> upper, siteCount -> array
    private void initParam(int stateCount, int siteCount) {
        // siteCount = num of codons, overwrite in CodonAlignment /= 3
        states = new int[siteCount];
        storedStates = new int[siteCount];

        // 0 - 59/60, ignore lowerValueInput upperValueInput
        lower = 0;
        upper = stateCount - 1; // not deal unknown codon

//        siteIsDirty = new boolean[siteCount];

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
        assert statesList.size() == states.length;
        // Java 8
        states = statesList.stream().mapToInt(i->i).toArray();
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
//    @Override
    public String getID() {
        return Integer.toString(getNodeNr());
    }

    /**
     * @return  null if not tip, otherwise return taxon name
     */
    public String getTipID() {
        return tipID;
    }

    public int getSiteCount() {
        return states.length;
    }

    /**
     * Get an internal node states. The state starts from 0.
     *
     * @return the codon states of the internal node sequence.
     */
    public int[] getStates() {
        // internal node index starts from getTipsCount();
        return states;
    }

    /**
     * Set the states to an internal node.
     * The node index has to convert to the array index before setValue.
     *
     * @param states int[]
     */
    public void setStates(final int[] states) {
        // internal node index starts from getTipsCount();
//        startEditing(null);
        System.arraycopy(states, 0, this.states, 0, states.length);
//        Arrays.fill(siteIsDirty, true);
    }

    /**
     * Set a section of states to an internal node.
     * @param startInclusive  start position (included) in the whole sequence
     * @param endExclusive    end position (excluded) in the whole sequence
     * @param states
     */
    public void setStates(int startInclusive, int endExclusive, final int[] states) {
        // internal node index starts from getTipsCount();
//        startEditing(null);
        System.arraycopy(states, 0, this.states, startInclusive, states.length);
//        Arrays.fill(siteIsDirty, startInclusive, endExclusive, true);
    }

    /**
     * Get a codon state from the site at the node.
     * The state starts from 0.
     *
     * @param codonNr the site index.
     * @return
     */
    public int getState(final int codonNr) {
        return states[codonNr];
    }


    /**
     * Set a codon state to the site at the internal node.
     * The state starts from 0.
     * <code>matrix[i,j] = values[i * minorDimension + j]</code>
     *
     * @param codonNr the site index.
     * @param state   new state to set.
     */
    public void setState(final int codonNr, final int state) {
//        startEditing(null);

        states[codonNr] = state;
//        siteIsDirty[codonNr] = true;
    }

    /**
     * @param ns  another NodeStates
     * @return    true if they have the same states, no matter siteIsDirty[].
     */
    public boolean hasSameStates(NodeStates ns) {
        int[] states2 = ns.getStates();
        if (states2.length != states.length)
            return false;
        for (int i = 0; i < states.length; i++) {
            if (states2[i] != states[i]) {
                Log.err.println("Internal node " + nodeNr + " site " + i + " has different states !");
                return false;
            }
        }
        return true;
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

    /**
     * cleans up and deallocates arrays.
     */
    @Override
    public void finalize() throws Throwable {
        states = null;
//        storedStates = null;
//        siteIsDirty = null;
    }


    public void store() {
//        if (tipID != null)
//            throw new UnsupportedOperationException("Can not sample tip states yet !");

        System.arraycopy(states, 0, storedStates, 0, states.length);
    }

    public void restore() {
//        if (tipID != null)
//            throw new UnsupportedOperationException("Can not sample tip states yet !");

        final int[] tmp = storedStates;
        storedStates = states;
        states = tmp;
//        hasStartedEditing = false;
//        Arrays.fill(siteIsDirty, false);
    }

    /**
     * Used in {@link NodeStatesArray#copy()}.
     * @return a deep copy of NodeStates
     */
    public NodeStates deepCopy() {
        NodeStates copy = null;
        try {
//            copy = new NodeStates(nodeNr, states, upper+1);
            copy = (NodeStates) this.clone();
            copy.setStates(this.states);
        } catch (Exception e) {
            e.printStackTrace();
        }
        return copy;
    }

    /**
     * Print <code>nodeNr+1</code> and states.
     * Note: BEAST tree log uses nodeNr+1 in Newick,
     * see {@link Node#toSortedNewick(int[])}.
     * @param out {@link PrintStream}
//     * @param delimiter  delimiter between states, use "" to make no delimiter case.
     */
    public void log(PrintStream out) {//, String delimiter) {
        out.print((nodeNr+1) + "\t"); // Note: BEAST tree log makes nodeNr+1 in Newick
        for (int sta : states) {
            // 00,01,02,...,59
            if (sta < 10)
                out.print("0" + sta);
            else
                out.print(sta);
            // delimiter between states
//            out.print(delimiter);
        }
    }

//    @Override
//    public boolean equals(Object obj) {
//        NodeStates ns = (NodeStates) obj;
//        // nodeNr
//        // tipID
//        // states
//        // siteIsDirty
//        // lower upper
//        return super.equals(obj);
//    }

//TODO

//    @Override
//    public void assignTo(StateNode other) {
//        @SuppressWarnings("unchecked")
//        final NodeStates copy = (NodeStates) other;
//        copy.setID(getID());
//        copy.index = index;
//
//        copy.states = states.clone();
//        //storedStates = copy.storedStates.clone();
//        copy.lower = lower;
//        copy.upper = upper;
//        copy.siteIsDirty = new boolean[getSiteCount()];
//    }
//
//    @Override
//    public void assignFrom(StateNode other) {
//        @SuppressWarnings("unchecked")
//        final NodeStates source = (NodeStates) other;
//        setID(source.getID());
//        states = source.states.clone();
////        storedStates = source.storedStates.clone();
//        lower = source.lower;
//        upper = source.upper;
//        siteIsDirty = new boolean[source.getSiteCount()];
//    }
//
//    @Override
//    public void assignFromFragile(StateNode other) {
//        @SuppressWarnings("unchecked")
//        final NodeStates source = (NodeStates) other;
//        System.arraycopy(source.states, 0, states, 0, source.states.length);
//        Arrays.fill(siteIsDirty, false);
//    }
//
//    @Override
//    public void fromXML(org.w3c.dom.Node xmlNode) {
//        throw new UnsupportedOperationException("in dev");
//    }
//
//    @Override
//    public int scale(double scale) {
//        throw new UnsupportedOperationException("in dev");
//    }

//    @Override
//    public int getDimension() { // use getInternalNodeCount() & getSiteCount()
//        throw new UnsupportedOperationException("Unsupported");
//    }
//
//    @Override
//    public double getArrayValue(int dim) { // use getASite(int, int)
//        throw new UnsupportedOperationException("Unsupported");
//    }

//    @Override
//    public void setEverythingDirty(boolean isDirty) {
//        setSomethingIsDirty(isDirty);
//        Arrays.fill(siteIsDirty, isDirty);
//    }
//

}
