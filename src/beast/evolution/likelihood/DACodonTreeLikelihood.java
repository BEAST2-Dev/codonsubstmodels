package beast.evolution.likelihood;

import beast.core.Input;
import beast.core.State;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.CodonAlignment;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.InternalNodeStates;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;

/**
 * Data augmentation to fast codon tree likelihood calculation.
 *
 * No pattern, use site count.
 *
 * TODO 1: make it working in MCMC
 * TODO 2: try only log the cell of trans prob matrix once
 */
public class DACodonTreeLikelihood extends GenericTreeLikelihood {

    final public Input<Boolean> m_useAmbiguities = new Input<>("useAmbiguities",
            "flag to indicate that sites containing ambiguous states should be handled instead of ignored (the default)",
            false);
    final public Input<Boolean> m_useTipLikelihoods = new Input<>("useTipLikelihoods",
            "flag to indicate that partial likelihoods are provided at the tips", false);
//    final public Input<String> implementationInput = new Input<>("implementation",
//            "name of class that implements this treelikelihood potentially more efficiently. " +
//            "This class will be tried first, with the TreeLikelihood as fallback implementation. " +
//            "When multi-threading, multiple objects can be created.",
//            "beast.evolution.likelihood.BeagleTreeLikelihood");

//    public static enum Scaling {none, always, _default};
    final public Input<TreeLikelihood.Scaling> scaling = new Input<>("scaling",
        "type of scaling to use, one of " + Arrays.toString(TreeLikelihood.Scaling.values()) +
                ". If not specified, the -beagle_scaling flag is used.",
        TreeLikelihood.Scaling._default, TreeLikelihood.Scaling.values());


    final public Input<InternalNodeStates> internalNodeStatesInput = new Input<>("internalNodeStates",
            "The large 2-d matrix to store internal node sequences.", Input.Validate.REQUIRED);


    /**
     * states in tips is stored by CodonAlignment List<List<Integer>> counts
     */
    CodonAlignment codonAlignment;
    /**
     * internal node sequences 2d matrix = [(internal) nodeNr , codonNr]
     */
    InternalNodeStates internalNodeStates;

    /**
     * calculation engine *
     */
    protected DAStatesLikelihoodCore daLdCore;
//    protected BeagleTreeLikelihood beagle;

    /**
     * BEASTObject associated with inputs. Since none of the inputs are StateNodes, it
     * is safe to link to them only once, during initAndValidate.
     */
    protected SubstitutionModel substitutionModel;
    protected SiteModel.Base siteModel;
    protected BranchRateModel.Base branchRateModel;

    /**
     * flag to indicate the
     * // when CLEAN=0, nothing needs to be recalculated for the node
     * // when DIRTY=1 indicates a node partial needs to be recalculated
     * // when FILTHY=2 indicates the indices for the node need to be recalculated
     * // (often not necessary while node partial recalculation is required)
     */
    protected int hasDirt;

    /**
     * Lengths of the branches in the tree associated with each of the nodes
     * in the tree through their node  numbers. By comparing whether the
     * current branch length differs from stored branch lengths, it is tested
     * whether a node is dirty and needs to be recomputed (there may be other
     * reasons as well...).
     * These lengths take branch rate models in account.
     */
    protected double[] branchLengths;
    protected double[] storedBranchLengths;

    /**
     * memory allocation for log likelihoods for each of node *
     */
    protected double[] nodeLogLikelihoods;
    /**
     * memory allocation for the root branch likelihood *
     */
//    protected double[] rootBranchLd;
    /**
     * memory allocation for probability tables obtained from the SiteModel *
     */
    protected double[] probabilities;

    protected int matrixSize;

    /**
     * dealing with proportion of site being invariant *
     */
    double proportionInvariant = 0;
    List<Integer> constantPattern = null;


    @Override
    public void initAndValidate() {

        internalNodeStates = internalNodeStatesInput.get();

        codonAlignment = CodonAlignment.toCodonAlignment(dataInput.get());
        // sanity check: alignment should have same #taxa as tree
        if (codonAlignment.getTaxonCount() != treeInput.get().getLeafNodeCount()) {
            throw new IllegalArgumentException("The number of nodes in the tree does not match the number of sequences");
        }
        Log.info.println("  " + codonAlignment.toString(true));
        // print startup messages via Log.print*

//        beagle = null;
//        beagle = new BeagleTreeLikelihood();
//        try {
//            beagle.initByName(
//                    "data", codonAlignment, "tree", treeInput.get(), "siteModel", siteModelInput.get(),
//                    "branchRateModel", branchRateModelInput.get(), "useAmbiguities", m_useAmbiguities.get(),
//                    "useTipLikelihoods", m_useTipLikelihoods.get(),"scaling", scaling.get().toString());
//            if (beagle.beagle != null) {
//                //a Beagle instance was found, so we use it
//                return;
//            }
//        } catch (Exception e) {
//            // ignore
//        }
//        // No Beagle instance was found, so we use the good old java likelihood core
//        beagle = null;

        int nodeCount = treeInput.get().getNodeCount();
        if (!(siteModelInput.get() instanceof SiteModel.Base)) {
            throw new IllegalArgumentException("siteModel input should be of type SiteModel.Base");
        }
        siteModel = (SiteModel.Base) siteModelInput.get();
        siteModel.setDataType(codonAlignment.getDataType());
        substitutionModel = siteModel.substModelInput.get();

        if (branchRateModelInput.get() != null) {
            branchRateModel = branchRateModelInput.get();
        } else {
            branchRateModel = new StrictClockModel();
        }
        branchLengths = new double[nodeCount];
        storedBranchLengths = new double[nodeCount];

        initCore(); // init DALikelihoodCore




        // TODO check
        proportionInvariant = siteModel.getProportionInvariant();
        siteModel.setPropInvariantIsCategory(false);
        // TODO no patterns, need to check
//        int patterns = codonAlignment.getPatternCount();
        if (proportionInvariant > 0) {
            throw new UnsupportedOperationException("proportionInvariant in dev");
//            calcConstantPatternIndices(patterns, stateCount);
        }

//        int siteCount = codonAlignment.getSiteCount();
//        rootBranchLd = new double[siteCount]; // improved from patterns * stateCount to siteCount

        if (codonAlignment.isAscertained) {
            throw new UnsupportedOperationException("Ascertainment correction is not available !");
        }
    }


    // no pattern, use codonAlignment.getSiteCount()
    protected void initCore() {
        final int stateCount = codonAlignment.getMaxStateCount();
        assert stateCount == 64;

        TreeInterface tree = treeInput.get();
        // init with states
        final int[][] tipStates = getTipStates(tree);

        daLdCore = new DAStatesLikelihoodCore(stateCount, tipStates, internalNodeStates, siteModel.getCategoryCount());
//        daLdCore.initialize();

        // transition probability matrix, P
//        matrixSize = stateCount * stateCount;
        matrixSize = daLdCore.getMatrixSize();
        probabilities = new double[matrixSize];
        Arrays.fill(probabilities, 1.0);

        if (m_useAmbiguities.get() || m_useTipLikelihoods.get()) {
            throw new UnsupportedOperationException("in development");
        }
        hasDirt = Tree.IS_FILTHY;

        final int nodeCount = tree.getNodeCount();
        // log likelihoods for each of node
        nodeLogLikelihoods = new double[nodeCount];

    }

    /**
     * This method samples the sequences based on the tree and site model.
     */
    @Override
    public void sample(State state, Random random) {
        throw new UnsupportedOperationException("Can't sample a fixed alignment!");
    }

    /**
     * get states from tips, row is nodeIndex, col is site
     * @param tree  TreeInterface
     * @return      int[][]
     * @see {@link DAStatesLikelihoodCore#getNodeStates(int)}
     * @see {@link InternalNodeStates#getNrStates(int)}
     */
    protected int[][] getTipStates(TreeInterface tree) {

        final int[][] states = new int[tree.getLeafNodeCount()][];

        for (Node node : tree.getExternalNodes()) {
            // set leaf node states
            int taxonIndex = getTaxonIndex(node.getID(), codonAlignment);
            // no patterns
            List<Integer> statesList = codonAlignment.getCounts().get(taxonIndex);

            assert statesList.size() == codonAlignment.getSiteCount();
            // Java 8
            states[taxonIndex] = statesList.stream().mapToInt(i->i).toArray();
        }
        return states;
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

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    double m_fScale = 1.01;
    int m_nScale = 0;
    int X = 100;

    @Override
    public double calculateLogP() {
        //TODO beagle
//        if (beagle != null) {
//            logP = beagle.calculateLogP();
//            return logP;
//        }
        logP = 0;

        final TreeInterface tree = treeInput.get();
        final double[] frequencies = substitutionModel.getFrequencies();

        try {
            if (updateNode(tree) != Tree.IS_CLEAN)
                logP = daLdCore.calculateLogLikelihoods(frequencies);
//            System.out.println("logP = " + logP);
        }
        catch (ArithmeticException e) {
            return Double.NEGATIVE_INFINITY;
        }

//        m_nScale++;
//        if (logP > 0 || (daLdCore.getUseScaling() && m_nScale > X)) {
////            System.err.println("Switch off scaling");
////            m_likelihoodCore.setUseScaling(1.0);
////            m_likelihoodCore.unstore();
////            m_nHasDirt = Tree.IS_FILTHY;
////            X *= 2;
////            updateNode(tree.getRoot());
////            calcLogP();
////            return logP;
//        } else if (logP == Double.NEGATIVE_INFINITY && m_fScale < 10 && !scaling.get().equals(TreeLikelihood.Scaling.none)) { // && !m_likelihoodCore.getUseScaling()) {
//            m_nScale = 0;
//            m_fScale *= 1.01;
//            Log.warning.println("Turning on scaling to prevent numeric instability " + m_fScale);
////TODO            daLdCore.setUseScaling(m_fScale);
//            daLdCore.unstore();
//            hasDirt = Tree.IS_FILTHY;
//            updateNode(tree.getRoot());
//
//            calcLogP();
//
//            return logP;
//        }
        return logP;
    }


    protected int updateNode(final TreeInterface tree) {

        //TODO strange code in SiteModel, node is not used
        final double[] proportions = siteModel.getCategoryProportions(tree.getRoot());

        int update = hasDirt;

        // exclude root node
        int rootIndex = tree.getNodeCount() - 1;
        if (rootIndex != tree.getRoot().getNr())
            throw new RuntimeException("Invalid root index " + tree.getRoot().getNr() + " != " + rootIndex);

//        for (Node node : tree.getNodesAsArray()) {
        for (int n = 0; n < rootIndex; n++) {
            final Node node = tree.getNode(n);
            final Node parent = node.getParent();
//            if (!node.isRoot()) {
            int nodeUpdate = node.isDirty();

            final int nodeIndex = node.getNr();
            final int parentNum = parent.getNr();

            final double branchRate = branchRateModel.getRateForBranch(node);
            final double branchTime = node.getLength() * branchRate;

//TODO deal with 0 branch length, such as SA
            if (branchTime < 1e-6)
                throw new UnsupportedOperationException(
//                System.err.println(
                        "Time from parent " + parentNum + " to node " + nodeIndex + " is 0 !  " +
                        "branch length = " + node.getLength() + ", branchRate = " + branchRate);

            //TODO how to distinguish branch len change and internal node seq change, when topology is same
            // ====== 1. update the transition probability matrix(ices) if the branch len changes ======
            if (nodeUpdate != Tree.IS_CLEAN || branchTime != branchLengths[nodeIndex]) {
                branchLengths[nodeIndex] = branchTime;

//                daLdCore.setNodeMatrixForUpdate(nodeIndex); // TODO implement updates
                for (int i = 0; i < siteModel.getCategoryCount(); i++) {
                    final double jointBranchRate = siteModel.getRateForCategory(i, node) * branchRate;
                    substitutionModel.getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), jointBranchRate, probabilities);
                    //System.out.println(node.getNr() + " " + Arrays.toString(probabilities));

                    daLdCore.setNodeMatrix(nodeIndex, i, probabilities); //TODO how to rm arraycopy
                }
                nodeUpdate |= Tree.IS_DIRTY; // TODO review
            }

            // ====== 2. recalculate likelihood if either child node wasn't clean ======
            if (nodeUpdate != Tree.IS_CLEAN) {
                // brLD is linked to the child node index down
//                daLdCore.setNodeBrLdForUpdate(parentNum); // TODO implement updates

                // populate branchLd[][excl. root], nodeIndex is child
                daLdCore.calculateNodeBrLdOverCategories(nodeIndex, parentNum, proportions);
            }

            update |= nodeUpdate;
//            }
        } // end i loop

        return update;
    }
//+++++++++++ time ++++++++++++
//    public int getNodeUpdate(Node node, Node parent, int nodeUpdate, int nodeIndex, double branchRate, double branchTime) {
//        if (nodeUpdate != Tree.IS_CLEAN || branchTime != branchLengths[nodeIndex]) {
//            branchLengths[nodeIndex] = branchTime;
//
//            long time1 = 0;
//            long time2 = 0;
////                daLdCore.setNodeMatrixForUpdate(nodeIndex); // TODO implement updates
//            for (int i = 0; i < siteModel.getCategoryCount(); i++) {
//                time1 += getPt(node, parent, branchRate, i);
//                //System.out.println(node.getNr() + " " + Arrays.toString(probabilities));
//
//                time2 += arraycopy(nodeIndex, i);
//            }
//            nodeUpdate |= Tree.IS_DIRTY; // TODO review
//
//            System.out.println("\ngetTransitionProbabilities time is " + time1 + " nanoseconds");
//            System.out.println("arraycopy time is " + time2 + " nanoseconds\n");
//        }
//        return nodeUpdate;
//    }
//
//    public long getPt(Node node, Node parent, double branchRate, int i) {
//        long start = System.nanoTime();
//        final double jointBranchRate = siteModel.getRateForCategory(i, node) * branchRate;
//        substitutionModel.getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), jointBranchRate, probabilities);
//        return System.nanoTime()-start;
//    }
//
//    public long arraycopy(int nodeIndex, int i) {
//        long start = System.nanoTime();
//        daLdCore.setNodeMatrix(nodeIndex, i, probabilities); //TODO how to rm arraycopy
//        return System.nanoTime()-start;
//    }


    // for testing
    public void getBranchPartials(int nodeIndex, double[] branchLdOut) {
        daLdCore.getNodeBranchLd(nodeIndex, branchLdOut);
    }

    /* return copy of node log likelihoods */
    public double [] getNodeLogLikelihoods() {
//        if (beagle != null) {
//            return beagle.getNodeLogLikelihoods();
//        }
        return nodeLogLikelihoods.clone();
    } // getNodeLogLikelihoods

    /** CalculationNode methods **/

    /**
     * check state for changed variables and update temp results if necessary *
     */
    @Override
    protected boolean requiresRecalculation() {
//        if (beagle != null) {
//            return beagle.requiresRecalculation();
//        }
        hasDirt = Tree.IS_CLEAN;

        if (dataInput.get().isDirtyCalculation()) {
            hasDirt = Tree.IS_FILTHY;
            return true;
        }
        if (siteModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_DIRTY;
            return true;
        }
        if (branchRateModel != null && branchRateModel.isDirtyCalculation()) {
            //m_nHasDirt = Tree.IS_DIRTY;
            return true;
        }
        return treeInput.get().somethingIsDirty();
    }

    @Override
    public void store() {
//        if (beagle != null) {
//            beagle.store();
//            super.store();
//            return;
//        }
        if (daLdCore != null) {
            daLdCore.store();
        }
        super.store();
        System.arraycopy(branchLengths, 0, storedBranchLengths, 0, branchLengths.length);
    }

    @Override
    public void restore() {
//        if (beagle != null) {
//            beagle.restore();
//            super.restore();
//            return;
//        }
        if (daLdCore != null) {
            daLdCore.restore();
        }
        super.restore();
        double[] tmp = branchLengths;
        branchLengths = storedBranchLengths;
        storedBranchLengths = tmp;
    }

    /**
     * @return a list of unique ids for the state nodes that form the argument
     */
    @Override
    public List<String> getArguments() {
        return Collections.singletonList(codonAlignment.getID());
    }

    /**
     * @return a list of unique ids for the state nodes that make up the conditions
     */
    @Override
    public List<String> getConditions() {
        return siteModel.getConditions();
    }



    /** // TODO
     * Determine indices of m_fRootProbabilities that need to be updates
     * // due to sites being invariant. If none of the sites are invariant,
     * // the 'site invariant' category does not contribute anything to the
     * // root probability. If the site IS invariant for a certain character,
     * // taking ambiguities in account, there is a contribution of 1 from
     * // the 'site invariant' category.
     */
//    void calcConstantPatternIndices(final int patterns, final int stateCount) {
//        constantPattern = new ArrayList<>();
//        for (int i = 0; i < patterns; i++) {
//            final int[] pattern = dataInput.get().getPattern(i);
//            final boolean[] isInvariant = new boolean[stateCount];
//            Arrays.fill(isInvariant, true);
//            for (final int state : pattern) {
//                final boolean[] isStateSet = dataInput.get().getStateSet(state);
//                if (m_useAmbiguities.get() || !dataInput.get().getDataType().isAmbiguousCode(state)) {
//                    for (int k = 0; k < stateCount; k++) {
//                        isInvariant[k] &= isStateSet[k];
//                    }
//                }
//            }
//            for (int k = 0; k < stateCount; k++) {
//                if (isInvariant[k]) {
//                    constantPattern.add(i * stateCount + k);
//                }
//            }
//        }
//    }

    /**
     * set all nodes states in likelihood core *
     */
//    protected void setStates(TreeInterface tree) {
//
//        for (Node node : tree.getNodesAsArray()) {
//            final int[] states;
//
//            if (node.isLeaf()) {
//                // set leaf node states
//                int taxonIndex = getTaxonIndex(node.getID(), codonAlignment);
//                // no patterns
//                List<Integer> statesList = codonAlignment.getCounts().get(taxonIndex);
//
//                assert statesList.size() == codonAlignment.getSiteCount();
//                // Java 8
//                states = statesList.stream().mapToInt(i->i).toArray();
//            } else {
//                // set internal node states
//                int Nr = node.getNr();
//                states = internalNodeStates.getNrStates(Nr);
//            }
//
//            daLdCore.setNodeStates(node.getNr(), states);
//        }
//    }

//    protected void setStates(Node node) {
//
//        if (node.isLeaf()) {
//            // set leaf node states
//            int taxonIndex = getTaxonIndex(node.getID(), codonAlignment);
//            // no patterns
//            List<Integer> statesList = codonAlignment.getCounts().get(taxonIndex);
//
//            assert statesList.size() == codonAlignment.getSiteCount();
//            // Java 8
//            int[] states = statesList.stream().mapToInt(i->i).toArray();
//
//            daLdCore.setNodeStates(node.getNr(), states);
//        } else {
//            // set internal node states
//            int Nr = node.getNr();
//            int[] states = internalNodeStates.getNrStates(Nr);
//
//            daLdCore.setNodeStates(node.getNr(), states);
//
//            // updateNode
//            setStates(node.getChild(0));
//            setStates(node.getChild(1));
//        }
//    }

} // class TreeLikelihood
