package beast.likelihood;

import beast.core.Input;
import beast.core.State;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.CodonAlignment;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.likelihood.*;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import beast.tree.InternalNodeStates;

import java.util.*;

/**
 * Data augmentation to fast codon tree likelihood calculation.
 *
 * No pattern, use site count.
 *
 *
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
    protected LikelihoodCore likelihoodCore;
//    protected BeagleTreeLikelihood beagle;

    /**
     * BEASTObject associated with inputs. Since none of the inputs are StateNodes, it
     * is safe to link to them only once, during initAndValidate.
     */
    protected SubstitutionModel substitutionModel;
    protected SiteModel.Base m_siteModel;
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
    protected double[] m_branchLengths;
    protected double[] storedBranchLengths;

    /**
     * memory allocation for likelihoods for each of the patterns *
     */
    protected double[] siteLogLikelihoods;
    /**
     * memory allocation for the root partials *
     */
    protected double[] m_fRootPartials;
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
        m_siteModel = (SiteModel.Base) siteModelInput.get();
        m_siteModel.setDataType(codonAlignment.getDataType());
        substitutionModel = m_siteModel.substModelInput.get();

        if (branchRateModelInput.get() != null) {
            branchRateModel = branchRateModelInput.get();
        } else {
            branchRateModel = new StrictClockModel();
        }
        m_branchLengths = new double[nodeCount];
        storedBranchLengths = new double[nodeCount];

        int stateCount = codonAlignment.getMaxStateCount();
        assert stateCount == 64;
        likelihoodCore = new DAStatesLikelihoodCore(stateCount);

        Log.info.println("  " + codonAlignment.toString(true));
        // print startup messages via Log.print*

        // TODO check
        proportionInvariant = m_siteModel.getProportionInvariant();
        m_siteModel.setPropInvariantIsCategory(false);
        // TODO no patterns, need to check
//        int patterns = codonAlignment.getPatternCount();
        if (proportionInvariant > 0) {
            throw new UnsupportedOperationException("proportionInvariant in dev");
//            calcConstantPatternIndices(patterns, stateCount);
        }

        initCore(); // init DALikelihoodCore

        int siteCount = codonAlignment.getSiteCount();

        siteLogLikelihoods = new double[siteCount];
        m_fRootPartials = new double[siteCount]; // improved from patterns * stateCount to siteCount
        // transition probability matrix, P
        matrixSize = stateCount * stateCount;
        probabilities = new double[matrixSize];
        Arrays.fill(probabilities, 1.0);

        if (codonAlignment.isAscertained) {
            throw new UnsupportedOperationException("Ascertainment correction is not available !");
        }
    }

    /** // TODO
     * Determine indices of m_fRootProbabilities that need to be updates
     * // due to sites being invariant. If none of the sites are invariant,
     * // the 'site invariant' category does not contribute anything to the
     * // root probability. If the site IS invariant for a certain character,
     * // taking ambiguities in account, there is a contribution of 1 from
     * // the 'site invariant' category.
     */
    void calcConstantPatternIndices(final int patterns, final int stateCount) {
        constantPattern = new ArrayList<>();
        for (int i = 0; i < patterns; i++) {
            final int[] pattern = dataInput.get().getPattern(i);
            final boolean[] isInvariant = new boolean[stateCount];
            Arrays.fill(isInvariant, true);
            for (final int state : pattern) {
                final boolean[] isStateSet = dataInput.get().getStateSet(state);
                if (m_useAmbiguities.get() || !dataInput.get().getDataType().isAmbiguousCode(state)) {
                    for (int k = 0; k < stateCount; k++) {
                        isInvariant[k] &= isStateSet[k];
                    }
                }
            }
            for (int k = 0; k < stateCount; k++) {
                if (isInvariant[k]) {
                    constantPattern.add(i * stateCount + k);
                }
            }
        }
    }

    // no pattern, use codonAlignment.getSiteCount()
    protected void initCore() {
        final int nodeCount = treeInput.get().getNodeCount();
        likelihoodCore.initialize(
                nodeCount,
                codonAlignment.getSiteCount(),
                m_siteModel.getCategoryCount(),
                true, m_useAmbiguities.get()
        );

        final int extNodeCount = nodeCount / 2 + 1;
        final int intNodeCount = nodeCount / 2;

        if (m_useAmbiguities.get() || m_useTipLikelihoods.get()) {
            throw new UnsupportedOperationException("in development");
//            setPartials(treeInput.get().getRoot(), codonAlignment.getSiteCount());
        } else {
//            setStates(treeInput.get().getRoot());
            setStates(treeInput.get());
        }
        hasDirt = Tree.IS_FILTHY;
        for (int i = 0; i < intNodeCount; i++) {
            likelihoodCore.createNodePartials(extNodeCount + i);
        }

    }

    /**
     * This method samples the sequences based on the tree and site model.
     */
    @Override
    public void sample(State state, Random random) {
        throw new UnsupportedOperationException("Can't sample a fixed alignment!");
    }

    /**
     * set all nodes states in likelihood core *
     */
    protected void setStates(TreeInterface tree) {

        for (Node node : tree.getNodesAsArray()) {
            final int[] states;

            if (node.isLeaf()) {
                // set leaf node states
                int taxonIndex = getTaxonIndex(node.getID(), codonAlignment);
                // no patterns
                List<Integer> statesList = codonAlignment.getCounts().get(taxonIndex);

                assert statesList.size() == codonAlignment.getSiteCount();
                // Java 8
                states = statesList.stream().mapToInt(i->i).toArray();
            } else {
                // set internal node states
                int Nr = node.getNr();
                states = internalNodeStates.getNrStates(Nr);
            }

            likelihoodCore.setNodeStates(node.getNr(), states);
        }
    }
    //TODO why speed is same?
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
//            likelihoodCore.setNodeStates(node.getNr(), states);
//        } else {
//            // set internal node states
//            int Nr = node.getNr();
//            int[] states = internalNodeStates.getNrStates(Nr);
//
//            likelihoodCore.setNodeStates(node.getNr(), states);
//
//            // traverse
//            setStates(node.getChild(0));
//            setStates(node.getChild(1));
//        }
//    }


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

        try {
            if (traverse(tree.getRoot()) != Tree.IS_CLEAN)
                calcLogP();
//            System.out.println("logP = " + logP);
        }
        catch (ArithmeticException e) {
            return Double.NEGATIVE_INFINITY;
        }

//        m_nScale++;
//        if (logP > 0 || (likelihoodCore.getUseScaling() && m_nScale > X)) {
////            System.err.println("Switch off scaling");
////            m_likelihoodCore.setUseScaling(1.0);
////            m_likelihoodCore.unstore();
////            m_nHasDirt = Tree.IS_FILTHY;
////            X *= 2;
////            traverse(tree.getRoot());
////            calcLogP();
////            return logP;
//        } else if (logP == Double.NEGATIVE_INFINITY && m_fScale < 10 && !scaling.get().equals(TreeLikelihood.Scaling.none)) { // && !m_likelihoodCore.getUseScaling()) {
//            m_nScale = 0;
//            m_fScale *= 1.01;
//            Log.warning.println("Turning on scaling to prevent numeric instability " + m_fScale);
////TODO            likelihoodCore.setUseScaling(m_fScale);
//            likelihoodCore.unstore();
//            hasDirt = Tree.IS_FILTHY;
//            traverse(tree.getRoot());
//
//            calcLogP();
//
//            return logP;
//        }
        return logP;
    }

    protected void calcLogP() {
        logP = 0.0;
        for (int i = 0; i < codonAlignment.getSiteCount(); i++) {
            //TODO no pattern, check here
//            logP += siteLogLikelihoods[i] * codonAlignment.getSiteWeight(i);
            logP += siteLogLikelihoods[i];
            System.out.println("i = " + i + " logP = " + logP + " siteLogLikelihoods = " + siteLogLikelihoods[i]);
        }
    }


    /* Assumes there IS a branch rate model as opposed to traverse() */
    int traverse(final Node node) {

        int update = (node.isDirty() | hasDirt);

        final int nodeIndex = node.getNr();

        final double branchRate = branchRateModel.getRateForBranch(node);
        final double branchTime = node.getLength() * branchRate;

        // First update the transition probability matrix(ices) for this branch
        //if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_StoredBranchLengths[nodeIndex])) {
        if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_branchLengths[nodeIndex])) {
            m_branchLengths[nodeIndex] = branchTime;
            final Node parent = node.getParent();
            likelihoodCore.setNodeMatrixForUpdate(nodeIndex);
            for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
                final double jointBranchRate = m_siteModel.getRateForCategory(i, node) * branchRate;
                substitutionModel.getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), jointBranchRate, probabilities);
                //System.out.println(node.getNr() + " " + Arrays.toString(m_fProbabilities));
                likelihoodCore.setNodeMatrix(nodeIndex, i, probabilities);
            }
            update |= Tree.IS_DIRTY;
        }

        // If the node is internal, update the partial likelihoods.
        if (!node.isLeaf()) {

            // Traverse down the two child nodes
            final Node child1 = node.getChild(0); //Two children
            final int update1 = traverse(child1);

            final Node child2 = node.getChild(1);
            final int update2 = traverse(child2);

            // If either child node was updated then update this node too
            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {

                final int childNum1 = child1.getNr();
                final int childNum2 = child2.getNr();

                likelihoodCore.setNodePartialsForUpdate(nodeIndex);
                update |= (update1 | update2);
                if (update >= Tree.IS_FILTHY) {
                    likelihoodCore.setNodeStatesForUpdate(nodeIndex);
                }

                if (m_siteModel.integrateAcrossCategories()) {
                    likelihoodCore.calculatePartials(childNum1, childNum2, nodeIndex);
                } else {
                    throw new RuntimeException("Error TreeLikelihood 201: Site categories not supported");
                    //m_pLikelihoodCore->calculatePartials(childNum1, childNum2, nodeNum, siteCategories);
                }

                // No parent this is the root of the beast.tree -
                // calculate the pattern likelihoods
                if (node.isRoot()) {
                    final double[] frequencies = substitutionModel.getFrequencies();

                    final double[] proportions = m_siteModel.getCategoryProportions(node);
                    // m_fRootPartials is integrated across categories, so length is siteCount
                    likelihoodCore.integratePartials(node.getNr(), proportions, m_fRootPartials);

                    //TODO in dev
                    if (constantPattern != null) { // && !SiteModel.g_bUseOriginal) {
                        throw new UnsupportedOperationException("proportionInvariant in dev");
//                        proportionInvariant = m_siteModel.getProportionInvariant();
//                        // some portion of sites is invariant, so adjust root partials for this
//                        for (final int i : constantPattern) {
//                            m_fRootPartials[i] += proportionInvariant;
//                        }
                    }

                    likelihoodCore.calculateLogLikelihoods(m_fRootPartials, frequencies, siteLogLikelihoods);

                    System.out.println("siteLogLikelihoods = " + Arrays.toString(siteLogLikelihoods));
                }

            }
        }
        return update;
    } // traverseWithBRM

    /* return copy of pattern log likelihoods for each of the patterns in the alignment */
    public double [] getSiteLogLikelihoods() {
//        if (beagle != null) {
//            return beagle.getSiteLogLikelihoods();
//        }
        return siteLogLikelihoods.clone();
    } // getSiteLogLikelihoods

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
        if (m_siteModel.isDirtyCalculation()) {
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
        if (likelihoodCore != null) {
            likelihoodCore.store();
        }
        super.store();
        System.arraycopy(m_branchLengths, 0, storedBranchLengths, 0, m_branchLengths.length);
    }

    @Override
    public void restore() {
//        if (beagle != null) {
//            beagle.restore();
//            super.restore();
//            return;
//        }
        if (likelihoodCore != null) {
            likelihoodCore.restore();
        }
        super.restore();
        double[] tmp = m_branchLengths;
        m_branchLengths = storedBranchLengths;
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
        return m_siteModel.getConditions();
    }

} // class TreeLikelihood
