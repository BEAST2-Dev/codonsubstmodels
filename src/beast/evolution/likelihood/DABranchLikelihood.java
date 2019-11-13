//package beast.evolution.likelihood;
//
//import beast.core.Distribution;
//import beast.core.State;
//import beast.evolution.branchratemodel.BranchRateModel;
//import beast.evolution.sitemodel.SiteModel;
//import beast.evolution.substitutionmodel.SubstitutionModel;
//import beast.evolution.tree.Node;
//import beast.evolution.tree.Tree;
//
//import java.util.Arrays;
//import java.util.Collections;
//import java.util.List;
//import java.util.Random;
//
///**
// * the branch is defined above the selected node.
// */
//@Deprecated
//public class DABranchLikelihood extends Distribution {
//
//    // exclude root
//    protected final Node node;
//    protected final Node parent;
//
//    /**
//     * calculation engine *
//     */
//    protected AbstrDABranchLikelihoodCore daBrLdCore;
////    protected BeagleTreeLikelihood beagle;
//
//    // data states in tip/internal nodes: 0-63.
//    // length is nrOfSites
//    protected int[] nodeStates;
//    protected final int[] parentNodeStates;
//
//    /**
//     * BEASTObject associated with inputs. Since none of the inputs are StateNodes, it
//     * is safe to link to them only once, during initAndValidate.
//     */
//    final protected SubstitutionModel substitutionModel;
//    final protected SiteModel.Base siteModel;
//    final protected BranchRateModel.Base branchRateModel;
//
//    /**
//     * flag to indicate the
//     * // when CLEAN=0, nothing needs to be recalculated for the node
//     * // when DIRTY=1 indicates a node partial needs to be recalculated
//     * // when FILTHY=2 indicates the indices for the node need to be recalculated
//     * // (often not necessary while node partial recalculation is required)
//     */
//    protected int hasDirt;
//
//    /**
//     * Lengths of the branches in the tree associated with each of the nodes
//     * in the tree through their node  numbers. By comparing whether the
//     * current branch length differs from stored branch lengths, it is tested
//     * whether a node is dirty and needs to be recomputed (there may be other
//     * reasons as well...).
//     * These lengths take branch rate models in account.
//     */
//    protected double branchLength;
//    protected double storedBranchLength;
//
//    /**
//     * memory allocation for log likelihoods for each of node *
//     */
//    protected double nodeLogLikelihoods;
//    /**
//     * memory allocation for probability tables obtained from the SiteModel *
//     */
//    protected double[] probabilities;
//
//    //TODO need a states class? to fix requiresRecalculation()
//    // stats[] are passing by reference
//    public DABranchLikelihood(Node node, int[] nodeStates,  int[] parentNodeStates,
//                              SubstitutionModel substitutionModel, SiteModel.Base siteModel,
//                              BranchRateModel.Base branchRateModel) {
//        this.node = node;
//        this.nodeStates = nodeStates;
//        this.parentNodeStates = parentNodeStates;
//        this.substitutionModel = substitutionModel;
//        this.siteModel = siteModel;
//        this.branchRateModel = branchRateModel;
//
//        assert !node.isRoot();
//        parent = node.getParent();
//    }
//
//    @Override
//    public void initAndValidate() {
//
//        final int nrOfStates = getNrOfStates();
////        assert nrOfStates == 64;
//        final int nrOfSites = getNrOfSites();
//        final int categoryCount = siteModel.getCategoryCount();
//
//        daBrLdCore = new DABranchLikelihoodCore(nrOfStates, nrOfSites, categoryCount);
//
//        // matrixSize = stateCount * stateCount;
//        final int matrixSize = daBrLdCore.getMatrixSize();
//        // transition probability matrix, P
//        probabilities = new double[matrixSize];
//        Arrays.fill(probabilities, 1.0);
//
////        if (m_useAmbiguities.get() || m_useTipLikelihoods.get()) {
////            throw new UnsupportedOperationException("in development");
////        }
//        hasDirt = Tree.IS_FILTHY;
//
//        //TODO more?
//    }
//
//    public int getNrOfSites() {
//        return nodeStates.length;
//    }
//
//    public int getNrOfStates() {
//        return substitutionModel.getStateCount();
//    }
//
//    public double getNodeLogLikelihoods() {
//        //TODO beagle.getNodeLogLikelihoods(); double[]
//        return nodeLogLikelihoods;
//    }
//
//    /**
//     * This method samples the sequences based on the tree and site model.
//     */
//    @Override
//    public void sample(State state, Random random) {
//        throw new UnsupportedOperationException("Can't sample a fixed alignment!");
//    }
//
//
//    /**
//     * Calculate the log likelihood of the current state.
//     *
//     * @return the log likelihood.
//     */
//    double m_fScale = 1.01;
//    int m_nScale = 0;
//    int X = 100;
//
//    /**
//     * exclude double[] frequencies,
//     * need to be added later in DATreeLikelihood
//     */
//    @Override
//    public double calculateLogP() {
//        //TODO beagle
////        if (beagle != null) {
////            logP = beagle.calculateLogP();
////            return logP;
////        }
//        logP = 0;
//
//        //TODO strange code in SiteModel, node is not used
//        final double[] proportions = siteModel.getCategoryProportions(node);
//
//        try {
//            if (updateNode(proportions) != Tree.IS_CLEAN)
//                logP = daBrLdCore.calculateLogLikelihoods();
////            System.out.println("logP = " + logP);
//        }
//        catch (ArithmeticException e) {
//            return Double.NEGATIVE_INFINITY;
//        }
//
////        m_nScale++;
////        if (logP > 0 || (daLdCore.getUseScaling() && m_nScale > X)) {
//////            System.err.println("Switch off scaling");
//////            m_likelihoodCore.setUseScaling(1.0);
//////            m_likelihoodCore.unstore();
//////            m_nHasDirt = Tree.IS_FILTHY;
//////            X *= 2;
//////            updateNode(tree.getRoot());
//////            calcLogP();
//////            return logP;
////        } else if (logP == Double.NEGATIVE_INFINITY && m_fScale < 10 && !scaling.get().equals(TreeLikelihood.Scaling.none)) { // && !m_likelihoodCore.getUseScaling()) {
////            m_nScale = 0;
////            m_fScale *= 1.01;
////            Log.warning.println("Turning on scaling to prevent numeric instability " + m_fScale);
//////TODO            daLdCore.setUseScaling(m_fScale);
////            daLdCore.unstore();
////            hasDirt = Tree.IS_FILTHY;
////            updateNode(tree.getRoot());
////
////            calcLogP();
////
////            return logP;
////        }
//        return logP;
//    }
//
//    /**
//     * flag to indicate the
//     * // when CLEAN=0, nothing needs to be recalculated for the node
//     * // when DIRTY=1 indicates a node partial needs to be recalculated
//     * // when FILTHY=2 indicates the indices for the node need to be recalculated
//     * // (often not necessary while node partial recalculation is required)
//     */
//    protected int updateNode(final double[] proportions) {
//
//        int update = (node.isDirty() | hasDirt);
//
//        final int nodeIndex = node.getNr();
//        final int parentNum = parent.getNr();
//
//        final double branchRate = branchRateModel.getRateForBranch(node);
//        final double branchTime = node.getLength() * branchRate;
//
////TODO deal with 0 branch length, such as SA
//        if (branchTime < 1e-6)
//            throw new UnsupportedOperationException(
////                System.err.println(
//                    "Time from parent " + parentNum + " to node " + nodeIndex + " is 0 !  " +
//                            "branch length = " + node.getLength() + ", branchRate = " + branchRate);
//
//        //TODO how to distinguish branch len change and internal node seq change, when topology is same
//
//        // ====== 1. update the transition probability matrix(ices) if the branch len changes ======
//        if (update != Tree.IS_CLEAN || branchTime != branchLength) {
//            branchLength = branchTime;
//
//            daBrLdCore.setNodeMatrixForUpdate();
//            for (int i = 0; i < siteModel.getCategoryCount(); i++) {
//                final double jointBranchRate = siteModel.getRateForCategory(i, node) * branchRate;
//                //TODO take > 90% time here, how to improve?
//                substitutionModel.getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), jointBranchRate, probabilities);
//                //System.out.println(node.getNr() + " " + Arrays.toString(probabilities));
//
//                daBrLdCore.setNodeMatrix(i, probabilities);
//            }
//            update |= Tree.IS_DIRTY; // TODO review
//        }
//
//        // ====== 2. recalculate likelihood if node wasn't clean ======
//        if (update != Tree.IS_CLEAN) {
//            // brLD is linked to the child node index down
//            daBrLdCore.setNodeBrLdForUpdate();
//
//            // populate branchLd[][excl. root], nodeIndex is child
//            daBrLdCore.calculateNodeBrLdOverCategories(nodeStates, parentNodeStates, proportions);
//        }
//
//        return update;
//    }
//
//
//    /** CalculationNode methods **/
//
//    /**
//     * check state for changed variables and update temp results if necessary *
//     */
//    @Override
//    protected boolean requiresRecalculation() {
////        if (beagle != null) {
////            return beagle.requiresRecalculation();
////        }
//        hasDirt = Tree.IS_CLEAN;
//
////TODO to fix        if (dataInput.get().isDirtyCalculation()) {
////            hasDirt = Tree.IS_FILTHY;
////            return true;
////        }
//        if (siteModel.isDirtyCalculation()) {
//            hasDirt = Tree.IS_DIRTY;
//            return true;
//        }
//        if (branchRateModel != null && branchRateModel.isDirtyCalculation()) {
//            //m_nHasDirt = Tree.IS_DIRTY;
//            return true;
//        }
//        return node.isDirty() != Tree.IS_CLEAN;
//    }
//
//    @Override
//    public void store() {
////        if (beagle != null) {
////            beagle.store();
////            super.store();
////            return;
////        }
//        if (daBrLdCore != null) {
//            daBrLdCore.store();
//        }
//        super.store();
//        storedBranchLength = branchLength;
//    }
//
//    @Override
//    public void restore() {
////        if (beagle != null) {
////            beagle.restore();
////            super.restore();
////            return;
////        }
//        if (daBrLdCore != null) {
//            daBrLdCore.restore();
//        }
//        super.restore();
//        double tmp = branchLength;
//        branchLength = storedBranchLength;
//        storedBranchLength = tmp;
//    }
//
////TODO how to cooperate with requiresRecalculation() ?
//    /**
//     * set node states which is passed by reference
//     */
//    public void setNodeStates(final int[] newStates) {
//        this.nodeStates = newStates;
//    }
//
//    public void setState(final int newState, final int site) {
//        this.nodeStates[site] = newState;
//    }
//
//
//    /**
//     * @return a list of unique ids for the state nodes that form the argument
//     */
//    @Override
//    public List<String> getArguments() {
//        return Collections.singletonList(node.getID());
//    }
//
//    /**
//     * @return a list of unique ids for the state nodes that make up the conditions
//     */
//    @Override
//    public List<String> getConditions() {
//        return siteModel.getConditions();
//    }
//
//
//} // class
