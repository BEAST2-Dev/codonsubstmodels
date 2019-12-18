package beast.evolution.likelihood;

import beast.app.BeastMCMC;
import beast.core.Input;
import beast.core.State;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.NodeStatesArray;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * Data augmentation to fast codon tree likelihood calculation.
 *
 * No pattern, use site count.
 *
 * TODO 1: make it working in MCMC
 */
public class DataAugTreeLikelihood extends GenericDATreeLikelihood {

//    final public Input<Boolean> m_useAmbiguities = new Input<>("useAmbiguities",
//            "flag to indicate that sites containing ambiguous states should be handled instead of ignored (the default)",
//            false);
//    final public Input<Boolean> m_useTipLikelihoods = new Input<>("useTipLikelihoods",
//            "flag to indicate that partial likelihoods are provided at the tips", false);
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

    final public Input<Integer> maxNrOfThreadsInput = new Input<>("threads",
            "maximum number of threads to use, if less than 1 the number of threads " +
                    "in BeastMCMC is used (default -1)", -1);


    /**
     * number of threads to use, changes when threading causes problems
     */
    private int threadCount;

    private ExecutorService executor = null;
    /**
     * multi-threading {@link DABranchLikelihoodCallable}
     */
    private final List<Callable<Double>> likelihoodCallers = new ArrayList<Callable<Double>>();


    /**
     * flag to indicate the
     * // when CLEAN=0, nothing needs to be recalculated for the node
     * // when DIRTY=1 indicates a node partial needs to be recalculated
     * // when FILTHY=2 indicates the indices for the node need to be recalculated
     * // (often not necessary while node partial recalculation is required)
     */
//    protected int hasDirt; // TODO not need?

    /****** caching ******/

    /**
     * Lengths of the branches in the tree associated with each of the nodes
     * in the tree through their node  numbers. By comparing whether the
     * current branch length differs from stored branch lengths, it is tested
     * whether a node is dirty and needs to be recomputed (there may be other
     * reasons as well...).
     * These lengths take branch rate models in account.
     */
    protected double[] branchLengths; // length = nodeCount-1
    protected double[] storedBranchLengths;

    /**
     * caching log likelihoods for each of branch
     * root index is used to store frequencies prior at root
     */
    protected double[] branchLogLikelihoods; // length = nodeCount
    protected double[] storedBranchLogLikelihoods;
    /**
     * caching probability tables obtained from substitutionModel,
     * size = stateCount * stateCount
     */
    protected double[] probabilities;

    /****** TODO rest ******/
    // dealing with proportion of site being invariant
    double proportionInvariant = 0;


    public DataAugTreeLikelihood(){}

    public DataAugTreeLikelihood(NodeStatesArray nodesStates, TreeInterface tree, SiteModel.Base siteModel,
                                 BranchRateModel.Base branchRateModel, int threadCount) {
        initByName("nodesStates", nodesStates, "tree", tree, "siteModel", siteModel,
                "branchRateModel", branchRateModel, "threads", threadCount);
//TODO
    }

    @Override
    public void initAndValidate() {
        // init threadCount
        threadCount = BeastMCMC.m_nThreads;
        // threadCount is overwritten by TreeLikelihood threads
        if (maxNrOfThreadsInput.get() > 0) {
            threadCount = Math.max(maxNrOfThreadsInput.get(), BeastMCMC.m_nThreads);
        }
        // threadCount is overwritten by instanceCount
        String instanceCount = System.getProperty("beast.instance.count");
        if (instanceCount != null && instanceCount.length() > 0) {
            threadCount = Integer.parseInt(instanceCount);
        }

        // data, tree and models
        super.initAndValidate();

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

        // init DALikelihoodCore
        initCore();

        // TODO
        proportionInvariant = siteModel.getProportionInvariant();
        siteModel.setPropInvariantIsCategory(false);
//        int patterns = codonAlignment.getPatternCount();
        if (proportionInvariant > 0) {
            throw new UnsupportedOperationException("proportionInvariant in dev");
//            calcConstantPatternIndices(patterns, stateCount);
        }

    }



    protected void initCore() {

        final int stateCount = nodesStates.getStateCount();
        // no pattern, use getSiteCount()
        final int siteCount = nodesStates.getSiteCount();
        final int nodeCount = nodesStates.getNodeCount();
        // caching branch lengths and log likelihoods for each of branch,
        branchLengths = new double[nodeCount-1];
        storedBranchLengths = new double[nodeCount-1];
        // root index is used to store frequencies prior at root
        branchLogLikelihoods = new double[nodeCount];
        storedBranchLogLikelihoods = new double[nodeCount];

        // init likelihood core using branch index (child node index below the branch)
        for (int n = 0; n < getRootIndex(); n++) {
            final Node node = tree.getNode(n);
            assert node.getNr() == n;
            // make every nodes dirty
            node.makeDirty(Tree.IS_FILTHY);
            // init by the node below the branch
            daBranchLdCores[n] = new DABranchLikelihoodCore(n, stateCount,
                    siteCount, siteModel.getCategoryCount());
        }
        tree.getRoot().makeDirty(Tree.IS_FILTHY);

        // multi-threading
        if (threadCount > 1) {
            executor = Executors.newFixedThreadPool(threadCount);
            for (int n = 0; n < getRootIndex(); n++) {
                // n is branchNr
                likelihoodCallers.add(new DABranchLikelihoodCallable(daBranchLdCores[n], n));
            }
        }

        final int matrixSize = stateCount * stateCount; // matrixSize = stateCount * stateCount;
        // transition probability matrix, P
        probabilities = new double[matrixSize];
        Arrays.fill(probabilities, 1.0);

//        if (m_useAmbiguities.get() || m_useTipLikelihoods.get()) {
//            throw new UnsupportedOperationException("in development");
//        }
//        hasDirt = Tree.IS_FILTHY;


    }

    /**
     * Be careful to use.
     * @see ExecutorService#shutdown()
     */
    public void shutdown() {
        executor.shutdown();
    }

    /**
     * This method samples the sequences based on the tree and site model.
     */
    @Override
    public void sample(State state, Random random) {
        throw new UnsupportedOperationException("Can't sample a fixed alignment!");
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

        // exclude root node, branches = nodes - 1
        final int rootIndex = getRootIndex();

        if (threadCount <= 1) {
            // branch likelihoods indexes excludes root index
            for (int n = 0; n < rootIndex; n++) {
                final Node node = tree.getNode(n);
                DABranchLikelihoodCore daBranchLdCore = daBranchLdCores[n];
                try {
                    // caching branchLogLikelihoods[nodeNr]
                    if (updateBranch(daBranchLdCore, node) != Tree.IS_CLEAN) {
                        branchLogLikelihoods[n] = daBranchLdCore.calculateBranchLogLikelihood();
                    }
//            System.out.println("logP = " + logP);
                } catch (ArithmeticException e) {
                    System.err.println(e.getMessage());
                    return Double.NEGATIVE_INFINITY;
                }
            } // end n loop

        } else {
            try {
                // likelihoodCallers.add(new DABranchLikelihoodCallable(daBranchLdCores[n]));
                executor.invokeAll(likelihoodCallers);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }

        } // end if
        // TODO how to check the dirty condition at root?
        // if (nodesStates.isNodeStatesDirty(rootIndex) || siteModel.isDirtyCalculation()) {
//        if (tree.getRoot().isDirty() != Tree.IS_CLEAN || siteModel.isDirtyCalculation()) {
            final double[] frequencies = substitutionModel.getFrequencies();
        // nodeLogLikelihoods[rootIndex] = log likelihood for frequencies prior at root
            branchLogLikelihoods[rootIndex] = calculateRootLogLikelihood(rootIndex, frequencies);
//        }

        // sum logP
        logP =0;
        // sum up all branch ld, plus frequencies at root
        for (int i = 0; i < branchLogLikelihoods.length; i++) {
            logP += branchLogLikelihoods[i];
        }

        //TODO Scaling

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
//            daLdCore.setUseScaling(m_fScale);
//            daLdCore.unstore();
//            hasDirt = Tree.IS_FILTHY;
//            updateNode(tree.getRoot());
//
//            calcLogP();
//
//            return logP;
//        }
//        System.out.println("tree logP = " + logP);
        return logP;
    }

    // log likelihood for frequencies prior at root
    protected double calculateRootLogLikelihood(int rootIndex, double[] frequencies) {
        int siteCount = nodesStates.getSiteCount();
        double[] rootPrior = new double[siteCount];
        for (int k = 0; k < siteCount; k++) {
            // hard code for root node
            int state = nodesStates.getState(rootIndex, k); // 0-63
            rootPrior[k] = frequencies[state];
        }

        return DABranchLikelihoodCore.integrateLogLikelihood(rootPrior,
                DABranchLikelihoodCore.scalingThreshold); //+ getLogScalingFactor(k); TODO
    }

    //    protected synchronized void sumLogP(double logPBr) {
//        logP += logPBr;
//    }
    //TODO calculateNodeBrLdOverCategories
    // cache per site to avoid recalculation, when only sequence at a site is changed
    protected int updateBranch(final DABranchLikelihoodCore daBranchLdCore, final Node node) {
        // the branch between node and parent
        // root is excluded from node when creating DABranchLikelihoodCore
        final Node parent = node.getParent();
        final int nodeIndex = node.getNr();
        final int parentNum = parent.getNr();

        // if tips, always false
        boolean seqUpdate = nodesStates.isNodeStatesDirty(nodeIndex) || nodesStates.isNodeStatesDirty(parentNum);

        int nodeUpdate = node.isDirty() | parent.isDirty(); // TODO need to review

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
        if (seqUpdate || nodeUpdate != Tree.IS_CLEAN || branchTime != branchLengths[nodeIndex]) {
            branchLengths[nodeIndex] = branchTime;
            daBranchLdCore.setNodeMatrixForUpdate(); // TODO review the index
            // rate category
            for (int i = 0; i < siteModel.getCategoryCount(); i++) {
                final double jointBranchRate = siteModel.getRateForCategory(i, node) * branchRate;
                substitutionModel.getTransitionProbabilities(node, parent.getHeight(), node.getHeight(),
                        jointBranchRate, probabilities);
                //System.out.println(node.getNr() + " " + Arrays.toString(probabilities));

                for (int j=0; j < probabilities.length; j++)
                    //TODO P(t) cannot be 0, but short branch causes numeric precision error.
                    if (probabilities[j] <= 0) {
//                        System.err.println(Arrays.toString(probabilities));
                        throw new ArithmeticException("Return -Inf as short branch causes P(t) = 0 ! " +
                                "node Nr = " + j + ", node Nr = " + nodeIndex + ", branchTime = " + branchTime);
                    }

                daBranchLdCore.setNodeMatrix(i, probabilities); //cannot rm arraycopy
            }
            nodeUpdate |= Tree.IS_DIRTY;
        }
// TODO only some sites are changed
//       else if (seqUpdate) { }

        // ====== 2. recalculate likelihood if either child node wasn't clean ======
        if (nodeUpdate != Tree.IS_CLEAN) {
            // code in SiteModel, node is not used
            final double[] proportions = siteModel.getCategoryProportions(node);
            final int[] nodeStates = nodesStates.getStates(nodeIndex);
            final int[] parentNodeStates = nodesStates.getStates(parentNum);

            // brLD is linked to the child node index down
            daBranchLdCore.setBranchLdForUpdate(); // TODO review the index
            // populate branchLd[][excl. root], nodeIndex is child
            daBranchLdCore.calculateBranchLd(parentNodeStates, nodeStates, proportions);
        }

        return nodeUpdate;
    }


    // for multi-threading
    class DABranchLikelihoodCallable implements Callable<Double> {
        private final DABranchLikelihoodCore brLDCore;
        private final int branchNr;

        // per branch
        public DABranchLikelihoodCallable(DABranchLikelihoodCore brLDCore, int branchNr) {
            this.brLDCore = brLDCore;
            this.branchNr = branchNr;
        }

        @Override
        public Double call() throws Exception {
            try {
                final Node node = tree.getNode(branchNr); //TODO synchronized?
                // caching branchLogLikelihoods[nodeNr]
                if (updateBranch(brLDCore, node) != Tree.IS_CLEAN)
                    branchLogLikelihoods[branchNr] = brLDCore.calculateBranchLogLikelihood();

            } catch (Exception e) {
                System.err.println("Something wrong to calculate branch likelihood above node " +
                        branchNr + " during multithreading !");
                e.printStackTrace();
                System.exit(0);
            }
//            System.out.println("Branch likelihood logP = " + branchLogLikelihoods[branchNr] + " above node " + branchNr);
            return branchLogLikelihoods[branchNr];
        }

    }


    /** CalculationNode methods **/

    /**
     * check state for changed variables and update temp results if necessary *
     */
    @Override
    protected boolean requiresRecalculation() {
//        if (beagle != null) {
//            return beagle.requiresRecalculation();
//        }
//        hasDirt = Tree.IS_CLEAN;
// TODO check isDirtyCalculation()?
        if (nodesStates.somethingIsDirty()) {
//            hasDirt = Tree.IS_FILTHY;
            return true;
        }
        if (siteModel.isDirtyCalculation()) {
//            hasDirt = Tree.IS_DIRTY;
            return true;
        }
        if (branchRateModel != null && branchRateModel.isDirtyCalculation()) {
            //m_nHasDirt = Tree.IS_DIRTY;
            return true;
        }
        return tree.somethingIsDirty();
    }

    @Override
    public void store() {
//        if (beagle != null) {
//            beagle.store();
//            super.store();
//            return;
//        }
        for (DABranchLikelihoodCore daBrLdCore : daBranchLdCores)
            daBrLdCore.store();

        super.store(); // storedLogP = logP; isDirty = false
        System.arraycopy(branchLengths, 0, storedBranchLengths, 0, getNrOfBranches());
        System.arraycopy(branchLogLikelihoods, 0, storedBranchLogLikelihoods, 0, getNrOfBranches()+1);
    }

    @Override
    public void restore() {
//        if (beagle != null) {
//            beagle.restore();
//            super.restore();
//            return;
//        }
        for (DABranchLikelihoodCore daBrLdCore : daBranchLdCores)
            daBrLdCore.restore();

        super.restore(); // logP = storedLogP; isDirty = false

        // pass reference in restore, but have to copy array in store.
        double[] tmp1 = branchLengths;
        branchLengths = storedBranchLengths;
        storedBranchLengths = tmp1;

        double[] tmp2 = branchLogLikelihoods;
        branchLogLikelihoods = storedBranchLogLikelihoods;
        storedBranchLogLikelihoods = tmp2;
    }

    // for testing, nodeIndex is the child node below this branch
    public void getBranchLdFromCore(int nodeIndex, double[] branchLdOut) {
        daBranchLdCores[nodeIndex].getBranchLikelihoods(branchLdOut);
    }

    public int getRootIndex() {
        int rootIndex = tree.getNodeCount() - 1;
        assert rootIndex == tree.getRoot().getNr();
        return rootIndex;
    }

    public int getNrOfBranches() {
        int branchCount = tree.getNodeCount() - 1;
        assert branchCount == branchLengths.length;
        return branchCount;
    }


} // class DATreeLikelihood
