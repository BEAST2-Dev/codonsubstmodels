package beast.evolution.likelihood;

import beast.core.Input;
import beast.core.State;
import beast.core.util.Log;
import beast.evolution.tree.Node;
import beast.evolution.tree.NodeStates;
import beast.evolution.tree.Tree;
import beast.util.ThreadHelper;

import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;

/**
 * Data augmentation to fast codon tree likelihood calculation.
 * No pattern, use site count.
 *
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

    final public Input<Integer> maxNrOfBranchesPerTaskInput = new Input<>("branchesPerTask",
            "maximum number of branches assigned to a task for thread pool, " +
                    "if less than 1 use Java default.", -1);


    /****** calculation engine ******/
//    protected BeagleTreeLikelihood beagle;
    /**
     * calculation engine for each branch, excl. root index, nrOfNodes-1
     */
    protected DABranchLikelihoodCore[] daBranchLdCores;
    protected DABranchLikelihoodCore daRootLdCores;
    /**
     * number of threads to use, changes when threading causes problems
     */
    private ThreadHelper threadHelper;
//    private int threadCount;
    private int maxBrPerTask;
    // New thread pool from BEAST MCMC.
//    private ExecutorService executor = null;
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

    /****** TODO rest ******/
    // dealing with proportion of site being invariant
    double proportionInvariant = 0;


    public DataAugTreeLikelihood(){}

//    public DataAugTreeLikelihood(NodeStatesArray nodesStates, TreeInterface tree, SiteModel.Base siteModel,
//                                 BranchRateModel.Base branchRateModel, int threadCount) {
//        initByName("nodesStates", nodesStates, "tree", tree, "siteModel", siteModel,
//                "branchRateModel", branchRateModel, "threads", threadCount);
////TODO
//    }

    @Override
    public void initAndValidate() {
        // data, tree and models
        super.initAndValidate();

        threadHelper = new ThreadHelper(maxNrOfThreadsInput.get(), null);
        if (threadHelper.getThreadCount() > 1) {
            // set ThreadHelper and add Callers in NodeStatesArray
            nodesStates.setThreadHelper(threadHelper);
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

        // tree has to be here, otherwise BEAST will treat NodeStatesArray as CalculationNode
        nodesStates.initAllNodesStates(tree);

        // init DALikelihoodCore
        initLikelihoodCore();

        // TODO
        proportionInvariant = siteModel.getProportionInvariant();
        siteModel.setPropInvariantIsCategory(false);
//        int patterns = codonAlignment.getPatternCount();
        if (proportionInvariant > 0) {
            throw new UnsupportedOperationException("proportionInvariant in dev");
//            calcConstantPatternIndices(patterns, stateCount);
        }

    }

    protected void initLikelihoodCore() {
        final int stateCount = nodesStates.getStateCount();
        // no pattern, use getSiteCount()
        final int siteCount = nodesStates.getSiteCount();
        // nodeCount = 2 * codonAlignment.getTaxonCount() - 1
        final int nodeCount = nodesStates.getNodeCount();
        // caching branch lengths and log likelihoods for each of branch,
        branchLengths = new double[nodeCount-1];
        storedBranchLengths = new double[nodeCount-1];
        // root index is used to store frequencies prior at root
        branchLogLikelihoods = new double[nodeCount];
        storedBranchLogLikelihoods = new double[nodeCount];

        // branch Likelihood excl. root index
        daBranchLdCores = new DABranchLikelihoodCore[nodeCount-1];
        // init likelihood core using branch index (child node index below the branch)
        for (int n = 0; n < getRootIndex(); n++) {
            final Node node = tree.getNode(n);
            assert node.getNr() == n;
            // make every nodes dirty
            node.makeDirty(Tree.IS_FILTHY);
            // init by the node below the branch
            daBranchLdCores[n] = new DABranchLikelihoodCore(n, stateCount, siteCount, siteModel.getCategoryCount());
        }
        tree.getRoot().makeDirty(Tree.IS_FILTHY);
        // root special
        daRootLdCores = new DABranchLikelihoodCore(getRootIndex(), stateCount, siteCount);

        // multi-threading
        int threadCount = threadHelper.getThreadCount();
        if (threadCount > 1) {
            Log.info("Data augmentation tree likelihood threads = " + threadCount);

            maxBrPerTask = maxNrOfBranchesPerTaskInput.get();
            if (maxBrPerTask > 1) {
                int tasks = (int) Math.floor((double)(nodeCount-1) / (double)maxBrPerTask);
                Log.info("Set " + maxBrPerTask + " branches per task, total tasks = " + tasks);
                if ( tasks < threadCount)
                    throw new IllegalArgumentException("The number of groups " + ((nodeCount-1) / maxBrPerTask) +
                            " cannot be less than number of threads " + threadCount + " !");

                // group the BrLd cores [0, nodeCount-2] by maxBrPerTask to reduce the lock/unlock
                for (int i = 0; i < tasks; i++) {
                    List<DABranchLikelihoodCore> cores = new ArrayList<>();
                    for (int n = 0; n < getRootIndex(); n++) { // n is branchNr
                        if (n % tasks == i)
                            cores.add(daBranchLdCores[n]);
                    }
                    // add root special calculation method
                    if (i == (tasks-1))
                        cores.add(daRootLdCores);

                    DABranchLikelihoodCallable group = new DABranchLikelihoodCallable(cores);
                    likelihoodCallers.add(group);
                }
            } else {
                // TODO if number of element is small, esp. 1, this is slower than without List.
                for (int n = 0; n < getRootIndex(); n++) { // n is branchNr
                    List<DABranchLikelihoodCore> cores = new ArrayList<>();
                    cores.add(daBranchLdCores[n]);
                    likelihoodCallers.add(new DABranchLikelihoodCallable(cores));
                }
                // add root special calculation method
                List<DABranchLikelihoodCore> cores = Collections.singletonList(daRootLdCores);
                likelihoodCallers.add(new DABranchLikelihoodCallable(cores));
            }

//            for (int i = 0; i < threadCount; i++) {
//                List<DABranchLikelihoodCore> cores = new ArrayList<>();
//                for (int n = 0; n < getRootIndex(); n++) { // n is branchNr
//                    if (n % threadCount == i)
//                        cores.add(daBranchLdCores[n]);
//                }
//                DABranchLikelihoodCallable group = new DABranchLikelihoodCallable(cores);
//                likelihoodCallers.add(group);
//            }

        }

//        final int matrixSize = stateCount * stateCount; // matrixSize = stateCount * stateCount;
//        // transition probability matrix, P
//        probabilities = new double[matrixSize];
//        Arrays.fill(probabilities, 1.0);

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
        threadHelper.shutdown();
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

        if (threadHelper.getThreadCount() <= 1) {
            // exclude root node, branches = nodes - 1
            final int rootIndex = getRootIndex();

            final double[] frequencies = substitutionModel.getFrequencies();
            NodeStates rootStates = nodesStates.getNodeStates(rootIndex);

            // branch likelihoods indexes excludes root index
            for (int n = 0; n < rootIndex; n++) {
                final Node node = tree.getNode(n);
                DABranchLikelihoodCore daBranchLdCore = daBranchLdCores[n];
                try {
                    // caching branchLogLikelihoods[nodeNr]
                    if (updateBranch(daBranchLdCore, node) != Tree.IS_CLEAN) {
                        this.branchLogLikelihoods[n] = daBranchLdCore.calculateBranchLogLikelihood();
                    }
//            System.out.println("logP = " + logP);
                } catch (ArithmeticException e) {
                    Log.err.println(e.getMessage());
                    return Double.NEGATIVE_INFINITY;
                }
            } // end n loop
            // root special
            this.branchLogLikelihoods[rootIndex] =
                    daRootLdCores.calculateRootLogLikelihood(rootIndex, rootStates, frequencies);

        } else {
            try {
                // include root special
                threadHelper.invokeAll(likelihoodCallers);
            } catch (ArithmeticException e) {
                Log.err.println(e.getMessage());
                return Double.NEGATIVE_INFINITY;
            } catch (InterruptedException e) {
                e.printStackTrace();
            }

        } // end if
        // TODO how to check the dirty condition at root?
        // if (nodesStates.isNodeStatesDirty(rootIndex) || siteModel.isDirtyCalculation()) {
//        if (tree.getRoot().isDirty() != Tree.IS_CLEAN || siteModel.isDirtyCalculation()) {
        // nodeLogLikelihoods[rootIndex] = log likelihood for frequencies prior at root
//        this.branchLogLikelihoods[rootIndex] = calculateRootLogLikelihood(rootIndex);
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

    // log likelihood at root given codon frequencies
//    protected double calculateRootLogLikelihood(int rootIndex) {
//        final double[] frequencies = substitutionModel.getFrequencies();
//
//        NodeStates rootStates = nodesStates.getNodeStates(rootIndex);
//        int siteCount = nodesStates.getSiteCount();
//        double[] siteLdAtRoot = new double[siteCount];
//        int state;
//        for (int k = 0; k < siteCount; k++) {
//            // hard code for root node
//            state = rootStates.getState(k); // 0-63
//            siteLdAtRoot[k] = frequencies[state];
//        }
//
//        return daRootLdCores.logIntegratedLikelihood(siteLdAtRoot); //+ getLogScalingFactor(k); TODO
//    }

    //    protected synchronized void sumLogP(double logPBr) {
//        logP += logPBr;
//    }
    //TODO calculateNodeBrLdOverCategories
    protected int updateBranch(final DABranchLikelihoodCore daBranchLdCore, final Node node) {
        // the branch between node and parent
        // root is excluded from node when creating DABranchLikelihoodCore
        final Node parent = node.getParent();
        final int nodeNr = node.getNr();
        final int parentNum = parent.getNr();

        // if tips, always false
        //TODO cache per site to avoid recalculation, when only sequence at a site is changed
        boolean seqUpdate = isSeqUpdate(nodeNr, parentNum);

        int nodeUpdate = node.isDirty() | parent.isDirty();

        final double branchRate = branchRateModel.getRateForBranch(node);
        // do not use getLength
        final double branchTime = (parent.getHeight() - node.getHeight()) * branchRate;

//TODO deal with 0 branch length, such as SA
        if (branchTime == 0)
            throw new UnsupportedOperationException("0 branch length, such as SA, not supported !");
        if (branchTime < 1e-10)
            throw new ArithmeticException("Reject proposal : " +
                    "time is 0 at the branch between parent node " + parentNum + " and node " + nodeNr +
                    " !\n" + "branch length = " + node.getLength() + ", branchRate = " + branchRate);

        //TODO how to distinguish branch len change and internal node seq change, when topology is same

        // ====== 1. update the transition probability matrix(ices) if the branch len changes ======
        if (seqUpdate || nodeUpdate != Tree.IS_CLEAN || branchTime != branchLengths[nodeNr]) {
            setTransProbMatrix(daBranchLdCore, node, parent, nodeNr, branchRate, branchTime);
            nodeUpdate |= Tree.IS_DIRTY;
        }
// TODO only some sites are changed
//       else if (seqUpdate) { }

        // ====== 2. recalculate likelihood if either child node wasn't clean ======
        if (nodeUpdate != Tree.IS_CLEAN) {
            calculateBranchLd(daBranchLdCore, node, nodeNr, parentNum);
        }

        return nodeUpdate;
    }

    public boolean isSeqUpdate(int nodeNr, int parentNum) {
        return nodesStates.isNodeStatesDirty(nodeNr) || nodesStates.isNodeStatesDirty(parentNum);
    }

    public void setTransProbMatrix(DABranchLikelihoodCore daBranchLdCore, Node node, Node parent, int nodeNr, double branchRate, double branchTime) {
        this.branchLengths[nodeNr] = branchTime;
        daBranchLdCore.setNodeMatrixForUpdate(); // TODO review the index
        // rate category
        for (int i = 0; i < siteModel.getCategoryCount(); i++) {
            final double jointBranchRate = siteModel.getRateForCategory(i, node) * branchRate;
            // pass the reference of array for caching probability,
            // Note: cannot move it in this class because of multithreading by nodes.
            double[] probabilities = daBranchLdCore.getProbRef();

            double[] iexp = daBranchLdCore.getIexpRef();
            // this new code is faster
            substitutionModel.getTransiProbs(parent.getHeight(), node.getHeight(),
                    jointBranchRate, iexp, probabilities);
            //System.out.println(node.getNr() + " " + Arrays.toString(probabilities));

//                for (int j=0; j < probabilities.length; j++)
//                    //TODO P(t) cannot be 0, but short branch causes numeric precision error.
//                    if (probabilities[j] <= 0) {
////                        System.err.println(Arrays.toString(probabilities));
//                        throw new ArithmeticException("Select " + probabilities[j] + " probability in P(t) matrix " +
//                                ", possibly caused by a short branch ! " + "matrix index = " + j +
//                                ", branch Nr = " + nodeNr + ", branchTime = " + branchTime);
//                    }

            daBranchLdCore.setNodeMatrix(i, probabilities); //cannot rm arraycopy
        }
    }

    public void calculateBranchLd(DABranchLikelihoodCore daBranchLdCore, Node node, int nodeNr, int parentNum) {
        // code in SiteModel, node is not used
        final double[] proportions = siteModel.getCategoryProportions(node);
        final int[] nodeStates = nodesStates.getStates(nodeNr);
        final int[] parentNodeStates = nodesStates.getStates(parentNum);

        // brLD is linked to the child node index down
        daBranchLdCore.setBranchLdForUpdate();
        // populate branchLd[][excl. root], nodeIndex is child
        daBranchLdCore.calculateBranchLd(parentNodeStates, nodeStates, proportions);
    }


    public ThreadHelper getThreadHelper() {
        return threadHelper;
    }


    // ArithmeticException if branch time < 1e-10
    class DABranchLikelihoodCallable implements Callable<Double> {
//        private final DABranchLikelihoodCore brLDCore;
//        private final int branchNr; // used to make thread safe
//
//        // per branch
//        public DABranchLikelihoodCallable(DABranchLikelihoodCore brLDCore, int branchNr) {
//            this.brLDCore = brLDCore;
//            this.branchNr = branchNr;
//        }

        private final List<DABranchLikelihoodCore> brLDCores;

        // by group
        public DABranchLikelihoodCallable(List<DABranchLikelihoodCore> brLDCores) {
            this.brLDCores = brLDCores;
        }

        @Override
        public Double call() throws Exception {
//            try {
            double logP = 0;
            for (DABranchLikelihoodCore core : brLDCores) {
                final int branchNr = core.getBranchNr();
                final Node node = tree.getNode(branchNr);

                final int rootIndex = getRootIndex();
                if (branchNr == rootIndex) {
                    final double[] freqs = substitutionModel.getFrequencies();
                    NodeStates rootStates = nodesStates.getNodeStates(rootIndex);

                    branchLogLikelihoods[rootIndex] = core.calculateRootLogLikelihood(rootIndex, rootStates, freqs);
                    logP += branchLogLikelihoods[rootIndex];
                } else {
                    // caching branchLogLikelihoods[nodeNr]
                    if (updateBranch(core, node) != Tree.IS_CLEAN)
                        branchLogLikelihoods[branchNr] = core.calculateBranchLogLikelihood();

                    logP += branchLogLikelihoods[branchNr];
                }
            }
//            } catch (Exception e) {
//                System.err.println("Something wrong to calculate branch likelihood above node " +
//                        branchNr + " during multithreading !");
//                e.printStackTrace();
//                System.exit(0);
//            }
//            System.out.println("Branch likelihood logP = " + branchLogLikelihoods[branchNr] + " above node " + branchNr);
            return logP;
        }

    }

//    class Store implements Runnable{
//        private final int branchNr;
//        public Store(int branchNr){
//            this.branchNr = branchNr;
//        }
//        public void run(){
////            try{
//                daBranchLdCores[branchNr].store();
//                storedBranchLengths[branchNr] = branchLengths[branchNr];
//                storedBranchLogLikelihoods[branchNr] = branchLogLikelihoods[branchNr];
////            }catch(Exception err){
////                err.printStackTrace();
////            }
//        }
//    }
//
//    class Restore implements Runnable{
//        private final int branchNr;
//        public Restore(int branchNr){
//            this.branchNr = branchNr;
//        }
//        public void run(){
////            try{
//            daBranchLdCores[branchNr].restore();
//            branchLengths[branchNr] = storedBranchLengths[branchNr];
//            branchLogLikelihoods[branchNr] = storedBranchLogLikelihoods[branchNr];
////            }catch(Exception err){
////                err.printStackTrace();
////            }
//        }
//    }


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
        super.store(); // storedLogP = logP; isDirty = false

//        if (threadCount > 1) { // multi-threading
//            // daBranchLdCores.length = nodeCount-1
//            for (int n = 0; n < daBranchLdCores.length; n++) {
//                executor.execute(new Store(n));
//            }
//            // root index
//            storedBranchLogLikelihoods[daBranchLdCores.length] = branchLogLikelihoods[daBranchLdCores.length];
//        } else {
            for (DABranchLikelihoodCore daBrLdCore : daBranchLdCores)
                daBrLdCore.store();

            System.arraycopy(branchLengths, 0, storedBranchLengths, 0, getNrOfBranches());
            System.arraycopy(branchLogLikelihoods, 0, storedBranchLogLikelihoods, 0, getNrOfBranches()+1);
//        }
    }

    @Override
    public void restore() {
//        if (beagle != null) {
//            beagle.restore();
//            super.restore();
//            return;
//        }
        super.restore(); // logP = storedLogP; isDirty = false

//        if (threadCount > 1) { // multi-threading
//            // daBranchLdCores.length = nodeCount-1
//            for (int n = 0; n < daBranchLdCores.length; n++) {
//                executor.execute(new Restore(n));
//            }
//            // root index
//            branchLogLikelihoods[daBranchLdCores.length] = storedBranchLogLikelihoods[daBranchLdCores.length];
//        } else {
            for (DABranchLikelihoodCore daBrLdCore : daBranchLdCores)
                daBrLdCore.restore();

            // pass reference in restore, but have to copy array in store.
            double[] tmp1 = branchLengths;
            branchLengths = storedBranchLengths;
            storedBranchLengths = tmp1;

            double[] tmp2 = branchLogLikelihoods;
            branchLogLikelihoods = storedBranchLogLikelihoods;
            storedBranchLogLikelihoods = tmp2;
//        }
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

    /**
     * @param branchNr the node Nr of child node below the branch
     * @return {@link DABranchLikelihoodCore} containing cached values.
     */
    public DABranchLikelihoodCore getDaBranchLdCores(int branchNr) {
        return daBranchLdCores[branchNr];
    }


} // class DATreeLikelihood
