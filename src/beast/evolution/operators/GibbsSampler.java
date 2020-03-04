package beast.evolution.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.util.Log;
import beast.evolution.likelihood.DABranchLikelihoodCore;
import beast.evolution.likelihood.DataAugTreeLikelihood;
import beast.evolution.tree.Node;
import beast.evolution.tree.NodeStates;
import beast.evolution.tree.NodeStatesArray;
import beast.evolution.tree.Tree;
import beast.util.RandomUtils;
import beast.util.Randomizer;
import beast.util.ThreadHelper;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.Future;

/**
 * Gibbs sampler to sample the internal node states.
 *
 * <p><code>w_i</code> is the ith state of total 60/61 states:
 * <p><code>w_i ~ P_z{w_i}(t) * P_{w_i}x(t) * P_{w_i}y(t)</code>
 * <p>At the root:
 * <p><code>w_i ~ P_{w_i}(t) * P_{w_i}x(t) * P_{w_i}y(t)</code>,
 * where P_{w_i}(t) is the equilibrium frequency.
 * <p>Renormalise all <code>w_i</code> and choose <code>w</code> from the distribution.
 *
 * @author Walter Xie
 */
@Description("Gibbs sampler to sample the internal node states.")
public class GibbsSampler extends Operator {

    final public Input<Tree> treeInput = new Input<>("tree",
            "beast.tree on which this operation is performed", Input.Validate.REQUIRED);
    final public Input<NodeStatesArray> nodesStatesInput = new Input<>("nodesStates",
            "States in all nodes for sampling with the beast.tree", Input.Validate.REQUIRED);
    // to get P(t)
    final public Input<DataAugTreeLikelihood> daTreeLdInput = new Input<>("DATreeLikelihood",
            "The data augmentation tree likelihood used for the given nodes states and tree",
            Input.Validate.REQUIRED);

    final public Input<Boolean> storeNodesStatesInput = new Input<>("storeNodesStates",
            "Gibbs sampler alone does not need to store states, " +
                    "set to 'false' to ignore store/restore sequences to make faster.",
            true, Input.Validate.OPTIONAL);

    final public Input<Integer> maxNrOfThreadsInput = new Input<>("threads",
            "maximum number of threads to use, if less than 1 the number of threads " +
                    "in BeastMCMC is used (default -1)", -1);


//    final public Input<String> selectionInput = new Input<>("select",
//            "The method how to select nodes during Gibbs sampling, " +
//            "including RandomN, TowardsRoot, AwayRoot", "RandomN");

    private Tree tree; // Gibbs operates on sequences but not on the tree
    private NodeStatesArray nodesStates;
    private DataAugTreeLikelihood daTreeLd; // DATreeLikelihood is not allowed to change here

    // control which node will be operated, if only working on a node
//    private int opNodeNr = -1;

    private ThreadHelper threadHelper;
    private final List<Callable<NodeStates>> callers = new ArrayList<>();
    private int threads;

    // avoid to create new T[]
    private int[] newStates; // proposed states
    private double[] cpd_w;  // cumulative probability distribution of states at w node

    private boolean storeNodesStates;

    private void validate() {
        nodesStates.validateTree(tree);
        nodesStates.validateNodeStates();

        // make sure using same object
        assert daTreeLd.getNodesStates() == nodesStates;
        assert daTreeLd.getTree() == tree;
    }

    @Override
    public void initAndValidate() {
        this.tree = treeInput.get();
        this.nodesStates = nodesStatesInput.get();
        this.daTreeLd = daTreeLdInput.get();
        validate();

//        if ("RandomN".equalsIgnoreCase(selectionInput.get()))
//            throw new UnsupportedOperationException("Only RandomN method is available at the moment ! ");

        storeNodesStates = storeNodesStatesInput.get();
        newStates = new int[nodesStates.getSiteCount()];
        cpd_w = new double[nodesStates.getStateCount()];

        this.threadHelper = daTreeLd.getThreadHelper();
        if (threadHelper == null) {
            threadHelper = new ThreadHelper(maxNrOfThreadsInput.get(), null);
        }
        threads = threadHelper.getThreadCount();

        // TODO multithreading by sites
        if (threads > 1) {
            int chunks = threads; //TODO more chunks per thread
            MyTask[] tasks = new MyTask[chunks];
            for (int i = 0; i < chunks; ++i) {
                //You'd better cast to double and round here
                tasks[i] = new MyTask(iterations / chunks * i, iterations / chunks * (i + 1));
            }


            callers.add(new GibbsSamplingByChunks(i, this));

        }

    }

    /**
     * Gibbs sampling at one internal node
     * @see Operator#proposal()
     */
    @Override
    public double proposal() {
        // nextInt between 0 and n-1
        final int nodeNr = tree.getLeafNodeCount() +
                Randomizer.nextInt(tree.getNodeCount()-tree.getLeafNodeCount());
        final Node node = tree.getNode(nodeNr);

        final NodeStatesArray nodesStates;
        if (storeNodesStates) // State makes a copy and register this operator
            nodesStates = nodesStatesInput.get(this);
        else {
            nodesStates = this.nodesStates;
            Log.info("Gibbs sampler is ignoring store/restore sequences to make faster.");
        }

        if (threads > 1) {

            List<Future<int[]>> resp = threadHelper.invokeAll(callers);
            Iterator<Future<Response>> respIt = resp.iterator();

            //You'll have to handle exceptions here
            Response bestResponse = respIt.next().get();

            while (respIt.hasNext()) {
                Response r = respIt.next().get();
                if (r.result > bestResponse.result) {
                    bestResponse = r;
                }
            }

            try {
                List<Future<int[]>> resp = threadHelper.invokeAll(callers);



            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        } else {
            // proposed states in newStates[]
            gibbsSampling(newStates, node, nodesStates, daTreeLd, 0, newStates.length, cpd_w);
            // array copy newStates to nodesStates, which is in the State that registers this operator
            nodesStates.setStates(nodeNr, newStates);
        }

        // always be accepted, change to 0.0 for debug store/restore
        return Double.POSITIVE_INFINITY;//0.0;//
    }

    /**
     * Gibbs sampling all sites for a given node, and set new states.
     * @param newStates   proposed states from Gibbs sampling
     * @param node        {@link Node}
     * @param nodesStates {@link NodeStatesArray}
     * @param daTreeLd    {@link DataAugTreeLikelihood}
     * @param startInclusive   the site (codon) index to start in loop
     * @param endExclusive     the site (codon) index to end in loop
     * @param cpd_w        avoid to create new double[]
     * @see #gibbsSampling(int, int, DABranchLikelihoodCore, DABranchLikelihoodCore,
     *                     double[], double[], double[])
     * @see #gibbsSampling(int, int, int, DABranchLikelihoodCore, DABranchLikelihoodCore,
     *                     DABranchLikelihoodCore, double[], double[])
     */
    public void gibbsSampling(int[] newStates, final Node node, final NodeStatesArray nodesStates,
                              final DataAugTreeLikelihood daTreeLd, int startInclusive,
                              int endExclusive, double[] cpd_w) {
        // get parameters
        final double[] frequencies = daTreeLd.getSubstitutionModel().getFrequencies();
        final double[] proportions = daTreeLd.getSiteModel().getCategoryProportions(node);

        final int nodeNr = node.getNr();
        final int ch1Nr = node.getChild(0).getNr();
        final int ch2Nr = node.getChild(1).getNr();

        // w-x branch
        final DABranchLikelihoodCore wxBranchLd = daTreeLd.getDaBranchLdCores(ch1Nr);
        // w-y branch
        final DABranchLikelihoodCore wyBranchLd = daTreeLd.getDaBranchLdCores(ch2Nr);

        final boolean isRoot = node.isRoot();
        // z-w branch
        int parentNr = -1;
        DABranchLikelihoodCore zwBranchLd = null;
        if (!isRoot) {
            // z-w branch
            zwBranchLd = daTreeLd.getDaBranchLdCores(nodeNr);
            parentNr = node.getParent().getNr();
        }

        // sampling all sites for a given node
//        int oldState;
//        int changes = 0;
        for (int k = startInclusive; k < endExclusive; k++) { // the site (codon) index
            // states at child nodes
            final int x = nodesStates.getState(ch1Nr, k);
            final int y = nodesStates.getState(ch2Nr, k);
            if (isRoot) {
                newStates[k] = gibbsSampling(x, y, wxBranchLd, wyBranchLd,
                        cpd_w, frequencies, proportions);
            } else {
                final int z = nodesStates.getState(parentNr, k);

                newStates[k] = gibbsSampling(x, y, z, wxBranchLd, wyBranchLd, zwBranchLd,
                        cpd_w, proportions);
            }

//            oldState = nodesStates.getState(nodeNr, k);
//            if (newStates[k] != oldState) {
//                System.out.println("Node " + nodeNr + " site " + k + " : state changes from " +
//                        oldState + " to " + newStates[k]);
//                changes++;
//            }
        }
//        System.out.println("Node " + nodeNr + " changed " + changes + " sites.");
    }

    /**
     * Gibbs sampling state at 1 site and the root.
     * As no parent, the equilibrium frequencies are used.
     * @param x            codon state at 1 site in child node x
     * @param y            codon state at 1 site in child node y
     * @param wxBranchLd   branch likelihood calculation core at branch w->x
     * @param wyBranchLd   branch likelihood calculation core at branch w->y
     * @param cpd_w        cumulative probability distribution of states at w node,
     *                     to avoid using <code>new double[]</code>.
     * @param frequencies  equilibrium frequencies
     * @param proportions  the array of proportions of sites in different categories.
     * @return the proposed state sampled prob from 3 branch likelihoods,
     *         or negative if something is wrong.
     */
    protected int gibbsSampling(final int x, final int y, final DABranchLikelihoodCore wxBranchLd,
                                final DABranchLikelihoodCore wyBranchLd, double[] cpd_w,
                                final double[] frequencies, final double[] proportions) {

        double pwx,pwy,pr_w = 0;
        // no z
        for (int w=0; w < cpd_w.length; w++) {
            // n = w + i * state + j
            pwx = wxBranchLd.calculateBranchLdAtSite(w, x, proportions);
            pwy = wyBranchLd.calculateBranchLdAtSite(w, y, proportions);
            // w_i ~ P_{w_i}(t) * P_{w_i}x(t) * P_{w_i}y(t)
            pr_w = frequencies[w] * pwx * pwy;
            // cumulate pr_w
            cumulatePr(cpd_w, pr_w, w);
        } // end w loop

        // choose final state w from the distribution
        double random = Randomizer.nextDouble() * cpd_w[cpd_w.length-1];

        return RandomUtils.binarySearchSampling(cpd_w, random);
    }

    /**
     * Gibbs sampling state at 1 site 1 node. The node is not the root.
     * @param x            codon state at 1 site in child node x
     * @param y            codon state at 1 site in child node y
     * @param z            codon state at 1 site in parent node x
     * @param wxBranchLd   branch likelihood calculation core at branch w->x
     * @param wyBranchLd   branch likelihood calculation core at branch w->y
     * @param zwBranchLd   branch likelihood calculation core at branch z->w
     * @param cpd_w        cumulative probability distribution of states at w,
     *                     to avoid using <code>new double[]</code>.
     * @param proportions  the array of proportions of sites in different categories.
     * @return the proposed state sampled prob from 3 branch likelihoods,
     *         or negative if something is wrong.
     */
    protected int gibbsSampling( final int x, final int y, final int z,
                                 final DABranchLikelihoodCore wxBranchLd,
                                 final DABranchLikelihoodCore wyBranchLd,
                                 final DABranchLikelihoodCore zwBranchLd,
                                 double[] cpd_w, final double[] proportions) {
        double pzw,pwx,pwy,pr_w = 0;
        for (int w=0; w < cpd_w.length; w++) {
            pzw = zwBranchLd.calculateBranchLdAtSite(z, w, proportions);
            pwx = wxBranchLd.calculateBranchLdAtSite(w, x, proportions);
            pwy = wyBranchLd.calculateBranchLdAtSite(w, y, proportions);
            // w_i ~ P_z{w_i}(t) * P_{w_i}x(t) * P_{w_i}y(t)
            pr_w = pzw * pwx * pwy;
            // cumulate pr_w
            cumulatePr(cpd_w, pr_w, w);
        } // end w loop

        // choose final state w from the distribution
        double random = Randomizer.nextDouble() * cpd_w[cpd_w.length-1];

        return RandomUtils.binarySearchSampling(cpd_w, random);
    }

    // cumulate pr in cpd[], w is the index of cpd[]
    private void cumulatePr(double[] cpd, double pr, int w) {
        if (w == 0)
            cpd[0] = pr;
        else
            cpd[w] = cpd[w - 1] + pr;
    }


    // multi-threading by chunks of sites
    class GibbsSamplingByChunks implements Callable<int[]> {
        private final int startInclusive;
        private final int endExclusive;
        private final NodeStatesArray nodesStates;
        private final DataAugTreeLikelihood daTreeLd;

        private Node internalNode;

        private int[] newStates;
        private double[] cpd_w;

        public GibbsSamplingByChunks(NodeStatesArray nodesStates, DataAugTreeLikelihood daTreeLd,
                                     int startInclusive, int endExclusive){
            this.nodesStates = nodesStates;
            this.daTreeLd = daTreeLd;
            this.startInclusive = startInclusive;
            this.endExclusive = endExclusive;

            newStates = new int[nodesStates.getSiteCount()];
            cpd_w = new double[nodesStates.getStateCount()];
        }

        public void setInternalNode(Node internalNode) {
            this.internalNode = internalNode;
        }

        public int[] call() throws Exception {
            // internal nodes only
            gibbsSampling(newStates, internalNode, nodesStates, daTreeLd,
                    startInclusive, endExclusive, cpd_w);
            return newStates;
        }
    }


    public DataAugTreeLikelihood getDaTreeLd() {
        return daTreeLd;
    }

    @Override
    public String getCitations() {
        return null;
    }

//    public int getCurrentNodeNr() {
//        return currentNodeNr;
//    }

    /**
     * @param opNodeNr the internal node Nr for Gibbs sampling
     */
//    public void setOpNodeNr(int opNodeNr) {
//        this.opNodeNr = opNodeNr;
//    }

}

/***********************
 //        if ("TowardsRoot".equalsIgnoreCase(selectionInput.get())) {
 //            Log.info.println("Gibbs Sampling states at all " + tree.getInternalNodeCount() + " internal nodes.");
 //
 //        } else if ("AwayRoot".equalsIgnoreCase(selectionInput.get())) {
 //            Log.info.println("Gibbs Sampling states at all " + tree.getInternalNodeCount() + " internal nodes.");
 //
 //        } else
 if ("RandomN".equalsIgnoreCase(selectionInput.get())) {
 //            randomN = randomNInput.get();
 //            int threads = threadHelper.getThreadCount();
 //            if (randomN < threads) randomN = threads;
 // avoid to access the same internal node
 //            assert randomN < tree.getInternalNodeCount();
 //            if (randomN > 0)
 //                Log.info.println("Gibbs Sampling states at random " + randomN + " internal node(s).");
 } else {
 throw new UnsupportedOperationException("Only RandomN method is available at the moment ! ");
 //            throw new IllegalArgumentException("No such method ! " + selectionInput.get());
 }

 //        if ("TowardsRoot".equalsIgnoreCase(selectionInput.get()))
 //            gibbsSamplingTowardsRoot(tree.getRoot(), this);// (node, this)
 //        else if ("AwayRoot".equalsIgnoreCase(selectionInput.get()))
 //            gibbsSamplingAwayRoot(tree.getRoot(), this);// (node, this)
 //        } else { // "Random"

 //
 public void gibbsRandomNode(final Tree tree, final Operator operator) {
 int threads = threadHelper.getThreadCount();
 // multithreading by sites
 if (false && threads > 1) {
 callers.clear();
 //            for (Integer i : ) {
 //                callers.add(new GibbsSamplingNode(i, this));
 //            }

 try {
 threadHelper.invokeAll(callers);
 } catch (InterruptedException e) {
 e.printStackTrace();
 }
 } else { // 1 thread 1 node
 // nextInt between 0 and n-1
 final int nodeNr = tree.getLeafNodeCount() +
 Randomizer.nextInt(tree.getNodeCount()-tree.getLeafNodeCount());
 final Node node = tree.getNode(nodeNr);
 gibbsSampling(node, operator);
 }

 //        final Node node;// internal node
 //        if (opNodeNr < tree.getLeafNodeCount()) {
 //            node = RandomUtils.getRandomInternalNode(tree);
 //            setOpNodeNr(node.getNr());
 //        } else {
 //            node = tree.getNode(opNodeNr);
 //        }
 //        System.out.println("GibbsSampler at node " + nodeNr);

 // sampling all sites for given node, and set new states
 //        gibbsSampling(node, null);// (node, this)

 // this makes selecting node random every time,
 // unless setCurrentNodeNr(Nr) is called before proposal()
 //        setOpNodeNr(-1);
 }

 @Deprecated
 public void gibbsSamplingTowardsRoot(final Node node, final Operator operator) {
 // Traverse down the two child nodes
 final Node child1 = node.getChild(0);
 final Node child2 = node.getChild(1);

 if (!child1.isLeaf())
 gibbsSamplingTowardsRoot(child1, operator);
 if (!child2.isLeaf())
 gibbsSamplingTowardsRoot(child2, operator);

 gibbsSampling(node, operator);
 }

 @Deprecated
 public void gibbsSamplingAwayRoot(final Node node, final Operator operator) {
 gibbsSampling(node, operator);

 // Traverse down the two child nodes
 final Node child1 = node.getChild(0);
 final Node child2 = node.getChild(1);

 if (!child1.isLeaf())
 gibbsSamplingAwayRoot(child1, operator);
 if (!child2.isLeaf())
 gibbsSamplingAwayRoot(child2, operator);
 }

 ******************/