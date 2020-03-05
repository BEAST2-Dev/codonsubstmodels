package beast.evolution.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.util.Log;
import beast.evolution.likelihood.DataAugTreeLikelihood;
import beast.evolution.tree.Node;
import beast.evolution.tree.NodeStatesArray;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import beast.util.ThreadHelper;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
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
public class GibbsSamplingOperator extends Operator {

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

    // Gibbs Sampler
    private GibbsSampler[] gibbsSamplers;

    private ThreadHelper threadHelper;
    private final List<Callable<GibbsSampler>> callers = new ArrayList<>();
    private int threads;

    // flag whether to ignore store/restore sequences to make faster
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
        if (!storeNodesStates)
            Log.info("Gibbs sampler is ignoring store/restore sequences to make faster.");

        this.threadHelper = daTreeLd.getThreadHelper();
        if (threadHelper == null) {
            threadHelper = new ThreadHelper(maxNrOfThreadsInput.get(), null);
        }
        threads = threadHelper.getThreadCount();

        int siteCount = nodesStates.getSiteCount();
        if (threads <= 1) {
            gibbsSamplers = new GibbsSampler[1];
            gibbsSamplers[0] = new GibbsSampler(nodesStates.getStateCount(),
                    0, siteCount);

        } else { // TODO multithreading by sites
            Log.info("Gibbs sampling at " + threads + " threads.");

            int chunks = threads; //TODO more chunks per thread
            gibbsSamplers = new GibbsSampler[chunks];
            String msg = "Split codons into " + chunks + " chunks (startInclusive - endExclusive) : ";

            int endExclusive;
            for (int i = 0; i < chunks; ++i) {
                int startInclusive = siteCount / chunks * i;
                endExclusive = siteCount / chunks * (i + 1);
                // in case chunks is the factor of siteCount
                if (i == chunks-1 && endExclusive != siteCount)
                    endExclusive = siteCount;

                msg += startInclusive + "-" + endExclusive + ", ";

                gibbsSamplers[i] = new GibbsSampler(nodesStates.getStateCount(),
                        startInclusive, endExclusive);
                callers.add(new GibbsSamplingByChunks(gibbsSamplers[i]));
            }
            msg = msg.substring(0, msg.length()-2);
            Log.info(msg);

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
        }

        if (threads <= 1) {
//            gibbsSamplers[0].update(node, nodesStates, daTreeLd);
            // proposed states in newStates[]
            int[] newStates = gibbsSamplers[0].gibbsSampling(node, nodesStates, daTreeLd);
            // array copy newStates to nodesStates, which is in the State that registers this operator
            nodesStates.setStates(nodeNr, newStates);

        } else { // multithreading by sites
            for (int i = 0; i < callers.size(); ++i) {
                GibbsSamplingByChunks gsByChunks = (GibbsSamplingByChunks) callers.get(i);
                // Note: Need to {@link #update(Node, NodeStatesArray)} before each call.
                gsByChunks.update(node, nodesStates, daTreeLd);
            }

            // Execute all tasks and get reference to Future objects
            List<Future<GibbsSampler>> resultList = null;
            try {
                resultList = threadHelper.invokeAll(callers);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }

            // collect all pieces of newStates[]
            for (int i = 0; i < resultList.size(); i++) {
                Future<GibbsSampler> future = resultList.get(i);
                try {
                    GibbsSampler gs = future.get();

                    int[] newStates = gs.getNewStates();
                    int startInclusive = gs.getStartInclusive();
                    int endExclusive = gs.getEndExclusive();

                    nodesStates.setStates(nodeNr, startInclusive, endExclusive, newStates);

                } catch (InterruptedException | ExecutionException e) {
                    e.printStackTrace();
                }
            }




        }

        // always be accepted, change to 0.0 for debug store/restore
        return Double.POSITIVE_INFINITY;//0.0;//
    }



    // multi-threading by chunks of sites
    class GibbsSamplingByChunks implements Callable<GibbsSampler> {
        private GibbsSampler gibbsSampler;

        public GibbsSamplingByChunks(GibbsSampler gibbsSampler){
            this.gibbsSampler = gibbsSampler;
        }

        public void update(Node node, NodeStatesArray nodesStates, DataAugTreeLikelihood daTreeLd) {
            this.gibbsSampler.update(node, nodesStates, daTreeLd);
        }

        /**
         * Need to {@link #update(Node, NodeStatesArray, DataAugTreeLikelihood)} before each call.
         * @return
         * @throws Exception
         */
        public GibbsSampler call() throws Exception {
            // internal nodes only
            gibbsSampler.gibbsSampling();
            return gibbsSampler;
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