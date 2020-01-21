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
import beast.util.ThreadHelper;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;

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

    final public Input<String> selectionInput = new Input<>("select",
            "The method how to select nodes during Gibbs sampling, " +
            "including RandomN, ByNr, TowardsRoot, AwayRoot", "RandomN");

    final public Input<Integer> randomNInput = new Input<>("randomN",
            "the number for RandomN", 0);

    final public Input<Integer> maxNrOfThreadsInput = new Input<>("threads",
            "maximum number of threads to use, if less than 1 the number of threads " +
                    "in BeastMCMC is used (default -1)", -1);


    private Tree tree; // Gibbs operates on sequences but not on the tree
    private NodeStatesArray nodesStates;
    private DataAugTreeLikelihood daTreeLd; // DATreeLikelihood is not allowed to change here

    // control which node will be operated, if only working on a node
//    private int opNodeNr = -1;

    private ThreadHelper threadHelper;
    private final List<Callable<NodeStates>> callers = new ArrayList<>();

    private int randomN = 0;

    public GibbsSampler() { }

    public GibbsSampler(Tree tree, NodeStatesArray nodesStates, DataAugTreeLikelihood daTreeLd) {
        this.tree = tree;
        this.nodesStates = nodesStates;
        this.daTreeLd = daTreeLd;
    }

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

        this.threadHelper = daTreeLd.getThreadHelper();
        if (threadHelper == null) {
            threadHelper = new ThreadHelper(maxNrOfThreadsInput.get(), null);
        }

        if ("ByNr".equalsIgnoreCase(selectionInput.get())) {
            // internal nodes only
            for (int i = tree.getLeafNodeCount(); i < tree.getNodeCount(); i++) {
                callers.add(new GibbsSamplingNode(i, this));
            }
            Log.info.println("Gibbs Sampling states at all " + tree.getInternalNodeCount() + " internal nodes.");

        } else if ("TowardsRoot".equalsIgnoreCase(selectionInput.get())) {
            Log.info.println("Gibbs Sampling states at all " + tree.getInternalNodeCount() + " internal nodes.");

        } else if ("AwayRoot".equalsIgnoreCase(selectionInput.get())) {
            Log.info.println("Gibbs Sampling states at all " + tree.getInternalNodeCount() + " internal nodes.");

        } else if ("RandomN".equalsIgnoreCase(selectionInput.get())) {
            randomN = randomNInput.get();
//            int threads = threadHelper.getThreadCount();
//            if (randomN < threads) randomN = threads;
            // avoid to access the same internal node
            assert randomN < tree.getInternalNodeCount();
            Log.info.println("Gibbs Sampling states at random " + randomN + " internal nodes.");

        } else {
            throw new IllegalArgumentException("No such method ! " + selectionInput.get());
        }

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

    /**
     * Gibbs sampling at one internal node
     * @see Operator#proposal()
     */
    @Override
    public double proposal() {
        if ("ByNr".equalsIgnoreCase(selectionInput.get())) {
            gibbsSamplingByNr(tree, this);// (node, this)

        } else if ("TowardsRoot".equalsIgnoreCase(selectionInput.get())) {
            gibbsSamplingTowardsRoot(tree.getRoot(), this);// (node, this)

        } else if ("AwayRoot".equalsIgnoreCase(selectionInput.get())) {
            gibbsSamplingAwayRoot(tree.getRoot(), this);// (node, this)

        } else { // "Random"
            gibbsRandomNNodes(randomN, tree, this);
        }

        // always be accepted, change to 0.0 for debug store/restore
        return Double.POSITIVE_INFINITY;//0.0;//
    }

    // One node each thread, if threads*2 < tree.getInternalNodeCount()
    public void gibbsRandomNNodes(final int n, final Tree tree, final Operator operator) {
        int threads = threadHelper.getThreadCount();
        if (threads > 1) {
            callers.clear();
            Set<Integer> selectedN =  RandomUtils.getRandomInt(tree.getLeafNodeCount(), tree.getNodeCount(), n);
            for (Integer i : selectedN) {
                callers.add(new GibbsSamplingNode(i, this));
            }

            try {
                threadHelper.invokeAll(callers);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        } else {
            // internal node Nr = [tipsCount, NodeCount-1]
            Set<Integer> selectedN =  RandomUtils.getRandomInt(tree.getLeafNodeCount(), tree.getNodeCount(), n);
            for (Integer i : selectedN) {
                final Node node = tree.getNode(i);
                gibbsSampling(node, operator);
            }
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

    // multi-threading
    class GibbsSamplingNode implements Callable<NodeStates> {
        private final int nodeNr;
        private final Operator operator;
        public GibbsSamplingNode(int nodeNr, Operator operator){
            this.nodeNr = nodeNr;
            this.operator = operator;
        }
        public NodeStates call() throws Exception {
            // internal nodes only
            Node internalNode = tree.getNode(nodeNr);
            gibbsSampling(internalNode, operator);// (node, this)
            return nodesStates.getNodeStates(nodeNr);
        }
    }


    public void gibbsSamplingByNr(final Tree tree, final Operator operator) {
        int threads = threadHelper.getThreadCount();
        if (threads > 1) {
            try {
                threadHelper.invokeAll(callers);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        } else {
            // all internal nodes ordered by Nr
            for (int i = tree.getLeafNodeCount(); i < tree.getNodeCount(); i++){
                Node internalNode = tree.getNode(i);
                gibbsSampling(internalNode, operator);
            }
        }

    }

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
    /**
     * Sampling all sites for given node, and set new states.
     * @param node     {@link Node}
     * @param operator if null, then direct set states without proposal.
     * @see #gibbsSampling(NodeStatesArray, Node, int, DataAugTreeLikelihood, double[])
     */
    public void gibbsSampling(final Node node, final Operator operator) {
        final NodeStatesArray nodesStates;
        if (operator == null)
            nodesStates = this.nodesStates;
        else // State makes a copy and register this operator
            nodesStates = nodesStatesInput.get(operator);

        gibbsSampling(node, nodesStates);
    }

    /**
     * Sampling all sites for given node, and set new states.
     * @param node          {@link Node}
     * @param nodesStates   {@link NodeStatesArray}
//     * @param daTreeLd      calculation {@link DataAugTreeLikelihood}
     */
    public void gibbsSampling(final Node node, final NodeStatesArray nodesStates){
//                              final DataAugTreeLikelihood daTreeLd) {
        final int stateCount = nodesStates.getStateCount();
        double[] pr_w = new double[stateCount];
        final int nodeNr = node.getNr();
        // sampling all sites for given node
        int newState;
//        int oldState;
//        int changes = 0;
        for (int k = 0; k < nodesStates.getSiteCount(); k++) {
//            oldState = nodesStates.getState(nodeNr, k);
            newState = gibbsSampling(nodesStates, node, k, daTreeLd, pr_w);
            // change state at k for nodeNr
            nodesStates.setState(nodeNr, k, newState); // TODO setStates(int, int[]) faster?
//            if (newState != oldState) {
//                System.out.println("Node " + nodeNr + " site " + k + " : state changes from " +
//                        oldState + " to " + newState);
//                changes++;
//            }
        }
//        System.out.println("Total " + changes + " sites are changed.");
    }

    /**
     * Gibbs sampling state at a given site and node.
     * If root, then use equilibrium frequencies.
     * @param nodesStates   {@link NodeStatesArray}
     * @param node     {@link Node}
     * @param siteNr   the site (codon) index
     * @param daTreeLd the cache to get P(t) in each branch {@link DABranchLikelihoodCore}.
     * @param pr_w     probability distribution of states, to avoid using <code>new</code>.
     * @return         the proposed state at the node, or -1 if node is null
     */
    protected int gibbsSampling(NodeStatesArray nodesStates, final Node node, final int siteNr,
                                final DataAugTreeLikelihood daTreeLd, double[] pr_w) {
        final int nodeNr = node.getNr();
        final int ch1Nr = node.getChild(0).getNr();
        final int ch2Nr = node.getChild(1).getNr();
        // states at child nodes
        final int x = nodesStates.getState(ch1Nr, siteNr);
        final int y = nodesStates.getState(ch2Nr, siteNr);

        final double[] proportions = daTreeLd.getSiteModel().getCategoryProportions(node);
        final double[] frequencies = daTreeLd.getSubstitutionModel().getFrequencies();

        double pzw,pwx,pwy,sum = 0;
        // w-x branch
        final DABranchLikelihoodCore wxBranchLd = daTreeLd.getDaBranchLdCores(ch1Nr);
        // w-y branch
        final DABranchLikelihoodCore wyBranchLd = daTreeLd.getDaBranchLdCores(ch2Nr);

        if (node.isRoot()) {
            // no z
            for (int w=0; w < pr_w.length; w++) {
                // n = w + i * state + j
                pwx = wxBranchLd.calculateBranchLdAtSite(w, x, proportions);
                pwy = wyBranchLd.calculateBranchLdAtSite(w, y, proportions);
                // w_i ~ P_{w_i}(t) * P_{w_i}x(t) * P_{w_i}y(t)
                pr_w[w] = frequencies[w] * pwx * pwy;
                sum += pr_w[w];
            } // end w loop

        } else {
            // z-w branch
            final DABranchLikelihoodCore zwBranchLd = daTreeLd.getDaBranchLdCores(nodeNr);

            final int parentNr = node.getParent().getNr();
            final int z = nodesStates.getState(parentNr, siteNr);

            for (int w=0; w < pr_w.length; w++) {
                pzw = zwBranchLd.calculateBranchLdAtSite(z, w, proportions);
                pwx = wxBranchLd.calculateBranchLdAtSite(w, x, proportions);
                pwy = wyBranchLd.calculateBranchLdAtSite(w, y, proportions);
                // w_i ~ P_z{w_i}(t) * P_{w_i}x(t) * P_{w_i}y(t)
                pr_w[w] = pzw * pwx * pwy;
                sum += pr_w[w];
            } // end w loop

        } // end if

        // Renormalise all w
        for (int w=0; w < pr_w.length; w++)
            pr_w[w] = pr_w[w] / sum;

        // choose final state w from the distribution
        int w = RandomUtils.randomIntegerFrom(pr_w, false);
        return w;
    }

    public DataAugTreeLikelihood getDaTreeLd() {
        return daTreeLd;
    }

    @Override
    public String getCitations() {
        return null;
    }
}
