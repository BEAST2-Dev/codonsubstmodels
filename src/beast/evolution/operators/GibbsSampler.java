package beast.evolution.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.evolution.likelihood.DABranchLikelihoodCore;
import beast.evolution.likelihood.GenericDATreeLikelihood;
import beast.evolution.tree.Node;
import beast.evolution.tree.NodeStatesArray;
import beast.evolution.tree.Tree;
import beast.util.RandomUtils;
import beast.util.Randomizer;

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
    final public Input<GenericDATreeLikelihood> daTreeLdInput = new Input<>("DATreeLikelihood",
            "The data augmentation tree likelihood used for the given nodes states and tree",
            Input.Validate.REQUIRED);

    public GibbsSampler() {
    }

    @Override
    public void initAndValidate() {

        final Tree tree = treeInput.get();
        final NodeStatesArray nodesStates = nodesStatesInput.get();

        nodesStates.validateTree(tree);
        nodesStates.validateNodeStates();

        // To change DATreeLikelihood is not allowed
        final GenericDATreeLikelihood daTreeLd = daTreeLdInput.get();
        assert daTreeLd.getNodesStates() == nodesStates;
        assert daTreeLd.getTree() == tree;

    }

    /**
     * @see Operator#proposal()
     */
    @Override
    public double proposal() {
        Node node = getRandomInternalNode();
        final int nodeNr = node.getNr();

        NodeStatesArray nodesStates = nodesStatesInput.get(this);
        final GenericDATreeLikelihood daTreeLd = daTreeLdInput.get();

        int state;
        for (int k = 0; k < nodesStates.getSiteCount(); k++) {
            state = gibbsSampling(node, k, nodesStates, daTreeLd);
            // change state at k for nodeNr
            nodesStates.setState(nodeNr, k, state);
        }
        // Gibbs operator should always be accepted
        return 0.0;//test store/restore //Double.POSITIVE_INFINITY;
    }

    /**
     * randomly select an internal node
     * @return {@link Node}, or null if no internal nodes
     */
    public Node getRandomInternalNode() {
        final Tree tree = treeInput.get(this);
        final int tipsCount = tree.getLeafNodeCount();
        // Abort if no internal nodes
        if (tipsCount < 2)
            throw new IllegalArgumentException("Tree must have at least 2 tips ! " + tipsCount);

        int nodeNr;
        do {
            // internal node Nr = [tipsCount, 2*tipsCount-2]
            nodeNr = tipsCount + Randomizer.nextInt(tipsCount-1);
        } while (nodeNr < tree.getLeafNodeCount() || nodeNr >= tree.getNodeCount());

        return tree.getNode(nodeNr);
    }

    /**
     * Gibbs sampling
     * @param node     {@link Node}
     * @param siteNr   the site (codon) index
     * @param nodesStates   {@link NodeStatesArray}
     * @return         the proposed state at the node, or -1 if node is null
     */
    protected int gibbsSampling(Node node, int siteNr, NodeStatesArray nodesStates,
                                final GenericDATreeLikelihood daTreeLd) {
        final int nodeNr = node.getNr();
        final int ch1Nr = node.getChild(0).getNr();
        final int ch2Nr = node.getChild(1).getNr();
        // states at child nodes
        final int x = nodesStates.getState(ch1Nr, siteNr);
        final int y = nodesStates.getState(ch2Nr, siteNr);

        final int stateCount = nodesStates.getStateCount();
        double[] pr_w = new double[stateCount];

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

    @Override
    public String getCitations() {
        return null;
    }
}
