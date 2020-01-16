package beast.evolution.operators;

import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

/**
 * run Gibbs sampler {@link GibbsSampler} after {@link Uniform}
 * tree operator proposal is accepted.
 *
 * @author Walter Xie
 */
public class UniformGibbs extends Uniform {

    final public Input<GibbsSampler> gibbsSamplerInput = new Input<>("gibbsSampler",
            "The Gibbs sampler will be used after this tree operator proposal is accepted.",
            Input.Validate.REQUIRED);

    private GibbsSampler gibbsSampler;

    // record which node is operated
    private int opNodeNr = -1;

    @Override
    public void initAndValidate() {
        gibbsSampler = gibbsSamplerInput.get();

        super.initAndValidate();
    }

    @Override
    public void accept() {
        super.accept();

        final Tree tree = treeInput.get();
        // Gibbs sampling at the operated node
        Node node = tree.getNode(opNodeNr);
        gibbsSampler.gibbsSampling(node, null);

        // run Gibbs sampler and directly set states
        //TODO follow lineage of operated node?
//        for (Node inNode : tree.getInternalNodes()) {
//            gibbsSampler.gibbsSampling(inNode, null);
//        }
        gibbsSampler.gibbsSamplingTowardsRoot(tree.getRoot(), null);
//        gibbsSampler.gibbsSamplingAwayRoot(tree.getRoot(), null);
    }

    /**
     * change the parameter and return the hastings ratio.
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {
        final Tree tree = treeInput.get(this);

        // randomly select internal node
        final int nodeCount = tree.getNodeCount();

        // Abort if no non-root internal nodes
        if (tree.getInternalNodeCount()==1)
            return Double.NEGATIVE_INFINITY;

        Node node;
        do {
            opNodeNr = nodeCount / 2 + 1 + Randomizer.nextInt(nodeCount / 2);
            node = tree.getNode(opNodeNr);
        } while (node.isRoot() || node.isLeaf());
        final double upper = node.getParent().getHeight();
        final double lower = Math.max(node.getLeft().getHeight(), node.getRight().getHeight());
        final double newValue = (Randomizer.nextDouble() * (upper - lower)) + lower;
        node.setHeight(newValue);

        return 0.0;
    }

}
