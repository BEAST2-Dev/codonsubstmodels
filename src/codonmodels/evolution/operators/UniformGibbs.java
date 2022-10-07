package codonmodels.evolution.operators;

import beast.base.core.Input;
import beast.base.evolution.operator.Uniform;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;

/**
 * run Gibbs sampler {@link GibbsSampler} after {@link Uniform}
 * tree operator proposal is accepted.
 *
 * @author Walter Xie
 */
public class UniformGibbs extends Uniform {

//    final public Input<NodeStatesArray> nodesStatesInput = new Input<>("nodesStates",
//            "States in all nodes for sampling with the beast.tree", Input.Validate.REQUIRED);

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

//        final Tree tree = treeInput.get();
        // Gibbs sampling at the operated node
//        Node node = tree.getNode(opNodeNr);
//        gibbsSampler.gibbsSampling(node, null);

        // run Gibbs sampler and directly set states
        //TODO follow lineage of operated node?
//        gibbsSampler.gibbsSamplingByNr(tree, null);
//        gibbsSampler.gibbsSamplingTowardsRoot(tree.getRoot(), null);
//        gibbsSampler.gibbsSamplingAwayRoot(tree.getRoot(), null);
    }

    /**
     * change the parameter and return the hastings ratio.
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {
        final Tree tree = treeInput.get();
        tree.startEditing(this);

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
        final double lower = Math.max(node.getChild(0).getHeight(), node.getChild(1).getHeight());
        final double newValue = (Randomizer.nextDouble() * (upper - lower)) + lower;
//        node.setHeight(newValue);
        node.setHeightDA(newValue);

//        NodeStatesArray nodesStates = nodesStatesInput.get(this);
//        int[] states = gibbsSampler.gibbsSampling(node, nodesStates);
//        nodesStates.setStates(opNodeNr, states);
//        gibbsSampler.gibbsSampling(node, this);

        return 0.0;
    }

}
