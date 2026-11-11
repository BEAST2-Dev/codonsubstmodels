package codonmodels.evolution.operators;

import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.spec.evolution.operator.ScaleTreeOperator;

/**
 * Run Gibbs sampler {@link GibbsSampler} after the tree scale operator proposal is accepted.
 *
 * BEAST 2.8: rebased from the legacy beast.base.evolution.operator.ScaleOperator (which scaled
 * either a parameter or a tree) onto beast.base.spec.evolution.operator.ScaleTreeOperator.
 * The parameter-scaling branch of the 2.7 class (commented out there) is gone with it.
 *
 * @author Walter Xie
 * @author Ugne Stolz (updated for BEAST 2.8)
 */
public class ScaleOperatorGibbs extends ScaleTreeOperator {

    final public Input<GibbsSampler> gibbsSamplerInput = new Input<>("gibbsSampler",
            "The Gibbs sampler will be used after this tree operator proposal is accepted.",
            Input.Validate.REQUIRED);

    private GibbsSampler gibbsSampler;
    // record which node is operated
//    private int opNodeNr = -1;

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
     * Same as {@link ScaleTreeOperator#proposal()} (Bactrian scale kernel; the Hastings ratio
     * is the log Jacobian), except that the root-only case uses {@link Node#setHeightDA(double)}
     * as the 2.7 class did, and the Gibbs sampling hooks are kept in place.
     */
    @Override
    public double proposal() {
        try {
            final Tree tree = treeInput.get();
            tree.startEditing(this);

            if (rootOnlyInput.get()) {
                final Node root = tree.getRoot();
                final double scale = getScaler(root.getNr(), root.getHeight());
                final double newHeight = root.getHeight() * scale;
                if (newHeight < Math.max(root.getLeft().getHeight(), root.getRight().getHeight())) {
                    return Double.NEGATIVE_INFINITY;
                }
                root.setHeightDA(newHeight);
//                gibbsSampler.gibbsSampling(root, this);
                return Math.log(scale);
            } else {
                // scale the beast.tree; tree.scale returns the log Jacobian (dof * log(scale))
                final double scale = getScaler(0, Double.NaN);
                final double logJacobian = tree.scale(scale);
//                gibbsSampler.gibbsSamplingAwayRoot(tree.getRoot(), this);
                return logJacobian;
            }
        } catch (Exception e) {
            // whatever went wrong, we want to abort this operation...
            return Double.NEGATIVE_INFINITY;
        }
    }
}
