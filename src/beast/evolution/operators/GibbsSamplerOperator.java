package beast.evolution.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.evolution.tree.NodesStates;
import beast.evolution.tree.Tree;

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
public class GibbsSamplerOperator extends Operator {

    final public Input<Tree> treeInput = new Input<>("tree",
            "beast.tree on which this operation is performed", Input.Validate.REQUIRED);
    final public Input<NodesStates> nodesStatesInput = new Input<>("nodesStates",
            "States in all nodes for sampling with the beast.tree", Input.Validate.REQUIRED);




    @Override
    public void initAndValidate() {

        final Tree tree = treeInput.get();
        final NodesStates nodesStates = nodesStatesInput.get();

        nodesStates.validateTree(tree);
        nodesStates.validateNodeStates();

    }





    @Override
    public double proposal() {
        return 0;
    }


    @Override
    public String getCitations() {
        return null;
    }
}
