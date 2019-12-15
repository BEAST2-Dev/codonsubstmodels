package beast.evolution.operators;

import beast.core.Operator;

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
public class GibbsSamplerOperator extends Operator {





    @Override
    public double proposal() {
        return 0;
    }

    @Override
    public void initAndValidate() {

    }

    @Override
    public String getDescription() {
        return null;
    }

    @Override
    public String getCitations() {
        return null;
    }
}
