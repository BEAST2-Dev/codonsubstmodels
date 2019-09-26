package beast.operators;

import beast.core.Operator;

/**
 * Gibbs sampler to work on the states of internal node sequences.
 *
 *
 * @author Walter Xie, Fabio Mendes
 */
public class GibbsInternalNodeSeqsOperator extends Operator {





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
