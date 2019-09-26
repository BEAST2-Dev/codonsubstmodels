package beast.likelihood;

import beast.core.Input;
import beast.core.State;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.tree.InternalNodeSeqs;

import java.io.PrintStream;
import java.util.Random;

/**
 * Data augmentation to fast codon tree likelihood calculation.
 *
 *
 *
 * @author Walter Xie, Fabio Mendes
 */
public class DACodonTreeLikelihood extends GenericTreeLikelihood {


    final public Input<InternalNodeSeqs> internalNodeSeqsInput = new Input<>("internalNodeSeqs",
            "The large 2-d matrix to store internal node sequences.", Input.Validate.OPTIONAL);


    @Override
    public void initAndValidate() {



    }

    @Override
    public void sample(State state, Random random) {
        super.sample(state, random);
    }

    @Override
    protected boolean requiresRecalculation() {
        return super.requiresRecalculation();
    }

    @Override
    protected void accept() {
        super.accept();
    }

    @Override
    public double calculateLogP() {
        return super.calculateLogP();
    }

    @Override
    public void store() {
        super.store();
    }

    @Override
    public void restore() {
        super.restore();
    }

    @Override
    public void log(long sample, PrintStream out) {
        super.log(sample, out);
    }
}
