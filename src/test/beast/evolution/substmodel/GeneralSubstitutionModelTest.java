package test.beast.evolution.substmodel;

import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.evolution.tree.Node;
import beast.util.StringUtils;
import org.junit.Before;
import org.junit.Test;

import java.util.Arrays;
import java.util.stream.DoubleStream;

import static org.junit.Assert.assertTrue;

/**
 * Instantaneous rate q_ij can be 0, but the transition prob p_ij(t) cannot.
 *
 * @author Walter Xie
 */
public class GeneralSubstitutionModelTest {
    GeneralSubstitutionModel geneSubstModel;

    @Before
    public void setUp() {
        RealParameter f = new RealParameter(new Double[]{0.3333333, 0.3333333, 0.3333333});
        Frequencies freqs = new Frequencies();
        freqs.initByName("frequencies", f, "estimate", false);

        // A -> B -> C, not A -> C
        // off-diagonal: NrOfStates * (NrOfStates - 1)
        Double[] r = new Double[]{0.1, 0.0, 0.1, 0.2, 0.0, 0.2};
        RealParameter rates = new RealParameter(r);
        geneSubstModel = new GeneralSubstitutionModel();
        geneSubstModel.initByName("frequencies", freqs, "rates", rates);
    }

    @Test
    public void getTransitionProbabilities() {
        double startTime = 1E-5;
        double endTime = 0;
        double rate = 1;

        System.out.println("freqs = \n" + Arrays.toString(geneSubstModel.getFrequencies()) + "\n");

        int len = geneSubstModel.getStateCount();
        double[] prob = new double[len*len];
        geneSubstModel.getTransitionProbabilities(new Node(), startTime, endTime, rate, prob, true);

        System.out.println("relative rates :\n" + Arrays.toString(geneSubstModel.getRelativeRates()) + "\n");
        System.out.println("renormalised rate matrix :\n" + StringUtils.get2DMatrixString(geneSubstModel.getRateMatrix(), null));
        System.out.println("transition prob :\n" + StringUtils.get2DMatrixString(prob, null));

        // P(t) row sum to 1
        for (int i=0; i < len; i++) {
            double[] row = new double[len];
            System.arraycopy(prob, i*len, row, 0, len);
            double sum = DoubleStream.of(row).sum();
            System.out.println("row " + i + " prob sum = " + sum);
        }

        for (int i=0; i < prob.length; i++)
            assertTrue(prob[i] > 0);
    }
}
