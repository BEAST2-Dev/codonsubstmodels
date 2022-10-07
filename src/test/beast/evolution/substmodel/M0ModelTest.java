package test.beast.evolution.substmodel;

import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Alignment;
import codonmodels.CodonFrequencies;
import codonmodels.M0Model;
import codonmodels.evolution.alignment.CodonAlignment;
import codonmodels.evolution.datatype.Codon;
import codonmodels.util.StringUtils;

import org.junit.Before;
import org.junit.Test;
import test.beast.evolution.CodonData;

import java.util.Arrays;
import java.util.stream.DoubleStream;

import static junit.framework.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Instantaneous rate q_ij can be 0, but the transition prob p_ij(t) cannot.
 *
 * @author Walter Xie
 */
public class M0ModelTest {
    M0Model m0Model;
    Codon codon;

    @Before
    public void setUp() {

        Alignment data = CodonData.getAlig6T333();
        // create Codon Alignment
        CodonAlignment codonAlignment = new CodonAlignment();
        codonAlignment.initByName("data", data, "dataType", "codon", "geneticCode", "vertebrateMitochondrial");
        codon = codonAlignment.getDataType();

        // F3x4
        CodonFrequencies codonFreq = new CodonFrequencies();
        codonFreq.initByName("pi", "F3X4", "data", codonAlignment, "verbose", true);

        // omega 0.1, kappa 1
        RealParameter omega = new RealParameter("0.1");
        RealParameter kappa = new RealParameter("1");

        m0Model = new M0Model();
        m0Model.initByName("omega", omega, "kappa", kappa, "frequencies", codonFreq, "verbose", true);
    }

    @Test
    public void getTransitionProbabilities() {
        double startTime = 1;
//        double startTime = 1E-5; // when genetic distance -> 1E-5, P(t) may has 0.
        double endTime = 0;
        double rate = 1;

        System.out.println("\nfreqs = \n" + Arrays.toString(m0Model.getFrequencies()) + "\n");

        int len = m0Model.getStateCount();
        double[] prob = new double[len*len];
        double[] iexp = new double[len*len];
        m0Model.getTransiProbs(startTime, endTime, rate, iexp, prob);

//        System.out.println("relative rates :\n" + Arrays.toString(m0Model.getRelativeRates()) + "\n");
        System.out.println("renormalised rate matrix :\n" + StringUtils.
                get2DMatrixString(m0Model.getRateMatrix(), codon));
        System.out.println("P(t) :\n" + StringUtils.get2DMatrixString(prob, codon));

        // row sum to 1
        for (int i=0; i < len; i++) {
            double[] row = new double[len];
            System.arraycopy(prob, i*len, row, 0, len);
            double sum = DoubleStream.of(row).sum();
            System.out.println("row " + i + " prob sum = " + sum);
            assertEquals(1, sum, 1e-12);
        }

        for (int i=0; i < prob.length; i++)
            assertTrue("Transition prob[" + i +"] = 0 !", prob[i] > 0);
    }
}
