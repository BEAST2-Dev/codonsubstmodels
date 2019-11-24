package test.beast.evolution.substmodel;

import beast.core.BEASTObject;
import beast.evolution.alignment.CodonAlignment;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Node;
import beast.util.DALikelihoodBenchmarking1;
import beast.util.StringUtils;
import beast.util.XMLParserException;
import codonmodels.M0Model;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.util.Arrays;
import java.util.stream.DoubleStream;

import static org.junit.Assert.assertTrue;

/**
 * Instantaneous rate q_ij can be 0, but the transition prob p_ij(t) cannot.
 *
 * @author Walter Xie
 */
public class M0ModelTest {
    M0Model m0Model;
    CodonAlignment codonAlignment;

    @Before
    public void setUp() throws IOException, XMLParserException {
//        RealParameter f = new RealParameter(new Double[]{.03495605, .02489258, .00004395, .00010254, .04003906, .07020996,
//                .00504883, .01969238, 0.0, .01495605, 0.0, .00504883, .03990723, .06005371, .0051123, .02483398, .005,
//                .03998535, 0.0, .00499023, .02514648, .03479004, .0000293, .01006348, .005, .01996094, 0.0, 0.0, .06370605,
//                .04036133, .01095703, .02984375, 0.0, .01499512, 0.0, .00000488, .01999512, .01997559, 0.0, .01006836, .005,
//                .01500977, 0.0, .01498535, .02500488, .0251123, .005, .00993164, 0.0, .04482422, 0.0, .0052002, .015,
//                .01503906, .00000977, .005, .01, .00503906, .005, 0.0, .03520508, .03976562, .00008789, .02001465});
//        assert f.getDimension() == 64;

        DALikelihoodBenchmarking1 util = new DALikelihoodBenchmarking1(
                System.getProperty("user.home") + "/WorkSpace/codonsubstmodels/evolver/1024T1K/1024T1K.xml");

        // omega 0.3, kappa 5
        BEASTObject[] models = util.initF3X4(1024, 200, true, true);
        codonAlignment = (CodonAlignment) models[0];
        SiteModel siteModel = (SiteModel) models[1];

        m0Model = (M0Model) siteModel.getSubstitutionModel();
    }

    @Test
    public void getTransitionProbabilities() {
//        double startTime = 1;
        double startTime = 1E-5;
        double endTime = 0;
        double rate = 1;

        System.out.println("freqs = \n" + Arrays.toString(m0Model.getFrequencies()) + "\n");

        int len = m0Model.getStateCount();
        double[] prob = new double[len*len];
        m0Model.getTransitionProbabilities(new Node(), startTime, endTime, rate, prob, true);

        System.out.println("relative rates :\n" + Arrays.toString(m0Model.getRelativeRates()) + "\n");
        System.out.println("renormalised rate matrix :\n" + StringUtils.
                get2DMatrixString(m0Model.getRateMatrix(), codonAlignment.getDataType()));
        System.out.println("P(t) :\n" + StringUtils.
                get2DMatrixString(prob, codonAlignment.getDataType()));

        // row sum to 1
        for (int i=0; i < len; i++) {
            double[] row = new double[len];
            System.arraycopy(prob, i*len, row, 0, len);
            double sum = DoubleStream.of(row).sum();
            System.out.println("row " + i + " prob sum = " + sum);
        }

        for (int i=0; i < prob.length; i++)
            assertTrue("Transition prob[" + i +"] = 0 !", prob[i] > 0);
    }
}
