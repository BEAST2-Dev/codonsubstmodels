package test.beast.evolution.likelihood;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.likelihood.TreeLikelihood;
import codonmodels.CodonFrequencies;
import codonmodels.evolution.alignment.CodonAlignment;

import org.junit.Before;
import org.junit.Test;
import test.beast.evolution.CodonData;

import static junit.framework.Assert.assertEquals;

/**
 * @author Walter Xie
 */
public class M0TreeLikelihoodTest {
    CodonAlignment codonAlignment;

    @Before
    public void setUp() {
        Alignment data = CodonData.getAlig6T333();
        // create Codon Alignment
        codonAlignment = new CodonAlignment();
        codonAlignment.initByName("data", data, "dataType", "codon", "geneticCode", "vertebrateMitochondrial");
    }

    @Test
    public void testEqualFreq() {
        CodonFrequencies codonFreq = new CodonFrequencies();
        codonFreq.initByName("pi", "equal", "data", codonAlignment);

        String newickTree = "((Human_Horai: 0.013346, Human_Arnason: 0.013346): 0.191383," +
                "(Chimp_Horai: 0.003004, Chimp_Arnason: 0.003004): 0.201725, (Gorilla_Horai: 0.003761, " +
                "Gorilla_Arnason: 0.003761): 0.200968);";

        double logP = getLogTreeLikelihoodM0("0.14593", "10.69924",
                codonFreq, codonAlignment, newickTree);
        System.out.println("logTreeLikelihood = " + logP);
        assertEquals(-2098.8642114431645, logP, 1e-6);
    }

    @Test
    public void testF1X4() {
        CodonFrequencies codonFreq = new CodonFrequencies();
        codonFreq.initByName("pi", "F1X4", "data", codonAlignment);

        String newickTree = "((Human_Horai: 0.014702, Human_Arnason: 0.014702): 0.251960, " +
                "(Chimp_Horai: 0.003339, Chimp_Arnason: 0.003339): 0.263323, (Gorilla_Horai: 0.004302, " +
                "Gorilla_Arnason: 0.004302): 0.262361);";

        double logP = getLogTreeLikelihoodM0("0.10195", "13.49811",
                codonFreq, codonAlignment, newickTree);
        System.out.println("logTreeLikelihood = " + logP);
        assertEquals(-2009.8343141884977, logP, 1e-6);
    }

    @Test
    public void testF3X4() {
        CodonFrequencies codonFreq = new CodonFrequencies();
        codonFreq.initByName("pi", "F3X4", "data", codonAlignment);

        String newickTree = "((Human_Horai: 0.012988, Human_Arnason: 0.012988): 0.245731, " +
                "(Chimp_Horai: 0.002959, Chimp_Arnason: 0.002959): 0.255761, (Gorilla_Horai: 0.003918, " +
                "Gorilla_Arnason: 0.003918): 0.254801);";

        double logP = getLogTreeLikelihoodM0("0.08000", "15.34858",
                codonFreq, codonAlignment, newickTree);
        System.out.println("logTreeLikelihood = " + logP);
        assertEquals(-1960.7661713238033, logP, 1e-6);
    }

    @Test
    public void testF60() {
        CodonFrequencies codonFreq = new CodonFrequencies();
        codonFreq.initByName("pi", "F60/F61", "data", codonAlignment);

        String newickTree = "((Human_Horai: 0.013483, Human_Arnason: 0.013483): 0.294406, " +
                "(Chimp_Horai: 0.003084, Chimp_Arnason: 0.003084): 0.304805, (Gorilla_Horai: 0.004137, " +
                "Gorilla_Arnason: 0.004137): 0.303752);";

        double logP = getLogTreeLikelihoodM0("0.08327", "15.54039",
                codonFreq, codonAlignment, newickTree);
        System.out.println("logTreeLikelihood = " + logP);
        assertEquals(-1915.2251567137826, logP, 1e-6);
    }

    @Test
    public void testLikelihoodGivenEqualFreq() {
        CodonFrequencies codonFreq = new CodonFrequencies();
        codonFreq.initByName("data", codonAlignment, "frequencies",
                "0.016666666666666666 0.016666666666666666 0.016666666666666666 0.016666666666666666 0.016666666666666666 " +
                        "0.016666666666666666 0.016666666666666666 0.016666666666666666 0.016666666666666666 " +
                        "0.016666666666666666 0.016666666666666666 0.016666666666666666 0.016666666666666666 " +
                        "0.016666666666666666 0.016666666666666666 0.016666666666666666 0.016666666666666666 " +
                        "0.016666666666666666 0.016666666666666666 0.016666666666666666 0.016666666666666666 " +
                        "0.016666666666666666 0.016666666666666666 0.016666666666666666 0.016666666666666666 " +
                        "0.016666666666666666 0.016666666666666666 0.016666666666666666 0.016666666666666666 " +
                        "0.016666666666666666 0.016666666666666666 0.016666666666666666 0.016666666666666666 " +
                        "0.016666666666666666 0.016666666666666666 0.016666666666666666 0.016666666666666666 " +
                        "0.016666666666666666 0.016666666666666666 0.016666666666666666 0.016666666666666666 " +
                        "0.016666666666666666 0.016666666666666666 0.016666666666666666 0.016666666666666666 " +
                        "0.016666666666666666 0.016666666666666666 0.016666666666666666 0.016666666666666666 " +
                        "0.016666666666666666 0.016666666666666666 0.016666666666666666 0.016666666666666666 " +
                        "0.016666666666666666 0.016666666666666666 0.016666666666666666 0.016666666666666666 " +
                        "0.016666666666666666 0.016666666666666666 0.016666666666666666");

        String newickTree = "((Human_Horai: 0.013346, Human_Arnason: 0.013346): 0.191383," +
                "(Chimp_Horai: 0.003004, Chimp_Arnason: 0.003004): 0.201725, (Gorilla_Horai: 0.003761, " +
                "Gorilla_Arnason: 0.003761): 0.200968);";

        double logP = getLogTreeLikelihoodM0("0.14593", "10.69924",
                codonFreq, codonAlignment, newickTree);
        System.out.println("logTreeLikelihood = " + logP);
        assertEquals(-2098.8642114431645, logP, 1e-6);
    }

    @Test
    public void testLikelihoodGivenF3X4Freqs() {
        CodonFrequencies codonFreq = new CodonFrequencies();
        codonFreq.initByName("data", codonAlignment, "frequencies",
                // copied from F3X4 after change decimal place to 20 in printCodonFrequencies()
                "0.02704482396143513 0.028943334782568986 0.0035462749300424878 0.012035842186810867 " +
                        "0.04634205909572173 0.04959521026403068 0.006076640861558214 0.020623750802864244 " +
                        "0.012974598350806787 0.005395377532018664 0.05781996846927511 0.06187885367307853 " +
                        "0.0075816912297460085 0.025731800537319785 0.022940460433412455 0.0245508503711222 " +
                        "0.003008086864778587 0.010209264510763691 0.03930911787794282 0.042068565887917624 " +
                        "0.005154440622405748 0.01749385908210436 0.010283654677046962 0.011005553614640989 " +
                        "0.0013484527324869526 0.004576566849652688 0.04904512230591628 0.052488024931364716 " +
                        "0.006431082262630083 0.02182670343680513 0.011580168525492549 0.012393081018010569 " +
                        "0.0015184591841374336 0.005153558443133107 0.01984294129301906 0.021235889489747552 " +
                        "0.0026019221033230294 0.008830765926429675 0.005191110028669074 0.005555519077039221 " +
                        "6.806885997857461E-4 0.00231021585381829 0.024757601675190966 0.026495552521263975 " +
                        "0.0032463610143627893 0.011017952533594923 0.012471518239643547 0.005186175901633951 " +
                        "0.01996852952905083 0.021370293853606715 0.0026183899647364665 0.008886656850014673 " +
                        "0.005223965155432802 0.005590680590185039 6.849967554805926E-4 0.002324837473146254 " +
                        "0.024914295356679515 0.026663245891651723 0.003266907603061288 0.011087686410389827");

        String newickTree = "((Human_Horai: 0.012988, Human_Arnason: 0.012988): 0.245731, " +
                "(Chimp_Horai: 0.002959, Chimp_Arnason: 0.002959): 0.255761, (Gorilla_Horai: 0.003918, " +
                "Gorilla_Arnason: 0.003918): 0.254801);";

        double logP = getLogTreeLikelihoodM0("0.08000", "15.34858",
                codonFreq, codonAlignment, newickTree);
        System.out.println("logTreeLikelihood = " + logP);
        assertEquals(-1960.7661713238033, logP, 1e-6);
    }

    // M0 get logP, 0.3, 5
    private static double getLogTreeLikelihoodM0(String omegaValue, String kappaValue, CodonFrequencies codonFreq,
                                                 Alignment data, String newickTree) {
        TreeLikelihood likelihood = CodonData.getTreeLikelihoodM0(omegaValue, kappaValue,
                codonFreq, data, newickTree, true, false);
        return likelihood.calculateLogP();
    }

}
