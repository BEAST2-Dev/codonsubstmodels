package test.beast.evolution.likelihood;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.CodonAlignment;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Tree;
import codonmodels.CodonFrequencies;
import codonmodels.M0Model;
import org.junit.Before;
import org.junit.Test;
import test.beast.evolution.CodonTestData;

import static junit.framework.Assert.assertEquals;
import static org.junit.Assert.assertArrayEquals;

/**
 * @author Walter Xie
 */
public class M0TreeLikelihoodTest {
    CodonAlignment codonAlignment;

    @Before
    public void setUp() {
        Alignment data = CodonTestData.getAlignment();
        // create Codon Alignment
        codonAlignment = new CodonAlignment();
        codonAlignment.initByName("data", data, "dataType", "codon", "geneticCode", "vertebrateMitochondrial");
    }

    @Test
    public void testCodonPositionBaseFrequencies(){
        double[][] freq = codonAlignment.getCodonPositionBaseFrequencies();
        // result from codeml
        // Codon position * base (3x4) table + overall
        // position 1 : A C G T
        assertArrayEquals(new double[]{.36937,.31331,.15816,.15916}, freq[0], 0.00001);
        // position 2 : A C G T
        assertArrayEquals(new double[]{.18869,.32332,.08458,.4034},freq[1], 0.00001);
        // position 3 : A C G T
        assertArrayEquals(new double[]{.37788,.4044,.04955,.16817}, freq[2], 0.00001);
        // average    : A C G T
        assertArrayEquals(new double[]{.31198,.34701,.09743,.24358}, freq[3], 0.00001);
    }

    @Test
    public void testCodonFrequencies(){
        int[][] usage = codonAlignment.getCodonUsage();
        CodonFrequencies codonFreq = new CodonFrequencies();
        codonFreq.initByName("data", codonAlignment);

        // 60 states, exclude stop codon in freq
        double[] freq = codonFreq.getCodonFrequenciesByUsage(usage);
        // result from codeml
        assertArrayEquals(
                new double[]{.03203203,.03653654,.00500501,.01501502,.03503504,.05705706,0,.02302302,
                        0,.01051051,0,0,.04604605,.05605606,.00800801,.04504505,
                        .02902903,.02702703,.003003,.0035035,.02152152,.04304304,.003003,.01201201,
                        .00900901,.00850851,0,.001001,.08058058,.03403403,.01601602,.02202202,
                        .00800801,.00900901,.002002,.002002,.01801802,.04904905,.001001,.00800801,
                        .01051051,.01601602,.0015015,.002002,.01401401,.00800801,.003003,.00600601,
                        0,.01301301,0,.0035035,.03003003,.01201201,.002002,.00850851,
                        .01851852,.003003,.003003,.001001,.02552553,.02152152,.002002,.01551552},
                freq, 0.00000001);
    }

    @Test
    public void testEqualFreq(){
        CodonFrequencies codonFreq = new CodonFrequencies();
        codonFreq.initByName("pi", "equal", "data", codonAlignment);

        String startingTree = "((Human_Horai: 0.013346, Human_Arnason: 0.013346): 0.191383," +
                "(Chimp_Horai: 0.003004, Chimp_Arnason: 0.003004): 0.201725, (Gorilla_Horai: 0.003761, " +
                "Gorilla_Arnason: 0.003761): 0.200968);";

        double logP = getLogTreeLikelihood("0.14593", "10.69924",
                codonFreq, codonAlignment, startingTree);
        System.out.println("logTreeLikelihood = " + logP);
        assertEquals(-2098.8642114431645, logP, 1e-6);
    }

    @Test
    public void testF1X4(){
        CodonFrequencies codonFreq = new CodonFrequencies();
        codonFreq.initByName("pi", "F1X4", "data", codonAlignment);

        String startingTree = "((Human_Horai: 0.014702, Human_Arnason: 0.014702): 0.251960, " +
                "(Chimp_Horai: 0.003339, Chimp_Arnason: 0.003339): 0.263323, (Gorilla_Horai: 0.004302, " +
                "Gorilla_Arnason: 0.004302): 0.262361);";

        double logP = getLogTreeLikelihood("0.10195", "13.49811",
                codonFreq, codonAlignment, startingTree);
        System.out.println("logTreeLikelihood = " + logP);
        assertEquals(-2009.8343141884977, logP, 1e-6);
    }

    @Test
    public void testF3X4(){
        CodonFrequencies codonFreq = new CodonFrequencies();
        codonFreq.initByName("pi", "F3X4", "data", codonAlignment);

        String startingTree = "((Human_Horai: 0.012988, Human_Arnason: 0.012988): 0.245731, " +
                "(Chimp_Horai: 0.002959, Chimp_Arnason: 0.002959): 0.255761, (Gorilla_Horai: 0.003918, " +
                "Gorilla_Arnason: 0.003918): 0.254801);";

        double logP = getLogTreeLikelihood("0.08000", "15.34858",
                codonFreq, codonAlignment, startingTree);
        System.out.println("logTreeLikelihood = " + logP);
        assertEquals(-1960.7661713238033, logP, 1e-6);
    }

    @Test
    public void testF60(){
        CodonFrequencies codonFreq = new CodonFrequencies();
        codonFreq.initByName("pi", "F60/F61", "data", codonAlignment);

        String startingTree = "((Human_Horai: 0.013483, Human_Arnason: 0.013483): 0.294406, " +
                "(Chimp_Horai: 0.003084, Chimp_Arnason: 0.003084): 0.304805, (Gorilla_Horai: 0.004137, " +
                "Gorilla_Arnason: 0.004137): 0.303752);";

        double logP = getLogTreeLikelihood("0.08327", "15.54039",
                codonFreq, codonAlignment, startingTree);
        System.out.println("logTreeLikelihood = " + logP);
        assertEquals(-1915.2251567137826, logP, 1e-6);
    }

    @Test
    public void testLikelihoodGivenEqualFreq(){
        CodonFrequencies codonFreq = new CodonFrequencies();
        codonFreq.initByName("data", codonAlignment, "frequencies",
                "0.016666666666666666 0.016666666666666666 0.016666666666666666 0.016666666666666666 0.016666666666666666 " +
                        "0.016666666666666666 0.016666666666666666 0.016666666666666666 0.0 0.016666666666666666 0.0 " +
                        "0.016666666666666666 0.016666666666666666 0.016666666666666666 0.016666666666666666 " +
                        "0.016666666666666666 0.016666666666666666 0.016666666666666666 0.016666666666666666 " +
                        "0.016666666666666666 0.016666666666666666 0.016666666666666666 0.016666666666666666 " +
                        "0.016666666666666666 0.016666666666666666 0.016666666666666666 0.016666666666666666 " +
                        "0.016666666666666666 0.016666666666666666 0.016666666666666666 0.016666666666666666 " +
                        "0.016666666666666666 0.016666666666666666 0.016666666666666666 0.016666666666666666 " +
                        "0.016666666666666666 0.016666666666666666 0.016666666666666666 0.016666666666666666 " +
                        "0.016666666666666666 0.016666666666666666 0.016666666666666666 0.016666666666666666 " +
                        "0.016666666666666666 0.016666666666666666 0.016666666666666666 0.016666666666666666 " +
                        "0.016666666666666666 0.0 0.016666666666666666 0.0 0.016666666666666666 0.016666666666666666 " +
                        "0.016666666666666666 0.016666666666666666 0.016666666666666666 0.016666666666666666 " +
                        "0.016666666666666666 0.016666666666666666 0.016666666666666666 0.016666666666666666 " +
                        "0.016666666666666666 0.016666666666666666 0.016666666666666666");

        String startingTree = "((Human_Horai: 0.013346, Human_Arnason: 0.013346): 0.191383," +
                "(Chimp_Horai: 0.003004, Chimp_Arnason: 0.003004): 0.201725, (Gorilla_Horai: 0.003761, " +
                "Gorilla_Arnason: 0.003761): 0.200968);";

        double logP = getLogTreeLikelihood("0.14593", "10.69924",
                codonFreq, codonAlignment, startingTree);
        System.out.println("logTreeLikelihood = " + logP);
        assertEquals(-2098.8642114431645, logP, 1e-6);
    }

    @Test
    public void testLikelihoodGivenF3X4Freqs(){
        CodonFrequencies codonFreq = new CodonFrequencies();
        codonFreq.initByName("data", codonAlignment, "frequencies",
                // copied from F3X4 but change codon frequencies decimal place to 20
                "0.02704482396143513 0.028943334782568986 0.0035462749300424878 0.012035842186810867 " +
                        "0.04634205909572173 0.04959521026403068 0.006076640861558214 0.020623750802864244 0.0 " +
                        "0.012974598350806787 0.0 0.005395377532018664 0.05781996846927511 0.06187885367307853 " +
                        "0.0075816912297460085 0.025731800537319785 0.022940460433412455 0.0245508503711222 " +
                        "0.003008086864778587 0.010209264510763691 0.03930911787794282 0.042068565887917624 " +
                        "0.005154440622405748 0.01749385908210436 0.010283654677046962 0.011005553614640989 " +
                        "0.0013484527324869526 0.004576566849652688 0.04904512230591628 0.052488024931364716 " +
                        "0.006431082262630083 0.02182670343680513 0.011580168525492549 0.012393081018010569 " +
                        "0.0015184591841374336 0.005153558443133107 0.01984294129301906 0.021235889489747552 " +
                        "0.0026019221033230294 0.008830765926429675 0.005191110028669074 0.005555519077039221 " +
                        "6.806885997857461E-4 0.00231021585381829 0.024757601675190966 0.026495552521263975 " +
                        "0.0032463610143627893 0.011017952533594923 0.0 0.012471518239643547 0.0 0.005186175901633951 " +
                        "0.01996852952905083 0.021370293853606715 0.0026183899647364665 0.008886656850014673 " +
                        "0.005223965155432802 0.005590680590185039 6.849967554805926E-4 0.002324837473146254 " +
                        "0.024914295356679515 0.026663245891651723 0.003266907603061288 0.011087686410389827");

        String startingTree = "((Human_Horai: 0.012988, Human_Arnason: 0.012988): 0.245731, " +
                "(Chimp_Horai: 0.002959, Chimp_Arnason: 0.002959): 0.255761, (Gorilla_Horai: 0.003918, " +
                "Gorilla_Arnason: 0.003918): 0.254801);";

        double logP = getLogTreeLikelihood("0.08000", "15.34858",
                codonFreq, codonAlignment, startingTree);
        System.out.println("logTreeLikelihood = " + logP);
        assertEquals(-1960.7661713238033, logP, 1e-6);
    }

    // get logP
    protected double getLogTreeLikelihood(String omegaValue, String kappaValue, CodonFrequencies codonFreq,
                                          Alignment data, String startingTree) {
        RealParameter omega = new RealParameter(omegaValue);
        RealParameter kappa = new RealParameter(kappaValue);

        M0Model m0 = new M0Model();
        m0.initByName("omega", omega, "kappa", kappa, "frequencies", codonFreq);

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("substModel", m0);

        Tree tree = CodonTestData.getStartingTree(data, startingTree);

        System.setProperty("java.only","true");
        TreeLikelihood likelihood = new TreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel);
        return likelihood.calculateLogP();
    }

}
