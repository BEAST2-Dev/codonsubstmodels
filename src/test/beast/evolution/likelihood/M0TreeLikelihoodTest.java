package test.beast.evolution.likelihood;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.CodonAlignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import codonmodels.CodonFrequencies;
import codonmodels.M0Model;
import org.junit.Before;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;
import static org.junit.Assert.assertArrayEquals;

/**
 * @author Walter Xie
 */
public class M0TreeLikelihoodTest {
    CodonAlignment codonAlignment;

    @Before
    public void setUp() {
        Alignment data = getAlignment();
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
        // 60 states, exclude stop codon in freq
        double[] freq = codonAlignment.getCodonFrequencies(false);
        // result from codeml
        assertArrayEquals(
                new double[]{.03203203,.03653654,.00500501,.01501502,.03503504,.05705706,0,.02302302,.01051051,0,
                .04604605,.05605606,.00800801,.04504505,.02902903,.02702703,.003003,.0035035,.02152152,.04304304,
                .003003,.01201201,.00900901,.00850851,0,.001001,.08058058,.03403403,.01601602,.02202202,
                .00800801,.00900901,.002002,.002002,.01801802,.04904905,.001001,.00800801,.01051051,.01601602,
                .0015015,.002002,.01401401,.00800801,.003003,.00600601,.01301301,.0035035,.03003003,.01201201,
                .002002,.00850851,.01851852,.003003,.003003,.001001,.02552553,.02152152,.002002,.01551552},
                freq, 0.00000001);
    }

    @Test
    public void testEqualFreq(){
        CodonFrequencies codonFreq = new CodonFrequencies();
        codonFreq.initByName("pi", "equal", "data", codonAlignment);

        double logP = getLogTreeLikelihood("0.14593", "10.69924", codonFreq, codonAlignment);
        System.out.println("logTreeLikelihood = " + logP);
        assertEquals(-2098.8642114431645, logP, 1e-6);
    }

    @Test
    public void testF1X4(){
        CodonFrequencies codonFreq = new CodonFrequencies();
        codonFreq.initByName("pi", "F1X4", "data", codonAlignment);

        double logP = getLogTreeLikelihood("0.10195", "13.49811", codonFreq, codonAlignment);
        System.out.println("logTreeLikelihood = " + logP);
        assertEquals(-2009.8343141884977, logP, 1e-6);
    }

    @Test
    public void testF3X4(){
        CodonFrequencies codonFreq = new CodonFrequencies();
        codonFreq.initByName("pi", "F3X4", "data", codonAlignment);

        double logP = getLogTreeLikelihood("0.08000", "15.34858", codonFreq, codonAlignment);
        System.out.println("logTreeLikelihood = " + logP);
        assertEquals(-1960.7661713238033, logP, 1e-6);
    }

    @Test
    public void testF60(){
        CodonFrequencies codonFreq = new CodonFrequencies();
        codonFreq.initByName("pi", "F60/F61", "data", codonAlignment);

        double logP = getLogTreeLikelihood("0.08327", "15.54039", codonFreq, codonAlignment);
        System.out.println("logTreeLikelihood = " + logP);
        assertEquals(-1915.2251567137826, logP, 1e-6);
    }

    // get logP
    protected double getLogTreeLikelihood(String omegaValue, String kappaValue,
                                          CodonFrequencies codonFreq, Alignment data) {
        RealParameter omega = new RealParameter(omegaValue);
        RealParameter kappa = new RealParameter(kappaValue);

        M0Model m0 = new M0Model();
        m0.initByName("omega", omega, "kappa", kappa, "frequencies", codonFreq);

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("substModel", m0);

        Tree tree = getStartingTree(data);

        System.setProperty("java.only","true");
        TreeLikelihood likelihood = new TreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel);
        return likelihood.calculateLogP();
    }

    public static Alignment getAlignment() {
        Sequence s1 = new Sequence("Human_Horai", "CTACCCGCCGCAGTACTGATCATTCTATTTCCCCCTCTATTGATCCCCACCTCCAAATATCTCATCAACAACCGACTAATTACCACCCAACAATGACTAATCAAACTAACCTCAAAACAAATGATAGCCATACACAACACTAAAGGACGAACCTGATCTCTTATACTAGTATCCTTAATCATTTTTATTGCCACAACTAACCTCCTCGGACTCCTGCCTCACTCATTTACACCAACCACCCAACTATCTATAAACCTAGCCATGGCCATCCCCTTATGAGCGGGCGCAGTGATTATAGGCTTTCGCTCTAAGATTAAAAATGCCCTAGCCCACTTCTTACCACAAGGCACACCTACACCCCTTATCCCCATACTAGTTATTATCGAAACCATCAGCCTACTCATTCAACCAATAGCCCTGGCCGTACGCCTAACCGCTAACATTACTGCAGGCCACCTACTCATGCATCTAATTGGAAGCGCCACCCTAGCAATATCAACCATTAACCTTCCCTCTACACTTATCATCTTCACAATTCTAATTCTACTGACTATCCTAGAAATCGCTGTCGCCTTAATCCAAGCCTACGTTTTCACACTTCTAGTAAGCCTCTACCTGCACGACAACACACCCCAACTAAATACTACCGTATGGCCCACCATAATTACCCCCATACTCCTTACACTATTTCTCATCACCCAACTAAAAATATTAAACACAAGCTACCACTTACCTCCCTCACCAAAGCCCATAAAAATAAAAAATTATAACAAACCCTGAGAACCAACCCCAATACGCAAAATTAACCCCCTAATAAAATTAATTAACCACTCATTCATCGACCTCCCCACCCCATCCAACATCTCCGCATGATGAAACTTCGGCTCACTCCTTGGCGCCTGCCTGATCCTCCAAATCACCACAGGACTATTCCTAGCCATGCACTACTCACCAGACGCCTCAACCGCCTTTTCATCAATCGCCCACATCACTCGAGAC");
        Sequence s2 = new Sequence("Human_Arnason", "CTACCCGCCGCAGTACTGATCATTCTATTTCCCCCTCTATTGATCCCCACCTCCAAATATCTCATCAACAACCGACTAATCACCACCCAACAATGACTAATCAAACTAACCTCAAAACAAATGATAACCATACACAACACTAAAGGACGAACCTGATCTCTTATACTAGTATCCTTAATCATTTTTATTGCCACAACTAACCTCCTCGGACTCCTGCCTCACTCATTTACACCAACCACCCAACTATCTATAAACCTAGCCATGGCCATCCCCTTATGAGCGGGCGCAGTGATTATAGGCTTTCGCTCTAAGATTAAAAATGCCCTAGCCCACTTCTTACCACAAGGCACACCTACACCCCTTATCCCCATACTAGTTATTATCGAAACCATCAGCCTACTCATTCAACCAATAGCCCTGGCCGTACGCCTAACCGCTAACATTACTGCAGGCCACCTACTCATGCACCTAATTGGAAGCGCCACCCTAGCAATATCAACCATTAACCTTCCCTCTACACTTATCATCTTCACAATTCTAATTCTACTGACTATCCTAGAAATCGCTGTCGCCTTAATCCAAGCCTACGTTTTCACACTTCTAGTAAGCCTCTACCTGCACGACAACACACCCCAACTAAATACTACCGTATGGCCCACCATAATTACCCCCATACTCCTTACACTATTCCTCATCACCCAACTAAAAATATTAAACACAAACTACCACCTACCTCCCTCACCAAAGCCCATAAAAATAAAAAATTATAACAAACCCTGAGAACCAACCCCAATACGCAAAATTAACCCCCTAATAAAATTAATTAACCGCTCATTCATCGACCTCCCCACCCCATCCAACATCTCCGCATGATGAAACTTCGGCTCACTCCTTGGCACCTGCCTGATCCTCCAAATCACCACAGGACTATTCCTAGCCATGCACTACTCACCAGACGCCTCAACCGCCTTTTCATCAATCGCCCACATCACTCGAGAC");
        Sequence s3 = new Sequence("Chimp_Horai", "TTACCCGCCGCAGTACTAATCATTCTATTCCCCCCTCTACTGGTCCCCACTTCTAAACATCTCATCAACAACCGACTAATTACCACCCAACAATGACTAATTCAACTGACCTCAAAACAAATAATAACTATACACAGCACTAAAGGACGAACCTGATCTCTCATACTAGTATCCTTAATCATTTTTATTACCACAACCAATCTTCTTGGGCTTCTACCCCACTCATTCACACCAACCACCCAACTATCTATAAACCTAGCCATGGCTATCCCCCTATGAGCAGGCGCAGTAGTCATAGGCTTTCGCTTTAAGACTAAAAATGCCCTAGCCCACTTCTTACCGCAAGGCACACCTACACCCCTTATCCCCATACTAGTTATCATCGAAACTATTAGCCTACTCATTCAACCAATAGCCTTAGCCGTACGTCTAACCGCTAACATTACTGCAGGCCACCTACTCATGCACCTAATTGGAAGCGCCACACTAGCATTATCAACTATCAATCTACCCTATGCACTCATTATCTTCACAATTCTAATCCTACTGACTATTCTAGAGATCGCCGTCGCCTTAATCCAAGCCTACGTTTTTACACTTCTAGTGAGCCTCTACCTGCACGACAACACACCCCAACTAAATACCGCCGTATGACCCACCATAATTACCCCCATACTCCTGACACTATTTCTCGTCACCCAACTAAAAATATTAAATTCAAATTACCATCTACCCCCCTCACCAAAACCCATAAAAATAAAAAACTACAATAAACCCTGAGAACCAACCCCGACACGCAAAATTAACCCACTAATAAAATTAATTAATCACTCATTTATCGACCTCCCCACCCCATCCAACATTTCCGCATGATGGAACTTCGGCTCACTTCTCGGCGCCTGCCTAATCCTTCAAATTACCACAGGATTATTCCTAGCTATACACTACTCACCAGACGCCTCAACCGCCTTCTCGTCGATCGCCCACATCACCCGAGAC");
        Sequence s4 = new Sequence("Chimp_Arnason", "TTACCCGCCGCAGTACTAATCATTCTATTCCCCCCTCTACTGGTCCCCACTTCTAAACATCTCATCAACAACCGACTAATTACCACCCAACAATGACTAATTCAACTGACCTCAAAACAAATAATAACTATACACAGCACTAAAGGACGAACCTGATCTCTCATACTAGTATCCTTAATCATTTTTATTACCACAACCAATCTTCTTGGACTTCTACCCCACTCATTCACACCAACCACCCAACTATCTATAAACCTAGCCATGGCTATCCCCCTATGAGCAGGCGCAGTAGTCATAGGCTTTCGCTTTAAGACTAAAAATGCCCTAGCCCACTTCTTACCGCAAGGCACACCTACACCCCTTATCCCCATACTAGTTATCATCGAAACTATTAGCCTACTCATTCAACCAATAGCCTTAGCCGTACGTCTAACCGCTAACATTACTGCAGGCCACCTACTCATGCACCTAATTGGAAGCGCCACACTAGCATTATCAACTATCAATCTACCCTCTGCACTCATTATCTTCACAATTCTAATCCTACTGACTATTCTAGAGATCGCCGTCGCCTTAATCCAAGCCTACGTTTTTACACTTCTAGTGAGCCTCTACCTGCACGACAACACACCCCAACTAAATACCGCCGTATGACCCACCATAATTACCCCCATACTCCTGACACTATTTCTCGTCACCCAACTAAAAATATTAAATTCAAATTACCATCTACCCCCCTCACCAAAACCCATAAAAATAAAAAACTACAATAAACCCTGAGAACCAACCCCGACACGCAAAATTAACCCACTAATAAAATTAATTAATCACTCATTTATCGACCTCCCCACCCCATCCAACATTTCCGCATGATGGAACTTCGGCTCACTTCTCGGCGCCTGCCTAATCCTTCAAATTACCACAGGATTATTCCTAGCTATACACTACTCACCAGACGCCTCAACCGCCTTCTCGTCGATCGCCCACATCACCCGAGAC");
        Sequence s5 = new Sequence("Gorilla_Horai", "TTACCCGCCGCAGTATTAATTATCCTACTTCCCCCTCTACTGATCCCCACCTCCAAATATCTCATCAACAACCGACTGATTGCCACCCAACAGTGACTAATCCAACTAACCTCAAAACAAATAATAACTATACATAACGCCAAGGGACGAACCTGATCCCTTATGTTAATATGATTAATTATTTTTATTGCCACAACCAACCTCCTCGGACTCTTGCCCCACTCATTCACACCAACTACCCAGCTATCTATAAACCTGGCCATAGCCATCCCCCTGTGAGCAGGCGCAGTAACTACAGGCTTTCGCTCCAAGACTAAAAATGCCCTAGCCCACCTACTACCACAAGGCACCCCTACACCCCTTATCCCTATACTAGTCATCATTGAAACCATCAGCCTATTCATCCAACCAATAGCCCTAGCTGTACGCCTAACCGCTAACATCACTGCAGGTCACCTACTTATGCACCTAATCGGAAGCGCCACACTAGCAATATCAACTACCAATCTTCCCTCAACACTCATTATCTTTACAGTCCTAATTTTATTAACTATACTAGAAATCGCTGTCGCCCTCATCCAAGCCTACGTTTTCACACTTTTAGTGAGCCTCTACCTGCACGAGAACACACCCCAGCTAAATACCACCGTATGGCCCACCATAATTGCCCCAATACTCCTCACACTATTTCTCATTACCCAACTAAAAGTTTTAAACACAAATTACCACCTACCCCCCTTACCAAAAACTATAAAAATAAAAAACTTCTGTAAACCCTGAGAACCAACCCCTATACGCAAAACTAACCCACTAGCAAAACTAATTAACCACTCATTCATTGACCTCCCTACCCCGTCCAACATCTCCACATGATGAAACTTCGGCTCACTCCTTGGTGCCTGCTTAATCCTTCAAATCACCACAGGGCTATTCCTAGCCATACACTACTCACCTGATGCCTCAACCGCCTTCTCATCAATTGCCCACATCACCCGAGAT");
        Sequence s6 = new Sequence("Gorilla_Arnason", "TTACCCGCCGCAGTATTAATTATCCTACTTCCCCCTCTACTGATCCCCACCTCCAAATATCTCATCAACAACCGACTGATTGCCACCCAACAGTGACTAATCCAACTAACCTCAAAACAAATAATAACTATACATAACGCCAAGGGACGAACCTGATCCCTTATGTTAATATCATTAATTATTTTTATTGCCACAACCAACCTCCTCGGACTCTTGCCCCACTCATTCACACCAACTACCCAGCTATCTATAAACCTGGCCATAGCCATCCCCCTGTGAGCAGGCGCAGTAACTACAGGCTTTCGCTCCAAGACTAAAAATGCCCTAGCCCACCTACTACCACAAGGCACCCCTACACCCCTTATCCCTATACTAGTCATCATTGAAACCATCAGCCTATTCATCCAACCAATAGCCCTAGCTGTACGCCTAACCGCTAACATCACTGCAGGTCACCTACTTATGCACCTAATCGGAAGCGCCACACTAGCAATATCAACTACCAATCTTCCCTCAACACTCATTATCTTTACAGTCCTAATTTTATTAACTATACTAGAAATCGCTGTCGCCCTCATCCAAGCCTACGTTTTCACACTTTTAGTGAGCCTCTACCTGCACGAGAACACACCCCAGCTAAATACCACCGTATGGCCCACCATAATTGCCCCAATACTCCTCACACTATTTCTCATTACCCAACTAAAAGTTTTAAACACAAATTACCACCTACCCCCCTTACCAAAAACTATAAAAATAAAAAACTTCTGTAAACCCTGAGAATCAACCCCTATACGCAAAACTAACCCACTAGCAAAACTAATTAACCACTCATTCATTGACCTCCCTACCCCGTCCAACATCTCCACATGATGAAACTTCGGCTCACTCCTTGGTGCCTGCTTAATCCTTCAAATCACCACAGGGCTATTCCTAGCCATACACTACTCACCTGATGCCTCAACCGCCTTCTCATCAATCGCCCACATCACCCGAGAT");

        Alignment data = new Alignment();
        data.initByName("sequence", s1, "sequence", s2, "sequence", s3, "sequence", s4, "sequence", s5, "sequence", s6,
                "dataType", "nucleotide"
        );
        return data;
    }

    public static Tree getStartingTree(Alignment data) {
        String startingTree = "((Human_Horai: 0.012988, Human_Arnason: 0.012988): " +
                "0.245731, (Chimp_Horai: 0.002959, Chimp_Arnason: 0.002959): 0.255761, " +
                "(Gorilla_Horai: 0.003918, Gorilla_Arnason: 0.003918): 0.254801);";

        TreeParser t = new TreeParser();
        t.initByName("taxa", data, "newick", startingTree, "IsLabelledNewick", true);
        return t;
    }
}
