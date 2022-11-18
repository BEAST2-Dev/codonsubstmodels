package codonmodels.evolution;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import beast.pkgmgmt.BEASTClassLoader;
import codonmodels.CodonFrequencies;
import codonmodels.M0Model;
import codonmodels.evolution.datatype.Codon;

import java.io.File;

/**
 * @author Walter Xie
 */
public class CodonData {

    public static final Codon codonUNIVERSAL = new Codon();

    public static void initDataTypes() {
        String PROJECT_DIR = System.getProperty("user.dir");
        System.out.println("user.dir = " + PROJECT_DIR);
        File versionXML = new File(PROJECT_DIR + File.separator + "version.xml");
        if (!versionXML.exists())
            throw new IllegalArgumentException("Cannot locate codonsubstmodels/version.xml !");
        BEASTClassLoader.initServices();
        BEASTClassLoader.addServices(versionXML.getAbsolutePath());
    }


    // create a Sequence obj using Codon.codeMap
    public static Sequence getSeqCodeMap() {
        return new Sequence("codemap", codonUNIVERSAL.getCodeMap());
    }

    public static int[] getIntegers() {
        return new int[]{0, 20, 40, 59, 60, 63, 64, 121};
    }


    public static Alignment getAlig6T333() {
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

    /**
     * Create BEAST tree from alignment and newick string
     */
    public static Tree getTree(Alignment data, String newickTree, boolean adjustTipHeights, boolean verbose) {
        if (verbose)
            System.out.println("Tree is " + newickTree + "\n");
        TreeParser t = new TreeParser();
        t.initByName("taxa", data, "newick", newickTree, "IsLabelledNewick", true, "adjustTipHeights", adjustTipHeights);
        return t;
    }

    public static SiteModel getSiteModel(String omegaValue, String kappaValue, CodonFrequencies codonFreq, boolean verbose) {
        RealParameter omega = new RealParameter(omegaValue);
        RealParameter kappa = new RealParameter(kappaValue);

        M0Model m0 = new M0Model();
        m0.initByName("omega", omega, "kappa", kappa, "frequencies", codonFreq, "verbose", verbose);

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("substModel", m0);

        return siteModel;
    }

    // M0 get logP, 0.3, 5
    public static TreeLikelihood getTreeLikelihoodM0(String omegaValue, String kappaValue, CodonFrequencies codonFreq,
                                                     Alignment data, String newickTree, boolean adjustTipHeights,
                                                     boolean verbose) {
        SiteModel siteModel = CodonData.getSiteModel(omegaValue, kappaValue, codonFreq, verbose);
        Tree tree = CodonData.getTree(data, newickTree, adjustTipHeights, verbose);

        System.setProperty("java.only","true");
        TreeLikelihood likelihood = new TreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel);
        return likelihood;
    }

}