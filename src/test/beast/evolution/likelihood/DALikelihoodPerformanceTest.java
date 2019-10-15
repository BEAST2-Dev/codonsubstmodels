package test.beast.evolution.likelihood;

import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.CodonAlignment;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Tree;
import beast.likelihood.DACodonTreeLikelihood;
import beast.tree.InternalNodeStates;
import beast.util.XMLParserException;
import codonmodels.CodonFrequencies;
import org.junit.Before;
import org.junit.Test;
import test.beast.evolution.CodonTestData;

import java.text.DecimalFormat;
import java.util.Arrays;


/**
 * Performance tests between standard TreeLikelihood and
 * DACodonTreeLikelihood.
 *
 * @author Walter Xie
 */
public class DALikelihoodPerformanceTest {
    CodonAlignment codonAlignment;
    SiteModel siteModel;
    Tree tree;

    @Before
    public void setUp() {
//        init6T333("F3X4");
        init16T10K("F3X4");
//        init32T10K("F3X4");
    }

    private void init6T333(String pi) {
        Alignment data = CodonTestData.getAlig6T333();
        // create Codon Alignment
        codonAlignment = new CodonAlignment();
        codonAlignment.initByName("data", data, "dataType", "codon",
                "geneticCode", "vertebrateMitochondrial", "verbose", false);

        if (pi=="F3X4") {
            String newickTree = "(((Human_Horai: 0.1, Human_Arnason: 0.1): 0.15, (Chimp_Horai: 0.15, " +
                    "Chimp_Arnason: 0.15): 0.1): 0.25, (Gorilla_Horai: 0.2, Gorilla_Arnason: 0.2): 0.3);";

            CodonFrequencies codonFreq = new CodonFrequencies();
            codonFreq.initByName("pi", pi, "data", codonAlignment, "verbose", false);
            siteModel = CodonTestData.getSiteModel("0.08000", "15.34858", codonFreq, false);
            tree = CodonTestData.getTree(codonAlignment, newickTree);

            System.setProperty("java.only","true");
        } else {
            throw new IllegalArgumentException("Invalid pi " + pi);
        }
    }

    private void init16T10K(String pi) {
        Alignment data = null;
        try {
            data = CodonTestData.getAlig16T10K();
        } catch (XMLParserException e) {
            e.printStackTrace();
        }
        // create Codon Alignment
        codonAlignment = new CodonAlignment();
        codonAlignment.initByName("data", data, "dataType", "codon",
                "geneticCode", "vertebrateMitochondrial", "verbose", false);

        if (pi=="F3X4") {
            String newickTree = "((t6:0.8328606829,t5:0.07987240003):0.8731191335,(((((t10:0.4489309546,t2:0.5118521685):0.6122735683,t12:0.3223016059):0.7295662321,(t15:0.6618216026,(t4:0.6256095201,t1:0.7471578235):0.8436422099):0.3587985209):0.07577025215,(((t11:0.375966636,(t14:0.01772094774,t9:0.5878796382):0.5867906068):0.2359019916,t8:0.2657493246):0.675317778,t16:0.8670396307):0.3450425724):0.1524739896,(t3:0.4887036914,(t13:0.4686036934,t7:0.1729974856):0.1101265033):0.1990886456):0.8794986582);";

            CodonFrequencies codonFreq = new CodonFrequencies();
            codonFreq.initByName("pi", pi, "data", codonAlignment, "verbose", false);
            siteModel = CodonTestData.getSiteModel("0.3", "5", codonFreq, false);
            tree = CodonTestData.getTree(codonAlignment, newickTree);

            System.setProperty("java.only","true");
        } else {
            throw new IllegalArgumentException("Invalid pi " + pi);
        }
    }

    private void init32T10K(String pi) {
        Alignment data = null;
        try {
            data = CodonTestData.getAlig32T10K();
        } catch (XMLParserException e) {
            e.printStackTrace();
        }
        // create Codon Alignment
        codonAlignment = new CodonAlignment();
        codonAlignment.initByName("data", data, "dataType", "codon",
                "geneticCode", "vertebrateMitochondrial", "verbose", false);

        if (pi=="F3X4") {
            String newickTree = "(((t15:0.9542181483,t19:0.3497522059):0.3365416529,t2:0.5435624667):0.9214519702,((((((t21:0.7185004468,t20:0.5407348312):0.7264935684,(t17:0.5117207766,(t8:0.3790637162,(t30:0.7958646573,(t14:0.4777974519,t18:0.9716053747):0.7753737369):0.6706928122):0.6051157189):0.890223223):0.1895852804,(t22:0.4900219755,t13:0.05010107416):0.04109163722):0.06207286357,((t6:0.5607458302,t16:0.9570369916):0.4422796173,t31:0.8856992035):0.4335002408):0.7788049644,(((t7:0.7984650929,t4:0.9124820076):0.224632408,(t1:0.544289869,t32:0.702123753):0.5266149947):0.8062962941,((t24:0.9591187276,t27:0.7903348505):0.9343056851,t10:0.924266933):0.8483990482):0.8787279353):0.4969263494,(((t5:0.596928519,t23:0.6392580203):0.3214438602,t9:0.1940213139):0.2096933059,((t11:0.8579683306,(t28:0.9736292204,t26:0.06027402054):0.07778765121):0.1267342451,((t29:0.8047343672,t25:0.236599003):0.0191930905,(t12:0.9822205401,t3:0.1412437484):0.1372304766):0.6479954512):0.6275273843):0.5978206636):0.5197886776);";

            CodonFrequencies codonFreq = new CodonFrequencies();
            codonFreq.initByName("pi", pi, "data", codonAlignment, "verbose", false);
            siteModel = CodonTestData.getSiteModel("0.3", "5", codonFreq, false);
            tree = CodonTestData.getTree(codonAlignment, newickTree);

            System.setProperty("java.only","true");
        } else {
            throw new IllegalArgumentException("Invalid pi " + pi);
        }
    }

    @Test
    public void benchmarkingF3X4(){

        // =============== Standard likelihood ===============
        int iteration = 1;
        long[] elapsedTimeMillis = new long[iteration];

        // Get current time
        long start = System.currentTimeMillis();

        TreeLikelihood likelihood = new TreeLikelihood();
        likelihood.initByName("data", codonAlignment, "tree", tree, "siteModel", siteModel);

        long standardInit = System.currentTimeMillis()-start;

        double timeStandard = 0;
        for (int i=0; i<iteration; i++) {
            start = System.currentTimeMillis();

            double logP = likelihood.calculateLogP();

            // Get elapsed time in milliseconds
            elapsedTimeMillis[i] = System.currentTimeMillis()-start;
            timeStandard += elapsedTimeMillis[i];

            System.out.println("i = " + i + " logP = " + logP);
        }
        double timeStandardStdev = 0;
        for (int i=0; i<iteration; i++) {
            timeStandardStdev += (elapsedTimeMillis[i]-timeStandard/iteration)*(elapsedTimeMillis[i]-timeStandard/iteration);
        }
        timeStandardStdev = Math.sqrt(timeStandardStdev / iteration);

        // =============== DA likelihood ===============

        // Get current time
        start = System.currentTimeMillis();

        int tipCount = tree.getLeafNodeCount();
        int internalNodeCount = tree.getInternalNodeCount();
        int siteCount = codonAlignment.getSiteCount();
        // vert MT stop codon states = 8, 10, 48, 50
        int[] stopCodons = codonAlignment.getGeneticCode().getStopCodonStates();
        System.out.println("Stop codon states " + Arrays.toString(stopCodons));

        InternalNodeStates internalNodeStates = new InternalNodeStates(internalNodeCount, siteCount);

        // internal nodes
        for (int i=tipCount; i<tree.getNodeCount(); i++) {
            int[] states = new int[siteCount];
            // 0 - 63
            for (int j=0; j < states.length; j++) {
                states[j] = (int)(Math.random() * 64);
                // not stop codon in vertebrateMitochondrial
                while(states[j] == 8 || states[j] == 10 || states[j] == 48 || states[j] == 50)
                    states[j] = (int)(Math.random() * 64);
            }
            internalNodeStates.setNrStates(i, states);
        }

        DACodonTreeLikelihood daLikelihood = new DACodonTreeLikelihood();

        daLikelihood.initByName("data", codonAlignment, "tree", tree, "siteModel", siteModel,
                "internalNodeStates", internalNodeStates);

        long daInit = System.currentTimeMillis()-start;

        long[] elapsedTimeMillis2 = new long[iteration];
        double timeDA = 0;
        for (int i=0; i<iteration; i++) {
            start = System.currentTimeMillis();

            double logPDA = daLikelihood.calculateLogP();

            // Get elapsed time in milliseconds
            elapsedTimeMillis2[i] = System.currentTimeMillis()-start;
            timeDA += elapsedTimeMillis2[i];

            System.out.println("i = " + i + " DA logP = " + logPDA);
        }
        double timeDAStdev = 0;
        for (int i=0; i<iteration; i++) {
            timeDAStdev += (elapsedTimeMillis2[i]-timeDA/iteration)*(elapsedTimeMillis2[i]-timeDA/iteration);
        }
        timeDAStdev = Math.sqrt(timeDAStdev / iteration);

        // =============== report ===============

        DecimalFormat df = new DecimalFormat("#.00");

        System.out.println("\n=============== Standard likelihood ===============\n");
        System.out.println("Init time " + standardInit + " milliseconds");
        System.out.println(iteration + " iteration(s) " + timeStandard + " milliseconds");
        timeStandard /= iteration;
        System.out.println("Calculation time = " + timeStandard + " +- " + df.format(timeStandardStdev) +
                " milliseconds per iteration in average.");


        System.out.println("\n=============== DA likelihood ===============\n");
        System.out.println("Init time " + daInit + " milliseconds");
        System.out.println(iteration + " iteration(s) " + timeDA + " milliseconds");
        timeDA /= iteration;
        System.out.println("Calculation time = " + timeDA + " +- " + df.format(timeDAStdev) +
                " milliseconds per iteration in average.");

        System.out.println("\n" + df.format(timeStandard/timeDA) + " times faster ");
    }


    @Test
    public void benchmarkingForLoop(){

        double[][] m = new double[64][64];
        double sum = 0;

        long start = System.currentTimeMillis();
        for (int n=0; n<1000000; n++) {
            for (int i = 0; i < m.length; i++) {
                for (int j = 0; j < m.length; j++) {
                    sum += m[i][j];
                }
            }
        }
        long end1 = System.currentTimeMillis()-start;
        System.out.println("Sum 64*64 time is " + end1 + " milliseconds");

        start = System.currentTimeMillis();
        for (int n=0; n<1000000; n++) {
//        for (int i=0; i<m.length; i++) {
            for (int j=0; j<m.length; j++) {
            sum += m[16][j];
            }
//        }
        }
        long end2 = System.currentTimeMillis()-start;
        System.out.println("Sum 64 time is " + end2 + " milliseconds");


        start = System.currentTimeMillis();
        for (int n=0; n<1000000; n++) {
//        for (int i=0; i<m.length; i++) {
//            for (int j=0; j<m.length; j++) {
            sum += m[16][63];
//            }
//        }
        }
        long end3 = System.currentTimeMillis()-start;
        System.out.println("Take 1 element " + end3 + " milliseconds");


        System.out.println("\nend1 / end3 = " + (end1 / end3) + " times");
        System.out.println("end2 / end3 = " + (end2 / end3) + " times");
    }

}
