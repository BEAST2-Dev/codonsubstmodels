package test.beast.evolution.likelihood;

import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.CodonAlignment;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Tree;
import beast.likelihood.DACodonTreeLikelihood;
import beast.tree.InternalNodeStates;
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
        Alignment data = CodonTestData.getAlignment();
        // create Codon Alignment
        codonAlignment = new CodonAlignment();
        codonAlignment.initByName("data", data, "dataType", "codon",
                "geneticCode", "vertebrateMitochondrial", "verbose", false);

    }

    private void init(String pi) {
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

    @Test
    public void benchmarkingF3X4(){
        init("F3X4");

        // =============== Standard likelihood ===============
        int iteration = 10;
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

//            System.out.println("i = " + i + " logP = " + logP);
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

        long[] elapsedTimeMillis2 = new long[iteration];

        DACodonTreeLikelihood daLikelihood = new DACodonTreeLikelihood();

        daLikelihood.initByName("data", codonAlignment, "tree", tree, "siteModel", siteModel,
                "internalNodeStates", internalNodeStates);

        long daInit = System.currentTimeMillis()-start;

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


}
