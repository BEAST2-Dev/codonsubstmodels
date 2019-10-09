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


/**
 * @author Walter Xie
 */
public class DACodonTreeLikelihoodTest {
    CodonAlignment codonAlignment;
    String newickTree;

    @Before
    public void setUp() {
        Alignment data = CodonTestData.getAlignment();
        // create Codon Alignment
        codonAlignment = new CodonAlignment();
        codonAlignment.initByName("data", data, "dataType", "codon",
                "geneticCode", "vertebrateMitochondrial", "verbose", false);

        newickTree = "((Human_Horai: 0.012988, Human_Arnason: 0.012988): 0.245731, " +
                "(Chimp_Horai: 0.002959, Chimp_Arnason: 0.002959): 0.255761, (Gorilla_Horai: 0.003918, " +
                "Gorilla_Arnason: 0.003918): 0.254801);";
    }


    @Test
    public void benchmarkingF3X4(){
        CodonFrequencies codonFreq = new CodonFrequencies();
        codonFreq.initByName("pi", "F3X4", "data", codonAlignment, "verbose", false);

        SiteModel siteModel = CodonTestData.getSiteModel("0.08000", "15.34858", codonFreq, false);

        Tree tree = CodonTestData.getTree(codonAlignment, newickTree);

        System.setProperty("java.only","true");

        // =============== Standard likelihood ===============
        int iteration = 1;
        long[] elapsedTimeMillis = new long[iteration];
        double timeStandard = 0;

        // Get current time
        long start = System.currentTimeMillis();

        TreeLikelihood likelihood = new TreeLikelihood();
        likelihood.initByName("data", codonAlignment, "tree", tree, "siteModel", siteModel);

        long standardInit = System.currentTimeMillis()-start;

        for (int i=0; i<iteration; i++) {
            start = System.currentTimeMillis();

            double logP = likelihood.calculateLogP();

            // Get elapsed time in milliseconds
            elapsedTimeMillis[i] = System.currentTimeMillis()-start;
            timeStandard += elapsedTimeMillis[i];

            System.out.println("i = " + i + " logP = " + logP);
        }

        // =============== DA likelihood ===============

        // Get current time
        start = System.currentTimeMillis();

        int tipCount = tree.getLeafNodeCount();
        int internalNodeCount = tree.getInternalNodeCount();
        int siteCount = codonAlignment.getSiteCount();

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

        elapsedTimeMillis = new long[iteration];
        double timeDA = 0;

        DACodonTreeLikelihood daLikelihood = new DACodonTreeLikelihood();

        daLikelihood.initByName("data", codonAlignment, "tree", tree, "siteModel", siteModel,
                "internalNodeStates", internalNodeStates);

        long daInit = System.currentTimeMillis()-start;


        for (int i=0; i<iteration; i++) {
            start = System.currentTimeMillis();

            double logPDA = daLikelihood.calculateLogP();

            // Get elapsed time in milliseconds
            elapsedTimeMillis[i] = System.currentTimeMillis()-start;
            timeDA += elapsedTimeMillis[i];

            System.out.println("i = " + i + " logPDA = " + logPDA);
        }

        // =============== report ===============

        System.out.println("\n=============== Standard likelihood ===============\n");
        System.out.println("Init time " + standardInit + " milliseconds");
        System.out.println(iteration + " iteration(s) " + timeStandard + " milliseconds");
        timeStandard /= iteration;
        System.out.println("Calculation time = " + timeStandard + " milliseconds per iteration in average.");


        System.out.println("\n=============== DA likelihood ===============\n");
        System.out.println("Init time " + daInit + " milliseconds");
        System.out.println(iteration + " iteration(s) " + timeDA + " milliseconds");
        timeDA /= iteration;
        System.out.println("Calculation time = " + timeDA + " milliseconds per iteration in average.");

    }


}
