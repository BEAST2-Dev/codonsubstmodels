package test.beast.evolution.likelihood;

import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.CodonAlignment;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Tree;
import beast.likelihood.DACodonTreeLikelihood;
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
        codonAlignment.initByName("data", data, "dataType", "codon", "geneticCode", "vertebrateMitochondrial");

        newickTree = "((Human_Horai: 0.012988, Human_Arnason: 0.012988): 0.245731, " +
                "(Chimp_Horai: 0.002959, Chimp_Arnason: 0.002959): 0.255761, (Gorilla_Horai: 0.003918, " +
                "Gorilla_Arnason: 0.003918): 0.254801);";


    }

    @Test
    public void testCodonPositionBaseFrequencies(){
    }

    @Test
    public void testCodonFrequencies(){
    }

    @Test
    public void testEqualFreq(){
    }

    @Test
    public void testF1X4(){
    }

    @Test
    public void benchmarkingF3X4(){
        CodonFrequencies codonFreq = new CodonFrequencies();
        codonFreq.initByName("pi", "F3X4", "data", codonAlignment);

        SiteModel siteModel = CodonTestData.getSiteModel("0.08000", "15.34858", codonFreq);

        Tree tree = CodonTestData.getTree(codonAlignment, newickTree);

        System.setProperty("java.only","true");

        // =============== Standard likelihood ===============
        // Get current time
        long start = System.currentTimeMillis();
        int iteration = 100;

        for (int i=0; i<iteration; i++) {
            TreeLikelihood likelihood = new TreeLikelihood();
            likelihood.initByName("data", codonAlignment, "tree", tree, "siteModel", siteModel);

            double logP = likelihood.calculateLogP();
        }

        // Get elapsed time in milliseconds
        long elapsedTimeMillis = System.currentTimeMillis()-start;
        double time = (double) elapsedTimeMillis/iteration;

        System.out.println("\n=============== Standard likelihood ===============\n");
        System.out.println(iteration + " iteration(s) " + elapsedTimeMillis + " milliseconds");
        System.out.println("Calculation time = " + time + " milliseconds per iteration in average.");

        // =============== DA likelihood ===============

        start = System.currentTimeMillis();

        for (int i=0; i<iteration; i++) {
            DACodonTreeLikelihood likelihood = new DACodonTreeLikelihood();
            likelihood.initByName("data", codonAlignment, "tree", tree, "siteModel", siteModel);

            double logP = likelihood.calculateLogP();
        }

        // Get elapsed time in milliseconds
        elapsedTimeMillis = System.currentTimeMillis()-start;
        time = (double) elapsedTimeMillis/iteration;

        System.out.println("\n=============== DA likelihood ===============\n");
        System.out.println(iteration + " iteration(s) " + elapsedTimeMillis + " milliseconds");
        System.out.println("Calculation time = " + time + " milliseconds per iteration in average.");

    }

    @Test
    public void testF60(){
    }

}
