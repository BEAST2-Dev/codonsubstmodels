package codonmodels.evolution.likelihood;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;
import codonmodels.CodonFrequencies;
import codonmodels.evolution.CodonData;
import codonmodels.evolution.alignment.CodonAlignment;
import codonmodels.evolution.tree.NodeStatesArray;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.assertEquals;


/**
 * test Multithreading DataAugTreeLikelihood.
 *
 * @author Walter Xie
 */
public class DALikelihoodMultithreadingTest {

    private int[] threads = new int[]{1,2,4};
    DataAugTreeLikelihood[] testLds;

//    DecimalFormat df = new DecimalFormat("#.00");

    @Before
    public void setUp() {
        Alignment data = CodonData.getAlig6T333();
        // create Codon Alignment
        CodonAlignment codonAlignment = new CodonAlignment();
        codonAlignment.initByName("data", data, "dataType", "codon",
                "geneticCode", "vertebrateMitochondrial", "verbose", false);

        CodonFrequencies codonFreq = new CodonFrequencies();
        codonFreq.initByName("pi", "F3X4", "data", codonAlignment, "verbose", false);
        SiteModel siteModel = CodonData.getSiteModel("0.3", "5", codonFreq, false);

        String newickTree = "(((Human_Horai: 1, Human_Arnason: 1.1): 1, " +
                "(Chimp_Horai: 1, Chimp_Arnason: 1.1): 0.8): 0.5, " +
                "(Gorilla_Horai: 1, Gorilla_Arnason: 1.1): 1.5);";
        Tree tree = CodonData.getTree(codonAlignment, newickTree, false, true);

        System.setProperty("java.only","true");


        // =============== DA tree likelihood ===============
        Randomizer.setSeed(777);
        testLds = new DataAugTreeLikelihood[threads.length];
        // set internal node states on the fly
        NodeStatesArray nodeStatesArray = new NodeStatesArray();
        nodeStatesArray.initByName("data", codonAlignment, "initINS", "parsimony");

        // different threads
        for (int i = 0; i < threads.length; i++) {
//            NodeStatesArray newNSA = (NodeStatesArray) nodeStatesArray.copy();
            int th = threads[i];
            testLds[i] = new DataAugTreeLikelihood();
            testLds[i].initByName("nodesStates", nodeStatesArray, "tree", tree,
                    "siteModel", siteModel, "threads", th);
        }

        NodeStatesArray nodesStates = testLds[0].getNodesStates();
        // test internal node states are same
        for (int i = 1; i < threads.length; i++) {
            NodeStatesArray nodesStates2 = testLds[i].getNodesStates();
            if (!nodesStates2.hasSameInternalNodeStates(nodesStates))
                throw new IllegalArgumentException("Internal node states cannot be different !");
        }
        System.out.println("\nSetup the same internal node states.");
    }

    @After
    public void close() {
        for (int i = 0; i < threads.length; i++)
            testLds[i].shutdown(); // shut down ExecutorService
    }

    @Test
    public void testMultithreading(){
        double[] logPs = new double[threads.length];
        long[] timeDAs = new long[threads.length];
        System.out.println("\n=============== DA tree likelihood multithreading test ===============\n");
//int i=1; // 2 threads
        long start = System.nanoTime();
        for (int i = 0; i < threads.length; i++){
            logPs[i] = testLds[i].calculateLogP();
            timeDAs[i] = System.nanoTime() - start;
            System.out.println("Thread " + threads[i] + " DA tree likelihood logP = " + logPs[i] +
                    ", time = " + timeDAs[i] + " nanoseconds.");
        }

        assertEquals(logPs[0], logPs[1], 1e-10);
        assertEquals(logPs[0], logPs[2], 1e-10);

        // not true for single core
//        assertTrue("Invalid test for single core", timeDAs[0] > timeDAs[2]);
    }



}
