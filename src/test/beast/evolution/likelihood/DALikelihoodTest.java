package test.beast.evolution.likelihood;

import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.CodonAlignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.likelihood.DataAugTreeLikelihood;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.NodeStatesArray;
import beast.evolution.tree.Tree;
import codonmodels.CodonFrequencies;
import org.junit.Before;
import org.junit.Test;
import test.beast.evolution.CodonData;

import java.util.Arrays;

import static junit.framework.Assert.assertEquals;


/**
 * Use BEAST 2 core TreeLikelihood result to test DACodonTreeLikelihood.
 *
 * @author Walter Xie
 */
public class DALikelihoodTest {

    TreeLikelihood treeLikelihood;
    DataAugTreeLikelihood daTreeLikelihood;

//    DecimalFormat df = new DecimalFormat("#.00");

    @Before
    public void setUp() {

//        String newickTree = "(t1:0.5, t2:0.0):0.0;"; // TODO bug for 0 branch length ?
        String newickTree = "(t1:0.5, t2:0.1):0.0;";
//        boolean adjustTipHeights = true;
        boolean adjustTipHeights = false;

        // =============== Standard tree likelihood ===============
        treeLikelihood = getTreeLikelihood("equal", newickTree, adjustTipHeights);

        // =============== DA tree likelihood ===============
        // need to set internal node states
        daTreeLikelihood = getDATreeLikelihood("equal", newickTree, adjustTipHeights);

//        rootNr = daTreeLikelihood.get.getRoot().getNr();
//        System.out.println("Root index = " + rootNr);
    }

    @Test
    public void testDALikelihood1Site(){
        System.out.println("\n=============== Standard tree likelihood ===============\n");

        double logP = treeLikelihood.calculateLogP();
        System.out.println(" logP = " + logP);

        double[] partial = treeLikelihood.getRootPartials(); // root
        System.out.println("Root partial without freqs = " + Arrays.toString(partial));

        System.out.println("\n=============== DA tree likelihood ===============");

        NodeStatesArray nodesStates = daTreeLikelihood.getNodesStates();
        int rootNr = daTreeLikelihood.getRootIndex();
        System.out.println("Root index = " + rootNr);

        // assume vertebrateMitochondrial
        final double freq = 1.0/60.0;
        double daLogP = 0;
        double daLogP2 = 0; // nodeLd * freq
        // 0 - 59
        for (int s=0; s < 60; s++) {
            // internal nodes
            nodesStates.setStates(rootNr, new int[]{s});
            System.out.println("\nSet internal node state to " + s);

            double logPDA = daTreeLikelihood.calculateLogP();
            System.out.println(" DA logP = " + logPDA);

            daLogP += Math.exp(logPDA);

            // [branch][site]
            double[][] branchLd = new double[2][1];
            daTreeLikelihood.getBranchLdFromCore(0, branchLd[0]); // child 1
            daTreeLikelihood.getBranchLdFromCore(1, branchLd[1]); // child 2
            System.out.println("Branch likelihood child1 child2 without freqs = " +
                    Arrays.toString(branchLd[0]) + ", " + Arrays.toString(branchLd[1]));

            double nodeLd = branchLd[0][0] * branchLd[1][0];
            System.out.println("Root likelihood without freqs = " + nodeLd);
            System.out.println("Root partial without freqs = " + partial[s]);

            assertEquals(partial[s], nodeLd, 1e-15);

            daLogP2 += nodeLd * freq; // frequencies prior on the root
        }
        daLogP = Math.log(daLogP);
        daLogP2 = Math.log(daLogP2);
        System.out.println("\nTotal :");
        System.out.println("Total DA tree likelihood logP = " + daLogP);
        System.out.println("Total DA tree likelihood logP2 = " + daLogP2);
        System.out.println("Standard tree likelihood logP = " + logP);
        //Total DA tree likelihood logP = -15.667929146385474
        //Standard tree likelihood logP = -15.667929146369291
        assertEquals(logP, daLogP, 1e-10);
        assertEquals(logP, daLogP2, 1e-10);
    }


    private CodonAlignment getAlig2Tips1Site() {
        Sequence s1 = new Sequence("t1", "AAA");
        Sequence s2 = new Sequence("t2", "CCC");

        Alignment data = new Alignment();
        data.initByName("sequence", s1, "sequence", s2, "dataType", "nucleotide" );

        // create Codon Alignment
        CodonAlignment codonAlignment = new CodonAlignment();
        codonAlignment.initByName("data", data, "dataType", "codon",
                "geneticCode", "vertebrateMitochondrial", "verbose", false);

        return codonAlignment;
    }

    // need to set internal node states after this
    private DataAugTreeLikelihood getDATreeLikelihood(String pi, String newickTree, boolean adjustTipHeights) {

        CodonAlignment codonAlignment = getAlig2Tips1Site();

        CodonFrequencies codonFreq = new CodonFrequencies();
        codonFreq.initByName("pi", pi, "data", codonAlignment, "verbose", true);
        SiteModel siteModel = CodonData.getSiteModel("0.3", "5", codonFreq, false);

        Tree tree = CodonData.getTree(codonAlignment, newickTree, adjustTipHeights);
        System.out.println("Tree is " + newickTree + "\n");

        System.setProperty("java.only","true");

        // =============== DA tree likelihood ===============
        // this only init NodeStates for tips
        NodeStatesArray nodesStates = new NodeStatesArray(codonAlignment, tree);
        // set internal node states on the fly

        DataAugTreeLikelihood daTreeLikelihood = new DataAugTreeLikelihood();
        daTreeLikelihood.initByName("nodesStates", nodesStates, "tree", tree, "siteModel", siteModel);
        return daTreeLikelihood;
    }

    private TreeLikelihood getTreeLikelihood(String pi, String newickTree, boolean adjustTipHeights) {

        CodonAlignment codonAlignment = getAlig2Tips1Site();

        CodonFrequencies codonFreq = new CodonFrequencies();
        codonFreq.initByName("pi", pi, "data", codonAlignment, "verbose", true);
        SiteModel siteModel = CodonData.getSiteModel("0.3", "5", codonFreq, false);

        Tree tree = CodonData.getTree(codonAlignment, newickTree, adjustTipHeights);
        System.out.println("Tree is " + newickTree + "\n");

        System.setProperty("java.only","true");

        // =============== Standard tree likelihood ===============
        TreeLikelihood treeLikelihood = new TreeLikelihood();
        treeLikelihood.initByName("data", codonAlignment, "tree", tree, "siteModel", siteModel);

        return treeLikelihood;
    }

}
