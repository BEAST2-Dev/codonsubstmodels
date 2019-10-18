package test.beast.evolution.likelihood;

import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.CodonAlignment;
import beast.evolution.alignment.Sequence;
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

import static junit.framework.Assert.assertEquals;


/**
 *
 *
 * @author Walter Xie
 */
public class DALikelihoodTest {
    TreeLikelihood treeLikelihood;
    DACodonTreeLikelihood daTreeLikelihood;
    int rootNr;
    InternalNodeStates internalNodeStates;

    DecimalFormat df = new DecimalFormat("#.00");

    @Before
    public void setUp() {

//        String newickTree = "(t1:0.5, t2:0.0):0.0;"; // TODO bug for 0 branch length ?
        String newickTree = "(t1:0.5, t2:0.1):0.0;";

        // =============== Standard tree likelihood ===============
        treeLikelihood = getTreeLikelihood("equal", newickTree, false);

        // =============== DA tree likelihood ===============
        // need to set internal node states
        daTreeLikelihood = getDATreeLikelihood("equal", newickTree, false);

        rootNr = daTreeLikelihood.treeInput.get().getRoot().getNr();
        System.out.println("Root index = " + rootNr);
    }

    @Test
    public void testDALikelihood1Site(){
        System.out.println("\n=============== Standard tree likelihood ===============\n");

        double logP = treeLikelihood.calculateLogP();
        System.out.println(" logP = " + logP);

        double[] partial = treeLikelihood.getRootPartials(); // root
        System.out.println("Root partial without freqs = " + Arrays.toString(partial));

        System.out.println("\n=============== DA tree likelihood ===============");
        final double freq = 1.0/60.0;
        double totalLogPDA = 0;
        // 0 - 63
        for (int s=0; s < 64; s++) {
            // stop codon states in vertebrateMitochondrial
            if (s != 8 && s != 10 && s != 48 && s != 50) {
                // internal nodes
                internalNodeStates.setNrStates(rootNr, new int[]{s});
                System.out.println("\nSet internal node state to " + s);

                double logPDA = daTreeLikelihood.calculateLogP();
                System.out.println(" DA logP = " + logPDA);

                // [branch][site]
                double[][] branchLd = new double[2][1];
                daTreeLikelihood.getBranchPartials(0, branchLd[0]); // child 1
                daTreeLikelihood.getBranchPartials(1, branchLd[1]); // child 2
                System.out.println("Branch likelihood child1 child2 without freqs = " +
                        Arrays.toString(branchLd[0]) + ", " + Arrays.toString(branchLd[1]));

                double nodeLd = branchLd[0][0] * branchLd[1][0];
                System.out.println("Root likelihood without freqs = " + nodeLd);
                System.out.println("Root partial without freqs = " + partial[s]);

                assertEquals(partial[s], nodeLd, 1e-16);

                totalLogPDA += nodeLd * freq;
            }
        }
        System.out.println("\nTotal :");
        System.out.println("Total DA tree likelihood logP = " + Math.log(totalLogPDA));
        System.out.println("Standard tree likelihood logP = " + logP);
        //Total DA tree likelihood logP = -15.667929146385474
        //Standard tree likelihood logP = -15.667929146369291
        assertEquals(logP, totalLogPDA, 1e-10);
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
    private DACodonTreeLikelihood getDATreeLikelihood(String pi, String newickTree, boolean adjustTipHeights) {

        CodonAlignment codonAlignment = getAlig2Tips1Site();

        CodonFrequencies codonFreq = new CodonFrequencies();
        codonFreq.initByName("pi", pi, "data", codonAlignment, "verbose", true);
        SiteModel siteModel = CodonTestData.getSiteModel("0.3", "5", codonFreq, false);

        Tree tree = CodonTestData.getTree(codonAlignment, newickTree, adjustTipHeights);
        System.out.println("Tree is " + newickTree + "\n");

        System.setProperty("java.only","true");

        // =============== DA tree likelihood ===============
        int internalNodeCount = tree.getInternalNodeCount();
        int siteCount = codonAlignment.getSiteCount();
        System.out.println(internalNodeCount + " internal nodes " + siteCount + " codon(s).");
        internalNodeStates = new InternalNodeStates(internalNodeCount, siteCount);

        DACodonTreeLikelihood daTreeLikelihood = new DACodonTreeLikelihood();
        daTreeLikelihood.initByName("data", codonAlignment, "tree", tree, "siteModel", siteModel,
                "internalNodeStates", internalNodeStates);
        return daTreeLikelihood;
    }

    private TreeLikelihood getTreeLikelihood(String pi, String newickTree, boolean adjustTipHeights) {

        CodonAlignment codonAlignment = getAlig2Tips1Site();

        CodonFrequencies codonFreq = new CodonFrequencies();
        codonFreq.initByName("pi", pi, "data", codonAlignment, "verbose", true);
        SiteModel siteModel = CodonTestData.getSiteModel("0.3", "5", codonFreq, false);

        Tree tree = CodonTestData.getTree(codonAlignment, newickTree, adjustTipHeights);
        System.out.println("Tree is " + newickTree + "\n");

        System.setProperty("java.only","true");

        // =============== Standard tree likelihood ===============
        TreeLikelihood treeLikelihood = new TreeLikelihood();
        treeLikelihood.initByName("data", codonAlignment, "tree", tree, "siteModel", siteModel);

        return treeLikelihood;
    }

}
