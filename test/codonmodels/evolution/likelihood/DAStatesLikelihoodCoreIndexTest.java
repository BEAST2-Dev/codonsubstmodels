package codonmodels.evolution.likelihood;


import org.junit.Before;
import org.junit.Test;

import java.util.Arrays;

import static org.junit.Assert.assertEquals;

/**
 * Build transition probability matrices under 3-node tree,
 * 2 tips + 1 internal node, to test if the index lookup is correct
 * during the data augmentation likelihood calculation.
 *
 * @author Walter Xie
 */
public class DAStatesLikelihoodCoreIndexTest {

    DABranchLikelihoodCore daLdCore1siteBr1;
    DABranchLikelihoodCore daLdCore1siteBr2;
    DABranchLikelihoodCore daLdCoreBr1;
    DABranchLikelihoodCore daLdCoreBr2;
    // flattened transition probability matrix
    double[] p0;
    double[] p1;

    @Before
    public void setUp() throws Exception {
        final int nrOfState = 61;

        // ======= transition probability matrix =======
        int matrixSize = nrOfState * nrOfState;
        p0 = new double[matrixSize];
        p1 = new double[matrixSize];
        for (int i=0; i < matrixSize; i++) {
            p0[i] = i * 0.1;
            p1[i] = i * 0.01;
        }

        // ======= 2 tips + 1 internal node
//        Tree tree = new Tree("((0:1.0,1:1.0)2:0.0);");

        // ======= 1 codon, 1 category =======
        int[] statesParentNode = new int[]{8};
        int[] statesChildNode1 = new int[]{1};
        int[] statesChildNode2 = new int[]{60};
        // branch 1
        daLdCore1siteBr1 = new DABranchLikelihoodCore(0, nrOfState, 1, 1);
        daLdCore1siteBr1.setNodeMatrix(0, p0);
        // set codon state to 1 site, state = [0, 63]
        daLdCore1siteBr1.calculateBranchLd(statesParentNode, statesChildNode1, new double[]{1});

        // branch 2
        daLdCore1siteBr2 = new DABranchLikelihoodCore(1, nrOfState, 1, 1);
        daLdCore1siteBr2.setNodeMatrix(0, p1);
        daLdCore1siteBr2.calculateBranchLd(statesParentNode, statesChildNode2, new double[]{1});

        // ======= 2 codons, 4 category =======
        statesParentNode = new int[]{8,9};
        statesChildNode1 = new int[]{1,2};
        statesChildNode2 = new int[]{60,59};
        double[] proportions = new double[]{0.1, 0.2, 0.3, 0.1};
        // branch 1
        daLdCoreBr1 = new DABranchLikelihoodCore(0, nrOfState, 2, 4);
        for (int i = 0; i < daLdCoreBr1.getNrOfCategories(); i++)
            daLdCoreBr1.setNodeMatrix(i, p0);
        daLdCoreBr1.calculateBranchLd(statesParentNode, statesChildNode1, proportions);

        // branch 2
        daLdCoreBr2 = new DABranchLikelihoodCore(1, nrOfState, 2, 4);
        for (int i = 0; i < daLdCoreBr2.getNrOfCategories(); i++)
            daLdCoreBr2.setNodeMatrix(i, p1);
        daLdCoreBr2.calculateBranchLd(statesParentNode, statesChildNode2, proportions);
    }

    @Test
    public void testInit() {
        assertEquals(61, daLdCore1siteBr1.getNrOfStates());
        assertEquals(1, daLdCore1siteBr1.getNrOfSites());
        assertEquals(1, daLdCore1siteBr1.getNrOfCategories()); // nrOfCategories
        assertEquals(3721, daLdCore1siteBr1.getMatrixSize()); // nrOfStates^2
        assertEquals(1, daLdCore1siteBr1.getBranchLdSize()); // = nrOfSites

        assertEquals(2, daLdCoreBr1.getNrOfSites());
        assertEquals(4, daLdCoreBr1.getNrOfCategories());
        assertEquals(3721, daLdCoreBr1.getMatrixSize());
        assertEquals(2, daLdCoreBr1.getBranchLdSize());
    }


    // DA intermediate likelihood per site at the branch
    @Test
    public void testNodeBranchLdIndices() {
        // test index
        double[] brLd = new double[daLdCore1siteBr1.getBranchLdSize()];
        daLdCore1siteBr1.getBranchLikelihoods(brLd);
        System.out.println("branchLd = " + Arrays.toString(brLd));
        // 1 site 8=>1 : p0_8_1 = 8*(61*0.1) + 1*0.1 = 48.9
        assertEquals(48.9, brLd[0], 1e-6);

        daLdCore1siteBr2.getBranchLikelihoods(brLd);
        System.out.println("branchLd = " + Arrays.toString(brLd));
        // 1 site 8=>60 : p1_8_60 = 8 * (61*0.01) + 60*0.01 = 5.48
        assertEquals(5.48, brLd[0], 1e-6);
    }

    @Test
    public void testIntegratedBrLd() {
        // test index category 1 :
        double[] integratedBrLd = new double[daLdCoreBr1.getBranchLdSize()];
        daLdCoreBr1.getBranchLikelihoods(integratedBrLd);
        System.out.println("branchLd = " + Arrays.toString(integratedBrLd));
        // 2 sites, proportions = 0.1 + 0.2 + 0.3 + 0.1 = 0.7
        // p0_8_1 = 48.9, p0_9_2 = 55.1
        assertEquals(48.9 * 0.7, integratedBrLd[0], 1e-6);
        assertEquals(55.1 * 0.7, integratedBrLd[1], 1e-6);

        daLdCoreBr2.getBranchLikelihoods(integratedBrLd);
        System.out.println("branchLd = " + Arrays.toString(integratedBrLd));
        // 2 sites, proportions = 0.1 + 0.2 + 0.3 + 0.1 = 0.7
        // p1_8_60 = 5.48, p1_9_59 = 6.08
        assertEquals(5.48 * 0.7, integratedBrLd[0], 1e-6);
        assertEquals(6.08 * 0.7, integratedBrLd[1], 1e-6);
    }

}
