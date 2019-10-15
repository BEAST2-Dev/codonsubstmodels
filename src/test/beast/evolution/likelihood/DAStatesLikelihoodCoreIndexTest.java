package test.beast.evolution.likelihood;


import beast.likelihood.DAStatesLikelihoodCore;
import beast.tree.InternalNodeStates;
import org.junit.Before;
import org.junit.Test;

import java.util.Arrays;

import static junit.framework.Assert.assertEquals;

/**
 * Build transition probability matrices under 3-node tree,
 * 2 tips + 1 internal node, to test if the index lookup is correct
 * during the data augmentation likelihood calculation.
 *
 * @author Walter Xie
 */
public class DAStatesLikelihoodCoreIndexTest {

    DAStatesLikelihoodCore daLDCore1site;
    DAStatesLikelihoodCore daLDCore;
    // flattened transition probability matrix
    double[] p0;
    double[] p1;

    @Before
    public void setUp() throws Exception {
        final int nrOfState = 64;

        // ======= transition probability matrix =======
        int matrixSize = nrOfState * nrOfState;
        p0 = new double[matrixSize];
        p1 = new double[matrixSize];
        for (int i=0; i < matrixSize; i++) {
            p0[i] = i * 0.1;
            p1[i] = i * 0.01;
        }

        // ======= 2 tips + 1 internal node, 1 codon, 1 category =======

        // set codon state to 1 site, state = [0, 63]
        InternalNodeStates internalNodeStates = new InternalNodeStates(nrOfState, new int[][]{{8}});
        daLDCore1site = new DAStatesLikelihoodCore(nrOfState);
        // 2 tips + 1 internal node, 1 codon, 1 category
        daLDCore1site.initialize( new int[][]{{1},{60}}, internalNodeStates, 1);

        daLDCore1site.setNodeMatrix(0, 0, p0);
        daLDCore1site.setNodeMatrix(1, 0, p1);

        // DA intermediate likelihood per site at parent node
        daLDCore1site.integrateNodeBranchLdOverCategories(0,1,2, new double[]{1});


        // ======= 2 tips + 1 internal node, 2 codons, 4 category =======

        // set codon state to 1 site, state = [0, 63]
        internalNodeStates = new InternalNodeStates(nrOfState, new int[][]{{8,9}});
        daLDCore = new DAStatesLikelihoodCore(nrOfState);
        // 2 tips + 1 internal node, 2 codons, 4 category
        daLDCore.initialize( new int[][]{{1,2},{60,61}}, internalNodeStates, 4);

        // set same p to 4 categories
        for (int i=0; i < daLDCore.getNrOfCategories(); i++) {
            daLDCore.setNodeMatrix(0, i, p0);
            daLDCore.setNodeMatrix(1, i, p1);
        }

        double[] proportions = new double[]{0.1, 0.2, 0.3, 0.1};
        // DA intermediate likelihood per site at parent node
        daLDCore.integrateNodeBranchLdOverCategories(0,1,2, proportions);

    }

    @Test
    public void testInit() {
        assertEquals(64, daLDCore1site.getNrOfStates());
        assertEquals(3, daLDCore1site.getNrOfNodes());
        assertEquals(1, daLDCore1site.getNrOfSites());
        assertEquals(1, daLDCore1site.getNrOfCategories()); // nrOfCategories
        assertEquals(4096, daLDCore1site.getMatrixSize()); // nrOfStates^2
        assertEquals(1, daLDCore1site.getBranchLdSize()); // = nrOfSites * nrOfCategories;

        assertEquals(2, daLDCore.getNrOfSites());
        assertEquals(4, daLDCore.getNrOfCategories());
        assertEquals(4096, daLDCore.getMatrixSize());
        assertEquals(2, daLDCore.getBranchLdSize());
    }


    @Test
    public void testNodeBranchLdIndices() {

        // test index: p0_1_8 = 6.4 + 0.8 = 7.2, p1_60_8 = 60 * 0.64 + 0.08 = 38.48
        double[] brLd = new double[daLDCore1site.getBranchLdSize()];
        daLDCore1site.getNodeBranchLd(2, brLd);
        System.out.println("branchLd = " + Arrays.toString(brLd));

        // 1 site
        assertEquals(7.2 * 38.48, brLd[0], 1e-6);
    }

    @Test
    public void testIntegratedBrLd() {

        // test index category 1 :
        // p0_1_8 = 6.4 + 0.8 = 7.2, p1_60_8 = 60 * 0.64 + 0.08 = 38.48
        // p0_2_9 = 6.4 * 2 + 0.9 = 13.7, p1_61_9 = 61 * 0.64 + 0.09 = 39.13
        double[] integratedBrLd = new double[daLDCore.getBranchLdSize()];
        daLDCore.getNodeBranchLd(2, integratedBrLd);
        System.out.println("branchLd = " + Arrays.toString(integratedBrLd));

        // 2 sites, proportions = 0.1 + 0.2 + 0.3 + 0.1 = 0.7
        assertEquals(7.2 * 38.48 * 0.7, integratedBrLd[0], 1e-6);
        assertEquals(13.7 * 39.13 * 0.7, integratedBrLd[1], 1e-6);
    }



}
