package test.beast.evolution.likelihood;


import beast.likelihood.DAStatesLikelihoodCore;
import org.junit.Before;
import org.junit.Test;

import java.util.Arrays;

import static junit.framework.Assert.assertEquals;

/**
 * @author Walter Xie
 */
public class DAStatesLikelihoodCoreTest {

    DAStatesLikelihoodCore daLDCore1site;
    DAStatesLikelihoodCore daLDCore;
    // flattened transition probability matrix
    double[] p0;
    double[] p1;

    @Before
    public void setUp() throws Exception {
        daLDCore1site = new DAStatesLikelihoodCore(64);
        // 2 tips + 1 internal node, 1 codon, 1 category
        daLDCore1site.initialize( 3, 1, 1,
                true, false   );

        daLDCore = new DAStatesLikelihoodCore(64);
        // 2 tips + 1 internal node, 2 codons, 4 category
        daLDCore.initialize( 3, 2, 4,
                true, false   );

        // transition probability matrix
        int matrixSize = 64 * 64;
        p0 = new double[matrixSize];
        p1 = new double[matrixSize];
        for (int i=0; i < matrixSize; i++) {
            p0[i] = i * 0.1;
            p1[i] = i * 0.01;
        }

    }

    @Test
    public void testInit() {
        assertEquals(64, daLDCore1site.getNrOfStates());
        assertEquals(3, daLDCore1site.getNrOfNodes());
        assertEquals(1, daLDCore1site.getNrOfSites());
        assertEquals(1, daLDCore1site.getNrOfCategories()); // nrOfMatrices
        assertEquals(4096, daLDCore1site.getMatrixSize()); // nrOfStates^2
        assertEquals(1, daLDCore1site.getPartialsSize()); // = nrOfSites * nrOfMatrices;

        assertEquals(2, daLDCore.getNrOfSites());
        assertEquals(4, daLDCore.getNrOfCategories());
        assertEquals(4096, daLDCore.getMatrixSize());
        assertEquals(8, daLDCore.getPartialsSize());
    }


    @Test
    public void calculatePartials() {

        daLDCore1site.createNodePartials(2);

        // set codon state to 1 site, state = [0, 63]
        daLDCore1site.setNodeStates(0, new int[]{1});
        daLDCore1site.setNodeStates(1, new int[]{60});
        daLDCore1site.setNodeStates(2, new int[]{8});

        daLDCore1site.setNodeMatrix(0, 0, p0);
        daLDCore1site.setNodeMatrix(1, 0, p1);

        // DA intermediate likelihood per site at parent node
        daLDCore1site.calculatePartials(0,1,2);

        // test index: p0_1_8 = 6.4 + 0.8 = 7.2, p1_60_8 = 60 * 0.64 + 0.08 = 38.48
        double[] partials = new double[daLDCore1site.getPartialsSize()];
        daLDCore1site.getNodePartials(2, partials);
        System.out.println("partials = " + Arrays.toString(partials));

        // 1 site
        assertEquals(7.2 * 38.48, partials[0], 1e-6);

        // ======= 2 tips + 1 internal node, 2 codons, 4 category =======
        daLDCore.createNodePartials(2);

        // set codon state to 1 site, state = [0, 63]
        daLDCore.setNodeStates(0, new int[]{1,2});
        daLDCore.setNodeStates(1, new int[]{60,61});
        daLDCore.setNodeStates(2, new int[]{8,9});

        daLDCore.setNodeMatrix(0, 0, p0);
        daLDCore.setNodeMatrix(1, 0, p1);
        daLDCore.setNodeMatrix(0, 1, p0);
        daLDCore.setNodeMatrix(1, 1, p1);
        daLDCore.setNodeMatrix(0, 2, p0);
        daLDCore.setNodeMatrix(1, 2, p1);
        daLDCore.setNodeMatrix(0, 3, p0);
        daLDCore.setNodeMatrix(1, 3, p1);

        // DA intermediate likelihood per site at parent node
        daLDCore.calculatePartials(0,1,2);

        // test index category 1 :
        // p0_1_8 = 6.4 + 0.8 = 7.2, p1_60_8 = 60 * 0.64 + 0.08 = 38.48
        // p0_2_9 = 6.4 * 2 + 0.9 = 13.7, p1_61_9 = 61 * 0.64 + 0.09 = 39.13
        partials = new double[daLDCore.getPartialsSize()];
        daLDCore.getNodePartials(2, partials);
        System.out.println("partials = " + Arrays.toString(partials));

        // 2 sites
        assertEquals(7.2 * 38.48, partials[0], 1e-6);
        assertEquals(13.7 * 39.13, partials[1], 1e-6);

        // TODO partials = [277.05600000000004, 536.0810000000001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
//        assertEquals(7.2 * 38.48, partials[2], 1e-6);
//        assertEquals(13.7 * 39.13, partials[3], 1e-6);
    }
}
