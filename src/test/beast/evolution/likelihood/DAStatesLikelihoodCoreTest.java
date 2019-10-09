package test.beast.evolution.likelihood;


import beast.likelihood.DAStatesLikelihoodCore;
import org.junit.Before;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;

/**
 * @author Walter Xie
 */
public class DAStatesLikelihoodCoreTest {

    DAStatesLikelihoodCore daLDCore1site;
    DAStatesLikelihoodCore daLDCore2sites;

    @Before
    public void setUp() throws Exception {
        daLDCore1site = new DAStatesLikelihoodCore(64);
        // 2 tips + 1 internal node, 1 codon, 1 category
        daLDCore1site.initialize( 3, 1, 1,
                true, false   );

        daLDCore2sites = new DAStatesLikelihoodCore(64);
        // 2 tips + 1 internal node, 1 codon, 1 category
        daLDCore2sites.initialize( 3, 2, 1,
                true, false   );
    }

    @Test
    public void testInit() {
        assertEquals(64, daLDCore1site.getNrOfStates());
        assertEquals(3, daLDCore1site.getNrOfNodes());
        assertEquals(1, daLDCore1site.getNrOfSites());
        assertEquals(1, daLDCore1site.getNrOfCategories());
        assertEquals(4096, daLDCore1site.getMatrixSize()); // nrOfStates^2
        assertEquals(1, daLDCore1site.getPartialsSize()); // = nrOfSites * nrOfMatrices;

        assertEquals(2, daLDCore2sites.getNrOfSites());
        assertEquals(2, daLDCore2sites.getPartialsSize());
    }


    @Test
    public void calculatePartials() {

        daLDCore1site.createNodePartials(2);

//        daLDCore1site.setNodeStates(0, );
//
//        daLDCore1site.setNodeMatrix(0, );

        daLDCore1site.calculatePartials(0,1,2);

        double[] partials = new double[daLDCore1site.getPartialsSize()];
        daLDCore1site.getNodePartials(2, partials);

    }
}
