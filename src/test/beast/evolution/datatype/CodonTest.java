package test.beast.evolution.datatype;

import beast.evolution.datatype.Codon;
import org.junit.Before;
import org.junit.Test;
import test.beast.evolution.CodonTestData;

import static junit.framework.Assert.assertEquals;

/**
 * @author Walter Xie
 */
public class CodonTest {
    protected Codon codon;

    @Before
    public void setUp() {
        codon = CodonTestData.codonUNIVERSAL;
    }


    @Test
    public void testGetTriplet() {
        int[] states = CodonTestData.getTestStates();
        for (int i = 0; i < states.length; i++) {
            int state = states[i];
            String triplet = codon.getTriplet(state);
            System.out.println(state + "  " + triplet);
            assertEquals("Triplet : ", Codon.CODON_TRIPLETS[state], triplet);
        }
    }

    @Test(expected = RuntimeException.class)
    public void testGetTripletException1() {
        String triplet = codon.getTriplet(-1);
    }

    @Test(expected = RuntimeException.class)
    public void testGetTripletException2() {
        String triplet = codon.getTriplet(66);
    }


}