package test.beast.evolution.alignment;

import beast.evolution.alignment.Sequence;
import org.junit.Test;
import test.beast.evolution.CodonTestData;

import java.util.List;

import static junit.framework.Assert.assertEquals;

/**
 * @author Walter Xie
 */
public class CodonAlignmentTest {

    @Test
    public void testTripletsStates() {
        Sequence seq1 = CodonTestData.getSeq1();
        System.out.println(seq1);

        List<Integer> codonStates = seq1.getSequence(CodonTestData.codonUNIVERSAL);
        System.out.println(codonStates);

        // define codon triplets states as indices of Codon.CODON_TRIPLETS
        for (int i=0; i < codonStates.size(); i++) {
            int tripletState = codonStates.get(i);
            assertEquals("Triplets states : ", i, tripletState);
        }
    }

    @Test
    public void convertCodonToState() {

    }
}