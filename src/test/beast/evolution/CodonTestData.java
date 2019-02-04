package test.beast.evolution;

import beast.evolution.alignment.Sequence;
import beast.evolution.datatype.Codon;

/**
 * @author Walter Xie
 */
public class CodonTestData {

    public static final Codon codonUNIVERSAL = new Codon();

    // concatenate Codon.CODON_TRIPLETS to create a Sequence obj
    public static Sequence getSeq1() {
        return new Sequence("seq1", codonUNIVERSAL.getCodeMap());
    }

    public static int[] getTestStates() {
        return new int[]{0, 20, 40, 60, 63, 64, 65};
    }
}
