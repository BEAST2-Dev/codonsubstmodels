package test.beast.evolution;

import beast.evolution.alignment.Sequence;
import beast.evolution.datatype.Codon;
import beast.util.StringUtils;

/**
 * @author Walter Xie
 */
public class CodonTestData {

    public static final Codon codonUNIVERSAL = new Codon();

    // concatenate Codon.CODON_TRIPLETS to a string
    public static String getTriplets() {
        return StringUtils.concatenateToString(Codon.CODON_TRIPLETS);
    }

    // create a Sequence obj using Codon.codeMap
    public static Sequence getSeqCodeMap() {
        return new Sequence("codemap", codonUNIVERSAL.getCodeMap());
    }

    public static int[] getTestStates() {
        return new int[]{0, 20, 40, 60, 63, 64, 65};
    }
}
