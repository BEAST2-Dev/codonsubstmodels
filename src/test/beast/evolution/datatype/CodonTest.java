package test.beast.evolution.datatype;

import beast.evolution.datatype.Codon;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.Nucleotide;
import org.junit.Before;
import org.junit.Test;
import test.beast.evolution.CodonTestData;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

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
    public void testCodeMap(){
        int stateCount = codon.getStateCount();
        System.out.println("stateCount = " + stateCount);
        assertEquals("Codon stateCount : ", 64, stateCount);

        int codeLength = codon.getCodeLength();
        System.out.println("codeLength = " + codeLength);
        assertEquals("Codon codeLength : ", 3, codeLength);

        String codeMap = codon.getCodeMap();
        System.out.println("codeMap = " + codeMap);
        assertEquals("Codon codeMap : ", CodonTestData.getTriplets(), codeMap);
    }

//    @Test
//    public void testGetStatesForCode() {
//        for (int i = 0; i < codon.getStateCountAmbiguous(); i++) {
//            int[] stateSet = codon.getStatesForCode(i);
//            System.out.println("mapCodeToStateSet " + i + " = " + Arrays.toString(stateSet));
//        }
//        assertEquals("Codon mapCodeToStateSet : ", 0, codon.getStatesForCode(0)[0]);
//        assertEquals("Codon mapCodeToStateSet : ", 10, codon.getStatesForCode(10)[0]);
//        assertEquals("Codon mapCodeToStateSet : ", 49, codon.getStatesForCode(48)[0]);
//        assertEquals("Codon mapCodeToStateSet : ", 51, codon.getStatesForCode(49)[0]);
//        assertEquals("Codon mapCodeToStateSet : ", 57, codon.getStatesForCode(54)[0]);
//        assertEquals("Codon mapCodeToStateSet : ", 63, codon.getStatesForCode(60)[0]);
//        // stop codon
//        assertEquals("Codon mapCodeToStateSet : ", 48, codon.getStatesForCode(61)[0]);
//        assertEquals("Codon mapCodeToStateSet : ", 50, codon.getStatesForCode(62)[0]);
//        assertEquals("Codon mapCodeToStateSet : ", 56, codon.getStatesForCode(63)[0]);
//    }

    @SafeVarargs
    private final <T> List<T> toList(T... a) {
        return new ArrayList<>(Arrays.asList(a));
    }

    @Test
    public void testStringToEncoding() {
        String data = CodonTestData.getTriplets();
        List<Integer> states = codon.stringToEncoding(data);
        System.out.println("StringToEncoding : " + data);
        System.out.println("StringToEncoding : " + states);

        data = "AAATTT";
        states = codon.stringToEncoding(data);
        assertEquals("StringToEncoding : ", toList(0, 63), states);

        // stop codon
        data = "TAATAGTGA";
        states = codon.stringToEncoding(data);
        assertEquals("StringToEncoding : ", toList(48, 50, 56), states);

        // ambiguous
        data = "???---";
        states = codon.stringToEncoding(data);
        assertEquals("StringToEncoding : ", toList(64, 65), states);
    }

    @Test
    public void testEncodingToString() {
        String data = codon.encodingToString(toList(0, 63));
        assertEquals("StringToEncoding : ", "AAATTT", data);
        data = codon.encodingToString(toList(48, 50, 56));
        assertEquals("StringToEncoding : ", "TAATAGTGA", data);
        data = codon.encodingToString(toList(64, 65));
        assertEquals("StringToEncoding : ", "???---", data);
    }

    @Test
    public void testNucleotideState() {
        Nucleotide nuc = new Nucleotide();
        int gapState = nuc.getCodeMap().indexOf(DataType.GAP_CHAR);
        assertEquals("Nucleotide '-' state : ", Codon.NUC_GAP_STATE, gapState);
        int missState = nuc.getCodeMap().indexOf(DataType.MISSING_CHAR);
        assertEquals("Nucleotide '?' state : ", Codon.NUC_MISSING_STATE, missState);
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

    @Test
    public void testStateToAminoAcid() {
        int[] states = new int[]{0, 8, 63};
        String aa = codon.stateToAminoAcid(states);
        System.out.println("StateToAminoAcid : 0, 8, 63 => " + aa);
        assertEquals("StateToAminoAcid : ", "KRF", aa);

        states = new int[]{48, 50, 56};
        aa = codon.stateToAminoAcid(states);
        System.out.println("StateToAminoAcid : 48, 50, 56 => " + aa);
        assertEquals("StateToAminoAcid : ", "***", aa);

        states = new int[]{64, 65};
        aa = codon.stateToAminoAcid(states);
        System.out.println("StateToAminoAcid : 64, 65 => " + aa);
        assertEquals("StateToAminoAcid : ", "?-", aa);
    }
}