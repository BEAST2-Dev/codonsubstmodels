package codonmodels.evolution.datatype;

import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.Nucleotide;
import codonmodels.evolution.CodonData;
import org.junit.Before;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static junit.framework.Assert.*;

/**
 * @author Walter Xie
 */
public class CodonTest {
    protected Codon codon;

    @Before
    public void setUp() {
        codon = CodonData.codonUNIVERSAL;
    }

    @Test
    public void testCodeMap(){
        int stateCount = codon.getStateCount();
        System.out.println("stateCount = " + stateCount);
        assertEquals("Codon stateCount : ", 61, stateCount);

        stateCount = codon.getStateCountAmbiguous();
        System.out.println("Max stateCount including ambiguous = " + stateCount);
        assertEquals("Max stateCount including ambiguous : ", 122, stateCount);

        int codeLength = codon.getCodeLength();
        System.out.println("codeLength = " + codeLength);
        assertEquals("Codon codeLength : ", 3, codeLength);

        String codeMap = codon.getCodeMap();
        System.out.println("codeMap = " + codeMap);
        assertEquals("Codon codeMap : ", codon.getTripletsNoStopCodon(), codeMap);
    }

    @SafeVarargs
    private final <T> List<T> toList(T... a) {
        return new ArrayList<>(Arrays.asList(a));
    }

    @Test
    public void testStringToEncoding() {
        String data = codon.getTripletsNoStopCodon();
        List<Integer> states = codon.stringToEncoding(data);
        System.out.println("StringToEncoding : " + data);
        System.out.println("StringToEncoding : " + states);

        // 1st and last
        data = "AAATTT";
        states = codon.stringToEncoding(data);
        assertEquals("StringToEncoding : ", toList(0, 60), states);

        // stop codon
        Throwable caught = null;
        try {
            data = "TAA";
            codon.stringToEncoding(data);
        } catch (Throwable t) {
            caught = t;
        }
        assertNotNull(caught);
        assertSame(IllegalArgumentException.class, caught.getClass());

        caught = null;
        try {
            data = "TAG";
            codon.stringToEncoding(data);
        } catch (Throwable t) {
            caught = t;
        }
        assertNotNull(caught);
        assertSame(IllegalArgumentException.class, caught.getClass());

        caught = null;
        try {
            data = "TGA";
            codon.stringToEncoding(data);
        } catch (Throwable t) {
            caught = t;
        }
        assertNotNull(caught);
        assertSame(IllegalArgumentException.class, caught.getClass());

        // ambiguous
        data = "AA-A-----";
        states = codon.stringToEncoding(data);
        assertEquals("StringToEncoding : ", toList(61, 109, 121), states);
    }

    @Test
    public void testEncodingToString() {
        String data = codon.encodingToString(toList(0, 60));
        assertEquals("StringToEncoding : ", "AAATTT", data);
        data = codon.encodingToString(toList(48, 50, 56));
        assertNotSame("StringToEncoding : ", "TAATAGTGA", data);
        data = codon.encodingToString(toList(61, 109, 121));
        assertEquals("StringToEncoding : ", "AA-A-----", data);
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
        int[] states = CodonData.getIntegers();
        System.out.println("Codon states : " + Arrays.toString(states));
        String[] triExpected = new String[]{"AAA","CCA","GGA","TTG","TTT","AG-","AT-","---"};
        System.out.println("Expecting triplets : " + Arrays.toString(triExpected));
        for (int i = 0; i < states.length; i++) {
            int state = states[i];
            String triplet = codon.stateToTriplet(state);
            System.out.println(state + "  " + triplet);
            assertEquals("Triplet : ", triExpected[i], triplet);
        }
    }

    @Test(expected = RuntimeException.class)
    public void testGetTripletException1() {
        codon.stateToTriplet(-1);
    }

    // last state = getStateCountAmbiguous()-1
    @Test(expected = RuntimeException.class)
    public void testGetTripletException2() {
        codon.stateToTriplet(codon.getStateCountAmbiguous());
    }

    @Test
    public void testStateToAminoAcid() {
        int[] states = new int[]{0, 8, 60};
        String aa = codon.stateToAminoAcid(states);
        System.out.println("StateToAminoAcid : 0, 8, 60 => " + aa);
        assertEquals("StateToAminoAcid : ", "KRF", aa);

        // no stop codons
        states = new int[]{48, 50, 56};
        aa = codon.stateToAminoAcid(states);
        System.out.println("StateToAminoAcid : 48, 50, 56 => " + aa);
        assertEquals("StateToAminoAcid : ", "YSC", aa);

        // ambiguous states
        states = new int[]{61, 121};
        aa = codon.stateToAminoAcid(states);
        System.out.println("StateToAminoAcid : 61, 121 => " + aa);
        assertEquals("StateToAminoAcid : ", "--", aa);

        states = new int[]{125};
        Throwable caught = null;
        try {
            codon.stateToAminoAcid(states);
        } catch (Throwable t) {
            caught = t;
        }
        assertNotNull(caught);
        assertSame(IllegalArgumentException.class, caught.getClass());
    }
}