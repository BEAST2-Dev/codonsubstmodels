package codonmodels.evolution.datatype;

import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.Nucleotide;
import codonmodels.evolution.CodonData;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

/**
 * @author Walter Xie
 */
public class CodonTest {
    protected Codon codon;

    @BeforeEach
    public void setUp() {
        codon = CodonData.codonUNIVERSAL;
    }

    @Test
    public void testCodeMap(){
        int stateCount = codon.getStateCount();
        System.out.println("stateCount = " + stateCount);
        assertEquals(61, stateCount, "Codon stateCount : ");

        stateCount = codon.getStateCountAmbiguous();
        System.out.println("Max stateCount including ambiguous = " + stateCount);
        assertEquals(122, stateCount, "Max stateCount including ambiguous : ");

        int codeLength = codon.getCodeLength();
        System.out.println("codeLength = " + codeLength);
        assertEquals(3, codeLength, "Codon codeLength : ");

        String codeMap = codon.getCodeMap();
        System.out.println("codeMap = " + codeMap);
        assertEquals(codon.getTripletsNoStopCodon(), codeMap, "Codon codeMap : ");
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
        assertEquals(toList(0, 60), states, "StringToEncoding : ");

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
        assertEquals(toList(61, 109, 121), states, "StringToEncoding : ");
    }

    @Test
    public void testEncodingToString() {
        String data = codon.encodingToString(toList(0, 60));
        assertEquals("AAATTT", data, "StringToEncoding : ");
        data = codon.encodingToString(toList(48, 50, 56));
        assertNotSame("TAATAGTGA", data, "StringToEncoding : ");
        data = codon.encodingToString(toList(61, 109, 121));
        assertEquals("AA-A-----", data, "StringToEncoding : ");
    }

    @Test
    public void testNucleotideState() {
        Nucleotide nuc = new Nucleotide();
        int gapState = nuc.getCodeMap().indexOf(DataType.GAP_CHAR);
        assertEquals(Codon.NUC_GAP_STATE, gapState, "Nucleotide '-' state : ");
        int missState = nuc.getCodeMap().indexOf(DataType.MISSING_CHAR);
        assertEquals(Codon.NUC_MISSING_STATE, missState, "Nucleotide '?' state : ");
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
            assertEquals(triExpected[i], triplet, "Triplet : ");
        }
    }

    @Test
    public void testGetTripletException1() {
        assertThrows(RuntimeException.class, () -> codon.stateToTriplet(-1));
    }

    // last state = getStateCountAmbiguous()-1
    @Test
    public void testGetTripletException2() {
        assertThrows(RuntimeException.class, () -> codon.stateToTriplet(codon.getStateCountAmbiguous()));
    }

    @Test
    public void testStateToAminoAcid() {
        int[] states = new int[]{0, 8, 60};
        String aa = codon.stateToAminoAcid(states);
        System.out.println("StateToAminoAcid : 0, 8, 60 => " + aa);
        assertEquals("KRF", aa, "StateToAminoAcid : ");

        // no stop codons
        states = new int[]{48, 50, 56};
        aa = codon.stateToAminoAcid(states);
        System.out.println("StateToAminoAcid : 48, 50, 56 => " + aa);
        assertEquals("YSC", aa, "StateToAminoAcid : ");

        // ambiguous states
        states = new int[]{61, 121};
        aa = codon.stateToAminoAcid(states);
        System.out.println("StateToAminoAcid : 61, 121 => " + aa);
        assertEquals("--", aa, "StateToAminoAcid : ");

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