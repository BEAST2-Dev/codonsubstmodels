package codonmodels.evolution.datatype;

import codonmodels.evolution.CodonData;
import org.junit.Before;
import org.junit.Test;

import java.util.Arrays;

import static junit.framework.Assert.assertEquals;
import static org.junit.Assert.*;

/**
 * @author Walter Xie
 */
public class GeneticCodeTest {
    protected GeneticCode geneticCode;

    @Before
    public void setUp() {
        geneticCode = GeneticCode.UNIVERSAL;
    }


    // Code Table index != state
    @Test
    public void testStopCodonIndexInCodeTable() {
        int[] indices = geneticCode.getStopCodonIndices();
        assertArrayEquals("Stop codon : ", new int[]{48, 50, 56}, indices);

        assertTrue(geneticCode.isStopCodonIndex(48));
        assertTrue(geneticCode.isStopCodonIndex(50));
        assertTrue(geneticCode.isStopCodonIndex(56));

        assertFalse(geneticCode.isStopCodonIndex(0));
        assertFalse(geneticCode.isStopCodonIndex(49));
        assertFalse(geneticCode.isStopCodonIndex(63));

        //treat ? and - as not stop-codon
//        assertFalse(geneticCode.isStopCodonIndex(64));
//        assertFalse(geneticCode.isStopCodonIndex(120));

        GeneticCode geneticCode2 = GeneticCode.VERTEBRATE_MT;
        indices = geneticCode2.getStopCodonIndices();
        assertArrayEquals("Stop codon : ", new int[]{8,10,48,50}, indices);

        assertTrue(geneticCode2.isStopCodonIndex(8));
        assertTrue(geneticCode2.isStopCodonIndex(10));
        assertTrue(geneticCode2.isStopCodonIndex(48));
        assertTrue(geneticCode2.isStopCodonIndex(50));
    }

    // last index is 124 if including ambiguous state
    @Test
    public void testGetAminoAcid() {
        // assume indices
        int[] indices = CodonData.getIntegers();
        System.out.println("Testing code table indices : " + Arrays.toString(indices));
        String[] triExpected = new String[]{"AAA","CCA","GGA","TGT","TTA","TTT","AA-","--C"};
        System.out.println("Expecting triplets : " + Arrays.toString(triExpected));
        Character[] aaExpected = new Character[]{'K','P','G','C','L','F','-','-'};
        System.out.println("Expecting amino acids : " + Arrays.toString(aaExpected));

        // test code table indexing
        System.out.println("idx\ttri\taa");
        for (int i = 0; i < indices.length; i++) {
            int index = indices[i];
            Character aa = geneticCode.getAminoAcidFromCodeTable(index);
            String tri = Codon.CODON_TRIPLETS[index];
            System.out.println(index + "\t" + tri + "\t" + aa);
            assertEquals("Triplet : ", triExpected[i], tri);
            assertEquals("Amino Acid : ", aaExpected[i], aa);
        }
    }

    // last state is 121 or 120 if including ambiguous state
    @Test
    public void testGetAminoAcidCodonState() {
        int[] states = CodonData.getIntegers();
        System.out.println("Codon states : " + Arrays.toString(states));
        String[] triExpected = new String[]{"AAA","CCA","GGA","TTG","TTT","AG-","AT-","---"};
        System.out.println("Expecting triplets : " + Arrays.toString(triExpected));
        Character[] aaExpected = new Character[]{'K','P','G','L','F','-','-','-'};
        System.out.println("Expecting amino acids : " + Arrays.toString(aaExpected));

        Codon codon = new Codon(geneticCode);
        // test state indexing
        System.out.println("sta\ttri\taaS\taa");
        for (int i = 0; i < states.length; i++) {
            int state = states[i];
            int aaState = geneticCode.getAminoAcidState(state);
            Character aa = GeneticCode.AMINOACID_CHARS[aaState];
            String tri = codon.encodingToString(new int[]{state});
            System.out.println(state + "\t" + tri + "\t" + aaState + "\t" + aa);
            assertEquals("Triplet : ", triExpected[i], tri);
            assertEquals("Amino Acid : ", aaExpected[i], aa);
        }
    }

}