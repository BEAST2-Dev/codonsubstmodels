package test.beast.evolution.datatype;

import beast.evolution.datatype.Codon;
import beast.evolution.datatype.GeneticCode;
import org.junit.Before;
import org.junit.Test;
import test.beast.evolution.CodonData;

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


    @Test
    public void testStopCodon() {
        int[] states = geneticCode.getStopCodonStates();
        assertArrayEquals("Stop codon : ", new int[]{48, 50, 56}, states);

        assertTrue(geneticCode.isStopCodon(48));
        assertTrue(geneticCode.isStopCodon(50));
        assertTrue(geneticCode.isStopCodon(56));

        assertFalse(geneticCode.isStopCodon(0));
        assertFalse(geneticCode.isStopCodon(49));
        assertFalse(geneticCode.isStopCodon(63));

        //treat ? and - as not stop-codon
        assertFalse(geneticCode.isStopCodon(64));
        assertFalse(geneticCode.isStopCodon(120));

        GeneticCode geneticCode2 = GeneticCode.VERTEBRATE_MT;
        states = geneticCode2.getStopCodonStates();
        assertArrayEquals("Stop codon : ", new int[]{8,10,48,50}, states);

        assertTrue(geneticCode2.isStopCodon(8));
        assertTrue(geneticCode2.isStopCodon(10));
        assertTrue(geneticCode2.isStopCodon(48));
        assertTrue(geneticCode2.isStopCodon(50));
    }

    @Test
    public void testGetAminoAcid() {
        int[] states = CodonData.getTestStates();
        Character[] aaExpected = new Character[]{'K','P','G','L','F','-','-'};
        for (int i = 0; i < states.length; i++) {
            int state = states[i];
            Character aa = geneticCode.getAminoAcid(state);
            System.out.println(state + "  " + Codon.CODON_TRIPLETS[state] + "  " + aa);
            assertEquals("Amino Acid : ", aaExpected[i], aa);
        }
    }


    @Test
    public void testGetAminoAcidCodonState() {
        int[] states = CodonData.getTestStates();
        Character[] aaExpected = new Character[]{'K','P','G','L','F','-','-'};
        for (int i = 0; i < states.length; i++) {
            int state = states[i];
            int aaState = geneticCode.getAminoAcidState(state);
            Character aa = GeneticCode.AMINOACID_CHARS[aaState];
            System.out.println(state + "  " + Codon.CODON_TRIPLETS[state] + "  " + aaState + "  " + aa);
            assertEquals("Amino Acid : ", aaExpected[i], aa);
        }
    }

}