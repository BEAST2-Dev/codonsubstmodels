package test.beast.evolution.datatype;

import beast.evolution.datatype.Codon;
import beast.evolution.datatype.GeneticCode;
import org.junit.Test;
import test.beast.evolution.CodonTestData;

import static junit.framework.Assert.assertEquals;

/**
 * @author Walter Xie
 */
public class GeneticCodeTest {


    @Test
    public void testStopCodon() {


    }

    @Test
    public void testGetAminoAcid() {
        GeneticCode geneticCode = GeneticCode.UNIVERSAL;
        int[] states = CodonTestData.getTestStates();
        Character[] aaExpected = new Character[]{'K','P','G','L','F','?','-'};
        for (int i = 0; i < states.length; i++) {
            int state = states[i];
            Character aa = geneticCode.getAminoAcid(state);
            System.out.println(state + "  " + Codon.CODON_TRIPLETS[state] + "  " + aa);
            assertEquals("Amino Acid : ", aaExpected[i], aa);
        }
    }


    @Test
    public void testGetAminoAcidCodonState() {
        GeneticCode geneticCode = GeneticCode.UNIVERSAL;
        int[] states = CodonTestData.getTestStates();
        Character[] aaExpected = new Character[]{'K','P','G','L','F','?','-'};
        for (int i = 0; i < states.length; i++) {
            int state = states[i];
            int aaState = geneticCode.getAminoAcidState(state);
            Character aa = GeneticCode.AMINOACID_CHARS[aaState];
            System.out.println(state + "  " + Codon.CODON_TRIPLETS[state] + "  " + aaState + "  " + aa);
            assertEquals("Amino Acid : ", aaExpected[i], aa);
        }
    }

}