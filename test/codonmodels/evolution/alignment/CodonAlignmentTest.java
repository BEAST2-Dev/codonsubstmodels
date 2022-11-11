package codonmodels.evolution.alignment;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.datatype.DataType;
import beast.pkgmgmt.BEASTClassLoader;
import codonmodels.evolution.CodonData;
import codonmodels.evolution.datatype.Codon;
import org.junit.Test;

import java.util.List;
import java.util.Set;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;


/**
 * @author Walter Xie
 */
public class CodonAlignmentTest {

    @Test
    public void testGetSequence() {
        Sequence seq1 = CodonData.getSeqCodeMap();
        System.out.println(seq1);

        List<Integer> codonStates = seq1.getSequence(CodonData.codonUNIVERSAL);
        System.out.println(codonStates);

        // define codon triplets states as indices of Codon.CODON_TRIPLETS
        for (int i=0; i < codonStates.size(); i++) {
            int tripletState = codonStates.get(i);
            assertEquals("Triplets states : ", i, tripletState);
        }
    }

    @Test
    public void testGettingDataType() {
        //TODO
        Set<String> dataTypes = BEASTClassLoader.loadService(DataType.class);
        System.out.println(dataTypes);
        assertTrue("Cannot load Codon data type !", dataTypes.contains(Codon.CODON));
    }

    @Test
    public void testCodonUsage() {
        Alignment data = CodonData.getAlig6T333();
        // create Codon Alignment
        CodonAlignment codonAlignment = new CodonAlignment();
        //TODO
        codonAlignment.initByName("data", data, "dataType", "codon", "geneticCode", "vertebrateMitochondrial");

        Codon codon = codonAlignment.getDataType();
        int[][] usage = codonAlignment.getCodonUsage();

        int state = codon.stringToEncoding("ATC").get(0);
        int[] expectedUsage = new int[]{21,22,16,16,18,19};
        for (int i=0; i < usage.length; i++)
            assertEquals("ATC frequency at taxon " + i, expectedUsage[i], usage[i][state]);

        state = codon.stringToEncoding("TTA").get(0);
        expectedUsage = new int[]{7,6,9,9,10,10};
        for (int i=0; i < usage.length; i++)
            assertEquals("ATC frequency at taxon " + i, expectedUsage[i], usage[i][state]);
    }
}