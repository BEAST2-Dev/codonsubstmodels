package codonmodels.evolution.alignment;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.datatype.DataType;
import beast.pkgmgmt.BEASTClassLoader;
import codonmodels.evolution.CodonData;
import codonmodels.evolution.datatype.Codon;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.List;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;


/**
 * @author Walter Xie
 */
public class CodonAlignmentTest {

    @BeforeEach
    public void setUp() throws Exception {
        CodonData.initDataTypes();
    }

    @Test
    public void testGetSequence() {
        Sequence seq1 = CodonData.getSeqCodeMap();
        System.out.println(seq1);

        List<Integer> codonStates = seq1.getSequence(CodonData.codonUNIVERSAL);
        System.out.println(codonStates);

        // define codon triplets states as indices of Codon.CODON_TRIPLETS
        for (int i=0; i < codonStates.size(); i++) {
            int tripletState = codonStates.get(i);
            assertEquals(i, tripletState, "Triplets states : ");
        }
    }

    @Test
    public void testGettingDataType() {
        Set<String> dataTypes = BEASTClassLoader.loadService(DataType.class);
        System.out.println(dataTypes);
        assertTrue(dataTypes.contains(Codon.class.getName()), "Cannot load Codon data type !");
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
            assertEquals(expectedUsage[i], usage[i][state], "ATC frequency at taxon " + i);

        state = codon.stringToEncoding("TTA").get(0);
        expectedUsage = new int[]{7,6,9,9,10,10};
        for (int i=0; i < usage.length; i++)
            assertEquals(expectedUsage[i], usage[i][state], "ATC frequency at taxon " + i);
    }
}