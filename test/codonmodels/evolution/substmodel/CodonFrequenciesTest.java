package codonmodels.evolution.substmodel;

import beast.base.evolution.alignment.Alignment;
import codonmodels.CodonFrequencies;
import codonmodels.evolution.CodonData;
import codonmodels.evolution.alignment.CodonAlignment;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.assertArrayEquals;

/**
 * @author Walter Xie
 */
public class CodonFrequenciesTest {

    CodonAlignment codonAlignment;

    @Before
    public void setUp() {
        CodonData.initDataTypes();

        Alignment data = CodonData.getAlig6T333();
        // create Codon Alignment
        codonAlignment = new CodonAlignment();
        codonAlignment.initByName("data", data, "dataType", "codon", "geneticCode", "vertebrateMitochondrial");
    }

    @Test
    public void testCodonPositionBaseFrequencies(){
        double[][] freq = codonAlignment.getObservedBaseFrequencies();
        // result from codeml
        // Codon position * base (3x4) table + overall
        // position 1 : A C G T
        assertArrayEquals(new double[]{.36937,.31331,.15816,.15916}, freq[0], 0.00001);
        // position 2 : A C G T
        assertArrayEquals(new double[]{.18869,.32332,.08458,.4034},freq[1], 0.00001);
        // position 3 : A C G T
        assertArrayEquals(new double[]{.37788,.4044,.04955,.16817}, freq[2], 0.00001);
        // average    : A C G T
        assertArrayEquals(new double[]{.31198,.34701,.09743,.24358}, freq[3], 0.00001);
    }

    @Test
    public void testCodonFrequenciesByUsage(){
        int[][] usage = codonAlignment.getCodonUsage();
        CodonFrequencies codonFreq = new CodonFrequencies(); // default F3X4
        codonFreq.initByName("data", codonAlignment);

        // 60 states, exclude stop codon in freq
        double[] freq = codonFreq.getCodonFrequenciesByUsage(usage);
        // result from codeml
        assertArrayEquals(
                new double[]{.03203203, .03653654, .00500501, .01501502, .03503504, .05705706, 0, .02302302,
                        .01051051, 0, .04604605, .05605606, .00800801, .04504505, .02902903, .02702703,
                        .003003, .0035035, .02152152, .04304304, .003003, .01201201, .00900901, .00850851,
                        0, .001001, .08058058, .03403403, .01601602, .02202202, .00800801, .00900901,
                        .002002, .002002, .01801802, .04904905, .001001, .00800801, .01051051, .01601602,
                        .0015015, .002002, .01401401, .00800801, .003003, .00600601, .01301301, .0035035,
                        .03003003, .01201201, .002002, .00850851, .01851852, .003003, .003003, .001001,
                        .02552553, .02152152, .002002, .01551552 },
                freq, 0.00000001);
    }

}