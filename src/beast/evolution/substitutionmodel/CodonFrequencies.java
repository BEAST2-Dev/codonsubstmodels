
package beast.evolution.substitutionmodel;

import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.CodonAlignment;
import beast.evolution.datatype.Codon;
import beast.evolution.datatype.GeneticCode;
import beast.util.StringUtils;

import java.util.Arrays;


@Description("The equilibrium codon frequencies PI, including equal, F1X4, F3X4, and F6n.")
@Citation(value="Yang, Z. (2006). Computational molecular evolution. Oxford University Press.", year = 2006)
public class CodonFrequencies extends Frequencies {

    final public Input<String> piInput = new Input<>("pi", "The assumption of equilibrium codon frequencies PI, " +
            "including equal, F1X4, F3X4 (default), and F6n", "F3X4");

    public CodonFrequencies() {
        dataInput.setRule(Input.Validate.REQUIRED); // only use internally
    }

    @Override
    public void initAndValidate() {
        estimateInput.setValue(false, this);

        CodonAlignment codonAlignment = getCodonAlignment();
        Codon codonDataType = getDataType(codonAlignment);
        int codonStateCount = codonDataType.getStateCount();
        Log.info.println("Codon alignment " + codonAlignment.getID() + " : number of states = " + codonStateCount);

        int[] indices = codonAlignment.getGeneticCode().getStopCodonIndices();
        if (indices.length >= codonStateCount || indices.length == 0)
            throw new IllegalArgumentException("Invalid stop codon indices length " + indices.length);

        double[] freqs;
        if ("equal".equals(piInput.get())) {
            freqs = new double[codonStateCount];
            double equalFreq =  1.0 / (double) (codonStateCount - indices.length);
            Arrays.fill(freqs, equalFreq);
            // replace to 0 for stop codon
            for (int id : indices)
                freqs[id] = 0.0;

        } else if ("F1X4".equals(piInput.get())) {

            double[][] freqsCPB = codonAlignment.getCodonPositionBaseFrequencies();
            freqs = estimateFrequencies(freqsCPB[3], freqsCPB[3], freqsCPB[3], codonDataType);

        } else if ("F3X4".equals(piInput.get())) {

            double[][] freqsCPB = codonAlignment.getCodonPositionBaseFrequencies();
            freqs = estimateFrequencies(freqsCPB[0], freqsCPB[1], freqsCPB[2], codonDataType);

        } else if ("F6n".equals(piInput.get())) {
            freqs = codonAlignment.getCodonFrequencies();

        } else {
            throw new IllegalArgumentException("Invalid input of pi = " + piInput.get());
        }
        // convert double[] to Double[]
        if (freqs.length != codonStateCount)
            throw new IllegalArgumentException("Invalid codon frequencies, dimension " +
                    freqs.length + " != " + codonStateCount);
        Double[] values = new Double[codonStateCount];
        for (int i = 0; i < freqs.length; i++)
            values[i] = freqs[i];

        Log.info.println("Use " + piInput.get() + " equilibrium codon frequencies. ");
        Log.info.println("Initial values to CodonFrequencies = " +
                Arrays.toString(StringUtils.roundDoubleArrays(values, 5)));

        RealParameter frequencies = new RealParameter(values);
        frequencies.lowerValueInput.setValue(0.0, frequencies);
        frequencies.upperValueInput.setValue(1.0, frequencies);
        frequenciesInput.setValue(frequencies, this);

        super.initAndValidate();

        Log.info.println("Set frequencies dimension = " + getFreqs().length);
    }

    // F1X4 (nucFreq1==nucFreq2=nucFreq3) or F3X4, nuc order: A,C,G,T
    protected double[] estimateFrequencies(double[] nucFreq1, double[] nucFreq2, double[] nucFreq3,
                                           Codon codonDataType) {
        GeneticCode geneticCode = codonDataType.getGeneticCode();
        int codonStateCount = codonDataType.getStateCount();
        double[] freqs = new double[codonStateCount];
        double sum = 0;
        for (int i = 0; i < codonStateCount; i++) {
            if (! geneticCode.isStopCodon(i)) {
                int[] tripletStates = codonDataType.getTripletStates(i);
                for (int t = 0; t < tripletStates.length; t++) {
                    if (tripletStates[t] > 3)
                        throw new IllegalArgumentException("Invalid nucleotide state ! " + tripletStates[t]);
                }

                freqs[i] = nucFreq1[tripletStates[0]] * nucFreq2[tripletStates[1]] * nucFreq3[tripletStates[2]];
                sum += freqs[i];
            }
        }
        // re-normalize
        for (int i = 0; i < codonStateCount; i++)
            freqs[i] = freqs[i] / sum;
        return freqs;
    }

    public CodonAlignment getCodonAlignment() {
        Alignment alignment = dataInput.get();
        if (! (alignment instanceof CodonAlignment) )
            throw new IllegalArgumentException("Codon alignment is required by CodonFrequencies !");
        return (CodonAlignment) alignment;
    }

    public Codon getDataType(CodonAlignment codonAlignment) {
        return codonAlignment.getDataType();
    }

} // class Frequencies
