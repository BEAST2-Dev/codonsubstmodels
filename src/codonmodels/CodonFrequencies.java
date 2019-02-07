
package codonmodels;

import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.CodonAlignment;
import beast.evolution.datatype.Codon;
import beast.evolution.datatype.GeneticCode;
import beast.evolution.substitutionmodel.Frequencies;
import beast.util.StringUtils;

import java.util.Arrays;


@Description("The equilibrium codon frequencies PI, including equal, F1X4, F3X4, and F60/F61.")
@Citation(value="Yang, Z. (2006). Computational molecular evolution. Oxford University Press.", year = 2006)
public class CodonFrequencies extends Frequencies {

    final public Input<String> piInput = new Input<>("pi", "The assumption of equilibrium codon frequencies PI, " +
            "including equal, F1X4, F3X4 (default), and F60/F61", "F3X4");

    public CodonFrequencies() {
        estimateInput.setValue(false, this);
//        estimateInput.setRule(Input.Validate.FORBIDDEN);
//        frequenciesInput.setRule(Input.Validate.OPTIONAL);
//        dataInput.setRule(Input.Validate.REQUIRED); // only use internally
    }

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        Log.info.println("Set frequencies dimension = " + getFreqs().length);
    }

    @Override
    protected void update() {
        if (frequenciesInput.get() != null) {

            // if user specified, parse frequencies from space delimited string
            freqs = new double[frequenciesInput.get().getDimension()];

            for (int i = 0; i < freqs.length; i++) {
                freqs[i] = frequenciesInput.get().getValue(i);
            }

        } else {
            CodonAlignment codonAlignment = getCodonAlignment();
            Codon codonDataType = getDataType(codonAlignment);
            int codonStateCount = codonDataType.getStateCount();
            int stopCodonCount = codonDataType.getGeneticCode().getStopCodonCount();
            Log.info.println("Codon alignment " + (codonAlignment.getID()==null?"":codonAlignment.getID()) +
                    " : number of states = " + codonStateCount + ", including " + stopCodonCount + " stop codon");

            // codon states exclude stop codon, freqs should be stateCount
            if ("equal".equals(piInput.get())) {
                freqs = new double[codonStateCount];
                double equalFreq = 1.0 / (double) (codonStateCount - stopCodonCount);
                Arrays.fill(freqs, equalFreq);
//                // replace to 0 for stop codon
//                for (int id : indices)
//                    freqs[id] = 0.0;

            } else if ("F1X4".equals(piInput.get())) {

                double[][] freqsCPB = codonAlignment.getCodonPositionBaseFrequencies();
                freqs = estimateFrequencies(freqsCPB[3], freqsCPB[3], freqsCPB[3], codonDataType);

            } else if ("F3X4".equals(piInput.get())) {

                double[][] freqsCPB = codonAlignment.getCodonPositionBaseFrequencies();
                freqs = estimateFrequencies(freqsCPB[0], freqsCPB[1], freqsCPB[2], codonDataType);

            } else if ("F60/F61".equals(piInput.get())) {
                freqs = codonAlignment.getCodonFrequencies();

            } else {
                throw new IllegalArgumentException("Invalid input of pi = " + piInput.get());
            }
            // convert double[] to Double[]
            if (freqs.length != codonStateCount)
                throw new IllegalArgumentException("Invalid codon frequencies, dimension " +
                        freqs.length + " != " + codonStateCount);
//        Double[] values = new Double[codonStateCount];
//        for (int i = 0; i < freqs.length; i++)
//            values[i] = freqs[i];

            Log.info.println("Use " + piInput.get() + " equilibrium codon frequencies. ");
            Log.info.println("Initial values to codon frequencies = " +
                    Arrays.toString(StringUtils.roundDoubleArrays(freqs, 5)));

//        RealParameter frequencies = new RealParameter(values);
//        frequencies.lowerValueInput.setValue(0.0, frequencies);
//        frequencies.upperValueInput.setValue(1.0, frequencies);
//        frequenciesInput.setValue(frequencies, this);
        }
        needsUpdate = false;
    }

    // F1X4 (nucFreq1==nucFreq2=nucFreq3) or F3X4, nuc order: A,C,G,T
    protected double[] estimateFrequencies(double[] nucFreq1, double[] nucFreq2, double[] nucFreq3,
                                           Codon codonDataType) {
        GeneticCode geneticCode = codonDataType.getGeneticCode();
        int stateCount = codonDataType.getStateCount(); // 64

        assert stateCount==64;

        double[] freqs = new double[stateCount];
        double sum = 0;
        for (int i = 0; i < geneticCode.getCodeTableLength(); i++) {
            if (! geneticCode.isStopCodon(i)) {
                // convert codon state i into 3 nuc states
                int[] necStates = codonDataType.getTripletNucStates(i);
                for (int t = 0; t < necStates.length; t++) {
                    if (necStates[t] > 3)
                        throw new IllegalArgumentException("Invalid nucleotide state ! " + necStates[t]);
                }

                freqs[i] = nucFreq1[necStates[0]] * nucFreq2[necStates[1]] * nucFreq3[necStates[2]];
                sum += freqs[i];
            }
        }
        // re-normalize
        for (int i = 0; i < freqs.length; i++)
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
