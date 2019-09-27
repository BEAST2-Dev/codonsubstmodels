
package codonmodels;

import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.alignment.CodonAlignment;
import beast.evolution.datatype.Codon;
import beast.evolution.datatype.GeneticCode;
import beast.evolution.substitutionmodel.Frequencies;

import java.text.DecimalFormat;


@Description("The equilibrium codon frequencies PI, including equal, F1X4, F3X4, and F60/F61.")
@Citation(value="Yang, Z. (2006). Computational molecular evolution. Oxford University Press.", year = 2006)
public class CodonFrequencies extends Frequencies {

    final public Input<String> piInput = new Input<>("pi", "The assumption of equilibrium codon frequencies PI, " +
            "including equal, F1X4, F3X4 (default), and F60/F61", "F3X4");

    final public Input<Boolean> verboseInput = new Input<>("verbose",
            "Print the codon frequencies, etc.",
            Boolean.TRUE);

    CodonAlignment codonAlignment;
    Codon codonDataType; // contain selected genetic code
    int codonStateCount;
    GeneticCode geneticCode;
    int stopCodonCount;

    public CodonFrequencies() {
        estimateInput.setValue(false, this);
//        estimateInput.setRule(Input.Validate.FORBIDDEN);
        // change freqs rule from XOR to OPT
        frequenciesInput.setRule(Input.Validate.OPTIONAL);
        dataInput.setRule(Input.Validate.REQUIRED); // only use internally
    }

    @Override
    public void initAndValidate() {
        this.codonAlignment = CodonAlignment.toCodonAlignment(dataInput.get());
        this.codonDataType = codonAlignment.getDataType();
        codonStateCount = codonDataType.getStateCount();
        assert codonStateCount == 64;

        geneticCode = codonDataType.getGeneticCode();
        stopCodonCount = geneticCode.getStopCodonCount();
        Log.info.println("Codon alignment " + (codonAlignment.getID()==null?"":codonAlignment.getID()) +
                " : number of states = " + codonStateCount + ", including " + stopCodonCount + " stop codon");

        //======  Observe data
        if (verboseInput.get()) {
            int[][] usage = codonAlignment.getCodonUsage();
            // freqs.length = 64
            double[] freqs = getCodonFrequenciesByUsage(usage);
            printCodonFrequencies(freqs, "Observed codon frequencies by usage");
            Log.info.println();
        }

        //====== this includes update()
        super.initAndValidate();

        if (verboseInput.get()) {
            freqs = getFreqs();
            Log.info.println();
            Log.info.println("Use " + piInput.get() + " equilibrium codon frequencies. ");
            Log.info.println("Set frequencies dimension = " + freqs.length);

//        Log.info.println("Initial values to codon frequencies = " +
//                Arrays.toString(StringUtils.roundDoubleArrays(freqs, 8)));
            printCodonFrequencies(freqs, "Codon frequencies passed to CodonSubstitutionModel");
            Log.info.println();
        }
    }

    @Override
    protected void update() {

        if (frequenciesInput.get() != null) {
            // 64 codon frequencies (separated by white space) are in fixed order AAA AAC AAG AAT ... TTA TTC TTG TTT
            RealParameter freqsInput = frequenciesInput.get();
            assert freqsInput.getDimension() == codonStateCount;

            freqs = new double[freqsInput.getDimension()];
            double sum = 0;
            for (int i = 0; i < freqs.length; i++) {
                freqs[i] = freqsInput.getValue(i);
                sum += freqs[i];
            }
            if (Math.round(sum * 1.0e6) / 1.0e6 != 1) // see also Frequencies.initAndValidate
                throw new IllegalArgumentException("The codon frequencies do not sum up to 1 ! " + sum);

        } else {

            // codon states exclude stop codon, freqs should be stateCount
            if ("equal".equals(piInput.get())) {
                freqs = new double[codonStateCount];
                double equalFreq = 1.0 / (double) (codonStateCount - stopCodonCount);
                //  Log.info.println("equalFreq = " + equalFreq);

                for (int i = 0; i < freqs.length; i++) {
                    if (! geneticCode.isStopCodon(i))
                        freqs[i] = equalFreq;
                }

            } else if ("F1X4".equals(piInput.get())) {

                double[][] freqsCPB = codonAlignment.getCodonPositionBaseFrequencies();
                freqs = estimateFrequencies(freqsCPB[3], freqsCPB[3], freqsCPB[3]);

            } else if ("F3X4".equals(piInput.get())) {

                double[][] freqsCPB = codonAlignment.getCodonPositionBaseFrequencies();
                freqs = estimateFrequencies(freqsCPB[0], freqsCPB[1], freqsCPB[2]);

            } else if ("F60/F61".equals(piInput.get())) {
                int[][] usage = codonAlignment.getCodonUsage();
                freqs = getCodonFrequenciesByUsage(usage);

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

//        RealParameter frequencies = new RealParameter(values);
//        frequencies.lowerValueInput.setValue(0.0, frequencies);
//        frequencies.upperValueInput.setValue(1.0, frequencies);
//        frequenciesInput.setValue(frequencies, this);
        } // end frequenciesInput.get() != null
        needsUpdate = false;
    }

    // F1X4 (nucFreq1==nucFreq2=nucFreq3) or F3X4, nuc order: A,C,G,T
    protected double[] estimateFrequencies(double[] nucFreq1, double[] nucFreq2, double[] nucFreq3) {
        GeneticCode geneticCode = codonDataType.getGeneticCode();

        double[] freqs = new double[codonStateCount];
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


    /**
     * Codon frequencies from codon usage (AAA AAC AAG AAT ... TTT), excluding stop codon.
     * @param usage i is taxon, j is state. Not include totals, and for j cols > 63 are ambiguous states count.
     * @return 1d frequency array freqs[64]
     */
    public double[] getCodonFrequenciesByUsage(int[][] usage) {
        final int stateMax = codonStateCount; // 64
//        assert stateMax == 64;

        // include stop codon, same as codeml
        double[] freqs = new double[stateMax];
        double sum = 0;
        // loop through each taxon
        for (int i = 0; i < usage.length; i++) {
            // j is state
            for (int j = 0; j < usage[0].length; j++) {
                // j = [0, 63] no ambiguous
                if (j < stateMax) {
                    freqs[j] += usage[i][j];

                } else if (usage[i][j] > 0) { // j > 63 are ambiguous
                    // all non-ambiguous states for this ambiguous state
                    int[] states = codonDataType.getStatesForCode(j);
                    // equally distribute ambiguous into the count of each of non-stop-codon state
                    for (int s : states)
                        freqs[s] += (double) usage[i][j] / (double) states.length;

                } // ignore 0 usage
                sum += usage[i][j];
            }
        }

        if (sum == 0)
            throw new IllegalArgumentException("Invalid codon usage, the total is 0 !");

        // re-normalize
        for (int j = 0; j < freqs.length; j++)
            freqs[j] = freqs[j] / sum;
        return freqs;
    }


    /**
     * Print codon frequencies in the fixed order (AAA AAC AAG AAT ... TTT)
     * @param frequencies
     * @param title
     */
    protected void printCodonFrequencies(double[] frequencies, String title) {
        Log.info.println("\n============ " + title + " (AAA AAC AAG AAT ... TTA TTC TTG TTT) ============");
        DecimalFormat df = new DecimalFormat("#");
        df.setMaximumFractionDigits(8);
        for (int i = 0; i < frequencies.length; i++) {
//            int state = getDataType().getStatesForCode(i)[0];
//            Log.info.println(i + "  " + state);
            if (i % 8 == 0) {
                Log.info.print("\n" + df.format(frequencies[i]));
            } else {
                Log.info.print("\t" + df.format(frequencies[i]));
            }
        }
        Log.info.println();
    }

} // class Frequencies
