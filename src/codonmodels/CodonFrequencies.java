
package codonmodels;

import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import beast.base.core.Log;
import codonmodels.evolution.alignment.CodonAlignment;
import codonmodels.evolution.datatype.Codon;
import codonmodels.evolution.datatype.GeneticCode;
import beast.base.evolution.substitutionmodel.Frequencies;

import java.text.DecimalFormat;
import java.util.Arrays;


@Description("The equilibrium codon frequencies PI, including equal, F1X4, F3X4, and F60/F61.")
@Citation(value="Yang, Z. (2006). Computational molecular evolution. Oxford University Press.", year = 2006)
public class CodonFrequencies extends Frequencies {

    final public Input<String> piInput = new Input<>("pi", "The assumption of equilibrium codon frequencies PI, " +
            "including equal, F1X4, F3X4 (default), and F60/F61", "F3X4");

    final public Input<Boolean> verboseInput = new Input<>("verbose",
            "Print the codon frequencies, etc.",
            Boolean.TRUE);

    protected CodonAlignment codonAlignment;
    protected Codon codonDataType; // contain selected genetic code
    protected int codonStateCount;
//    protected GeneticCode geneticCode;
//    int stopCodonCount;

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
//        assert codonStateCount == 64;

        GeneticCode geneticCode = codonDataType.getGeneticCode();
//        stopCodonCount = geneticCode.getStopCodonCount();
        Log.info.println("Codon alignment " + (codonAlignment.getID()==null?"":codonAlignment.getID()) +
                " : number of states = " + codonStateCount + ", excluding " +
                geneticCode.getStopCodonCount() + " stop codon");

        //======  Observe data
        if (verboseInput.get()) {
            int[][] usage = codonAlignment.getCodonUsage();
            // freqs.length = 64
            double[] freqs = getCodonFrequenciesByUsage(usage);
            printCodonFrequencies(freqs, "Observed codon frequencies in alignment");
            Log.info.println();
        }

        //====== this includes update()
        super.initAndValidate();

        freqs = getFreqs();
        Log.info.println("Use " + piInput.get() + " equilibrium codon frequencies. ");
        Log.info.println("Set frequencies dimension = " + freqs.length);
//      Log.info.println("Initial values to codon frequencies = " +
//                Arrays.toString(StringUtils.roundDoubleArrays(freqs, 8)));

        // non-stop code freq has to > 0.
        for (int s=0; s < freqs.length; s++) {
            if ( freqs[s] <= 0 )
// TODO should throw new IllegalArgumentException?
            Log.warning.println("Zero frequency at non-stop code " +
                    codonDataType.encodingToString(new int[]{s}) + " ! freqs[" + s + "]=" + freqs[s]);
        }

        if (verboseInput.get())
            printCodonFrequencies(freqs, "Codon frequencies used in CodonSubstitutionModel");
        Log.info.println();
    }

    @Override
    protected void update() {

        if (frequenciesInput.get() != null) {
            // 60/61 codon frequencies (separated by white space) are in fixed order AAA AAC AAG AAT ... TTA TTC TTG TTT
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
                double equalFreq = 1.0 / (double) codonStateCount;
                //  Log.info.println("equalFreq = " + equalFreq);
                Arrays.fill(freqs, equalFreq);

            } else if ("F1X4".equals(piInput.get())) {
                // observed freq from data
                double[][] freqsCPB = codonAlignment.getObservedBaseFrequencies();
                // last row freqs[3] is used to estimate
                freqs = estimateFrequencies(freqsCPB[3], freqsCPB[3], freqsCPB[3]);

            } else if ("F3X4".equals(piInput.get())) { // default
                // observed freq from data
                double[][] freqsCPB = codonAlignment.getObservedBaseFrequencies();
                freqs = estimateFrequencies(freqsCPB[0], freqsCPB[1], freqsCPB[2]);

            } else if ("F60/F61".equals(piInput.get())) {
                int[][] usage = codonAlignment.getCodonUsage();
                freqs = getCodonFrequenciesByUsage(usage);

            } else {
                throw new IllegalArgumentException("Invalid input of pi = " + piInput.get());
            }

        } // end frequenciesInput.get() != null
        if (freqs.length != codonStateCount)
            throw new IllegalArgumentException("Invalid codon frequencies, dimension " +
                    freqs.length + " != " + codonStateCount);
        needsUpdate = false;
    }

    // F1X4 (nucFreq1==nucFreq2=nucFreq3) or F3X4, nuc order: A,C,G,T
    protected double[] estimateFrequencies(double[] nucFreq1, double[] nucFreq2, double[] nucFreq3) {
//        GeneticCode geneticCode = codonDataType.getGeneticCode();

        double[] freqs = new double[codonStateCount];
        double sum = 0;
        for (int i = 0; i < codonStateCount; i++) {
//            if (! geneticCode.isStopCodonIndex(i)) {
                // convert codon state i into 3 nuc states
                int[] necStates = codonDataType.getTripletNucStates(i);
                for (int t = 0; t < necStates.length; t++) {
                    if (necStates[t] > 3)
                        throw new IllegalArgumentException("Invalid nucleotide state ! " + necStates[t]);
                }

                freqs[i] = nucFreq1[necStates[0]] * nucFreq2[necStates[1]] * nucFreq3[necStates[2]];
                sum += freqs[i];
//            }
        }
        // re-normalize
        for (int i = 0; i < freqs.length; i++)
            freqs[i] = freqs[i] / sum;
        return freqs;
    }


    /**
     * Codon frequencies from codon usage (AAA AAC AAG AAT ... TTT), excluding stop codon.
     * @param usage i is taxon, j is state. Not include totals.
     *              If ambiguous states exit, then their counts are equally distributed to
     *              possible unambiguous states using {@link Codon#getStatesForCode(int)}.
     * @return 1d frequency array freqs[60]
     */
    public double[] getCodonFrequenciesByUsage(int[][] usage) {
        final int stateMax = codonStateCount; // 60/61
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
