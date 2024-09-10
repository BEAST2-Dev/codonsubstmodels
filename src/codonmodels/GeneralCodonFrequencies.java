
package codonmodels;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.inference.parameter.RealParameter;
import codonmodels.evolution.alignment.CodonAlignment;
import codonmodels.evolution.datatype.Codon;

import java.util.ArrayList;
import java.util.List;


@Description("Frequencies for codon models that cna be estimated")
public class GeneralCodonFrequencies extends Frequencies {

	final public Input<List<RealParameter>> freqInput = new Input<>("freq", 
			"determines base frequencies:"
			+ "single freq with dim=4: F1x4 "
			+ "three freq with dim=3: F3x4 "
			+ "single freq with dim=60 or 61: F60/F61 "
			+ "", new ArrayList<>()); 

    protected Codon codonDataType; // contain selected genetic code
    protected int codonStateCount;

    @Override
    public void initAndValidate() {
    	codonDataType = CodonAlignment.toCodonAlignment(dataInput.get()).getDataType();
        codonStateCount = codonDataType.getStateCount();
   }

    @Override
    protected void update() {
        // 60/61 codon frequencies (separated by white space) are in fixed order AAA AAC AAG AAT ... TTA TTC TTG TTT
        RealParameter freqsInput = frequenciesInput.get();

        if (freqsInput != null) {
        	if (freqsInput.getDimension() == 4) {
        		double [] freqs0 = freqsInput.getDoubleValues();
                freqs = estimateFrequencies(freqs0, freqs0, freqs0);
        	} else {
	            if (freqsInput.getDimension() != codonStateCount)
	                throw new IllegalArgumentException("frequencies input has " + freqsInput.getDimension() +
	                        " entries, but state count = " + codonStateCount + " !" );
	
	            freqs = new double[freqsInput.getDimension()];
	            double sum = 0;
	            for (int i = 0; i < freqs.length; i++) {
	                freqs[i] = freqsInput.getValue(i);
	                sum += freqs[i];
	            }
	            if (Math.round(sum * 1.0e6) / 1.0e6 != 1) // see also Frequencies.initAndValidate
	                throw new IllegalArgumentException("The codon frequencies do not sum up to 1 ! " + sum);
        	}

        } else {
        	List<RealParameter> freqList = freqInput.get();
        	switch (freqList.size()) {
        	case 1:
        		double [] freqs0 = freqList.get(0).getDoubleValues();
                freqs = estimateFrequencies(freqs0, freqs0, freqs0);
                break;
        	case 3:
        		double [] freqs0_ = freqList.get(0).getDoubleValues();
        		double [] freqs1_ = freqList.get(1).getDoubleValues();
        		double [] freqs2_ = freqList.get(2).getDoubleValues();
                freqs = estimateFrequencies(freqs0_, freqs1_, freqs2_);
                break;
        	default:
        		throw new IllegalArgumentException("Expected 1 or 4 frequencies, not " + freqList.size());
        	}
        }
        	
        if (freqs.length != codonStateCount)
            throw new IllegalArgumentException("Invalid codon frequencies, dimension " +
                    freqs.length + " != " + codonStateCount);
        needsUpdate = false;
    }

    // F1X4 (nucFreq1==nucFreq2=nucFreq3) or F3X4, nuc order: A,C,G,T
    protected double[] estimateFrequencies(double[] nucFreq1, double[] nucFreq2, double[] nucFreq3) {

        double[] freqs = new double[codonStateCount];
        double sum = 0;
        for (int i = 0; i < codonStateCount; i++) {
                int[] nucStates = codonDataType.getTripletNucStates(i);
                for (int t = 0; t < nucStates.length; t++) {
                    if (nucStates[t] > 3) {
                        throw new IllegalArgumentException("Invalid nucleotide state ! " + nucStates[t]);
                    }
                }

                freqs[i] = nucFreq1[nucStates[0]] * nucFreq2[nucStates[1]] * nucFreq3[nucStates[2]];
                sum += freqs[i];
        }
        for (int i = 0; i < freqs.length; i++) {
            freqs[i] = freqs[i] / sum;
        }
        return freqs;
    }

} // class Frequencies
