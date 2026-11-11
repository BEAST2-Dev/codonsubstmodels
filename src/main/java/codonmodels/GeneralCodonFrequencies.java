
package codonmodels;

import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.alignment.Alignment;
import beast.base.spec.evolution.likelihood.GenericTreeLikelihood;
import beast.base.spec.evolution.substitutionmodel.Frequencies;
import beast.base.spec.type.Simplex;
import beast.base.inference.CalculationNode;
import beast.base.inference.StateNode;
import beast.base.inference.util.InputUtil;
import codonmodels.evolution.alignment.CodonAlignment;
import codonmodels.evolution.datatype.Codon;

import java.util.ArrayList;
import java.util.List;


@Description("Frequencies for codon models that cna be estimated")
public class GeneralCodonFrequencies extends Frequencies {

	final public Input<List<Simplex>> freqInput = new Input<>("freq", 
			"determines base frequencies:"
			+ "single freq with dim=4: F1x4 "
			+ "three freq with dim=3: F3x4 "
			+ "single freq with dim=60 or 61: F60/F61 "
			+ "", new ArrayList<>()); 

    protected Codon codonDataType; // contain selected genetic code
    protected int codonStateCount;

    
    public GeneralCodonFrequencies() {
        frequenciesInput.setRule(Input.Validate.OPTIONAL);
        dataInput.setRule(Input.Validate.OPTIONAL); 
    }
    
    @Override
    public void initAndValidate() {
        Alignment data = null;
        if (dataInput.get() != null && dataInput.get() instanceof CodonAlignment) {
        	data = dataInput.get();
        } else {
        	data = getData(this);
        }
        
     	codonDataType = CodonAlignment.toCodonAlignment(data).getDataType();
        codonStateCount = codonDataType.getStateCount();
        needsUpdate = true;
   }

    private Alignment getData(BEASTInterface o2) {
		if (o2 instanceof GenericTreeLikelihood) {
			GenericTreeLikelihood tl = (GenericTreeLikelihood) o2;
			return tl.dataInput.get();
		}
		
    	for (BEASTInterface o : o2.getOutputs()) {
    		Alignment data = getData(o);
    		if (data != null) {
    			return data;
    		}
    	}
    	return null;
	}

	@Override
    protected void update() {
        // 60/61 codon frequencies (separated by white space) are in fixed order AAA AAC AAG AAT ... TTA TTC TTG TTT
        Simplex freqsInput = frequenciesInput.get();

        if (freqsInput != null) {
        	if (freqsInput.size() == 4) {
        		double [] freqs0 = toDoubleArray(freqsInput);
                freqs = estimateFrequencies(freqs0, freqs0, freqs0);
        	} else {
	            if (freqsInput.size() != codonStateCount)
	                throw new IllegalArgumentException("frequencies input has " + freqsInput.size() +
	                        " entries, but state count = " + codonStateCount + " !" );
	
	            freqs = new double[freqsInput.size()];
	            double sum = 0;
	            for (int i = 0; i < freqs.length; i++) {
	                freqs[i] = freqsInput.get(i);
	                sum += freqs[i];
	            }
	            if (Math.round(sum * 1.0e6) / 1.0e6 != 1) // see also Frequencies.initAndValidate
	                throw new IllegalArgumentException("The codon frequencies do not sum up to 1 ! " + sum);
        	}

        } else {
        	List<Simplex> freqList = freqInput.get();
        	switch (freqList.size()) {
        	case 1:
        		double [] freqs0 = toDoubleArray(freqList.get(0));
                freqs = estimateFrequencies(freqs0, freqs0, freqs0);
                break;
        	case 3:
        		double [] freqs0_ = toDoubleArray(freqList.get(0));
        		double [] freqs1_ = toDoubleArray(freqList.get(1));
        		double [] freqs2_ = toDoubleArray(freqList.get(2));
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

    private static double[] toDoubleArray(Simplex simplex) {
        double[] values = new double[simplex.size()];
        for (int i = 0; i < values.length; i++) {
            values[i] = simplex.get(i);
        }
        return values;
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
    
    @Override
    protected boolean requiresRecalculation() {
        if (InputUtil.isDirty(frequenciesInput)) {
            needsUpdate = true;
        }
    	if (InputUtil.isDirty(freqInput)) {
			needsUpdate = true;
    	}
    	return needsUpdate;
    }

    public Codon getDataType() {
    	return codonDataType;
    }
} // class Frequencies
