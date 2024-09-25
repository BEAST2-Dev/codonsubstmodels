package codonmodels;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.datatype.Aminoacid;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;

@Description("General codon model that satisfies rate matrix from aminoacid model")
public class EmpiricalAACodonModel extends GeneralCodonSubstitutionModel {

	final public Input<GeneralSubstitutionModel> aaModelInput = new Input<>("aaModel", "amino acid substitution model (with frequencies) used to constraint codon rates", Validate.REQUIRED);

	private GeneralSubstitutionModel aaModel;
	private double [] codonRates;
	private Map<String, List<Integer>> map;
	private String aminoAcidCodeMap;
	private boolean needsToConvertAArates;
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();
		
		aaModel = aaModelInput.get();
		codonRates = new double[category.length];
		
		map = new HashMap<>();
		aminoAcidCodeMap = new Aminoacid().getCodeMap();

		for (int i = 0; i < fromAA.length; i++) {
			String a1 = aminoAcidCodeMap.charAt(fromAA[i]) + "" +  aminoAcidCodeMap.charAt(toAA[i]);
			if (!map.containsKey(a1)) {
				map.put(a1, new ArrayList<>());
			}
			map.get(a1).add(i);
		}
		
		needsToConvertAArates = true;
	}
	
	// ensure 
	private void convertAAratesToCodonRates() {
		double [] frequencies = getFrequencies();
		aaModel.setupRelativeRates();
		aaModel.setupRateMatrix();
		double [][] aaRates = aaModel.getRateMatrix();
		for (int i = 0; i < 20; i++) {
			char from = aminoAcidCodeMap.charAt(i);
			for (int j = 0; j < 20; j++) {
				if (i != j) {
					char to = aminoAcidCodeMap.charAt(j);
					List<Integer> range = map.get(from + "" + to);
					
					double codonFreqs = 0;
					for (int k : range) {
						// TODO: make sure this is the correct frequencies index
						codonFreqs += frequencies[k % nrOfStates];
					}
					
					double codonRate = aaRates[i][j] / codonFreqs;
					for (int k : range) {
						codonRates[k] = codonRate;
					}
				}
			}
		}
		needsToConvertAArates = false;
	}
	
    @Override
	public void setupRelativeRates() {
    	super.setupRelativeRates();
    	
    	if (needsToConvertAArates) {
    		convertAAratesToCodonRates();
    	}
    	
    	for (int i = 0; i < codonRates.length; i++) {
    		relativeRates[i] *= codonRates[i];
    	}
    }
	
	
	@Override
	protected boolean requiresRecalculation() {
		if (aaModel.isDirtyCalculation()) {
			needsToConvertAArates = true;
		}
		return super.requiresRecalculation();
	}
	
}
