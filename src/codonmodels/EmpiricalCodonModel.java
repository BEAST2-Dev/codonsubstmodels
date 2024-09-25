/*
 * EmpiricalCodonModel.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package codonmodels;

import java.io.File;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.inference.parameter.RealParameter;

/**
 * Empirical model of codon evolution
 *
 * Adapted from BEAST-X
 * @author Stefan Zoller
 */
@Description("Empirical model of codon evolution")
public class EmpiricalCodonModel extends M0Model {
	
	final public Input<RealParameter> mnmPenaltyInput = new Input<>("mnm-penalty", "multi-nucleotide mutation penalty. "
			+ "If specified, when 2 or 3 mutations are involved the rate is penalised by this parameter. "
			+ "Lower values mean lower likelihood of MNMs. "
			+ "If not specified, all mnm rates are set to 0.");
	
	final public Input<File> matrixDirInput = new Input<>("matrixDir", "director containing files with "
			+ "empirical rate matrix and empirical frequencies. Some examples can be found here: "
			+ "https://github.com/beast-dev/beast-mcmc/tree/master/examples/CodonModels/EmpiricalCodonModels/codon-data/ecmdata", Validate.REQUIRED);
	final public Input<String> rateFileInput = new Input<>("rateFile", "name of comma separated file containing "
			+ "empirical rate matrix", "rates.csv");
	final public Input<String> freqFileInput = new Input<>("freqFile", "name of comma separated file containing "
			+ "empirical frequencies", "freqs.csv");
		
	private RealParameter mnmPenalty;
	private EmpiricalRateMatrixReader rateMat;
	
	private int modelType;
	private final int ECM_OMEGA_2K = 2;
	private final int ECM_OMEGA_9K = 3;
	private final int ECM_OMEGA_NU = 4;
	private final int ECM_OMEGA = 1;


	public EmpiricalCodonModel() {
		kappaInput.setRule(Validate.OPTIONAL);
		omegaInput.setRule(Validate.OPTIONAL);
	}
	
	/**
     * constructors
     *
     * @param codonDataType				Data type as Codons.UNIVERSAL
     * @param omegaParam				Parameter: Omega
     * @param kappaParam				Parameter: Kappa (multidimensional)
     * @param mntParam					Parameter: Multi-nt
     * @param rMat						Initial rate matrix and frequencies
     * @param freqModel					Frequency model
     */	

	@Override
	public void initAndValidate() {
		super.initAndValidate();

		this.geneticCode = codonDataType.getGeneticCode();
		this.mnmPenalty = mnmPenaltyInput.get();
				
		EmpiricalRateMatrixReader rMat = new EmpiricalRateMatrixReader(codonDataType, matrixDirInput.get().getPath(), freqFileInput.get(), rateFileInput.get());
		this.rateMat = rMat;

		constructRateMap();
		
		checkForModelType();
		
		RealParameter freqs = frequencies.frequenciesInput.get();
		if (freqs == null) {
			throw new IllegalArgumentException("Expected a RealParameter in frequencies input, but found null. "
					+ "Please specify a frequencies parameter so it can be initialised from the empirical rates.");
		}
		freqs.setDimension(nrOfStates);
		for (int i = 0; i < nrOfStates; i++) {
			freqs.setValue(rateMat.frequencies[i]);
		}
		
	}
	
	// decide which model to use: ECM_OMEGA_2K, ECM_OMEGA_9K, ECM_OMEGA_NU or ECM_OMEGA
	private void checkForModelType() {
		this.modelType = 0;
		if(kappaInput.get() != null) {
			if(kappaInput.get().getDimension() == 2) {
				this.modelType = ECM_OMEGA_2K;
				Log.info("Using model ECM+omega+2k");
			} else {
				this.modelType = ECM_OMEGA_9K;
				Log.info("Using model ECM+omega+9k");
			}
		}
		if(mnmPenalty != null){
			this.modelType = ECM_OMEGA_NU;
			Log.info("Using model ECM+omega+nu");
		}
		if(kappaInput.get() == null && mnmPenalty == null) {
			this.modelType = ECM_OMEGA;
			Log.info("Using model ECM+omega");
		}
	}
    	
	// setup substitution matrix depending on model type
	public void setupRelativeRates(double[] rates) {
		switch(modelType) {
		case ECM_OMEGA:
			setupRelativeRatesECMOmega(rates);
			break;
		case ECM_OMEGA_2K:
			setupRelativeRatesECMOmega2k(rates);
			break;
		case ECM_OMEGA_9K:
			setupRelativeRatesECMOmega9k(rates);
			break;
		case ECM_OMEGA_NU:
			setupRelativeRatesECMOmegaNu(rates);
			break;
		}
	}
	
	// actual setup routines for different models
	private void setupRelativeRatesECMOmega(double[] rates) {
		double[] initRateMatrix = rateMat.getRates();
		double omega = getOmega();
		for (int i = 0; i < rateCount; i++) {
			switch (rateMap[i]) {
			case 1:											// 1ts, 0tv, syn
			case 3:											// 0ts, 1tv, syn
			case 5:											// 2ts, 0tv, syn
			case 7:											// 1ts, 1tv, syn
			case 9:											// 0ts, 2tv, syn
			case 11:										// 3ts, 0tv, syn
			case 13:										// 2ts, 1tv, syn
			case 15:										// 1ts, 2tv, syn
			case 17: rates[i] = initRateMatrix[i]; break;	// 0ts, 3tv, syn
			
			case 2:													// 1ts, 0tv, nonsyn
			case 4:													// 0ts, 1tv, nonsyn
			case 6:													// 2ts, 0tv, nonsyn
			case 8:													// 1ts, 1tv, nonsyn
			case 10:												// 0ts, 2tv, nonsyn
			case 12:												// 3ts, 0tv, nonsyn
			case 14:												// 2ts, 1tv, nonsyn
			case 16:												// 1ts, 2tv, nonsyn
			case 18: rates[i] = initRateMatrix[i] * omega; break;	// 0ts, 3tv, nonsyn
			}
		}
	}
	
	private void setupRelativeRatesECMOmega2k(double[] rates) {
		double[] initRateMatrix = rateMat.getRates();
		double omega = getOmega();
		double kts = getKappaTs();
		double ktv = getKappaTv();
		for (int i = 0; i < rateCount; i++) {
			switch (rateMap[i]) {
				case 1: rates[i] = initRateMatrix[i] * kts; break;						// 1ts, 0tv, syn
				case 2: rates[i] = initRateMatrix[i] * kts * omega; break;				// 1ts, 0tv, nonsyn
				case 3: rates[i] = initRateMatrix[i] * ktv; break;						// 0ts, 1tv, syn
				case 4: rates[i] = initRateMatrix[i] * ktv * omega; break;				// 0ts, 1tv, nonsyn
				case 5: rates[i] = initRateMatrix[i] * kts * kts; break;				// 2ts, 0tv, syn
				case 6: rates[i] = initRateMatrix[i] * kts * kts * omega; break;		// 2ts, 0tv, nonsyn
				case 7: rates[i] = initRateMatrix[i] * kts * ktv; break;				// 1ts, 1tv, syn
				case 8: rates[i] = initRateMatrix[i] * kts * ktv * omega; break;		// 1ts, 1tv, nonsyn
				case 9: rates[i] = initRateMatrix[i] * ktv * ktv; break;				// 0ts, 2tv, syn
				case 10: rates[i] = initRateMatrix[i] * ktv * ktv * omega; break;		// 0ts, 2tv, nonsyn
				case 11: rates[i] = initRateMatrix[i] * kts * kts * kts; break;			// 3ts, 0tv, syn
				case 12: rates[i] = initRateMatrix[i] * kts * kts * kts * omega; break;	// 3ts, 0tv, nonsyn
				case 13: rates[i] = initRateMatrix[i] * kts * kts * ktv; break;			// 2ts, 1tv, syn
				case 14: rates[i] = initRateMatrix[i] * kts * kts * ktv * omega; break;	// 2ts, 1tv, nonsyn
				case 15: rates[i] = initRateMatrix[i] * kts * ktv * ktv; break;			// 1ts, 2tv, syn
				case 16: rates[i] = initRateMatrix[i] * kts * ktv * ktv * omega; break;	// 1ts, 2tv, nonsyn
				case 17: rates[i] = initRateMatrix[i] * ktv * ktv * ktv; break;			// 0ts, 3tv, syn
				case 18: rates[i] = initRateMatrix[i] * ktv * ktv * ktv * omega; break;	// 0ts, 3tv, nonsyn
			}
		}
	}
	
	private void setupRelativeRatesECMOmega9k(double[] rates) {
		double[] initRateMatrix = rateMat.getRates();
		double omega = getOmega();
		double[] kappa = getKappa();
		for (int i = 0; i < rateCount; i++) {
			switch (rateMap[i]) {
				case 1: rates[i] = initRateMatrix[i] * kappa[0]; break;						// 1ts, 0tv, syn
				case 2: rates[i] = initRateMatrix[i] * kappa[0] * omega; break;				// 1ts, 0tv, nonsyn
				case 3: rates[i] = initRateMatrix[i] * kappa[1]; break;						// 0ts, 1tv, syn
				case 4: rates[i] = initRateMatrix[i] * kappa[1] * omega; break;				// 0ts, 1tv, nonsyn
				case 5: rates[i] = initRateMatrix[i] * kappa[2]; break;				// 2ts, 0tv, syn
				case 6: rates[i] = initRateMatrix[i] * kappa[2] * omega; break;		// 2ts, 0tv, nonsyn
				case 7: rates[i] = initRateMatrix[i] * kappa[3]; break;				// 1ts, 1tv, syn
				case 8: rates[i] = initRateMatrix[i] * kappa[3] * omega; break;		// 1ts, 1tv, nonsyn
				case 9: rates[i] = initRateMatrix[i] * kappa[4]; break;				// 0ts, 2tv, syn
				case 10: rates[i] = initRateMatrix[i] * kappa[4] * omega; break;		// 0ts, 2tv, nonsyn
				case 11: rates[i] = initRateMatrix[i] * kappa[5]; break;			// 3ts, 0tv, syn
				case 12: rates[i] = initRateMatrix[i] * kappa[5] * omega; break;	// 3ts, 0tv, nonsyn
				case 13: rates[i] = initRateMatrix[i] * kappa[6]; break;			// 2ts, 1tv, syn
				case 14: rates[i] = initRateMatrix[i] * kappa[6] * omega; break;	// 2ts, 1tv, nonsyn
				case 15: rates[i] = initRateMatrix[i] * kappa[7]; break;			// 1ts, 2tv, syn
				case 16: rates[i] = initRateMatrix[i] * kappa[7] * omega; break;	// 1ts, 2tv, nonsyn
				case 17: rates[i] = initRateMatrix[i] * kappa[8]; break;			// 0ts, 3tv, syn
				case 18: rates[i] = initRateMatrix[i] * kappa[8] * omega; break;	// 0ts, 3tv, nonsyn
			}
		}
	}
	
	private void setupRelativeRatesECMOmegaNu(double[] rates) {
		double[] initRateMatrix = rateMat.getRates();
		double omega = getOmega();
		double mnt = getMultiNt();
		for (int i = 0; i < rateCount; i++) {
			switch (rateMap[i]) {
			case 1:														// 1ts, 0tv, syn
			case 3: rates[i] = initRateMatrix[i]; break;				// 0ts, 1tv, syn
			case 2: 													// 1ts, 0tv, nonsyn
			case 4: rates[i] = initRateMatrix[i] * omega; break;		// 0ts, 1tv, nonsyn
			
			case 5:														// 2ts, 0tv, syn
			case 7: 													// 1ts, 1tv, syn
			case 9: 													// 0ts, 2tv, syn
			case 11: 													// 3ts, 0tv, syn
			case 13: 													// 2ts, 1tv, syn
			case 15: 													// 1ts, 2tv, syn
			case 17: rates[i] = initRateMatrix[i] * mnt; break;			// 0ts, 3tv, syn
			
			case 6: 													// 2ts, 0tv, nonsyn
			case 8: 													// 1ts, 1tv, nonsyn
			case 10: 													// 0ts, 2tv, nonsyn
			case 12: 													// 3ts, 0tv, nonsyn
			case 14: 													// 2ts, 1tv, nonsyn
			case 16: 													// 1ts, 2tv, nonsyn
			case 18: rates[i] = initRateMatrix[i] * mnt * omega; break;	// 0ts, 3tv, nonsyn
			}
		}
	}
	
    // getters for parameters

    public double getOmega() { return omegaInput.get().getValue(0); }
    
       
    public double getKappaTs() { 
    	if(kappaInput.get() != null) {
    		return kappaInput.get().getValue(0);
    	} else {
    		return 0.0;
    	}
    }
    
    public double getKappaTv() { 
    	if(kappaInput.get() != null) {
    		return kappaInput.get().getValue(1);
    	} else {
    		return 0.0;
    	}
    }
    
    public double[] getKappa() {
    	if(kappaInput.get() != null) {
    		return kappaInput.get().getDoubleValues();
    	} else {
    		return new double[9];
    	}
    }
    

    public double getMultiNt() { 
    	if(mnmPenalty != null) {
    		return mnmPenalty.getValue(0);
    	} else {
    		return 0.0;
    	}
	}
    
    
    /**
	 * Construct a map of the rate classes in the rate matrix using the current
	 * genetic code. Classes are:
	 * 	1-2: 1ts, 0tv (syn/nonsyn)
	 *  3-4: 0ts, 1tv
	 *  5-6: 2ts, 0tv
	 *  7-8: 1ts, 1tv
	 *  9-10: 0ts, 2tv
	 *  11-12: 3ts, 0tv
	 *  13-14: 2ts, 1tv
	 *  15-16: 1ts, 2tv
	 *  17-18: 0ts, 3tv
	 */
	protected void constructRateMap()
	{
		int u, v, i1, j1, k1, i2, j2, k2, ts, tv, non;
		byte rateClass;
		int[] codon;
		// int cs1, cs2;
		int aa1, aa2;

		int i = 0;

		rateMap = new short[rateCount];

		for (u = 0; u < nrOfStates; u++) {

			codon = getCodonStatesForRateClass(u, codonDataType);
			i1 = codon[0];
			j1 = codon[1];
			k1 = codon[2];

			aa1 = codon[3]; // geneticCode.getAminoAcidState(codonDataType.getCanonicalState(cs1));

			for (v = u + 1; v < nrOfStates; v++) {
				
				ts = 0;
				tv = 0;
				non = 0;
				rateClass = -1;

				codon = getCodonStatesForRateClass(v, codonDataType);
				i2 = codon[0];
				j2 = codon[1];
				k2 = codon[2];

				aa2 = codon[3]; 

				if (i1 != i2) {
					if ( (i1 == 0 && i2 == 2) || (i1 == 2 && i2 == 0) || // A <-> G
						 (i1 == 1 && i2 == 3) || (i1 == 3 && i2 == 1) ) { // C <-> T
						ts += 1; // Transition at position 1
					} else {
						tv += 1; // Transversion at position 1
					}
				}
				if (j1 != j2) {
					if ( (j1 == 0 && j2 == 2) || (j1 == 2 && j2 == 0) || // A <-> G
						(j1 == 1 && j2 == 3) || (j1 == 3 && j2 == 1) ) { // C <-> T
						ts += 1; // Transition
					} else {
						tv += 1; // Transversion
					}
				}
				if (k1 != k2) {
					if ( (k1 == 0 && k2 == 2) || (k1 == 2 && k2 == 0) || // A <-> G
						(k1 == 1 && k2 == 3) || (k1 == 3 && k2 == 1) ) { // C <-> T
						ts += 1; // Transition
					} else {
						tv += 1; // Transversion
					}
				}

	 			if (aa1 != aa2) {
					non = 1; // Is a non-synonymous change
				}

	 			// decide for rateClass
	 			switch(ts) {
	 				case 0:
	 					switch(tv) {
	 						case 1: rateClass = 3; break;	// 0ts, 1tv
							case 2: rateClass = 9; break;	// 0ts, 2tv
							case 3: rateClass = 17; break;	// 0ts, 3tv
							default: break;
	 					}
	 					break;
	 				case 1:
	 					switch(tv) {
	 						case 0: rateClass = 1; break;	// 1ts, 0tv
							case 1: rateClass = 7; break;	// 1ts, 1tv
							case 2: rateClass = 15; break;	// 1ts, 2tv
							default: break;
	 					}
	 					break;
	 				case 2:
	 					switch(tv) {
 							case 0: rateClass = 5; break;	// 2ts, 0tv
 							case 1: rateClass = 13; break;	// 2ts, 1tv
 							default: break;
	 					}
	 					break;
	 				case 3:
	 					rateClass = 11; break;	// 3ts, 0tv
	 				default: break;
	 			}
	 			
	 			if(non == 1) {
	 				rateClass += 1;
	 			}
				rateMap[i] = rateClass;
				i++;
			}
		}
	}

}
