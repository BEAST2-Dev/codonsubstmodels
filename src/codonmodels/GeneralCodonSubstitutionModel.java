
package codonmodels;


import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.datatype.Aminoacid;
import beast.base.evolution.datatype.DataType;
import codonmodels.evolution.datatype.Codon;
import codonmodels.evolution.datatype.GeneticCode;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.inference.parameter.RealParameter;

import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.List;

@Description("General codon substitution model "
		+ "based on nucleotide substitution models + "
		+ "non/synonymous model")
public class GeneralCodonSubstitutionModel extends GeneralSubstitutionModel {
	
	final public Input<List<GeneralSubstitutionModel>> substmodelInput = new Input<>("substModel", 
			"nucleotide subsitution model, either 1 shared among codon positions, "
			+ "or 3 seperate for each codon position", new ArrayList<>());
	
	final public Input<RealParameter> mnmPenaltyInput = new Input<>("mnmPenalty", "multi-nucleotide mutation pentalty (multiplied with rate) "
			+ "Set to zero if no more than 1 mutation is allowed (as in the M0 model). "
			+ "if dimension=2: multiplier penalty for 2 mutations + multiplier penalty for 3 mutations."
			+ "if dimension=1: multiplier penalties for both 2 and 3 mutations", Validate.REQUIRED);
	
	final public Input<RealParameter> substRateInput = new Input<>("substitutionRate","substitution rate for the three codon positions, if specified");
	
	final public Input<RealParameter> omegaInput = new Input<>("omega",
            "omega parameter to represent the nonsynonymous-synonymous rate ratio", Input.Validate.REQUIRED);

	final public Input<String> categoriesInput = new Input<>("categories",
            "comma separated list of amino acids falling in the same category. "
            + "If not specified, assume all amino acids are in their own category. ");

	
	protected int [] mutationCount;
	// nucleotide mutation for rate i, for codon position j. If negative, no mutation, otherwise index into relative rates matrix 
	protected int [] nucMutation1, nucMutation2, nucMutation3;
	protected int [] fromAA, toAA;
	protected int [] category;
	
	protected List<GeneralSubstitutionModel> substModels;
	
    protected GeneticCode geneticCode;
    protected RealParameter substRate;

    public GeneralCodonSubstitutionModel() {
        ratesInput.setRule(Input.Validate.FORBIDDEN); // only use internally
    }

    protected String [] categories;

    @Override
    public void initAndValidate() {
        this.frequencies = frequenciesInput.get();
        if (!(frequencies instanceof GeneralCodonFrequencies)) {
            throw new IllegalArgumentException("GeneralCodonFrequencies is required by GeneralCodonSubstitutionModel");
        }
        Codon codonDataType = ((GeneralCodonFrequencies) frequencies).codonDataType;
        geneticCode = codonDataType.getGeneticCode();

        //====== init states and rates ======
        updateMatrix = true;
        double[] freqs = frequencies.getFreqs();
        nrOfStates = freqs.length;
        if (nrOfStates != geneticCode.getStateCount()) {
        	throw new IllegalArgumentException("state count of genetic code and frequencies do not match");
        }

        try {
            eigenSystem = createEigenSystem();
        } catch (SecurityException | ClassNotFoundException | InstantiationException |
                IllegalAccessException | IllegalArgumentException | InvocationTargetException e) {
            throw new IllegalArgumentException(e.getMessage());
        }

        rateMatrix = new double[nrOfStates][nrOfStates];

        int rateCount = nrOfStates * (nrOfStates - 1);
        if (rateCount < 1) {
            throw new RuntimeException("Incorrect rate count = " + rateCount);
        }
        relativeRates = new double[rateCount];
        storedRelativeRates = new double[rateCount];

        if (categoriesInput.get() != null) {
        	categories = categoriesInput.get().split(",");
        } else {
        	categories = new String[] {};
        }

        constructStateMap(rateCount, nrOfStates, codonDataType);
        
        substModels = substmodelInput.get();
        
        substRate = substRateInput.get();
        if (substRate == null) {
        	substRate = new RealParameter("1.0");
        }
    
        if (omegaInput.get().getDimension() != categories.length + 1) {
        	omegaInput.get().setDimension(categories.length + 1);
        }
    }

    /**
     * Faster code to replace
     * {@link GeneralSubstitutionModel#getTransitionProbabilities(Node, double, double, double, double[], boolean)}.
     *
     * @param startTime parent.getHeight()
     * @param endTime   node.getHeight()
     * @param rate      joint rate = rate for a site category * mean branch rate.
     * @param iexp      iexp, without creating a new matrix each call.
     * @param matrix    P(t), without creating a new matrix each call.
     */
    public void getTransiProbs(double startTime, double endTime, double rate,
                               double[] iexp, double[] matrix) {//, boolean normalized) {
        double distance = (startTime - endTime) * rate;

        getTransiProbs(distance, iexp, matrix);
    }

    /**
     * Faster code to replace
     * {@link GeneralSubstitutionModel#getTransitionProbabilities(Node, double, double, double, double[], boolean)}.
     * @param distance  distance = (startTime - endTime) * mean branch rate * rate for a site category.
     * @param iexp      iexp, without creating a new matrix each call.
     * @param matrix    P(t), without creating a new matrix each call.
     */
    public void getTransiProbs(double distance, double[] iexp, double[] matrix) {
        int i, j, k;
        double temp;
        // this must be synchronized to avoid being called simultaneously by
        // two different likelihood threads - AJD
        synchronized (this) {
            if (updateMatrix) {
                setupRelativeRates();
                setupRateMatrix();
                eigenDecomposition = eigenSystem.decomposeMatrix(rateMatrix);
                updateMatrix = false;
            }
        }

        double[] Evec = eigenDecomposition.getEigenVectors();
        // inverse Eigen vectors
        double[] Ievc = eigenDecomposition.getInverseEigenVectors();
        // Eigen values
        double[] Eval = eigenDecomposition.getEigenValues();

        // faster computation, reviewed by AJD
        int x = 0;
        for (i = 0; i < nrOfStates; i++) {
            temp = Math.exp(distance * Eval[i]);
            for (j = 0; j < nrOfStates; j++) {
                // iexp[i * nrOfStates + j] = Ievc[i * nrOfStates + j] * temp;
                iexp[x] = Ievc[x] * temp;
                x++; // save time, plus once
            }
        }

        x = 0;
        int y;
        int u = 0;
        for (i = 0; i < nrOfStates; i++) {
            for (j = 0; j < nrOfStates; j++) {
                y = j;
                temp = 0.0;
                for (k = 0; k < nrOfStates; k++) {
                    temp += Evec[x + k] * iexp[y];
                    y += nrOfStates; // y = k + nrOfStates
                }

                matrix[u] = Math.abs(temp);
                u++; // u = u + j
            }
            x += nrOfStates; // x = i + nrOfStates
        }
    }


   
    protected void constructStateMap(int rateCount, int stateCount, Codon codon)	{
        mutationCount = new int[rateCount];
        nucMutation1 = new int[rateCount]; 
        nucMutation2 = new int[rateCount]; 
        nucMutation3 = new int[rateCount]; 
        fromAA = new int[rateCount];
        toAA = new int[rateCount];
        category = new int[rateCount];
        

        // this needs to match rateMatrix[i][j] <= relativeRates[] in setupRateMatrix()
        // i j is codonState
        for (int i = 0; i < stateCount; i++) {

            // i1, j1, k1, aa1
            int[] ids1 = getCodonStatesForRateClass(i, codon);

            int offset = 0;
            for (int j = 0; j < stateCount; j++) {
            	if (i != j) {
	                // i2, j2, k2, aa2
	                int[] ids2 = getCodonStatesForRateClass(j, codon);
	
	                int i1 = ids1[0];
	                int j1 = ids1[1];
	                int k1 = ids1[2];
	                
	                int i2 = ids2[0];
	                int j2 = ids2[1];
	                int k2 = ids2[2];
	                
	                int mutations = (i1!=i2?1:0) + (j1!=j2?1:0) + (k1!=k2?1:0); 
	
	                int index = i * (stateCount - 1) + j - offset;
	                mutationCount[index] = mutations;
	                nucMutation1[index] = nucMutation(i1, i2);
	                nucMutation2[index] = nucMutation(j1, j2);
	                nucMutation3[index] = nucMutation(k1, k2);

	                
//	                fromNucMutation[index][0] = i1;
//	                toNucMutation[index][0] = i2;
//	                fromNucMutation[index][1] = j1;
//	                toNucMutation[index][1] = j2;
//	                fromNucMutation[index][2] = k1;
//	                toNucMutation[index][2] = k2;
	                
	                fromAA[index] = ids1[3];
	                toAA[index] = ids2[3];
	                
	                category[index] = getCategory(ids1[3], ids2[3]);
            	} else {
            		offset = 1;
            	}

            }
        }
    }

    private int getCategory(int i, int j) {
		String aa =  new Aminoacid().getCodeMap();
		String a1 = aa.charAt(i) + "";
		String a2 = aa.charAt(j) + "";
        if (i == j) { 
        	return -1;
        } else {
    		for (int k = 0; k < categories.length; k++) {
    			if (categories[k].contains(a1) && categories[k].contains(a2)) {
    				return k;
    			}
    		}
        }
		return categories.length;
	}

	private int nucMutation(int i1, int i2) {
        if (i1 == i2) {
        	return -1;
        } else {
        	return i1 * 3 + i2 + (i2 > i1 ? -1 : 0);
        }
	}

	// return i1, j1, k1, aa1, given codonState
    public int[] getCodonStatesForRateClass(int codonState, Codon codon) {
        int i1, j1, k1, cs1, aa1;
        int[] nucStates  = codon.getTripletNucStates(codonState);
        i1 = nucStates[0];
        j1 = nucStates[1];
        k1 = nucStates[2];

        cs1 = codon.getCodonState(i1, j1, k1);

        // GeneticCode geneticCode = codon.getGeneticCode();
        aa1 = getGeneticCode().getAminoAcidState(cs1);
        return new int[]{i1, j1, k1, aa1};
    }



    @Override
    public boolean canHandleDataType(DataType dataType) {
        return dataType instanceof Codon;
    }

	public GeneticCode getGeneticCode() {
		return geneticCode;
	}


    @Override
	public void setupRelativeRates() {
    	int rateCount = mutationCount.length;
        RealParameter omega = omegaInput.get();
        Function mnmPenalty = mnmPenaltyInput.get();
        double mnmPenalty2 = mnmPenalty == null ? 0.0 : mnmPenalty.getArrayValue();
        double mnmPenalty3 = mnmPenalty == null ? 0.0 : 
        	(mnmPenalty.getDimension() > 1 ? mnmPenalty.getArrayValue(1) : mnmPenalty2);

        substModels.get(0).setupRelativeRates();
        if (substModels.size() > 1) {
            substModels.get(1).setupRelativeRates();
            substModels.get(2).setupRelativeRates();
        }
        
        double [] r1 = substModels.get(0).getRelativeRates();
        double [] r2 = substModels.size() > 1 ? substModels.get(1).getRelativeRates() : r1;
        double [] r3 = substModels.size() > 1 ? substModels.get(2).getRelativeRates() : r1;
        
        if (substRate.getDimension() == 3) {
	        double s = substRate.getValue(0);
	        for (int i = 0; i < r1.length; i++) {
	        	r1[i] *= s;
	        }
	        s = substRate.getValue(1);
	        for (int i = 0; i < r2.length; i++) {
	        	r2[i] *= s;
	        }
	        s = substRate.getValue(2);
	        for (int i = 0; i < r3.length; i++) {
	        	r3[i] *= s;
	        }
        }
        
        for (int i = 0; i < rateCount; i++) {
        	double rate = aaRate(i);

        	if (mutationCount[i] == 2) {
        		rate *= mnmPenalty2;
        	} else if (mutationCount[i] == 3) {
        		rate *= mnmPenalty3;
        	}
        	
        	if (nucMutation1[i] >= 0) {
        		rate *= r1[nucMutation1[i]];
        	}
        	if (nucMutation2[i] >= 0) {
        		rate *= r2[nucMutation2[i]];
        	}
        	if (nucMutation3[i] >= 0) {
        		rate *= r3[nucMutation3[i]];
        	}
//        	if (fromNucMutation[i][0] != toNucMutation[i][0]) {
//        		rate *= r1[fromNucMutation[i][0]][toNucMutation[i][0]];
//        	}
//        	if (fromNucMutation[i][1] != toNucMutation[i][1]) {
//        		rate *= r2[fromNucMutation[i][1]][toNucMutation[i][1]];
//        	}
//        	if (fromNucMutation[i][2] != toNucMutation[i][2]) {
//        		rate *= r3[fromNucMutation[i][2]][toNucMutation[i][2]];
//        	}
        	
            relativeRates[i] = rate;
        }
    }

    /** contribution of amino acid mutation to rate **/
	protected double aaRate(int i) {
		int index = category[i];
		if (index < 0) {
			return 1.0;
		}
		return omegaInput.get().getArrayValue(index);
	}
}