package codonmodels;
import java.util.Arrays;

import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.datatype.DataType;
import beast.evolution.substitutionmodel.EigenDecomposition;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;

@Description("Substitution model that approximates transition probability matrix of another model through splines")
public class ApproximateSubstModel extends SubstitutionModel.Base {
	final public Input<SubstitutionModel.Base> substModelInput = new Input<>("substModel", "substitution model we want to approximate", Validate.REQUIRED);
	
	boolean needsUpdate = true;
	SubstitutionModel.Base substModel;

	// temporary arrays for SplineInterpolator
	double [] x = new double[] {0.0, 1.0E-4, 0.10010000000000001, 0.30010000000000003, 0.5001, 0.7001, 0.9000999999999999, 1.2001000000000002, 1.5001000000000004, 1.9001000000000008, 2.400100000000001, 3.1001000000000016, 4.100100000000002, 5.600099999999997, 7.60009999999999, 9.900099999999982};
	double [] y;
	SplineInterpolator splineInterpolator;
	PolynomialFunction [][] functions;
	double [][] matrix;
	
	int stateCount;

	public ApproximateSubstModel() {
		frequenciesInput.setRule(Validate.OPTIONAL);
	}
	
	@Override
	public void initAndValidate() {
		substModel = substModelInput.get();
		frequenciesInput.setValue(substModel.frequenciesInput.get(), this);
		super.initAndValidate();
		needsUpdate = true;
		y = new double[x.length];
		splineInterpolator = new SplineInterpolator();
		stateCount = substModel.getStateCount();
		matrix = new double[x.length][stateCount * stateCount];
		functions = new PolynomialFunction[x.length][stateCount * stateCount];
	}
	
	@Override
	public void getTransitionProbabilities(Node node, double startTime, double endTime, double rate, double[] matrix) {
		if (needsUpdate) {
			update();
		}
		
        double distance = (startTime - endTime) * rate;
        
        if (distance > x[x.length-1]) {
        	// range of interpolation exceeded, take last matrix (containing base frequencies presumably)
        	System.arraycopy(this.matrix[x.length-1], 0, matrix, 0, stateCount * stateCount);
        }
        
        int k = Arrays.binarySearch(x, distance);
        if (k < 0) {
            k = -k - 2;
        }
        
		for (int i = 0; i < stateCount; i++) {
			for (int j = 0; j < stateCount; j++) {
				matrix[i * stateCount + j] = functions[k][i * stateCount + j].value(distance);
			}
		}

		// TODO: normalise matrix for numerical issues?
		
	}

	private void update() {
		// gather prob matrices
		for (int k = 0; k < x.length; k++) {
			substModel.getTransitionProbabilities(null, x[k], 0, 1.0, matrix[k]);
		}
		// TODO: make sure last matrix contains base frequencies?
		
		for (int i = 0; i < stateCount; i++) {
			for (int j = 0; j < stateCount; j++) {
				for (int k = 0; k < x.length; k++) {
					y[k] = matrix[k][i * stateCount + j];
				}
				PolynomialSplineFunction f = splineInterpolator.interpolate(x, y);
				PolynomialFunction [] fs = f.getPolynomials();
				for (int k = 0; k < x.length; k++) {
					functions[k][i * stateCount + j] = fs[k];
				}
			}			
		}
		needsUpdate = false;		
	}

	@Override
	public EigenDecomposition getEigenDecomposition(Node node) {
		return null;
	}

	@Override
	public boolean canHandleDataType(DataType dataType) {
		return substModel.canHandleDataType(dataType);
	}


	
	@Override
	protected void store() {
		super.store();
	}
	
	@Override
	protected void restore() {
		needsUpdate = true;
	}
	
	@Override
	protected boolean requiresRecalculation() {
		if (substModel.isDirtyCalculation()) {
			needsUpdate = true;
		}
		return super.requiresRecalculation();
	}
}
