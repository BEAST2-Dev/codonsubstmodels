package codonmodels;

import java.util.Arrays;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.datatype.DataType;
import beast.evolution.substitutionmodel.EigenDecomposition;
import beast.evolution.tree.Node;

@Description("Substitution model that approximates transition probability matrix of another model through splines")
public class ApproximateSubstModel extends CodonSubstitutionModel {
	final public Input<CodonSubstitutionModel> substModelInput = new Input<>("substModel", "substitution model we want to approximate", Validate.REQUIRED);
	
	boolean needsUpdate = true;
	CodonSubstitutionModel substModel;
	double [] iexp;

	// temporary arrays for SplineInterpolator
	// knots
	double [] x = new double[] {0.0, 1.0E-4, 0.10010000000000001, 0.30010000000000003, 0.5001, 0.7001, 0.9000999999999999, 1.2001000000000002, 1.5001000000000004, 1.9001000000000008, 2.400100000000001, 3.1001000000000016, 4.100100000000002, 5.600099999999997, 7.60009999999999, 9.900099999999982};
	// function values
	double [] y;
    // Differences between knot points
	double [] h;
    double mu[];
    double z[];
    // cubic spline coefficients --  b is linear, c quadratic, d is cubic (original y's are constants)
    double b[];
    double c[];
    double d[];
    // Number of intervals.  The number of data points is n + 1.
	int n;
	
	//PolynomialFunction [][] functions;
	double [][][] coefficients;
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
		stateCount = substModel.getStateCount();

        n = x.length - 1;
		matrix = new double[n+1][stateCount * stateCount];
		//functions = new PolynomialFunction[n+1][stateCount * stateCount];
		coefficients = new double[n+1][stateCount * stateCount][4];
        h = new double[n];
        for (int i = 0; i < n; i++) {
            h[i] = x[i + 1] - x[i];
        }
        mu = new double[n];
        z = new double[n + 1];
        b = new double[n];
        c = new double[n + 1];
        d = new double[n];
        
        iexp = new double[stateCount * stateCount];
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
        
        double [][] coeffs = coefficients[k];
		for (int i = 0; i < stateCount; i++) {
			for (int j = 0; j < stateCount; j++) {
				matrix[i * stateCount + j] = evaluate(coeffs[i * stateCount + j], distance);
			}
		}

		// TODO: normalise matrix for numerical issues?
		
	}

	private void update() {
		// gather prob matrices
		for (int k = 0; k < x.length; k++) {
//			substModel.getTransitionProbabilities(null, x[k], 0, 1.0, matrix[k]);
			substModel.getTransiProbs(x[k], 0, 1.0, iexp, matrix[k]);
		}
		// TODO: make sure last matrix contains base frequencies?
		
		for (int i = 0; i < stateCount; i++) {
			for (int j = 0; j < stateCount; j++) {
				for (int k = 0; k < x.length; k++) {
					y[k] = matrix[k][i * stateCount + j];
				}
				//PolynomialSplineFunction f = 
						interpolate(x, y, h, i, j);
				//PolynomialFunction [] fs = f.getPolynomials();
				//for (int k = 0; k < x.length; k++) {
				//	functions[k][i * stateCount + j] = fs[k];
				//}
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
	public void store() {
		super.store();
	}
	
	@Override
	public void restore() {
		needsUpdate = true;
	}
	
	@Override
	protected boolean requiresRecalculation() {
		if (substModel.isDirtyCalculation()) {
			needsUpdate = true;
		}        // Number of intervals.  The number of data points is n + 1.
        final int n = x.length - 1;

		return super.requiresRecalculation();
	}

    /**
     * adapted from apache commons library
     * 
     * Computes an interpolating function for the data set.
     * @param x the arguments for the interpolation points
     * @param y the values for the interpolation points
     * @param h intervals sizes x[i+1]-x[i]
     * @param i1, i2 indices for storing coefficients
     */
    public void interpolate(double x[], double y[], double [] h, int i1, int i2) {
//        if (x.length != y.length) {
//            throw new DimensionMismatchException(x.length, y.length);
//        }
//
//        if (x.length < 3) {
//            throw new NumberIsTooSmallException(LocalizedFormats.NUMBER_OF_POINTS,
//                                                x.length, 3, true);
//        }

//        MathArrays.checkOrder(x);

    	// resetting to 0 might not be necessary
        Arrays.fill(mu, 0);
        Arrays.fill(z, 0);
        Arrays.fill(b, 0);
        Arrays.fill(c, 0);
        Arrays.fill(d, 0);
        
        
        mu[0] = 0d;
        z[0] = 0d;
        double g = 0;
        for (int i = 1; i < n; i++) {
            g = 2d * (x[i+1]  - x[i - 1]) - h[i - 1] * mu[i -1];
            mu[i] = h[i] / g;
            z[i] = (3d * (y[i + 1] * h[i - 1] - y[i] * (x[i + 1] - x[i - 1])+ y[i - 1] * h[i]) /
                    (h[i - 1] * h[i]) - h[i - 1] * z[i - 1]) / g;
        }


        z[n] = 0d;
        c[n] = 0d;

        for (int j = n -1; j >=0; j--) {
            c[j] = z[j] - mu[j] * c[j + 1];
            b[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (c[j + 1] + 2d * c[j]) / 3d;
            d[j] = (c[j + 1] - c[j]) / (3d * h[j]);
        }

        //final PolynomialFunction polynomials[] = new PolynomialFunction[n];
        for (int i = 0; i < n; i++) {
            final double coefficients[] = this.coefficients[i][i1 * stateCount + i2];
            coefficients[0] = y[i];
            coefficients[1] = b[i];
            coefficients[2] = c[i];
            coefficients[3] = d[i];
            //polynomials[i] = new PolynomialFunction(coefficients);
        }

        //return new PolynomialSplineFunction(x, polynomials);
    }
    
    /**
     * @param coefficients Coefficients of the polynomial to evaluate.
     * @param argument Input value.
     * @return the value of the polynomial.
     */
    private static double evaluate(double[] coefficients, double argument) {
        double result = ((coefficients[3] * argument  + coefficients[2]) * argument + coefficients[1]) * argument + coefficients[0];
        return result;
    }


}
