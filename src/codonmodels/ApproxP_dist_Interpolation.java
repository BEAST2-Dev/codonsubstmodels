package codonmodels;

import beast.core.util.Log;
import beast.util.RandomUtils;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.interpolation.UnivariateInterpolator;

/**
 * Interpolation
 * @author Walter Xie
 */
public class ApproxP_dist_Interpolation extends ApproxP_dist_Piecewise {

    // knots
    final double [] X = new double[] {0.0, 1.0E-4, 0.1, 0.3, 0.5, 0.7, 0.9, 1.2,
            1.5, 1.9, 2.4, 3.1, 4.1, 5.6, 7.6, 9.9,
            13, 17, 22, 27, 35, 45, 55, 60};

    // 3600
    UnivariateFunction[] functions; // replace double[][] p_d_

    @Override
    public void initAndValidate() {
        // learn the knots
//        super.initAndValidate();

//        double[][] xyMinMax = getXYMinMax();
//        printXY(xyMinMax);
//        printXYMax(xyMinMax);

        double[] freq = initCodonSubstInput();

        double xMax = knots[knots.length-1];
        if (X[X.length-1] < xMax) {
            knots = new double[X.length+8];
            System.arraycopy(X, 0 , knots, 0, X.length);
            for (int i = 0; i < knots.length-X.length; i++)
                knots[X.length+i] = X[X.length-1] + (xMax - X[X.length-1]) / 8.0 * (i + 1);
        } else {
            knots = X;
        }

        p_d_ = new double[knots.length][prob.length];
        assert p_d_[0].length == nrOfStates * nrOfStates;

        functions = new UnivariateFunction[prob.length];
        createUnivariateFunctions();

        Log.info.println("Spline approximation of codon substitution model " + codonSubstModel.getID() +
                ", creating " + knots.length + " data points.\n");

    }

    @Override
    public void getTransiProbs(double distance, double[] matrix) {
        // > biggest distance
        int i = knots.length-1;
        if (distance >= knots[i]) {
            System.arraycopy(p_d_[i], 0 , matrix, 0, matrix.length);
            return;
        }

        // intervals[i-1] <= distance <= intervals[i]
        i = RandomUtils.binarySearchSampling(knots, distance);
        if (distance == knots[i]) {
            System.arraycopy(p_d_[i], 0 , matrix, 0, matrix.length);
        } else {
            // approximation
            for (int j = 0; j < p_d_[i].length; j++) {
                matrix[j] = functions[j].value(distance);
            }
        }

    }


    protected void createUnivariateFunctions() {

//        List<Double> intervalList = new ArrayList<>();
        UnivariateInterpolator interpolator = new SplineInterpolator();

        for (int i = 0; i < knots.length; i++) {
            // true y
            codonSubstModel.getTransiProbs(knots[i], iexp, prob);
            // collect y
            System.arraycopy(prob, 0 , p_d_[i], 0, prob.length);
        }

        double[] y = new double[knots.length];
        for (int j = 0; j < prob.length; j++) {
            for (int i = 0; i < knots.length; i++)
               y[i] = p_d_[i][j];

            functions[j] = interpolator.interpolate(knots, y);
        }


    }



    // use plotPt.R to observe data space
    public static void main(final String[] args) {
        String omega = "0.08", kappa = "15";
//        String[] omegas = new String[]{"0.01","0.08","0.1","1"};
//        String[] kappas = new String[]{"1","5","15","30"};
//        for (String omega : omegas) {
//            for (String kappa : kappas) {
        M0Model m0Model = getM0Model(omega, kappa);

        ApproxP_dist_Interpolation pd = new ApproxP_dist_Interpolation();
//                pd.initByName("substModel", m0Model, "file", "p_d_" + omega + "_" + kappa + ".txt");
        pd.initByName("substModel", m0Model);
//        pd.printP_d_();
        pd.testAccuracy(450);
//            }
//        }
    }


}
