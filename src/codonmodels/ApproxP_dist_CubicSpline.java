package codonmodels;

import beast.base.core.Log;
import codonmodels.util.RandomUtils;

import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Interpolation
 * @author Walter Xie
 */
public class ApproxP_dist_CubicSpline extends ApproxP_dist_Piecewise {

    // knots TODO Natural Spline
    final double [] X = new double[] {0.0, 1.0E-4, 0.1, 0.3, 0.5, 0.7, 0.9, 1.2,
            1.5, 1.9, 2.4, 3.1, 4.1, 5.6, 7.6, 9.9,   13, 17, 22, 27, 35, 45, 55, 60};

    // 3600
    PolynomialSplineFunction[] cubicSpline; // replace double[][] p_d_

    @Override
    public void initAndValidate() {
        // learn the knots
//        super.initAndValidate();

//        double[][] xyMinMax = getXYMinMax();
//        printXY(xyMinMax);
//        printXYMax(xyMinMax);

        double[] freq = initCodonSubstInput();

        double maxDist = getMaxDist(freq, 1E-3);

        createKnots(maxDist);

        // cache y, TODO use UnivariateFunction.value(double)
        p_d_ = new double[knots.length][prob.length];
        assert p_d_[0].length == nrOfStates * nrOfStates;

        cubicSpline = new PolynomialSplineFunction[prob.length];
        createCubicSplineFunctions();

        Log.info.println("Spline approximation of codon substitution model " + codonSubstModel.getID() +
                ", creating " + knots.length + " data points.");
        Log.info.println(Arrays.toString(knots));

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
            int a,b,j;
            double sum;
            for (a = 0; a < nrOfStates; a++) {
                sum = 0;
                for (b = 0; b < nrOfStates; b++) {
                    j = a * nrOfStates + b;
                    matrix[j] = cubicSpline[j].value(distance);
                    sum += matrix[j];
                }
                // normalise TODO it seems not improve much.
                for (b = 0; b < nrOfStates; b++) {
                    matrix[a * nrOfStates + b] /= sum;
                }
            }
        }

    }

//TODO different knots of each curve?
    protected void createKnots(double maxDist) {
        if (X[X.length-1] < maxDist) { // add extra 8 points uniformly distributed
            int add = 8;
            knots = new double[X.length + add];
            System.arraycopy(X, 0 , knots, 0, X.length);
            for (int i = 0; i < knots.length-X.length; i++)
                knots[X.length+i] = X[X.length-1] + (maxDist - X[X.length-1]) / add * (i + 1);

        } else { // add X[i] until > maxDist
            List<Double> knotPointsList = new ArrayList<>();
            int i = 0;
            while (X[i] <= maxDist) {
                knotPointsList.add(X[i]);
                i++;
            }
            knots = knotPointsList.stream().mapToDouble(j -> j).toArray();
        }
    }

    protected void createCubicSplineFunctions() {

//        List<Double> intervalList = new ArrayList<>();
        SplineInterpolator interpolator = new SplineInterpolator();

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

            cubicSpline[j] = interpolator.interpolate(knots, y);

//            Log.info.println(Arrays.toString(cubicSpline[j].getPolynomials()));
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

        ApproxP_dist_CubicSpline pd = new ApproxP_dist_CubicSpline();
//                pd.initByName("substModel", m0Model, "file", "p_d_" + omega + "_" + kappa + ".txt");
        pd.initByName("substModel", m0Model);
//        pd.printP_d_();
        pd.testAccuracy(450, 0.01); //STEP
//            }
//        }
    }


    protected double getMaxDist(double[] freq, final double diff) {
        final double step = 10.0;
        double dist = 10.0;
        // backwards search to 10
        while (dist < MaxDistance) {
            // eigen decomp to get points, startTime > endTime
            codonSubstModel.getTransiProbs(dist, iexp, prob);

            boolean isDiff = false;
            for (int c = 0; c < prob.length; c++) {
                int i = 0;
                if (Math.abs(prob[c] - freq[i]) > diff) {
                    isDiff = true;
                    break;
                }
                i++;
            }
            if (!isDiff)
                return dist;

            dist += step;
        }

        return dist;
    }

}
