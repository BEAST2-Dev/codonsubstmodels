package codonmodels;

import beast.core.util.Log;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.fitting.PolynomialCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoint;
import org.apache.commons.math3.fitting.WeightedObservedPoints;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Interpolation
 * @author Walter Xie
 */
public class ApproxP_dist_Polynomial extends ApproxP_dist_Piecewise {


    // 3600, 2nd[] is quantiles
    PolynomialFunction[][] polynomials; // replace double[][] p_d_
    private final int degree = 3;

    final int breaks = 8;
    double[] frequencies;

    @Override
    public void initAndValidate() {
//        super.initAndValidate();

//        double[][] xyMinMax = getXYMinMax();
//        printXY(xyMinMax);
//        printXYMax(xyMinMax);

        frequencies = initCodonSubstInput();
        // learn the knots
        super.createKnotsPd();
        Log.info.println("Select " + knots.length + " data points from linear approximation.\n");

//        double[][] xyMinMax = getXYMinMax();

        //reassign knots by breaks
        double[] x = new double[knots.length];
        System.arraycopy(knots, 0, x, 0, knots.length);
        createKnots(x);
        // 1st knot is 0, so knots = breaks + 1
        assert knots.length == breaks + 1;
        assert knots[0] == 0.0;

        polynomials = new PolynomialFunction[prob.length][breaks];

        double[] y = new double[x.length];
        for (int i = 0; i < prob.length; i++) {
            for (int j = 0; j < y.length; j++)
                y[j] = p_d_[j][i];
            polynomials[i] = createSegmentedPolynomials(degree, x, y);
        }

        Log.info.println("Polynomial fit " + breaks + " breaks for " +
                "approximation of codon substitution model " + codonSubstModel.getID());
        for (int i = 0; i < 5; i++)
            System.out.println(i + "\t" + Arrays.toString(polynomials[i]));

        // reassign p_d_?
        double[][] p_d_new = new double[knots.length][prob.length];
        int nrOfPoints = x.length / breaks;
        p_d_new[0] = p_d_[0];
        p_d_new[p_d_new.length-1] = p_d_[p_d_.length-1];
        for (int i = 1; i < p_d_new.length-1; i++) {
            p_d_new[i] = p_d_[i * nrOfPoints];
        }
        p_d_ = p_d_new;
    }

    @Override
    public void getTransiProbs(double distance, double[] matrix) {
        // > biggest distance
        int b = knots.length-1;
        if (distance >= knots[b]) { // TODO replace p_d_[knots.length-1] by freq[]?
            System.arraycopy(p_d_[b], 0 , matrix, 0, matrix.length);
            return;
        }
        if (distance == 0.0) {
            System.arraycopy(p_d_[0], 0 , matrix, 0, matrix.length);
            return;
        }

        // intervals[i-1] <= distance <= intervals[i]
        for (b = 1; b < knots.length; b++) {
//        i = RandomUtils.binarySearchSampling(knots, distance);
            if (distance == knots[b]) {
                System.arraycopy(p_d_[b], 0, matrix, 0, matrix.length);
                return;
            } else if (distance < knots[b] && distance > knots[b-1]) {
                // approximation
                for (int j = 0; j < matrix.length; j++) {
                    // polynomials[j].length == breaks == knots.length-1
                    matrix[j] = polynomials[j][b-1].value(distance);
                }
                return;
            }
        }
        throw new IllegalArgumentException("Invalid distance " + distance +
                " not in knot[]" + Arrays.toString(knots));
    }


    // nrOfPoints = x.length / breaks; midIndex[i] = i * nrOfPoints;
    protected void createKnots(double[] x) {
        assert x.length > 2 * breaks;
        // 1st knot is 0, so knots = breaks + 1
        knots = new double[breaks+1];
        int nrOfPoints = x.length / breaks;

        knots[0] = x[0];
        knots[knots.length-1] = x[x.length-1];
        for (int i = 1; i < knots.length-1; i++) {
            knots[i] = x[i * nrOfPoints];
        }
    }


    // create segmented polynomials at one curve from points learnt in linear piecewise
    protected PolynomialFunction[] createSegmentedPolynomials(int degree, double[] xObs, double[] yObs) {
        assert knots.length == breaks + 1;

        PolynomialFunction[] polynFuns = new PolynomialFunction[breaks];
        WeightedObservedPoints obs = new WeightedObservedPoints();
        PolynomialCurveFitter fitter = PolynomialCurveFitter.create(degree);
        double[] coeff;
        List<WeightedObservedPoint> obsList;
        int i = 0;
        for (int b = 0; b < polynFuns.length; b++) {
            // allow yObs.length < xObs.length during debug
            while (i < yObs.length) {
                obs.add(xObs[i], yObs[i]);
                i++;
                if (xObs[i] >= knots[b+1])
                    break; // start next segment
            }
            obsList = obs.toList();
            assert obsList.size() > 1;

            coeff = fitter.fit(obsList);
            polynFuns[b] = new PolynomialFunction(coeff);
            obs.clear();
        }

        return polynFuns;
    }




    // use plotPt.R to observe data space
    public static void main(final String[] args) {
        String omega = "0.08", kappa = "15";
//        String[] omegas = new String[]{"0.01","0.08","0.1","1"};
//        String[] kappas = new String[]{"1","5","15","30"};
//        for (String omega : omegas) {
//            for (String kappa : kappas) {
        M0Model m0Model = getM0Model(omega, kappa);

        ApproxP_dist_Polynomial pd = new ApproxP_dist_Polynomial();
//                pd.initByName("substModel", m0Model, "file", "p_d_" + omega + "_" + kappa + ".txt");
        pd.initByName("substModel", m0Model);
        pd.printP_d_("Polynomial fit");
//        pd.testAccuracy(0.16, 0.001); //STEP
        pd.testAccuracy(pd.getLastKnot()+50, 0.01); //STEP
//            }
//        }
    }

    // test by breaks
    public void testAccuracy(final double maxX, double step) {
        List<Double> intervalList = new ArrayList<>();
        List<Double> stepList = new ArrayList<>();
        List<Double> uptoList = new ArrayList<>();

        initTest(maxX, step, intervalList, stepList, uptoList);

        assert breaks == knots.length-1;
        System.out.println("\nTesting " + breaks + " breaks ...");
        List<Double> xList = new ArrayList<>();
        int i=0;
        double x;
        for (int b = 0; b < knots.length-1; b++) {
            while (i < intervalList.size()) {
                x = intervalList.get(i);
                xList.add(x);
                i++;
                if (x >= knots[b+1] && b < knots.length-2)
                    break; // start next segment
            }

            printTestResult(xList);
            System.out.println("Break " + (b+1) + " completed.");
            xList.clear();
        }

    }


}
