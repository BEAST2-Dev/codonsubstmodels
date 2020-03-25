package codonmodels;

import beast.core.util.Log;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.fitting.PolynomialCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;

/**
 * Interpolation
 * @author Walter Xie
 */
public class ApproxP_dist_Polynomial extends ApproxP_dist_Piecewise {


    // 3600
    PolynomialFunction[] polynomials; // replace double[][] p_d_

    @Override
    public void initAndValidate() {
//        super.initAndValidate();

//        double[][] xyMinMax = getXYMinMax();
//        printXY(xyMinMax);
//        printXYMax(xyMinMax);

        double[] freq = initCodonSubstInput();
        // learn the knots
        super.createKnotsPd();
        Log.info.println("Create " + knots.length + " data points using linear approximation.\n");

        double[][] xyMinMax = getXYMinMax();

        int ylen = 16;//knots.length;
        polynomials = new PolynomialFunction[prob.length];
        double[] y = new double[ylen];
        for (int i = 0; i < prob.length; i++) {
            for (int j = 0; j < y.length; j++)
                y[j] = p_d_[j][i];
            polynomials[i] = createPolynomialByQuantiles(10, knots, y, xyMinMax[i]);
        }

        Log.info.println("Polynomial fit approximation of codon substitution model " + codonSubstModel.getID());
        for (int i = 0; i < prob.length; i++) {
            System.out.println(i + "\t" + polynomials[i]);
        }

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
//        i = RandomUtils.binarySearchSampling(knots, distance);
//        if (distance == knots[i]) {
//            System.arraycopy(p_d_[i], 0 , matrix, 0, matrix.length);
//        } else {
            // approximation
            for (int j = 0; j < polynomials.length; j++) {
                matrix[j] = polynomials[j].value(distance);
            }
//        }

    }


    protected void createKnots() {
// per PolynomialCurveFitter?

    }


    // create polynomial from points learnt in linear piecewise
    protected PolynomialFunction createPolynomialByQuantiles(
            int nrOfQua, double[] xObs, double[] yObs, double[] minMaxLast) {

        double[] yQ = new double[nrOfQua];


        final WeightedObservedPoints obs = new WeightedObservedPoints();
        // if yObs.length < xObs.length
        for (int i = 0; i < yObs.length; i++) {
            obs.add(xObs[i], yObs[i]);
        }

        final PolynomialCurveFitter fitter = PolynomialCurveFitter.create(2);
        final double[] coeff = fitter.fit(obs.toList());

        return new PolynomialFunction(coeff);
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
//        pd.printP_d_();
        pd.testAccuracy(0.16, 0.001); //STEP
//            }
//        }
    }




}
