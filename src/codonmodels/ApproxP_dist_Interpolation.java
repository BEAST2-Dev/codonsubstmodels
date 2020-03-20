package codonmodels;

import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.interpolation.UnivariateInterpolator;

import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

/**
 * Interpolation
 * @author Walter Xie
 */
public class ApproxP_dist_Interpolation extends ApproxP_dist_Piecewise {

    // knots
    final double [] x = new double[] {0.0, 1.0E-4, 0.10010000000000001, 0.30010000000000003, 0.5001, 0.7001,
            0.9000999999999999, 1.2001000000000002, 1.5001000000000004, 1.9001000000000008, 2.400100000000001,
            3.1001000000000016, 4.100100000000002, 5.600099999999997, 7.60009999999999, 9.900099999999982};


    UnivariateInterpolator[][] interpolators; // replace double[][] p_d_

    @Override
    public void initAndValidate() {
//        initCodonSubstInput();

        // learn the knots
        super.initAndValidate();

        double[][] xyMinMax = getXYMinMax();
        printXY(xyMinMax);
        printXYMax(xyMinMax);


//        interpolators = new UnivariateInterpolator[][nrOfStates * nrOfStates];
//        // 3600
//        for (UnivariateInterpolator interpolator : interpolators)
//            interpolator = new SplineInterpolator();

//        createApproximations();
    }

    @Override
    public void getTransiProbs(double distance, double[] matrix) {




    }

    // get 3 x and y, where y is min, max, last
    protected double[][] getXYMinMax() {
        // 0 x min, 1 y min, 2 x max, ...
        double[][] xy = new double[nrOfStates * nrOfStates][6];
        for (int c = 0; c < p_d_[0].length; c++)
            xy[c][1] = 1; // y min
        for (int x = 0; x < p_d_.length; x++) {
            for (int c = 0; c < p_d_[x].length; c++) {
                if (p_d_[x][c] < xy[c][1]) {
                    xy[c][0] = x;
                    xy[c][1] = p_d_[x][c]; // y min
                }
                if (p_d_[x][c] > xy[c][3]) {
                    xy[c][2] = x;
                    xy[c][3] = p_d_[x][c]; // y max
                }
                if (x == p_d_.length - 1) {
                    xy[c][4] = x;
                    xy[c][5] = p_d_[x][c]; // last
                }
            }
        }
        return xy;
    }

//    @Override
//    protected void createApproximations() {
//        List<Double> intervalList = new ArrayList<>();
//        List<SplineInterpolator[]> splineList = new ArrayList<>();
//
////        double t=0,rate = 1.0;
//        double x1, xPre=0, d = 0;
//        int k,r=0;
//        double[] mid = new double[prob.length];
//        double[] tru = new double[prob.length];
//        double[] yPre = new double[prob.length];
//        double[] y1Array;
//
//
//
//
//
////        // add 1st point 0
////        codonSubstModel.getTransiProbs(d, iexp, prob);
////        add(intervalList, splineList, 0, prob);
////        // 2nd point
////        d = 1E-4;
////        codonSubstModel.getTransiProbs(d, iexp, prob);
////        add(intervalList, splineList, d, prob);
//
//
//
//
//    }

    private void add(List<Double> intervalList, List<SplineInterpolator[]> splineList,
                     double d, SplineInterpolator[] spli) {
        intervalList.add(d);
        splineList.add(spli);
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
//        pd.printStd();
//            }
//        }
    }


    private void printXY(double[][] xy) {
        for (int x = 0; x < xy.length; x++) {
            System.out.print(x);
            for (int i = 0; i < xy[x].length; i++) {
                System.out.print("\t");
                if (i % 2 == 0)
                    System.out.print((int) xy[x][i]);
                else
                    System.out.print(xy[x][i]);
            }
            System.out.println();
        }
    }

    private void printXYMax(double[][] xy) {
        Set<Integer> xMax = new TreeSet<>();
        for (int x = 0; x < xy.length; x++) {
            double max = Math.round(xy[x][3] * 1000.0) / 1000.0;
            double last = Math.round(xy[x][5] * 1000.0) / 1000.0;
            if (max != last)
              xMax.add((int) xy[x][2]);
        }
        System.out.println(xMax.size() + " x : " + Arrays.toString(xMax.toArray()));
    }
}
