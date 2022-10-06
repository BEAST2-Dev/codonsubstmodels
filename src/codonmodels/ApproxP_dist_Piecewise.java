package codonmodels;

import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.evolution.alignment.CodonAlignment;
import beast.base.evolution.alignment.Sequence;
import beast.util.RandomUtils;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

/**
 * Find best approximation of P(dist) function by piecewise regression.
 * <code>dist = time * rate</code>
 *
 * @author Walter Xie
 */
public class ApproxP_dist_Piecewise extends CodonSubstitutionModel {
    final public Input<CodonSubstitutionModel> substModelInput = new Input<>("substModel",
            "substitution model we want to approximate", Input.Validate.REQUIRED);

    final public Input<String> filePathInput = new Input<>("file",
            "file path to save all data points used in the linear approximation");

    protected CodonSubstitutionModel codonSubstModel; // such as MO

    protected final double MaxDistance = 1000.0;
    protected final double DIFF = 1E-5;
    protected final double STEP = 0.001;
    protected final int multiply = 10;

    protected double[] knots; // sorted
    // each p_d_[] is a 60*60 flattened P(t).
    protected double[][] p_d_; // 1st[] time interval, 2nd[] codon index

    double[] prob;
    double[] iexp;

    public ApproxP_dist_Piecewise() {
        frequenciesInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {
        double[] freq = initCodonSubstInput();

        // core method to create approx model
        createKnotsPd(); // need double[] freq ?

        Log.info.println("Linear approximation of codon substitution model " + codonSubstModel.getID() +
                ", creating " + knots.length + " data points.\n");

        assert p_d_[0].length == nrOfStates * nrOfStates;

        // if file is given, then write to file
        if (filePathInput.get() != null) { // "p_d_.txt"
            Path path = Paths.get(filePathInput.get());
            try {
                writeP_d_(path);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    protected double[] initCodonSubstInput() {
        codonSubstModel = substModelInput.get();
        frequenciesInput.setValue(codonSubstModel.frequenciesInput.get(), this);
        verboseInput.setValue(false, this); // avoid to print rate matrix twice
        super.initAndValidate();
        eigenDecomposition = null;

        if (codonSubstModel instanceof M0Model)
            System.out.println("\nomega = " + ((M0Model) codonSubstModel).omegaInput.get() +
                    ", kappa = " + ((M0Model) codonSubstModel).kappaInput.get());

        // for caching
        prob = new double[nrOfStates * nrOfStates];
        iexp = new double[nrOfStates * nrOfStates];

        double[] freq = getFrequencies();
//        System.out.println("\nfreqs = \n" + Arrays.toString(freq) + "\n");
        for (int i = 0; i < freq.length; i++) {
            if (freq[i] == 0)
                throw new IllegalArgumentException("Illegal 0 frequency at codon index " + i + " !");
            if (freq[i] < 1E-5)
                System.err.println("Small frequency (< 1E-5) at codon index " + i);
        }

        return freq;
    }

    /**
     * Approximate P(t) by caching the list of P(t) matrices in time intervals.
     */
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
                matrix[j] = linearApproximation(distance, knots[i - 1], knots[i],
                        p_d_[i-1][j], p_d_[i][j]);
            }
            // no need to normalise, sum is very close to 1
        }
    }

    /**
     * Approximate P(t) by caching the list of P(t) matrices in time intervals.
     */
    public void getTransiProbs(double startTime, double endTime, double rate, double[] iexp, double[] matrix) {
        // distance = (startTime - endTime) * mean branch rate * rate for a site category.
        double distance = (startTime - endTime) * rate;

        getTransiProbs(distance, matrix);
    }

    // return linear approximation, distance = (startTime - endTime) * rate,
    private double linearApproximation(double distance, double x1, double x2, double y1, double y2) {
        // y = (x-x1) * (y2-y1) / (x2-x1) + y1, where x1 < x < x2, y1 < y < y2
        return (distance - x1) * (y2 - y1) / (x2 - x1) + y1;
    }

    /**
     * Create approx which can be used by {@link #getTransiProbs(double, double[])}.
     * Here is <code>double[][] p_d_<code/> and <code>double[] intervals<code/>.
     * Simple method to minimise the difference between mean of two points
     * and true value in the curve.
     * Less computation, but need many points.
     */
    protected void createKnotsPd(){
        List<Double> knotPointsList = new ArrayList<>();
        List<double[]> p_d_List = new ArrayList<>();

//        double t=0,rate = 1.0;
        double x1, xPre=0, d = 0;
        int k,r=0;
        double[] mid = new double[prob.length];
        double[] tru = new double[prob.length];
        double[] yPre = new double[prob.length];
        double[] y1Array;

        // add 1st point 0
        codonSubstModel.getTransiProbs(d, iexp, prob);
        add(knotPointsList, p_d_List, 0, prob);
        // 2nd point
        d = 1E-4;
        codonSubstModel.getTransiProbs(d, iexp, prob);
        add(knotPointsList, p_d_List, d, prob);

        double step = STEP;
        // 3rd point
        d = step;
        // 3-point proposal: last point, new proposed point, and the middle point between them
        while (d < MaxDistance) {
            // propose a new point, d ~ prob
            d = roundPrecision(d); // round precision error
            // eigen decomp to get points, startTime > endTime
            codonSubstModel.getTransiProbs(d, iexp, prob);

            // retrieve last point
            k = knotPointsList.size() - 1;
            x1 = knotPointsList.get(k);
            y1Array = p_d_List.get(k);
//            System.arraycopy(p_d_List.get(k), 0, y1Array, 0, prob.length);

            // middle point
            for (int j = 0; j < prob.length; j++)
                mid[j] = linearApproximation((d+x1)/2, x1, d, y1Array[j], prob[j]);
            codonSubstModel.getTransiProbs((d+x1)/2, iexp, tru);

            // test diff in the middle point
            if (largeDiff(mid, tru)) {
                if (xPre == x1) {
                    // no previous proposed point(s) previously
                    add(knotPointsList, p_d_List, d, prob);
                    // step too big
                    step /= multiply;
                } else
                    // add previous proposed point
                    add(knotPointsList, p_d_List, xPre, yPre);
                r = 0;
            } else
                r++;

            // save it for next ite
            xPre = d;
            System.arraycopy(prob, 0, yPre, 0, prob.length);

            if (r > 5 && step <= 10) {
                step *= multiply; // jumping 10 grids, then increase step
                r = 0;
            }
            d += step;

        }

        knots = knotPointsList.stream().mapToDouble(i -> i).toArray();
        p_d_ = getP_dist_(p_d_List);

    }

    private void add(List<Double> knotList, List<double[]> p_d_List, double d, double[] prob) {
        knotList.add(d);
        double[] tmp = new double[prob.length];
        System.arraycopy(prob, 0, tmp, 0, prob.length);
        p_d_List.add(tmp);
    }

    // use plotPt.R to observe data space
    public static void main(final String[] args) {
        String omega = "0.08", kappa = "15";
//        String[] omegas = new String[]{"0.01","0.08","0.1","1"};
//        String[] kappas = new String[]{"1","5","15","30"};
//        for (String omega : omegas) {
//            for (String kappa : kappas) {
                M0Model m0Model = getM0Model(omega, kappa);

                ApproxP_dist_Piecewise pd = new ApproxP_dist_Piecewise();
//                pd.initByName("substModel", m0Model, "file", "p_d_" + omega + "_" + kappa + ".txt");
        pd.initByName("substModel", m0Model);
        pd.printP_d_("Piecewise linear model");
        pd.testAccuracy(pd.getLastKnot()+50, 0.01); //STEP
//            }
//        }
    }

    protected static M0Model getM0Model(String omega, String kappa) {
        // Not use sequence, but have to init
        Sequence s1 = new Sequence("aaa", "AAA");
        Sequence s2 = new Sequence("aac", "AAA");

        Alignment data = new Alignment();
        data.initByName("sequence", s1, "sequence", s2, "dataType", "nucleotide");
        // create Codon Alignment
        CodonAlignment codonAlignment = new CodonAlignment();
        codonAlignment.initByName("data", data, "dataType", "codon",
                "geneticCode", "vertebrateMitochondrial");

        // equal
        CodonFrequencies codonFreq = new CodonFrequencies();
        codonFreq.initByName("pi", "equal", "data", codonAlignment, "verbose", true);

        RealParameter omegaPara = new RealParameter(omega);
        RealParameter kappaPara = new RealParameter(kappa);

        M0Model m0Model = new M0Model();
        m0Model.initByName("omega", omegaPara, "kappa", kappaPara,
                "frequencies", codonFreq, "verbose", true);
        return m0Model;
    }

    public double getLastKnot() {
        return knots[knots.length-1];
    }

    // the difference between estimated val[] and true values trueVal[].
    protected boolean largeDiff(double[] val, double[] trueVal) {

        for (int i = 0; i < val.length; i++) {
            double diff = Math.abs(val[i] - trueVal[i]);
            if (diff > DIFF)
                return true;
        }
        return false;
    }


    // 1st[] is time, 2nd[] is flattened array of P(t)
    protected double[][] getP_dist_(List<double[]> pdList) {
        double[][] p_d_ = new double[pdList.size()][];

        for (int i = 0; i < pdList.size(); i++) {
            double[] probs = pdList.get(i);
            p_d_[i] = new double[probs.length];
            for (int j = 0; j < probs.length; j++) {
                p_d_[i][j] = probs[j];
            }
        }
        return p_d_;
    }


    protected void writeP_d_(Path path) throws IOException {

        assert knots.length == p_d_.length;

        try (BufferedWriter writer = Files.newBufferedWriter(path)) {
            for (int i = 0; i < knots.length; i++) {
                writer.write(knots[i] + "\t");
                double[] probs = p_d_[i];
                for (double pr : probs) {
                    writer.write("\t" + pr);
                }
                writer.newLine();
            }
        }

    }

    public void printP_d_(String model) {

        assert knots.length == p_d_.length;

        for (int i = 0; i < knots.length; i++) {//intervals.length
            System.out.print(knots[i] + "\t");
            double[] probs = p_d_[i];
            for (double pr : probs) {
                System.out.print("\t" + pr);
            }
            System.out.println();
        }
//        System.out.println("\n");
//        for (int i = intervals.length-5; i < intervals.length; i++) {//
//            System.out.print(intervals[i] + "\t");
//            double[] probs = p_d_[i];
//            for (double pr : probs) {
//                System.out.print("\t" + pr);
//            }
//            System.out.println();
//        }

        System.out.println("\n" + model + " " + knots.length +
                " data points, last dist = " + getLastKnot() + ".");
    }

    // round precision error
    protected double roundPrecision(double d) {
        return Math.round(d * 100000.0) / 100000.0;
    }

    public void testAccuracy(final double maxX, double step) {
        List<Double> intervalList = new ArrayList<>();
        List<Double> stepList = new ArrayList<>();
        List<Double> uptoList = new ArrayList<>();

        initTest(maxX, step, intervalList, stepList, uptoList);

        printTestResult(intervalList);
    }

    protected void printTestResult(List<Double> intervalList) {
        double[] trueP_d_ = new double[nrOfStates * nrOfStates];
        double[] approxP_d_ = new double[nrOfStates * nrOfStates];
        double[] std = new double[nrOfStates * nrOfStates];
        double diff,minSd=1,maxSd=0;

        double dist;
        for (int i = 0; i < intervalList.size(); i++) {
            dist = intervalList.get(i);
            // true
            codonSubstModel.getTransiProbs(dist, iexp, trueP_d_);
            // approx
            this.getTransiProbs(dist, approxP_d_);

            for (int j = 0; j < std.length; j++) {
                diff = approxP_d_[j] - trueP_d_[j];
                std[j] += diff * diff;
            }
        }
        for (int j = 0; j < std.length; j++) {
            std[j] = Math.sqrt(std[j] / std.length);
//            System.out.println(j + "\t" + std[j]);

            if (minSd > std[j])  minSd = std[j];
            if (maxSd < std[j])  maxSd = std[j];
        }
        System.out.println("\nTesting accuracy at " + intervalList.size() + " points, dist = [" +
                intervalList.get(0) + ", " + intervalList.get(intervalList.size()-1) + "].");
        System.out.println("Max std = " + maxSd + ", min std = " + minSd);
    }

    protected void initTest(double maxX, double step, List<Double> intervalList,
                            List<Double> stepList, List<Double> uptoList) {
        intervalList.add(0.0);
        intervalList.add(1E-5);
        intervalList.add(1E-4);
        if (step > 1E-3)  intervalList.add(1E-3);

        double r = 10;
        double dist = step;
        stepList.add(step);
        while (dist <= maxX) {
            dist = roundPrecision(dist);
            intervalList.add(dist);
            dist += step;
            if (dist > r && step < 10) {
                uptoList.add(dist);
                step *= multiply;
                stepList.add(step);
                r *= 5;
            }
        }
        uptoList.add(dist-step);

        System.out.println("\nInitialise test at " + intervalList.size() + " points, dist = [" +
                intervalList.get(0) + ", " + intervalList.get(intervalList.size()-1) + "].");
        System.out.println("Steps used         " + stepList);
        System.out.println("Step increased at  " + uptoList);
    }

    /*** analyse the curve to plot 3 x and y: min, max, last of P(d) ***/

    // 1st[] is codon index 3600, 2nd[] is min max
    protected double[][] getXYMinMax() {
        // 0 x min x index, 1 x min, 2 y min, 3 x max x index, ...
        double[][] xy = new double[nrOfStates * nrOfStates][9];
        for (int c = 0; c < p_d_[0].length; c++){
            xy[c][1] = 100000; // x min
            xy[c][2] = 2; // y min
        }
        for (int x = 0; x < p_d_.length; x++) {
            for (int c = 0; c < p_d_[x].length; c++) {
                if (p_d_[x][c] < xy[c][2]) {
                    xy[c][0] = x;
                    xy[c][1] = knots[x]; // x min
                    xy[c][2] = p_d_[x][c]; // y min
                }
                if (p_d_[x][c] > xy[c][5]) {
                    xy[c][3] = x;
                    xy[c][4] = knots[x]; // x max
                    xy[c][5] = p_d_[x][c]; // y max
                }
                if (x == p_d_.length - 1) {
                    xy[c][6] = x;
                    xy[c][7] = knots[x]; // x last
                    xy[c][8] = p_d_[x][c]; // y last
                }
            }
        }
        return xy;
    }

    // 3 x max x index, 4 x max, 5 y max, 8 y last
    protected void printXYMax(double[][] xy) {
        Set<Integer> xMaxId = new TreeSet<>();
        Set<Double> xMax = new TreeSet<>();
        for (int x = 0; x < xy.length; x++) {
            double yMax = Math.round(xy[x][5] * 1000.0) / 1000.0;
            double last = Math.round(xy[x][8] * 1000.0) / 1000.0;
            if (yMax != last) {
                xMaxId.add((int) xy[x][3]);
                xMax.add(xy[x][4]);
            }
        }
        // unique mapping
        System.out.println(xMax.size() + " x = " + Arrays.toString(xMax.toArray()));
        System.out.println(xMaxId.size() + " x index = " + Arrays.toString(xMaxId.toArray()));
    }

    protected void printXY(double[][] xy) {
        for (int x = 0; x < xy.length; x++) {
            System.out.print(x);
            for (int i = 0; i < xy[x].length; i++) {
                System.out.print("\t");
                if (i % 3 == 0)
                    System.out.print((int) xy[x][i]);
                else
                    System.out.print(xy[x][i]);
            }
            System.out.println();
        }
    }
}


