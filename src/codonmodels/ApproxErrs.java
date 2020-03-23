package codonmodels;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.CodonAlignment;
import beast.evolution.alignment.Sequence;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

/**
 * @author Walter Xie
 */
public class ApproxErrs {

    double[][] p_d_;
    double[] intervals;

    double[] iexp;

    public ApproxErrs(M0Model m0Model) {

//        ApproximateSubstModel approxSubstModel = new ApproximateSubstModel();
//        approxSubstModel.initByName("substModel", m0Model);

        int stateCount = m0Model.getStateCount();
        int len = stateCount*stateCount;

        iexp = new double[len];
        double [] transprobs = new double[len];

        double t = 0;
        int N = 101;
        double [][] ratesThroughTime = new double[N][len];
        double [] times = new double[N];

        int k = 0;
        while (t < 10) {
            m0Model.getTransiProbs(t, 0, 1.0, iexp, transprobs);
            System.out.print(t+ "\t");
            System.arraycopy(transprobs, 0, ratesThroughTime[k], 0, len);
            for (int i = 0; i < transprobs.length; i++) {
                System.out.print(transprobs[i] + "\t");
            }
            System.out.println();
            times[k] = t;

            if (t == 0) {
                t += 1e-4;
            } else {
                t += 0.1;
            }
            k++;
        }

        double minErr = 1000000;
        double maxErr = 0;
        for (int i = 0; i < len; i++) {
            double error = approximate(ratesThroughTime, times, i);
            System.out.println(i  + "\t" + error + "\t");
            if (minErr > error)
                minErr = error;
            if (maxErr < error)
                maxErr = error;
        }
        System.out.println("\nMax variance = " + maxErr + ", min variance = " + minErr);

    }

    private double approximate(double[][] ratesThroughTime, double[] times, int index) {
        SplineInterpolator interpolator = new SplineInterpolator();
        double [] x = new double[16]; // hard code with set(...)
        double [] y = new double[16];
        int i = 0;
        set(i++, x, y, 0, times, ratesThroughTime, index);
        set(i++, x, y, 1, times, ratesThroughTime, index);
        set(i++, x, y, 2, times, ratesThroughTime, index);
        set(i++, x, y, 4, times, ratesThroughTime, index);
        set(i++, x, y, 6, times, ratesThroughTime, index);

        set(i++, x, y, 8, times, ratesThroughTime, index);
        set(i++, x, y, 10, times, ratesThroughTime, index);
        set(i++, x, y, 13, times, ratesThroughTime, index);
        set(i++, x, y, 16, times, ratesThroughTime, index);
        set(i++, x, y, 20, times, ratesThroughTime, index);

        set(i++, x, y, 25, times, ratesThroughTime, index);
        set(i++, x, y, 32, times, ratesThroughTime, index);
        set(i++, x, y, 42, times, ratesThroughTime, index);
        set(i++, x, y, 57, times, ratesThroughTime, index);
        set(i++, x, y, 77, times, ratesThroughTime, index);
        set(i++, x, y, 100, times, ratesThroughTime, index);


        PolynomialSplineFunction f = interpolator.interpolate(x, y);
        double variance = calcError(f, times, ratesThroughTime, index);
        return variance;
//        System.out.println(index  + "\t" + error + "\t");
        //System.out.println(Arrays.toString(diffs));

    }

    private double calcError(PolynomialSplineFunction poly, double[] times, double[][] ratesThroughTime, int index) {
        double error1 = 0;
        double error2 = 0;
        int N = times.length;
        double [] diffs = new double[N];
        double maxDiff = 0;
        for (int i = 0; i < N; i++) {
            double f = poly.value(times[i]);
            double diff = (f - ratesThroughTime[i][index]);
            diffs[i] = diff;
            error2 += diff * diff;
            error1 += Math.abs(diff);
            maxDiff = Math.max(maxDiff, Math.abs(diff));
        }
        // return maxDiff;
        // return error1 / N;
        return error2 / N;
    }


    private void set(int i, double[] x, double[] y, int j, double[] times, double[][] ratesThroughTime, int index) {
        x[i] = times[j];
        y[i] = ratesThroughTime[j][index];
    }


    public static void main(final String[] args) {

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

        RealParameter omegaPara = new RealParameter("0.08");
        RealParameter kappaPara = new RealParameter("15");

        M0Model m0Model = new M0Model();
        m0Model.initByName("omega", omegaPara, "kappa", kappaPara,
                "frequencies", codonFreq, "verbose", true);

//        ApproxP_dist_Linear pd = new ApproxP_dist_Linear();
//        pd.initByName("substModel", m0Model, "file", "p_d_1.txt");

        ApproxErrs approxErrs = new ApproxErrs(m0Model);



    }


}
