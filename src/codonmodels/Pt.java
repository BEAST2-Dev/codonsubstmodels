package codonmodels;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.CodonAlignment;
import beast.evolution.alignment.Sequence;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Create points to plot P(t) by fixing other parameters.
 *
 * @author Walter Xie
 */
public class Pt {

    private double omega = 0.08;
    private double kappa = 15;

    private double[] freq;

    M0Model m0Model;

    double[] prob;
    double[] iexp;

    List<double[]> ptList = new ArrayList<>();

    public Pt() {
        CodonAlignment codonAlignment = initCodonAlignment();

        // equal
        CodonFrequencies codonFreq = new CodonFrequencies();
        codonFreq.initByName("pi", "equal", "data", codonAlignment, "verbose", true);

        initM0(codonFreq);
    }

    public Pt(CodonFrequencies codonFreq) {
        initM0(codonFreq);
    }

    private void initM0(CodonFrequencies codonFreq) {
        // omega 0.1, kappa 1
        RealParameter omegaPara = new RealParameter(Double.toString(omega));
        RealParameter kappaPara = new RealParameter(Double.toString(kappa));

        m0Model = new M0Model();
        m0Model.initByName("omega", omegaPara, "kappa", kappaPara,
                "frequencies", codonFreq, "verbose", true);

        freq = m0Model.getFrequencies();
        System.out.println("\nfreqs = \n" + Arrays.toString(freq) + "\n");
        for (int i = 0; i < freq.length; i++) {
            if (freq[i] == 0)
                throw new IllegalArgumentException("Frequency at codon index " + i + " is 0 !");
            if (freq[i] < 1E-5)
                System.err.println("Small frequency (< 1E-5) found at codon index " + i);
        }

        int len = m0Model.getStateCount();
        prob = new double[len*len];
        iexp = new double[len*len];
    }

    private CodonAlignment initCodonAlignment() {
        // Not use sequence, but have to init
        Sequence s1 = new Sequence("aaa", "AAA");
        Sequence s2 = new Sequence("aac", "AAA");

        Alignment data = new Alignment();
        data.initByName("sequence", s1, "sequence", s2, "dataType", "nucleotide");
        // create Codon Alignment
        CodonAlignment codonAlignment = new CodonAlignment();
        codonAlignment.initByName("data", data, "dataType", "codon",
                "geneticCode", "vertebrateMitochondrial");
        return codonAlignment;
    }


    public static void main(final String[] args) {
        Pt pt = new Pt();

        // rate = 1
        double jointRate = 1.0;

        double maxTime = pt.getMaxTime(jointRate);
        double[] intervals = pt.createTimeIntervals(maxTime);
        pt.getTransiProbsByTime(intervals,jointRate);

        Path path = Paths.get("p_t_.txt");
        try {
            pt.write(path, intervals);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public double getMaxTime(double rate) {

        int t;
        for (t = 10; t < 1000; t+=10) {
            // startTime > endTime, intervals[0] > 0
            m0Model.getTransiProbs(t, 0, rate, iexp, prob);

            boolean all = true;
            for (int i = 0; i < freq.length; i++) {
               if (prob[i]-freq[i] > 1E-4) // TODO check i->j
                   all = false;
            }

            if (all)
                return t;
        }
        throw new IllegalArgumentException("Cannot find the max time !");
    }


    // intervals[0] > 0
    public double[] createTimeIntervals(double maxTime){
        System.out.println("Max time = " + maxTime);

        int first = 3;
        int second = 30;
        int third = (int) maxTime / 4;
        List<Double> intervals = new ArrayList<>();

        int i = -1; // list index
        double end = 0;
        while (end < first) {
            i++;
            end += 0.05;
            end = Math.round(end*100) / 100.0;
            intervals.add(end);
        }
        System.out.println("Create " + intervals.size() + " intervals from " + intervals.get(0) +
                " to " + intervals.get(i));

        int j = i;
        while (end < second) {
            j++;
            end += 0.5;
            end = Math.round(end*10) / 10.0;
            intervals.add(end);
        }
        System.out.println("Create " + (intervals.size()-i+1) + " intervals from " + intervals.get(i+1) +
                " to " + intervals.get(j));

        int k = j;
        while (end < third) {
            k++;
            end += 2;
            end = Math.round(end);
            intervals.add(end);
        }
        System.out.println("Create " + (intervals.size()-j+1) + " intervals from " + intervals.get(j+1) +
                " to " + intervals.get(k));

        while (end < maxTime) {
            end += 20;
            end = Math.round(end);
            intervals.add(end);
        }
        System.out.println("Create last intervals " + (intervals.size()-k+1) +
                " intervals from " + intervals.get(k+1) + " to " + intervals.get(intervals.size()-1));

        return intervals.stream().mapToDouble(d -> d).toArray();
    }

    public void getTransiProbsByTime(double[] intervals, double rate) {

        for (int i = 0; i < intervals.length; i++) {

            // startTime > endTime, intervals[0] > 0
            m0Model.getTransiProbs(intervals[i], 0, rate, iexp, prob);
//            m0Model.getTransitionProbabilities(null, startTime, endTime, rate, prob);
            double[] tmp = new double[prob.length];
            System.arraycopy(prob, 0, tmp, 0, prob.length);
            ptList.add(tmp);

        }

    }

    public void getTransiProbsByOmega(double maxOmega, double maxTime, double rate) {
        final double step = 0.01;
        final double endTime = 0; // startTime > endTime
        double startTime = endTime + step;

        while (startTime < maxTime) {

            m0Model.getTransiProbs(startTime, endTime, rate, iexp, prob);
//            m0Model.getTransitionProbabilities(null, startTime, endTime, rate, prob);
            double[] tmp = new double[prob.length];
            System.arraycopy(prob, 0, tmp, 0, prob.length);
            ptList.add(tmp);

            startTime += step;

        }

    }

    // 1st[] is time, 2nd[] is flattened array of P(t)
    public double[][] getPt() {
        double[][] p_t_ = new double[ptList.size()][];

        for (int i = 0; i < ptList.size(); i++) {
            double[] probs = ptList.get(i);
            p_t_[i] = new double[probs.length];
            for (int j = 0; j < probs.length; j++) {
                p_t_[i][j] = probs[j];
            }
        }
        return p_t_;
    }


    public void write(Path path, double[] intervals) throws IOException {

        assert intervals.length == ptList.size();

        try (BufferedWriter writer = Files.newBufferedWriter(path)) {
            for (int i = 0; i < intervals.length; i++) {
                writer.write(intervals[i] + "\t");
                double[] probs = ptList.get(i);
                for (double pr : probs) {
                    writer.write("\t" + pr);
                }
                writer.newLine();
            }
        }


    }

}
