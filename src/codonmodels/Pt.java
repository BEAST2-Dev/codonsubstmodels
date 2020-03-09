package codonmodels;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.CodonAlignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.datatype.Codon;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Walter Xie
 */
public class Pt {

    M0Model m0Model;

    double[] prob;
    double[] iexp;

    List<double[]> ptList = new ArrayList<>();

    public Pt() {
        Sequence s1 = new Sequence("aaa", "AAA");
        Sequence s2 = new Sequence("aac", "AAA");

        Alignment data = new Alignment();
        data.initByName("sequence", s1, "sequence", s2, "dataType", "nucleotide");
        // create Codon Alignment
        CodonAlignment codonAlignment = new CodonAlignment();
        codonAlignment.initByName("data", data, "dataType", "codon",
                "geneticCode", "vertebrateMitochondrial");
        Codon codon = codonAlignment.getDataType();

        // F3x4
        CodonFrequencies codonFreq = new CodonFrequencies();
        codonFreq.initByName("pi", "equal", "data", codonAlignment, "verbose", true);

        // omega 0.1, kappa 1
        RealParameter omega = new RealParameter("0.1");
        RealParameter kappa = new RealParameter("1");

        m0Model = new M0Model();
        m0Model.initByName("omega", omega, "kappa", kappa, "frequencies", codonFreq, "verbose", true);

        System.out.println("\nfreqs = \n" + Arrays.toString(m0Model.getFrequencies()) + "\n");

        int len = m0Model.getStateCount();
        prob = new double[len*len];
        iexp = new double[len*len];

    }


    public static void main(final String[] args) {
        Pt pt = new Pt();

        // rate = 1
        pt.getTransiProbs(100,1);

        Path path = Paths.get("pt.txt");
        try {
            pt.write(path);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void getTransiProbs(double time, double rate) {
        final double step = 0.01;
        final double endTime = 0; // startTime > endTime
        double startTime = endTime + step;

        while (startTime < time) {

            m0Model.getTransiProbs(startTime, endTime, rate, iexp, prob);
//            m0Model.getTransitionProbabilities(null, startTime, endTime, rate, prob);
            double[] tmp = new double[prob.length];
            System.arraycopy(prob, 0, tmp, 0, prob.length);
            ptList.add(tmp);

            startTime += step;

        }

    }

    public void write(Path path) throws IOException {

        try (BufferedWriter writer = Files.newBufferedWriter(path)) {
            for (double[] probs : ptList) {
                for (double pr : probs) {
                    writer.write(pr + "\t");
                }
                writer.newLine();
            }
        }


    }

}
