package test.beast.util;

import beast.util.RandomUtils;
import beast.util.Randomizer;
import org.junit.Before;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static org.junit.Assert.assertArrayEquals;

/**
 * @author Walter Xie
 */
public class BinarySearchSamplingTest {
    final int nrOfStates = 60;
    final int ite = 100000000; // 100 million

    double[] prob;
    double[] cpd;
    @Before
    public void setUp() throws Exception {
        // unnormalized probabilities
        List<Double> list = new ArrayList<>();
        for (int i = 0; i < nrOfStates; i++)
            list.add((double) (i+1));

        //*** shuffle ***
        Collections.shuffle(list);
        double[] pr = list.stream().mapToDouble(d -> d).toArray();

        // normalized probabilities
        prob = new double[nrOfStates];
        // compute distribution
        computeDistribution(pr, prob);
        System.out.println("\nTrue probability distribution : " + Arrays.toString(prob) + "\n");

        // unnormalized cumulative probabilities
        cpd = new double[nrOfStates];
        cpd[0] = pr[0];
        for (int i = 1; i < pr.length; i++) {
            cpd[i] = cpd[i-1] + pr[i];
        }

    }

    @Test
    public void binarySearch() {
        //++++++  Linear time sampling ++++++//
        int w;
        int[] freq1 = new int[nrOfStates];
        long start = System.currentTimeMillis();
        for (int i = 0; i < ite; i++) {
            w = RandomUtils.linearTimeSampling(prob, -1, false);
            freq1[w]++;
        }
        long timeLinear = System.currentTimeMillis() - start;
        System.out.println("Linear sampling time : " + timeLinear + " milliseconds.");
        System.out.println("Freq : " + Arrays.toString(freq1) + "\n");

        double[] prob1 = new double[nrOfStates];
        computeDistribution(freq1, prob1);
        System.out.println("Normalized probability : " + Arrays.toString(prob1) + "\n");

        assertArrayEquals(prob, prob1, 1E-4);

        //++++++  Binary search sampling ++++++//
        double random;
        int[] freq2 = new int[nrOfStates];
        start = System.currentTimeMillis();
        for (int i = 0; i < ite; i++) {
            random = Randomizer.nextDouble() * cpd[cpd.length-1];
            w = RandomUtils.binarySearchSampling(cpd, random);
            freq2[w]++;
        }
        long timeBiSearch = System.currentTimeMillis() - start;
        System.out.println("Binary search sampling time : " + timeBiSearch + " milliseconds.");
        System.out.println("Freq : " + Arrays.toString(freq2) + "\n");

        double[] prob2 = new double[nrOfStates];
        computeDistribution(freq2, prob2);
        System.out.println("Normalized probability : " + Arrays.toString(prob2) + "\n");

        assertArrayEquals(prob, prob2, 1E-4);
    }

    private void computeDistribution(double[] pr, double[] distr) {
        double sum = 0;
        for (int i = 0; i < pr.length; i++)
            sum += pr[i];
        for (int i = 0; i < pr.length; i++)
            distr[i] = pr[i] / sum;
    }

    private void computeDistribution(int[] freq, double[] distr) {
        double sum = 0;
        for (int i = 0; i < freq.length; i++)
            sum += freq[i];
        for (int i = 0; i < freq.length; i++)
            distr[i] = freq[i] / sum;
    }


}