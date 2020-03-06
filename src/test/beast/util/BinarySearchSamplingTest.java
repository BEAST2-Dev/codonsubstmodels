package test.beast.util;

import beast.util.DistributionUtils;
import beast.util.RandomUtils;
import beast.util.Randomizer;
import org.junit.Before;
import org.junit.Test;

import java.util.Arrays;

import static org.junit.Assert.assertArrayEquals;

/**
 * @author Walter Xie
 */
public class BinarySearchSamplingTest {

    final int nrOfStates = 60;
    final int ite = 100000000; // 100 million


    double[] prob;
    double[] cf;

    @Before
    public void setUp() {
        // unnormalized probabilities
        double[] freq = new double[nrOfStates];
        for (int i = 0; i < nrOfStates; i++)
            freq[i] = Randomizer.nextDouble(); // may have 0

        // normalized probabilities
        prob = new double[nrOfStates];
        // compute distribution
        DistributionUtils.computeDistribution(freq, prob);
        System.out.println("\nTrue probability distribution : " + Arrays.toString(prob) + "\n");

        // unnormalized cumulative probabilities
        cf = new double[nrOfStates];
        cf[0] = freq[0];
        for (int i = 1; i < freq.length; i++) {
            cf[i] = cf[i-1] + freq[i];
        }

    }

    @Test
    public void linearSampling() {
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
        DistributionUtils.computeDistribution(freq1, prob1);
        System.out.println("Normalized probability : " + Arrays.toString(prob1) + "\n");

        assertArrayEquals(prob, prob1, 1E-4);

    }

    @Test
    public void linearSampling2() {

        //++++++  Linear time sampling ++++++//
        int w;
        int[] freq1 = new int[nrOfStates];

        long start = System.currentTimeMillis();
        double[] cpd = new double[cf.length];
        for (int i = 0; i < ite; i++) {
            // count the time to normalise
            for (int j = 0; j < cf.length; j++)
                cpd[j] = cf[j] / cf[cf.length-1];

            w = Randomizer.randomChoice(cpd);
            freq1[w]++;
        }
        long timeLinear = System.currentTimeMillis() - start;
        System.out.println("Linear sampling time : " + timeLinear + " milliseconds.");

        System.out.println("Freq : " + Arrays.toString(freq1) + "\n");

        double[] prob1 = new double[nrOfStates];
        DistributionUtils.computeDistribution(freq1, prob1);
        System.out.println("Normalized probability : " + Arrays.toString(prob1) + "\n");

        assertArrayEquals(prob, prob1, 1E-4);

    }

    @Test
    public void binarySearch() {
        int w;
        //++++++  Binary search sampling ++++++//
        double random;
        int[] freq2 = new int[nrOfStates];

        long start = System.currentTimeMillis();
        for (int i = 0; i < ite; i++) {
            random = Randomizer.nextDouble() * cf[cf.length-1];
            w = RandomUtils.binarySearchSampling(cf, random);
            freq2[w]++;
        }
        long timeBiSearch = System.currentTimeMillis() - start;
        System.out.println("Binary search sampling time : " + timeBiSearch + " milliseconds.");

        System.out.println("Freq : " + Arrays.toString(freq2) + "\n");

        double[] prob2 = new double[nrOfStates];
        DistributionUtils.computeDistribution(freq2, prob2);
        System.out.println("Normalized probability : " + Arrays.toString(prob2) + "\n");

        assertArrayEquals(prob, prob2, 1E-4);
    }



}