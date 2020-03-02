package test.beast.util;

import beast.util.RandomUtils;
import beast.util.Randomizer;
import org.junit.Test;

import java.util.Arrays;

import static org.junit.Assert.assertArrayEquals;

/**
 * @author Walter Xie
 */
public class BinarySearchSamplingTest {

    @Test
    public void binarySearch() {
        int nrOfStates = 60;
        // unnormalized probabilities
        double[] pr = new double[nrOfStates];
        // normalized probabilities
        double[] prob = new double[nrOfStates];
        double sum = 0;
        for (int i = 0; i < pr.length; i++) {
            pr[i] = i+1;
            sum += pr[i];
        }
        for (int i = 0; i < pr.length; i++)
            prob[i] = pr[i] / sum;
        System.out.println("Normalized probability distribution : " + Arrays.toString(prob) + "\n");

        int ite = 100000000; // 100 million
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

        sum = 0;
        double[] prob1 = new double[nrOfStates];
        for (int i = 0; i < freq1.length; i++)
            sum += freq1[i];
        for (int i = 0; i < freq1.length; i++)
            prob1[i] = freq1[i] / sum;
        System.out.println("Normalized probability : " + Arrays.toString(prob1) + "\n");

        assertArrayEquals(prob, prob1, 1E-4);

        //++++++  Binary search sampling ++++++//
        // unnormalized cumulative probabilities
        double[] cpd = new double[nrOfStates];
        cpd[0] = pr[0];
        for (int i = 1; i < pr.length; i++) {
            cpd[i] = cpd[i-1] + pr[i];
        }

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

        sum = 0;
        double[] prob2 = new double[nrOfStates];
        for (int i = 0; i < freq2.length; i++)
            sum += freq2[i];
        for (int i = 0; i < freq2.length; i++)
            prob2[i] = freq2[i] / sum;
        System.out.println("Normalized probability : " + Arrays.toString(prob2) + "\n");

        assertArrayEquals(prob, prob2, 1E-4);
    }


}