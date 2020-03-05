package test.beast.util;

import beast.util.RandomUtils;
import beast.util.Randomizer;
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

    double[] prob;
    double[] cf;

    public void setUp(int nrOfStates) {
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
        cf = new double[nrOfStates];
        cf[0] = pr[0];
        for (int i = 1; i < pr.length; i++) {
            cf[i] = cf[i-1] + pr[i];
        }

    }

    @Test
    public void binarySearch() {
        final int nrOfStates = 60;
        final int ite = 100000000; // 100 million

        setUp(nrOfStates);

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
            random = Randomizer.nextDouble() * cf[cf.length-1];
            w = RandomUtils.binarySearchSampling(cf, random);
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

    @Test
    public void binarySearchOptimalTest() {
        final int ite = 100000000; // 100 million
        final int maxState = 60;

        long[][] time = new long[2][maxState];

        for (int s = 4; s < maxState; s++) {
            System.out.println("\nMax states = " + s + "\n");
            setUp(s);

            //++++++  Linear time sampling ++++++//
            int w;
            int[] freq1 = new int[s];
            double[] cpd = new double[s];
            cpd[0] = prob[0];
            for (int i = 1; i < prob.length; i++) {
                cpd[i] = cpd[i-1] + prob[i];
            }

            long start = System.currentTimeMillis();
            for (int i = 0; i < ite; i++) {
                // linear time
                w = Randomizer.randomChoice(cpd); // require CPD
                freq1[w]++;
            }
            time[0][s] = System.currentTimeMillis() - start;
            System.out.println("Linear sampling time : " + time[0][s] + " milliseconds.");
            System.out.println("Freq : " + Arrays.toString(freq1) + "\n");

            double[] prob1 = new double[s];
            computeDistribution(freq1, prob1);
            System.out.println("Normalized probability : " + Arrays.toString(prob1) + "\n");

            assertArrayEquals(prob, prob1, 1E-4);


            //++++++  Binary search sampling ++++++//
            int[] freq2 = new int[s];
            start = System.currentTimeMillis();
            for (int i = 0; i < ite; i++) {
                double random = Randomizer.nextDouble() * cf[cf.length-1];
                // Binary search
                w = RandomUtils.binarySearchSampling(cf, random);
                freq2[w]++;
            }
            time[1][s] = System.currentTimeMillis() - start;
            System.out.println("Binary search sampling time : " + time[1][s] + " milliseconds.");
            System.out.println("Freq : " + Arrays.toString(freq2) + "\n");

            double[] prob2 = new double[s];
            computeDistribution(freq2, prob2);
            System.out.println("Normalized probability : " + Arrays.toString(prob2) + "\n");

            assertArrayEquals(prob, prob2, 1E-4);

        } // end s loop

        System.out.println("\nLinear sampling time : \n");
        for (int s = 4; s < maxState; s++) {
            System.out.println(s + " :\t" + time[0][s]);
        }
        System.out.println("\nBinary search sampling time : \n");
        for (int s = 4; s < maxState; s++) {
            System.out.println(s + " :\t" + time[1][s]);
        }
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