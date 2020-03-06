package beast.util;

import java.util.Arrays;

import static org.junit.Assert.assertArrayEquals;

/**
 * linear vs binary search.
 *
 * @author Walter Xie
 */
public class BinarySearchBenchmarking {

    static final int LINEAR = 0;
    static final int BINARY_SEARCH = 1;

    static final int TESTS = 50;
    final int ITE = 20000000; // 20 million

    private static int[] STATES = new int[]{4,21,40,60,80,100,128,256};

    public BinarySearchBenchmarking() {
//        for (int n = 0; n < states.length; n++)
//            states[n] = (int) Math.pow(2, (n+1));
    }

    public static void main(final String[] args) {

        BinarySearchBenchmarking benchmarking = new BinarySearchBenchmarking();

        int m = BINARY_SEARCH;//LINEAR;//

        long[][] time = new long[STATES.length][TESTS];

        for (int n = 0; n < STATES.length; n++) {
            for (int t = 0; t < TESTS; t++) {
                int state = STATES[n];
                System.out.println("\nNumber of states = " + state + ", test = " + (t + 1));
                // each test will generate random true probs in setUp(s)
                double[] trueProbs = benchmarking.setUpProbs(state);
                double[] cpd = benchmarking.getCPD(trueProbs);

//                //++++++  Binary search sampling ++++++//
//                benchmarking.benchmark(BINARY_SEARCH, n, t, cpd, trueProbs);
                //++++++  Linear time sampling ++++++//
                time[n][t] = benchmarking.benchmark(m, n, t, cpd, trueProbs);

            }
        } // end n loop

        benchmarking.printTime(m, time);
//        benchmarking.printTime(BINARY_SEARCH);

    }

    public long benchmark(int method, int stateIndex, int testIndex, double[] cpd, double[] trueProbs) {
        int state = STATES[stateIndex];
        int w;
        int[] freq1 = new int[state];

        long start = System.currentTimeMillis();
        for (int i = 0; i < ITE; i++) {
            if (method == LINEAR)
                w = Randomizer.randomChoice(cpd); // linear, require CPD
            else
                w = Randomizer.binarySearchSampling(cpd); // Binary search
            freq1[w]++;
        }
        long time = System.currentTimeMillis() - start;

        String msg = "Binary search sampling total time : ";
        if (method == LINEAR)
            msg = "Linear sampling total time :        ";
        System.out.println(msg + time + " milliseconds.");
//        System.out.println("\nFreq : " + Arrays.toString(freq1));

        double[] prob1 = new double[state];
        DistributionUtils.computeDistribution(freq1, prob1);
//        System.out.println("Normalized probability : " + Arrays.toString(prob1) + "\n");

        assertArrayEquals(trueProbs, prob1, 1E-3);
        return time;
    }


    public void printTime(int method, long[][] time) {
        long totTime = 0;
        if (method == LINEAR)
            System.out.println("\n#Linear sampling time :");
        else
            System.out.println("\n#Binary search sampling time :");
        for (int n = 0; n < STATES.length; n++) {
            int s = STATES[n];
            System.out.print(s + "\t");
            for (int t = 0; t < time[n].length; t++) {
                System.out.print(time[n][t] + "\t");
                totTime += time[n][t];
            }
            System.out.println(totTime);
        }
    }


    public double[] setUpProbs(int nrOfStates) {
        // unnormalized probabilities
        double[] freq = new double[nrOfStates];
        for (int i = 0; i < nrOfStates; i++)
            freq[i] = Randomizer.nextDouble(); // may have 0

        // normalized probabilities
        double[] prob = new double[nrOfStates];
        // compute distribution
        DistributionUtils.computeDistribution(freq, prob);
        System.out.println("\nTrue probability distribution : " + Arrays.toString(prob) + "\n");

        return prob;
    }

    public double[] getCPD(double[] prob) {
        // unnormalized cumulative probabilities
        double[] cpd = new double[prob.length];
//        long start = System.currentTimeMillis();
        cpd[0] = prob[0];
        for (int i = 1; i < prob.length; i++) {
            cpd[i] = cpd[i - 1] + prob[i];
        }
//        long time = System.currentTimeMillis() - start;
//        System.out.println("Renormalising time : " + time + " milliseconds.");
        return cpd;
    }

}
