package codonmodels.util;

import java.util.Set;

import beast.base.util.Randomizer;

/**
 * @author Walter Xie
 */
public class RandomUtils {

    /**
     * Randomly select an element from the set of {@link Integer}.
     * Use {@link Randomizer#setSeed(long)} to set seed.
     * @param candidates   a set of states to choose
     * @return
     */
    public static int randomSelect(Set<Integer> candidates) {
        int size = candidates.size();
        if (size < 1)
            throw new IllegalArgumentException("Set has no element !");
        // Randomizer.setSeed(long)
        int item = 0;
        if (size > 1) item = Randomizer.nextInt(size);
        int s = 0;
        for (Integer state : candidates) {
            if (s == item)
                return state;
            s++;
        }
        return -1;
    }

    /**
     * Binary search to sample an integer given an unnormalized cumulative probabilities.
     * Modified from {@link java.util.Arrays#binarySearch(double[], double)}.
     * @param cpd       unnormalized or normalized cumulative probabilities,
     *                  which have not to be a probability distribution.
     *                  If not a distribution (sum > 1), then <code>random = [0, sum)</code>.
     * @param randDoub  the random double between [0, sum), where sum is <code>cpd[cpd.length-1]</code>.
     *                  <code>double random = Randomizer.nextDouble() * cpd[cpd.length-1];</code>.
     * @return          a sample (index of <code>cpd[]</code>) according to
     *                  an unnormalized cumulative probabilities.
     */
    public static int binarySearchSampling(double[] cpd, double randDoub) {
        if (randDoub <= cpd[0])
            return 0;

        int mid,low = 0;
        int high = cpd.length - 1;
        double midVal;
        while (low <= high) {
            mid = (low + high) >>> 1;
            midVal = cpd[mid];

//            if (randDoub <= cpd[mid] && randDoub > cpd[mid - 1])
//                return mid; // take i where cpd[i - 1] < randDoub <= cpd[i]

            if (midVal < randDoub)
                low = mid + 1;
            else if (midVal > randDoub) {
                if (cpd[mid - 1] < randDoub)
                    return mid;
                high = mid - 1;
            } else
                return mid; // cpd == random

        }
        return -(low + 1);  // cpd not found.
    }

    /**
     * Chooses one category if a cumulative probability distribution is given.
     * Modified from {@link Randomizer#randomChoice(double[])}.
     * @param cf
     * @return
     */
    public static int randomChoice(double[] cf, double randDoub) {
        int s;
        if (randDoub <= cf[0]) {
            s = 0;
        } else {
            for (s = 1; s < cf.length; s++) {
                if (randDoub <= cf[s] && randDoub > cf[s - 1]) {
                    break;
                }
            }
        }
        return s;
    }

    /**
     * Replaced by {@link #binarySearchSampling(double[], double)}.
     * Generating a random integer with non-uniform distribution,
     * https://stackoverflow.com/questions/42456115/generating-a-random-integer-with-non-uniform-distribution.
     * Use {@link Randomizer#setSeed(long)} to set seed.
     * @param probs     non-uniform distribution, sum to 1.
     * @param validate  if validate probs[] sums to 1.
     * @return          an integer between 0 and the length of probs[] minus 1.
     */
    @Deprecated
    public static int linearTimeSampling(double[] probs, double sum, boolean validate) {
        if (sum > 0) {
            // Renormalise all w
            for (int w=0; w < probs.length; w++)
                probs[w] = probs[w] / sum;
        }

        if (validate) {
            double sumpr=0;
            for (double prob : probs)
                sumpr += prob;
            if (sumpr != 1)
                throw new RuntimeException("probabilities should sum to 1 " + sumpr);
        }
        double r = Randomizer.nextDouble(); // [0,1)
        for (int i = 0; i < probs.length; i++) {
            if (r < probs[i]) {
                return i;
            }
            r -= probs[i];
        }
        throw new RuntimeException("probabilities should sum to 1");
    }



    /**
     * randomly select n integers between min and max-min.
     * @param min
     * @param max
     * @param n
     * @param integerSet   result set, used to avoid <code>new HashSet<>()</code>.
     */
    public static void getRandomInt(final int min, final int max,
                                    final int n, Set<Integer> integerSet) {
        integerSet.clear();
        do {
            // Samples an int uniformly from between 0 and n-1.
            int nr = min + Randomizer.nextInt(max-min);
            integerSet.add(nr);
        } while (integerSet.size() < n);
    }


}
