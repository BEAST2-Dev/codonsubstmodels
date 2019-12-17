package beast.util;

import java.util.Set;

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
     * Generating a random integer with non-uniform distribution,
     * https://stackoverflow.com/questions/42456115/generating-a-random-integer-with-non-uniform-distribution.
     * Use {@link Randomizer#setSeed(long)} to set seed.
     * @param probs     non-uniform distribution, sum to 1.
     * @param validate  if validate probs[] sums to 1.
     * @return          an integer between 0 and the length of probs[] minus 1.
     */
    public static int randomIntegerFrom(double[] probs, boolean validate) {
        if (validate) {
            double sum=0;
            for (double prob : probs)
                sum += prob;
            if (sum != 1)
                throw new RuntimeException("probabilities should sum to 1 " + sum);
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


}
