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

}
