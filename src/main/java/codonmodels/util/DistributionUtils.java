package codonmodels.util;

/**
 * @author Walter Xie
 */
public class DistributionUtils {

    public static void computeDistribution(double[] pr, double[] distr) {
        double sum = 0;
        for (int i = 0; i < pr.length; i++)
            sum += pr[i];
        for (int i = 0; i < pr.length; i++)
            distr[i] = pr[i] / sum;
    }

    public static void computeDistribution(int[] freq, double[] distr) {
        double sum = 0;
        for (int i = 0; i < freq.length; i++)
            sum += freq[i];
        for (int i = 0; i < freq.length; i++)
            distr[i] = freq[i] / sum;
    }

}
