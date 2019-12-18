package beast.util;

import beast.evolution.alignment.CodonAlignment;
import beast.evolution.likelihood.DataAugTreeLikelihood;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.NodeStatesArray;
import beast.evolution.tree.Tree;

import java.io.IOException;


/**
 * Performance benchmarking between standard TreeLikelihood and DACodonTreeLikelihood.
 * <p>
 * <b>Benchmarking1: brute testing</b>
 * <p>
 * Record the time of 100 iterations to call either {@link TreeLikelihood#calculateLogP()}
 * or {@link DataAugTreeLikelihood#calculateLogP()}.
 * <p>
 * The tests cover the following of scenarios:
 * <ol>
 *     <li>Increasing the number of taxa</li>
 *     <li>Increasing the number of codon</li>
 *     <li>Symmetric or extremely asymmetric tree topology in 4 and 8 taxa</li>
 *     <li>TODO Increasing the number of threads</li>
 * </ol>
 * @author Walter Xie
 */
public class DALikelihoodBenchmarking1 extends BenchmarkingSetup {

    final int iteration;

//    public DALikelihoodBenchmarking1(String fileName, int iteration) throws IOException, XMLParserException {
//        super(fileName);
//        this.iteration = iteration;
//        this.states = null;
//    }

    public DALikelihoodBenchmarking1(String fileName, int iteration, int MAX_TIPS, int MAX_CODONS)
            throws IOException, XMLParserException {
        super(fileName, MAX_TIPS, MAX_CODONS);
        this.iteration = iteration;
    }


    @Override
    public long[][] run(boolean isDA, boolean verbose){
        if (isDA)
            System.out.println("\n=============== DA tree likelihood ===============\n");
        else
            System.out.println("\n=============== Standard tree likelihood ===============\n");

        return super.run(isDA, verbose);
    }

    @Override
    protected long[] timeTest(boolean isDA, CodonAlignment codonAlignment, SiteModel siteModel, Tree tree) {
        long[] tm;
        if (isDA) {
            tm = benchmarkingDA(codonAlignment, siteModel, tree);
        } else {
            tm = benchmarkingStandard(codonAlignment, siteModel, tree);
        }

        System.out.println("Init time " + tm[0] + " nanoseconds");
        System.out.println(iteration + " iteration(s) " + tm[1] + " nanoseconds");
        System.out.println("Calculation time = " +  (double) tm[1] / (double)iteration +
                " nanoseconds per iteration in average.");
        System.out.println();

        return tm;
    }

    public static void main(final String[] args) {
        System.out.print("Program arguments: ");
        for (String arg : args) {
            System.out.print("  " + arg);
        }
        System.out.print("\n");

        // 1st arg, if not "BEAST" then use DA tree likelihood
        boolean isDA = !"BEAST".equalsIgnoreCase(args[0]);

        int iteration = 1;
        if (args.length > 1) {
            try {
                iteration = Integer.parseInt(args[1]);
            } catch (NumberFormatException e) {
                System.err.println("2nd argument iteration has to be integer, but find " + args[1]);
            }
        }

        DALikelihoodBenchmarking1 likelihoodBenchmarking = null;

        try {
            if (args.length > 2)
                likelihoodBenchmarking = new DALikelihoodBenchmarking1(args[2], iteration, 1024, 1000);
            else
                likelihoodBenchmarking = new DALikelihoodBenchmarking1(null, iteration, 128, 500);

        } catch (IOException | XMLParserException e) {
            e.printStackTrace();
        }

        boolean verbose = false;
        long[][] time = likelihoodBenchmarking.run(isDA, verbose);
        likelihoodBenchmarking.report(time, isDA);
    }


    /**
     * Codon500 1 iteration(s) :
     * test 2 taxa : DA likelihood 25.00 times faster
     * test 4 taxa symmetric : DA likelihood 6.33 times faster
     * test 4 taxa asymmetric : DA likelihood 7.00 times faster
     * test 8 taxa symmetric : DA likelihood 2.60 times faster
     * test 8 taxa asymmetric : DA likelihood 2.60 times faster
     * test 16 taxa : DA likelihood 2.10 times faster
     * test 32 taxa : DA likelihood 2.50 times faster
     * test 64 taxa : DA likelihood 2.42 times faster
     * test 128 taxa : DA likelihood 2.80 times faster
     */
    // =============== report ===============
    @Override
    public void report(final long[][] time, boolean isDA) {

        System.out.print("\n\n=============== Summary of ");
        if (isDA)
            System.out.print("DA tree likelihood " + threads + " thread(s) ");
        else
            System.out.print("standard tree likelihood ");
        System.out.print(iteration + " iteration(s) ===============\n\n");

        super.printTime(time, "nanoseconds");
    }

    // =============== Standard likelihood ===============
    // return nano seconds
    @Override
    protected long[] benchmarkingStandard(CodonAlignment codonAlignment, SiteModel siteModel, Tree tree) {
        // init time
        long start = System.nanoTime();

        TreeLikelihood likelihood = new TreeLikelihood();
        likelihood.initByName("data", codonAlignment, "tree", tree, "siteModel", siteModel);

        long standardInit = System.nanoTime()-start;

        double logP = 0;
        start = System.nanoTime();
        // more iterations to avoid 0 millisecond
        for (int i=0; i<iteration; i++) {
            logP = likelihood.calculateLogP();
        }
        // Get elapsed time in nanoseconds
        long timeStandard = System.nanoTime()-start;
        System.out.println("Standard tree likelihood: logP = " + logP);

        return new long[]{standardInit, timeStandard};
    }

    // =============== DA likelihood ===============
    // return nano seconds
    @Override
    protected long[] benchmarkingDA(CodonAlignment codonAlignment, SiteModel siteModel, Tree tree) {
        // init time
        long start = System.nanoTime();

        NodeStatesArray nodesStates = new NodeStatesArray(codonAlignment, tree, "parsimony");

        DataAugTreeLikelihood daLikelihood = new DataAugTreeLikelihood();
        daLikelihood.initByName("nodesStates", nodesStates, "tree", tree, "siteModel", siteModel, "threads", threads);

        long daInit = System.nanoTime() - start;

        double logPDA = 0;
        start = System.nanoTime();
        // more iterations to avoid 0 millisecond
        for (int i = 0; i < iteration; i++) {
            logPDA = daLikelihood.calculateLogP();
        }
        // Get elapsed time in nanoseconds
        long timeDA = System.nanoTime() - start;
        System.out.println("DA tree likelihood: DA logP = " + logPDA);

        daLikelihood.shutdown(); // shut down ExecutorService

        return new long[]{daInit, timeDA};
    }

}
