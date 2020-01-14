package beast.util;

import beast.app.BeastMCMC;
import beast.core.util.Log;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * @author Walter Xie
 */
public class ThreadHelper {
    // New thread pool from BEAST MCMC.
    private ExecutorService executor = null;

    // number of threads to use
    private int threadCount;
    // maximum number of branches assigned to a task for thread pool
//    private int maxBrPerTask;

    /**
     * Initiate the number of threads from BEAST command line,
     * or XML threads=, or instance count.
     * @param maxNrOfThreads  the number of threads from XML.
     */
    public ThreadHelper(int maxNrOfThreads) {
        // init threadCount
        threadCount = BeastMCMC.m_nThreads;
        // overwritten by max(maxNrOfThreadsInput, BEAST thread in command line)
        if (maxNrOfThreads > 0) {
            threadCount = Math.max(maxNrOfThreads, BeastMCMC.m_nThreads);
        }
        // threadCount is overwritten by instanceCount
        String instanceCount = System.getProperty("beast.instance.count");
        if (instanceCount != null && instanceCount.length() > 0) {
            threadCount = Integer.parseInt(instanceCount);
        }
        Log.info("Data augmentation tree likelihood thread = " + threadCount);
    }

    /**
     * Initiate the number of threads and {@link ExecutorService}.
     * @param maxNrOfThreads   the number of threads from XML.
     * @param executor         {@link ExecutorService}
     * if null, then create a fixed thread pool when threadCount > 1.
     */
    public ThreadHelper(int maxNrOfThreads, ExecutorService executor) {
        // init threadCount
        this(maxNrOfThreads);

        if (threadCount > 1) {
            if (executor == null)
                this.executor = Executors.newFixedThreadPool(threadCount);
            else
                this.executor = executor;
        }
    }

    /**
     * Be careful to use.
     * @see ExecutorService#shutdown()
     */
    public void shutdown() {
        if (executor != null) executor.shutdown();
    }

    public int getThreadCount() {
        return threadCount;
    }

    public ExecutorService getExecutor() {
        return executor;
    }

//    public int getMaxBrPerTask() {
//        return maxBrPerTask;
//    }

}
