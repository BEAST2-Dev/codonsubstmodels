package beast.util;

import beast.app.BeastMCMC;

import java.util.Collection;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

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

//    private final List<Callable<T>> callers = new ArrayList<>();

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
     * @see ExecutorService#invokeAll(Collection)
     * @param tasks
     * @param <T>
     * @return
     * @throws InterruptedException
     */
    public <T> List<Future<T>> invokeAll(Collection<? extends Callable<T>> tasks)
            throws InterruptedException {
//        if (this.executor == null)
//            throw new IllegalArgumentException("ExecutorService cannot be null !");
        return this.executor.invokeAll(tasks);
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
