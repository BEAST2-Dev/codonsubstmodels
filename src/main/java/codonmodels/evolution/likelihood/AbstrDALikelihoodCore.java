

package codonmodels.evolution.likelihood;

/**
 * framework of general data augmentation likelihood core
 */

abstract public class AbstrDALikelihoodCore {

    /**
     * reserve memory for branchLd, indices and other
     * data structures required by the core *
     * call it inside constructor
     */
    protected abstract void initialize();
    //int[][] tipStates, InternalNodeStates internalNodeStates, int categoryCount

    /**
     * clean up after last likelihood calculation, if at all required *
     */
    @Override
    public abstract void finalize() throws Throwable;


    /**
     * flag to indicate whether scaling should be used in the
     * likelihood calculation. Scaling can help in dealing with
     * numeric issues (underflow).
     */
    boolean m_bUseScaling = false;

    abstract public void setUseScaling(double scale);

    public boolean getUseScaling() {
        return m_bUseScaling;
    }
    /**
     * return the cumulative scaling effect. Should be zero if no scaling is used *
     */
    abstract public double getLogScalingFactor(int patternIndex_);


    /**
     * store current state *
     */
    abstract public void store();

    /**
     * restore state *
     */
    abstract public void restore();

    /**
     * reset current state to stored state, only used when switching from non-scaled to scaled or vice versa *
     */
    abstract public void unstore();

}
