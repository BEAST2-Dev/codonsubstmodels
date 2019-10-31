

package beast.evolution.likelihood;

/**
 * data augmentation likelihood core
 */

abstract public class DABranchLikelihoodCore {

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
     * reset current state to stored state, only used when switching from non-scaled to scaled or vice versa *
     */
    abstract public void unstore();

    /**
     * restore state *
     */
    abstract public void restore();

//    /**
//     * reserve memory for branchLd for node with number nodeIndex *
//     */
//    public abstract void createNodeBranchLd(int nodeIndex);


    /**
     * indicate that the probability transition matrix for node
     * nodeIndex is about the be changed, that is, that the stored
     * state for node nodeIndex cannot be reused *
     */
    public abstract void setNodeMatrixForUpdate();

    /**
     * assign values of states for probability transition matrix for node with number nodeIndex *
     */
    public abstract void setNodeMatrix(int categoryIndex, double[] matrix);
//    public abstract void setNodeMatrix(double[] matrix);

//    public abstract void calculateNodeBranchLd(int childNum1, int childNum2, int nodeIndex);
//    public abstract void integrateBrLdOverCategories(double[] proportions, double[] integratedBrLd);

//    public abstract void calculateNodeBrLdOverCategories(int nodeIndex1, int nodeIndex2, int nodeIndex3, double[] proportions);

    public abstract void calculateNodeBrLdOverCategories(int[] childNodeStates, int[] parentNodeStates, double[] proportions);

    //    public abstract void calculateLogLikelihoods(double[] rootBranchLd, double[] frequencies, double[] siteLogLikelihoods);
    public abstract double calculateLogLikelihoods(boolean isRoot, double frequencies);


    //    @Override
//    public void setNodeBranchLd(double[] branchLd) {
//
//        if (this.branchLd[0][nodeIndex] == null) {
//            createNodeBranchLd(nodeIndex);
//        }
//        if (branchLd.length < branchLdSize) {
//            int k = 0;
//            for (int i = 0; i < nrOfCategories; i++) {
//                System.arraycopy(branchLd, 0, this.branchLd[0][nodeIndex], k, branchLd.length);
//                k += branchLd.length;
//            }
//        } else {
//            System.arraycopy(branchLd, 0, this.branchLd[0][nodeIndex], 0, branchLd.length);
//        }
//    }

    // suppose only used by unit test
    public abstract void getNodeBranchLd(double[] branchLdOut);

    public abstract void setNodeBrLdForUpdate();

//    public abstract int[] getNodeStates();
//
//    public abstract void setInternalNodeStates(int[] states);
}
