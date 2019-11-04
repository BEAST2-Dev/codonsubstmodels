

package beast.evolution.likelihood;

import beast.evolution.tree.InternalNodeStates;

/**
 * framework of data augmentation tree likelihood core
 */

public abstract class AbstrDATreeLikelihoodCore extends AbstrDALikelihoodCore {
    // to store branch likelihood calculation per site:
    // 1st dimension is matrix index (current, stored),
    // 2nd is node index,
    // 3rd is getNrOfSites(), Ld across categories are integrated
    protected double[][][] branchLd;

    // transition probability matrix(ices), P
    // 1st dimension is matrix index (current, stored),
    // 2nd is node index,
    // 3rd is nrOfCategories * matrixSize
    protected double[][][] matrices;
    // store the matrix index, instead of different matrices
    protected int[] currentMatrixIndex; // node count
    protected int[] storedMatrixIndex;
    protected int[] currentBrLdIndex;
    protected int[] storedBrLdIndex;


    // states in tip/internal nodes: 0-63
    // 1st is node index, 2nd is site index
    protected int[][] tipStates;
    protected InternalNodeStates internalNodeStates;


    /**
     * need to call initialize() after this
     * @param nrOfStates
     * @param tipStates
     * @param internalNodeStates
     * @param nrOfCategories
     */
    public AbstrDATreeLikelihoodCore(int nrOfStates, int[][] tipStates, InternalNodeStates internalNodeStates,
                                     int nrOfCategories) {
        super(nrOfStates, nrOfCategories);

        this.tipStates = tipStates;
        this.internalNodeStates = internalNodeStates;

    }


    /**
     * Store current state
     */
    @Override
    public void restore() {
        // Rather than copying the stored stuff back, just swap the pointers...
        int[] tmp1 = currentMatrixIndex;
        currentMatrixIndex = storedMatrixIndex;
        storedMatrixIndex = tmp1;

        int[] tmp2 = currentBrLdIndex;
        currentBrLdIndex = storedBrLdIndex;
        storedBrLdIndex = tmp2;
    }

    @Override
    public void unstore() {
        System.arraycopy(storedMatrixIndex, 0, currentMatrixIndex, 0, getNrOfNodes());
        System.arraycopy(storedBrLdIndex, 0, currentBrLdIndex, 0, getNrOfNodes());
    }

    /**
     * Restore the stored state
     */
    @Override
    public void store() {
        System.arraycopy(currentMatrixIndex, 0, storedMatrixIndex, 0, getNrOfNodes());
        System.arraycopy(currentBrLdIndex, 0, storedBrLdIndex, 0, getNrOfNodes());
    }

    public abstract int getNrOfNodes();

//    /**
//     * reserve memory for branchLd for node with number nodeIndex *
//     */
//    public abstract void createNodeBranchLd(int nodeIndex);


    /**
     * indicate that the probability transition matrix for node
     * nodeIndex is about the be changed, that is, that the stored
     * state for node nodeIndex cannot be reused *
     */
    public abstract void setNodeMatrixForUpdate(int nodeIndex);

    /**
     * assign values of states for probability transition matrix for node with number nodeIndex *
     */
    public abstract void setNodeMatrix(int nodeIndex, int categoryIndex, double[] matrix);
//    public abstract void setNodeMatrix(int nodeIndex, double[] matrix);

//    public abstract void calculateNodeBranchLd(int childNum1, int childNum2, int nodeIndex);
//    public abstract void integrateBrLdOverCategories(int nodeIndex, double[] proportions, double[] integratedBrLd);

//    public abstract void calculateNodeBrLdOverCategories(int nodeIndex1, int nodeIndex2, int nodeIndex3, double[] proportions);

    public abstract void calculateNodeBrLdOverCategories(int nodeIndex1, int nodeIndex3, double[] proportions);

    //    public abstract void calculateLogLikelihoods(double[] rootBranchLd, double[] frequencies, double[] siteLogLikelihoods);
    public abstract double calculateLogLikelihoods(double[] frequencies);


    //    @Override
//    public void setNodeBranchLd(int nodeIndex, double[] branchLd) {
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
    public abstract void getNodeBranchLd(int nodeIndex, double[] branchLdOut);

    public abstract void setNodeBrLdForUpdate(int nodeIndex);

    public abstract int[] getNodeStates(int nodeIndex);

    public abstract void setInternalNodeStates(int nodeIndex, int[] states);
}
