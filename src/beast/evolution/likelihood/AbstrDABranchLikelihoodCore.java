

package beast.evolution.likelihood;

/**
 * framework of data augmentation branch likelihood core
 */

public abstract class AbstrDABranchLikelihoodCore extends AbstrDALikelihoodCore {
//    private final Node node; //TODO bug : destroy store/restore
    private final int branchNr; // the node Nr of child node below the branch

    final protected int nrOfSites; // e.g. number of codons

    // to store branch likelihood calculation per site:
    // 1st dimension is matrix index (current, stored),
    // 2nd is getNrOfSites(), Ld across categories are integrated
    protected double[][] branchLd;

    // transition probability matrix(ices), P
    // 1st dimension is matrix index (current, stored),
    // 2nd is nrOfCategories * matrixSize
    protected double[][] matrices;
    // store the matrix index, instead of different matrices
    protected int currentMatrixIndex = 0;
    protected int storedMatrixIndex = 0;
    protected int currentBrLdIndex = 0;
    protected int storedBrLdIndex = 0;


    /**
     * data augmentation likelihood core based on a branch,
     * called initialize() inside
     * @param branchNr       the node Nr of child node below the branch
     * @param nrOfStates     number of states in the data, 64 for codon
     * @param nrOfSites      number of sites (codon)
     * @param nrOfCategories number of categories in the site model
     */
    public AbstrDABranchLikelihoodCore(int branchNr, int nrOfStates, int nrOfSites, int nrOfCategories) {
        super(nrOfStates, nrOfCategories);
        this.nrOfSites = nrOfSites;
        this.branchNr = branchNr;

        initialize();
    }


    /**
     * Restore the stored state
     */
    @Override
    public void store() {
        storedMatrixIndex = currentMatrixIndex;
        storedBrLdIndex = currentBrLdIndex;
    }

    /**
     * Store current state
     */
    @Override
    public void restore() {
        // Rather than copying the stored stuff back, just swap the pointers...
        int tmp1 = currentMatrixIndex;
        currentMatrixIndex = storedMatrixIndex;
        storedMatrixIndex = tmp1;

        int tmp2 = currentBrLdIndex;
        currentBrLdIndex = storedBrLdIndex;
        storedBrLdIndex = tmp2;
    }

    @Override
    public void unstore() {
        currentMatrixIndex = storedMatrixIndex;
        currentBrLdIndex = storedBrLdIndex;
    }

    /**
     * @return the current transition probability matrix P(t),
     *         which is flattened to n = w + i * state + j,
     *         where n is index of flattened transition probability matrix,
     *         i is child state, j is parent state, w is the category index.
     */
    public double[] getCurrentMatrix() {
        return matrices[currentMatrixIndex];
    }

    /**
     * branch likelihood calculation with site model categories.
     * @param stateParentNode   state in parent node, j.
     * @param stateChildNode    state in child node, i.
     * @param proportions       the proportions of sites in each category. length = nrOfCategories.
     * @return   the branch likelihood at one site
     */
    public double getBranchLdAtSite(final int stateParentNode, final int stateChildNode,
                                    final double[] proportions) {
        final double[] matrices = getCurrentMatrix();

        //branchLd is defined as the branch above the child node 1
        double branchLdSite = matrices[stateChildNode * nrOfStates + stateParentNode] * proportions[0];

        for (int l = 1; l < nrOfCategories; l++) {
            int w = l * matrixSize;
            //n = w + i * state + j
            branchLdSite += matrices[w + stateChildNode * nrOfStates + stateParentNode] * proportions[l];
        } // end l nrOfCategories
        return branchLdSite;
    }


    /**
     * indicate that the probability transition matrix for node
     * nodeIndex is about the be changed, that is, that the stored
     * state for node nodeIndex cannot be reused *
     */
    public abstract void setNodeMatrixForUpdate();

    /**
     * assign values of states for probability transition matrix for node *
     */
    public abstract void setNodeMatrix(int categoryIndex, double[] matrix);
//    public abstract void setNodeMatrix(double[] matrix);

//    public abstract void calculateNodeBranchLd(int childNum1, int childNum2, int nodeIndex);
//    public abstract void integrateBrLdOverCategories(double[] proportions, double[] integratedBrLd);

//    public abstract void calculateNodeBrLdOverCategories(int nodeIndex1, int nodeIndex2, int nodeIndex3, double[] proportions);

    public abstract void calculateBranchLdOverCategories(int[] parentNodeStates, int[] childNodeStates, double[] proportions);

    //    public abstract void calculateLogLikelihoods(double[] rootBranchLd, double[] frequencies, double[] siteLogLikelihoods);
    public abstract double calculateBranchLogLikelihood();


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
    public abstract void getBranchLikelihoods(double[] branchLdOut);

    public abstract void setBranchLdForUpdate();


    public int getNrOfSites() {
        return nrOfSites;
    }

    public int getBranchNr() {
        return branchNr;
    }
}
