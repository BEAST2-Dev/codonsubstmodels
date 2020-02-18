package beast.evolution.likelihood;


import beast.evolution.tree.NodeStates;

import java.util.Arrays;

/**
 * data augmentation likelihood core based on a branch for multithreading.
 * the branch is defined above the selected child node.
 *
 */
public class DABranchLikelihoodCore extends AbstrDALikelihoodCore {
    final private int branchNr; // the node Nr of child node below the branch

    final protected int nrOfStates; // e.g. 64
    final protected int matrixSize; // nrOfStates^2
    final protected int nrOfSites; // e.g. number of codons
    protected int nrOfCategories; // number of categories

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
     * caching probability tables obtained from substitutionModel,
     * size = stateCount * stateCount
     */
    protected double[] probabilities;

    protected double[] iexp;

    protected boolean useScaling = false;

    protected double[][] scalingFactors; //TODO

    static final public double SCALING_THRESHOLD = 1.0E-150; // MAX_VALUE 1.7*10^308
    double SCALE = 2;

    /**
     * no initialization, for calculating site likelihoods at the root.
     * @param branchNr       the node Nr of child node below the branch
     * @param nrOfStates     number of states in the data, 64 for codon
     * @param nrOfSites      number of sites (codon)
     */
    public DABranchLikelihoodCore(int branchNr, int nrOfStates, int nrOfSites) {
        this.branchNr = branchNr;
        this.nrOfStates = nrOfStates;
        this.matrixSize = nrOfStates * nrOfStates;
        this.nrOfSites = nrOfSites;
    }

    /**
     * data augmentation likelihood core based on a branch,
     * initialize() is called inside.
     * @param branchNr       the node Nr of child node below the branch
     * @param nrOfStates     number of states in the data, 64 for codon
     * @param nrOfSites      number of sites (codon)
     * @param nrOfCategories number of categories in the site model
     */
    public DABranchLikelihoodCore(int branchNr, int nrOfStates, int nrOfSites, int nrOfCategories) {
        this(branchNr, nrOfStates, nrOfSites);
        this.nrOfCategories = nrOfCategories;

        initialize();
    } // called initialize()


    /**
     * initializes states, likelihood arrays.
     */
    @Override
    protected void initialize() {
        // the branch above the node
        // merged likelihood for all categories
        branchLd = new double[2][nrOfSites];

        // matrixSize = stateCount * stateCount;
        // transition probability matrix P, [2] are current and stored
        matrices = new double[2][nrOfCategories * matrixSize];

        // used to cache transition probability matrix P
        probabilities = new double[matrixSize];
        Arrays.fill(probabilities, 1.0);

        iexp = new double[nrOfStates * nrOfStates];
    }

    /**
     * @return the reference of array caching probability,
     *         which has to be here to make thread safe.
     */
    public double[] getProbRef() {
        return this.probabilities;
    }

    public double[] getIexpRef() {
        return this.iexp;
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
     * cleans up and deallocates arrays.
     */
    @Override
    public void finalize() throws Throwable {
        nrOfCategories = 0; // TODO review

        branchLd = null;
        currentBrLdIndex = -1;
        storedBrLdIndex = -1;
//        states = null;
//        internalNodeStates = null;
        matrices = null;
        currentMatrixIndex = 0;
        storedMatrixIndex = 0;

        scalingFactors = null;
    }

    @Override
    public void setUseScaling(double scale) {
        useScaling = (scale != 1.0);

        if (useScaling) {
            scalingFactors = new double[2][nrOfSites];
        }
    }

    //============ states for a node ============

//    /**
//     * get the states at the node, if tips, then use <code>tipStates[][]</code>,
//     * if internal nodes, then use {@link InternalNodeStates}.
//     * @param nodeIndex
//     * @return
//     * @see {@link InternalNodeStates#getNrStates(int)}
//     */
//    @Override
//    public int[] getNodeStates(int nodeIndex) {
//        throw new RuntimeException("Invalid access !");
//    }
//
//    @Override
//    public void setInternalNodeStates(int nodeIndex, int[] states) {
//        throw new RuntimeException("Invalid access !");
//    }


    //============ transition probability matrix ============

    /**
     * @return the current transition probability matrix P(t),
     *         which is flattened to n = w + i * state + j,
     *         where n is index of flattened transition probability matrix,
     *         normally i is parent state, j is child state, w is the category index.
     */
    public double[] getCurrentMatrix() {
        return matrices[currentMatrixIndex];
    }

    /**
     * use before {@link #setNodeMatrix(int, double[])}
     */
    public void setNodeMatrixForUpdate() {
        currentMatrixIndex = 1 - currentMatrixIndex; // 0 or 1
    }


    /**
     * Sets probability matrix for a node
     */
    public void setNodeMatrix(int categoryIndex, double[] matrix) {
        System.arraycopy(matrix, 0, matrices[currentMatrixIndex],
                categoryIndex * matrixSize, matrixSize);
    }


    /**
     * Gets probability matrix for a node
     */
//    @Override
//	public void getNodeMatrix(int nodeIndex, int categoryIndex, double[] matrix) {
//        System.arraycopy(matrices[currentMatrixIndex[nodeIndex]][nodeIndex],
//                categoryIndex * matrixSize, matrix, 0, matrixSize);
//    }

    //============ branch likelihood ============

    /**
     * use before {@link #calculateBranchLd(int[], int[], double[])}
     */
    public void setBranchLdForUpdate() {
        currentBrLdIndex = 1 - currentBrLdIndex; // 0 or 1
    }


    /**
     * Branch likelihood calculation with site model categories from state i to j.
     * The flattened transition probability matrix P(t) index : n = w + i * state + j.
     * @param i   state i in parent node.
     * @param j    state j in child node.
     * @param proportions       the proportions of sites in each category. length = nrOfCategories.
     * @return   the branch likelihood at one site
     */
    public double calculateBranchLdAtSite(final int i, final int j, final double[] proportions) {
        final double[] matrices = getCurrentMatrix();

        //branchLd is defined as the branch above the child node 1
        double branchLdSite = matrices[i * nrOfStates + j] * proportions[0];

        for (int l = 1; l < nrOfCategories; l++) {
            int w = l * matrixSize;
            //n = w + i * state + j
            branchLdSite += matrices[w + i * nrOfStates + j] * proportions[l];
        } // end l nrOfCategories
        return branchLdSite;
    }

    /**
     * It calculates branch likelihoods at a branch,
     * and integrates branch likelihood across categories.
     *  @param statesParentNode  states at 'parent' node
     * @param statesChildNode   states at one 'child' node
     * @param proportions the proportions of sites in each category. length = nrOfCategories.
     */
    //TODO need to cache per site to avoid recalculation, when only the sequence at a site is changed
    public void calculateBranchLd(final int[] statesParentNode, final int[] statesChildNode,
                                  double[] proportions) {

        /**
         * Calculate DA branch likelihood per site per node.
         * matrix P(t) is flattened to n = w + i * state + j,
         * where n is index of flattened transition probability matrix (double[] matrices?),
         * i is parent state, j is child state, w is the category index.
         */

        // branch likelihoods for calculating
        double[] branchLd = this.branchLd[currentBrLdIndex];

        int i; // parent state
        int j; // child state
        // integrate over categories
        for (int k = 0; k < getNrOfSites(); k++) {

            i = statesParentNode[k]; // parent
            j = statesChildNode[k];

            if (j < nrOfStates) { // && state3 < nrOfStates && state2 < nrOfStates

                // sum_l(P_s3_s1(t) * proportions[l]), from i to j
                branchLd[k] = calculateBranchLdAtSite(i, j, proportions);

//                //branchLd[] is linked to the branch above the child node 1
//                branchLd[k] = matrices[state1 * nrOfStates + state3] * proportions[0];
//
//                for (int l = 1; l < nrOfCategories; l++) {
//                    int w = l * matrixSize;
//                    //n = w + i * state + j
//                    branchLd[k] += matrices[w + state1 * nrOfStates + state3] * proportions[l];
//                } // end l nrOfCategories

            } else {
                throw new UnsupportedOperationException("State out of range at site " + k +
                        ", child node state = " + j + " parent node state = " + i);

                //TODO child 1 has a gap or unknown state so treat it as unknown

//                    for (int j = 0; j < nrOfStates; j++) {
//                        branchLd3[v] = 1.0;
//                        v++;
//                    }
            }

            if (branchLd[k] == 0) {
                // transition probability matrix
                double[] matrices = getCurrentMatrix();
                for (int l = 0; l < nrOfCategories; l++)
                    System.err.println("Category " + l +  ": transition probability = " +
                            matrices[l * matrixSize + j * nrOfStates + i]);
//                System.err.println
                throw new RuntimeException("\nBranch above node " + getBranchNr() + " likelihood = 0 !\n" +
                        "At site " + k + ", child node state = " + j + ", parent node state = " + i);
            }

        } // end k  getNrOfSites()


//        } else {
//            throw new IllegalArgumentException("Every node must has states and length of " + getNrOfSites() + " ! \n" +
//                    "child 1 " + childNodeStates.length + ", parent " + node3States.length);
//        }


//        final double[] tmp = branchLd[currentBrLdIndex[nodeChild1]][nodeChild1];
//        for (int i=0; i < tmp.length; i++) {
//            if (tmp[i] == 0) {
//                // branchLd.length = getNrOfSites() * nrOfCategories, states.length = getNrOfSites()
//                int site = i % getNrOfSites();
//                System.err.println("i = " + i + " site = " + site + " : " +
//                        "child 1 id = " + nodeChild1 + " state = " + stateIndex1[site] +
//                        //", child 2 id = " + nodeIndex2 + " state = " + stateIndex2[site] +
//                        ", parent id = " + nodeParent + " state = " + stateIndex3[site] +
//                        ", branchLd[" + i + "] = " + tmp[i]);
//            }
//        }

        if (useScaling) {
            throw new UnsupportedOperationException("in dev");
//            scaleBranchLds(nodeParent);
        }

//
//        int k =0;
//        for (int i = 0; i < patternCount; i++) {
//            double f = 0.0;
//
//            for (int j = 0; j < stateCount; j++) {
//                f += branchLd[currentPartialsIndices[nodeParent]][nodeParent][k];
//                k++;
//            }
//            if (f == 0.0) {
//                Logger.getLogger("error").severe("A branch likelihood (node index = " + nodeParent + ", pattern = "+ i +") is zero for all states.");
//            }
//        }
    }

    /**
     * Log likelihood calculation engine.
     * Multiple given site likelihoods, and then log the product.
     * The likelihood input has been integrated across categories.
     * It can be also used for the site likelihoods at the root.
     * @param siteLd     likelihoods (not logged) by sites (codons).
//     * @param scalingThreshold   threshold for the product to call {@link Math#log(double)}.
     * @return           logged likelihood
     */
    public double logIntegratedLikelihood(double[] siteLd) {

        double product = 1.0;
        double logP = 0;

        // siteLd.length == getNrOfSites()
        for (int k = 0; k < siteLd.length; k++) {
            // multiple (not logged) site likelihoods
            product *= siteLd[k];

            // hard code to log when product is too small, Double.MAX_VALUE 1.79...e+308
            if (product < SCALING_THRESHOLD || siteLd[k] < SCALING_THRESHOLD) {
                // important check before implement log scaling
                if (product == 0)
                    throw new RuntimeException("Likelihood product -Inf ! " +
                            "\nlikelihood = " + siteLd[k] + " at site " + k);

                logP += Math.log(product); //+ getLogScalingFactor(k); TODO
                product = 1.0;
            } // if

        } // end k

        // log the rest
        if (product < 1.0)
            logP += Math.log(product); //+ getLogScalingFactor(k); TODO

        return logP;
    }


    /**
     * Calculates log likelihood at this branch, excluding root frequency prior.
     * The input branch likelihoods here have been integrated across categories.
     * frequency[] need to be added later in DATreeLikelihood
     */
    public double calculateBranchLogLikelihood() {
        return logIntegratedLikelihood(branchLd[currentBrLdIndex]);
    }

    // log likelihood at root given codon frequencies
    public double calculateRootLogLikelihood(int rootIndex, NodeStates rootStates, double[] frequencies) {
        int siteCount = rootStates.getSiteCount();
        double[] siteLdAtRoot = new double[siteCount];
        int state;
        for (int k = 0; k < siteCount; k++) {
            // hard code for root node
            state = rootStates.getState(k); // 0-63
            siteLdAtRoot[k] = frequencies[state];
        }

        return logIntegratedLikelihood(siteLdAtRoot); //+ getLogScalingFactor(k); TODO
    }


    /**
     * Scale the branch likelihoods at a given node. This uses a scaling suggested by Ziheng Yang in
     * Yang (2000) J. Mol. Evol. 51: 423-432
     * <p/>
     * This function looks over the branch likelihoods for each state at each pattern
     * and finds the largest. If this is less than the scalingThreshold (currently set
     * to 1E-40) then it rescales the branch likelihoods for that pattern by dividing by this number
     * (i.e., normalizing to between 0, 1). It then stores the log of this scaling.
     * This is called for every internal node after the branch likelihoods are calculated so provides
     * most of the performance hit. Ziheng suggests only doing this on a proportion of nodes
     * but this sounded like a headache to organize (and he doesn't use the threshold idea
     * which improves the performance quite a bit).
     *
     * @param nodeIndex
     */
//    protected void scaleBranchLds(int nodeIndex) {
//
//        int u = 0;
//
//        for (int i = 0; i < getNrOfSites(); i++) {
////TODO validate index
//            double scaleFactor = 0.0;
//            int v = u;
//            for (int k = 0; k < nrOfCategories; k++) {
////                for (int j = 0; j < nrOfStates; j++) {
//                    if (branchLd[currentBrLdIndex[nodeIndex]][nodeIndex][v] > scaleFactor) {
//                        scaleFactor = branchLd[currentBrLdIndex[nodeIndex]][nodeIndex][v];
//                    }
//                    v++;
////                }
////                v += (getNrOfSites() - 1) * nrOfStates;
//                v += (getNrOfSites() - 1);
//            }
//
//            if (scaleFactor < scalingThreshold) {
//
//                v = u;
//                for (int k = 0; k < nrOfCategories; k++) {
////                    for (int j = 0; j < nrOfStates; j++) {
//                        branchLd[currentBrLdIndex[nodeIndex]][nodeIndex][v] /= scaleFactor;
//                        v++;
////                    }
////                    v += (getNrOfSites() - 1) * nrOfStates;
//                    v += (getNrOfSites() - 1);
//                }
//                scalingFactors[currentBrLdIndex[nodeIndex]][nodeIndex][i] = Math.log(scaleFactor);
//
//            } else {
//                scalingFactors[currentBrLdIndex[nodeIndex]][nodeIndex][i] = 0.0;
//            }
////            u += nrOfStates;
//
//
//        }
//    }

    /**
     * This function returns the scaling factor for that pattern by summing over
     * the log scalings used at each node. If scaling is off then this just returns
     * a 0.
     *
     * @return the log scaling factor
     */
    @Override
    public double getLogScalingFactor(int patternIndex_) {
//    	if (m_bUseScaling) {
//    		return -(m_nNodeCount/2) * Math.log(SCALE);
//    	} else {
//    		return 0;
//    	}        
        double logScalingFactor = 0.0;
        if (useScaling) {
            throw new UnsupportedOperationException("in development");
//            for (int i = 0; i < nrOfNodes; i++) {
//                logScalingFactor += scalingFactors[currentBrLdIndex[i]][i][patternIndex_];
//            }
        }
        return logScalingFactor;
    }


    // suppose only used by unit test
    public void getBranchLikelihoods(double[] branchLdOut) {
        System.arraycopy(branchLd[currentBrLdIndex], 0, branchLdOut, 0, branchLdOut.length);
    }

    // ======= getters =======

    public int getNrOfSites() {
        return nrOfSites;
    }

    public int getBranchNr() {
        return branchNr;
    }

    // = getNrOfSites(), for unit test
    public int getBranchLdSize() {
        // branchLd[][root] == null
        return branchLd[0].length;
    }

    public int getNrOfStates() {
        return nrOfStates;
    }

    // transition probability matrix size = nrOfStates^2
    public int getMatrixSize() {
        return matrixSize;
    }

    // nrOfCategories
    public int getNrOfCategories() {
        return nrOfCategories;
    }
} // class
