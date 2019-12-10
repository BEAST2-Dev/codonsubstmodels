package beast.evolution.likelihood;


import beast.evolution.tree.Node;

/**
 * data augmentation likelihood core based on a branch for multithreading.
 * the branch is defined above the selected child node.
 *
 */
public class DABranchLikelihoodCore extends AbstrDABranchLikelihoodCore {

    protected boolean useScaling = false;

    protected double[][] scalingFactors; //TODO

    static final public double scalingThreshold = 1.0E-150; // MAX_VALUE 1.7*10^308
    double SCALE = 2;

    /**
     * data augmentation likelihood core based on a branch
     * @param node          the child {@link Node} below the branch
     * @param nrOfStates    64 for codon
     * @param nrOfSites     codon sites
     * @param categoryCount site model gamma category count
     */
    public DABranchLikelihoodCore(Node node, int nrOfStates, int nrOfSites, int categoryCount) {
        super(node, nrOfStates, nrOfSites, categoryCount);
    } // called initialize()


    /**
     * initializes states, likelihood arrays.
     */
    @Override
    protected void initialize() {
        // the branch above the node
        // merged likelihood for all categories
        branchLd = new double[2][nrOfSites];

        // 2 means current and stored
        matrices = new double[2][nrOfCategories * matrixSize];
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

    @Override
    public void setNodeMatrixForUpdate() {
        currentMatrixIndex = 1 - currentMatrixIndex; // 0 or 1
    }


    /**
     * Sets probability matrix for a node
     */
    @Override
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
     * Allocates branch likelihood for a node
     */
//    @Override
//    public void createNodeBranchLd(int nodeIndex) {
//        this.branchLd[0][nodeIndex] = new double[branchLdSize];
//        this.branchLd[1][nodeIndex] = new double[branchLdSize];
//    }

    /**
     * Sets branch likelihoods for a node
     */
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
    @Override
    public void getBranchLikelihoods(double[] branchLdOut) {
        System.arraycopy(branchLd[currentBrLdIndex], 0, branchLdOut, 0, branchLdOut.length);
    }

    @Override
    public void setBranchLdForUpdate() {
        currentBrLdIndex = 1 - currentBrLdIndex; // 0 or 1
    }

    /**
     * It calculates branch likelihoods at a branch,
     * and integrates branch likelihood across categories.
     *
     * @param childNodeStates   states at one 'child' node
     * @param parentNodeStates  states at 'parent' node
     * @param proportions the proportions of sites in each category. length = nrOfCategories.
     */
    //TODO need to cache per site to avoid recalculation, when only the sequence at a site is changed
    @Override
    public void calculateBranchLdOverCategories(final int[] childNodeStates, final int[] parentNodeStates,
                                                double[] proportions) {

        /**
         * Calculate DA branch likelihood per site per node.
         * matrix P(t) is flattened to n = w + i * state + j,
         * where n is index of flattened transition probability matrix (double[] matrices?),
         * i is child state, j is parent state, w is the category index.
         */

        // transition probability matrix
        double[] matrices = this.matrices[currentMatrixIndex];
        // branch likelihoods
        double[] branchLd = this.branchLd[currentBrLdIndex];

        int state1;
        int state3;
        // integrate over categories
        for (int k = 0; k < getNrOfSites(); k++) {

            state1 = childNodeStates[k];
            state3 = parentNodeStates[k];

            if (state1 < nrOfStates) { // && state3 < nrOfStates && state2 < nrOfStates

                //branchLd[] is linked to the branch above the child node 1
                branchLd[k] = matrices[state1 * nrOfStates + state3] * proportions[0];

                for (int l = 1; l < nrOfCategories; l++) {
                    int w = l * matrixSize;
                    //n = w + i * state + j
                    branchLd[k] += matrices[w + state1 * nrOfStates + state3] * proportions[l];
                } // end l nrOfCategories

            } else {
                throw new UnsupportedOperationException("State out of range at site " + k +
                        ", child node state = " + state1 + " parent node state = " + state3);

                //TODO child 1 has a gap or unknown state so treat it as unknown

//                    for (int j = 0; j < nrOfStates; j++) {
//                        branchLd3[v] = 1.0;
//                        v++;
//                    }
            }

            if (branchLd[k] == 0) {
                for (int l = 0; l < nrOfCategories; l++)
                    System.err.println("Category " + l +  ": transition probability = " +
                            matrices[l * matrixSize + state1 * nrOfStates + state3]);
//                System.err.println
                throw new RuntimeException("\nBranch above node " + getNodeNr() + " likelihood = 0 !\n" +
                        "At site " + k + ", child node state = " + state1 + ", parent node state = " + state3);
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
     * static calculation engine.
     * Multiple site likelihoods at this branch, if at root, then multiple frequencies.
     * Log the product.
     * The input likelihoods here have been integrated across categories.
     * @param siteLd     likelihoods (not logged), which can be either site likelihoods or frequencies at root.
     * @return           logged likelihood
     */
    public static double integrateLogLikelihood(double[] siteLd, double scalingThreshold) {

        double product = 1.0;
        double logP = 0;

        // exclude root frequency prior
//TODO review
        for (int k = 0; k < siteLd.length; k++) { // siteLd.length == getNrOfSites()
            // internal nodes, excl root
            product *= siteLd[k];

            // hard code to log when product is too small, Double.MAX_VALUE 1.79...e+308
            if (product < scalingThreshold || siteLd[k] < scalingThreshold) {
                // important check before implement log scaling
                if (product == 0)
                    throw new RuntimeException("Likelihood product -Inf ! " +
                            "\nlikelihood = " + siteLd[k] + " at site " + k);

                logP += Math.log(product); //+ getLogScalingFactor(k); TODO
                product = 1.0;
            } // if

        } // end k

        // the rest
        if (product < 1.0)
            logP += Math.log(product); //+ getLogScalingFactor(k); TODO

        return logP;
    }


    /**
     * Calculates log likelihood at this branch, excluding root frequency prior.
     * The input branch likelihoods here have been integrated across categories.
     * frequency[] need to be added later in DATreeLikelihood
     */
    @Override
    public double calculateBranchLogLikelihood() {
        return DABranchLikelihoodCore.
                integrateLogLikelihood(branchLd[currentBrLdIndex], scalingThreshold);
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


    // ======= getters for unit tests =======

    // for test, = getNrOfSites();
    public int getBranchLdSize() {
        // branchLd[][root] == null
        return branchLd[0].length;
    }

} // class
