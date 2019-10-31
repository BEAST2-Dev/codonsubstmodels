package beast.evolution.likelihood;


/**
 * general data augmentation likelihood core based on a branch
 * for multithreading
 *
 */
public class DAStatesBranchLikelihoodCore extends DABranchLikelihoodCore {
    final protected int nrOfStates; // e.g. 64
    final protected int matrixSize; // nrOfStates^2
//    final protected int nodeNr; // TODO need nodeNr?
    final protected boolean isTip;
    final protected int nrOfSites; // e.g. number of codons
    protected int nrOfCategories; // number of categories
//    protected int branchLdSize; // = getNrOfSites() * nrOfCategories;

    // to store branch likelihood calculation per site:
    // 1st dimension is matrix index (current, stored),
    // 2nd is getNrOfSites(), Ld across categories are integrated
    protected double[][] branchLd;

    // states in tip/internal nodes: 0-63
    // 1st is node index, 2nd is site index
//    protected int[][] tipStates;
//    protected InternalNodeStates internalNodeStates;

    // transition probability matrix(ices), P
    // 1st dimension is matrix index (current, stored),
    // 2nd is nrOfCategories * matrixSize
    protected double[][] matrices;
    // store the matrix index, instead of different matrices
    protected int currentMatrixIndex = -1;
    protected int storedMatrixIndex = -1;
    protected int currentBrLdIndex = -1;
    protected int storedBrLdIndex = -1;

    protected boolean useScaling = false;

    protected double[][] scalingFactors;

    private double scalingThreshold = 1.0E-100;
    double SCALE = 2;

    public DAStatesBranchLikelihoodCore(boolean isTip, int nrOfStates, int nrOfSites, int categoryCount) {
        this.isTip = isTip;
        this.nrOfStates = nrOfStates;
        matrixSize = nrOfStates * nrOfStates;
        this.nrOfSites = nrOfSites;
        this.nrOfCategories = categoryCount;

        initialize();
    } // c'tor


    /**
     * initializes states, likelihood arrays.
     */
    @Override
	protected void initialize() {
//        this.tipStates = tipStates;
//        this.internalNodeStates = internalNodeStates;
//        this.nrOfCategories = categoryCount;

//        if (getNrOfInterNodes() != getNrOfTips() - 1)
//            throw new IllegalArgumentException("Internal nodes " + getNrOfInterNodes() + " != tips " + getNrOfTips() + " - 1 !");
//        if (getNrOfSites() != internalNodeStates.getSiteCount())
//            throw new IllegalArgumentException("Tip site count " + getNrOfSites() +
//                    " != internal node site count " + internalNodeStates.getSiteCount() + " !");


//        this.nrOfNodes = nodeCount;
//        this.nrOfSites = siteCount;
//        this.integrateCategories = integrateCategories;
//        if (integrateCategories) { // always True
//            branchLdSize = siteCount * nrOfStates * categoryCount;
//        } else {
//            branchLdSize = siteCount * nrOfStates;
//        }
        // 1 likelihood for all sites
//        branchLd = new double[2][];

//        currentMatrixIndex = new int;
//        storedMatrixIndex = new int;
//
//        currentBrLdIndex = new int;
//        storedBrLdIndex = new int;

//        states = new int[];

        // init branchLd[][][]
//        for (int i = 0; i < getNrOfTips(); i++) {
//            // tips
//            branchLd[0][i] = null;
//            branchLd[1][i] = null;
//        }
//        final int branchLdSize = categoryCount * getNrOfSites();
//        for (int i = 0; i < getNrOfNodes()-1; i++) {
//        }
        // 1 likelihood for all sites
        if (isTip) {
            branchLd = new double[2][];
            branchLd[0] = null;
            branchLd[1] = null;
        } else {
            // internal nodes
            branchLd = new double[2][getNrOfSites()];
        }

        // 2 means current and stored
        matrices = new double[2][nrOfCategories * matrixSize];
    }

    /**
     * cleans up and deallocates arrays.
     */
    @Override
	public void finalize() throws Throwable {
        nrOfCategories = 0;

        branchLd = null;
        currentBrLdIndex = -1;
        storedBrLdIndex = -1;
//        states = null;
//        internalNodeStates = null;
        matrices = null;
        currentMatrixIndex = -1;
        storedMatrixIndex = -1;

        scalingFactors = null;
    }

    @Override
    public void setUseScaling(double scale) {
        useScaling = (scale != 1.0);

        if (useScaling) {
            scalingFactors = new double[2][getNrOfSites()];
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
    public void getNodeBranchLd(double[] branchLdOut) {
        System.arraycopy(branchLd[currentBrLdIndex], 0, branchLdOut, 0, branchLdOut.length);
    }

    @Override
    public void setNodeBrLdForUpdate() {
        currentBrLdIndex = 1 - currentBrLdIndex; // 0 or 1
    }

    /**
     * It calculates branch likelihoods at a node,
     * and integrates branch likelihood across categories.
     *
//     * @param nodeChild1  the 'child 1' node
//     * @param nodeParent  the 'parent' node
     * @param proportions the proportions of sites in each category. length = nrOfCategories.
     */
    // branchLd[] is linked to the branch above the child node 1
    @Override
    public void calculateNodeBrLdOverCategories(final int[] childNodeStates, final int[] parentNodeStates,
                                                double[] proportions) {
//        final int[] node1States = getNodeStates(nodeChild1);
//        final int[] node3States = getNodeStates(nodeParent);

//        if (node1States.length == getNrOfSites() && node3States.length == getNrOfSites()) {

            /**
             * Calculate DA branch likelihood per site per node.
             * matrix P(t) is flattened to n = w + i * state + j,
             * where n is index of flattened transition probability matrix (double[] matrices?),
             * i is child state, j is parent state, w is the category index.
             */

//            double[] matrices1 = matrices[currentMatrixIndex[nodeChild1]][nodeChild1];
        double[] matrices1 = matrices[currentMatrixIndex];
            //branchLd[][] is linked to the branch above the child node 1
//            double[] branchLd1 = branchLd[currentBrLdIndex[nodeChild1]][nodeChild1];
        double[] branchLd1 = branchLd[currentBrLdIndex];

            int state1;
            int state3;
            // integrate over categories
            for (int k = 0; k < getNrOfSites(); k++) {

                state1 = childNodeStates[k];
                state3 = parentNodeStates[k];

                if (state1 < nrOfStates) { // && state3 < nrOfStates && state2 < nrOfStates

                    //branchLd[] is linked to the branch above the child node 1
                    branchLd1[k] = matrices1[state1 * nrOfStates + state3] * proportions[0];

                    for (int l = 1; l < nrOfCategories; l++) {
                        int w = l * matrixSize;
//                    siteLd3 = matrices1[w + state1 * nrOfStates + state3] * matrices2[w + state2 * nrOfStates + state3];
                        branchLd1[k] += matrices1[w + state1 * nrOfStates + state3] * proportions[l];

                    } // end l nrOfCategories

                } else {
                    System.err.println("k = " + k + " state1 = " + state1 + " state3 = " + state3);
                    throw new UnsupportedOperationException("in dev");

                    //TODO child 1 has a gap or unknown state so treat it as unknown

//                    for (int j = 0; j < nrOfStates; j++) {
//                        branchLd3[v] = 1.0;
//                        v++;
//                    }
                }

                if (branchLd1[k] == 0) {
                    System.err.println("Site = " + k + " : state1 = " + state1 + " state3 = " + state3);
                    for (int l = 0; l < nrOfCategories; l++)
                        System.err.println("Product at category " + l + " = " +
                                matrices1[l * matrixSize + state1 * nrOfStates + state3]);
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
     * Calculates site log likelihoods at root node.
     * The input branch likelihoods here have been integrated across categories.
     *
     */
    @Override
    public double calculateLogLikelihoods(boolean isRoot, double frequency) {

//        int trunk = 1;
        double product = 1.0;
        if (isRoot) product = frequency;
        double logP = 0;

        // exclude root node
//        final int rootIndex = getNrOfNodes()-1;
//TODO review
        for (int k = 0; k < getNrOfSites(); k++) {
//            for (int i = 0; i < rootIndex; i++) {
                // internal nodes, excl root
                product *= branchLd[currentBrLdIndex][k];

                // log when product is too small, Double.MAX_VALUE 1.79...e+308
                if (product < 1e-200) {
                    if (product == 0)
                        throw new RuntimeException("Likelihood -Inf at site " + k + " ! " +
                                "\nbranch likelihood = " + branchLd[currentBrLdIndex][k]);

                    logP += Math.log(product); //+ getLogScalingFactor(k); TODO
                    product = 1.0;
//                trunk++;
//                }
            } // end i

            // hard code for root node
//            int state = internalNodeStates.getASite(rootIndex, k); // 0-63

//            if (frequencies[state] == 0)
//                throw new RuntimeException("frequencies[" + state + "] == 0 refers to stop codon, check the state or frequencies !");

        } // end k

        // the rest
        if (product < 1)
            logP += Math.log(product); //+ getLogScalingFactor(k); TODO

        return logP;
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
//        System.arraycopy(storedMatrixIndex, 0, currentMatrixIndex, 0, getNrOfNodes());
//        System.arraycopy(storedBrLdIndex, 0, currentBrLdIndex, 0, getNrOfNodes());
    }

    /**
     * Restore the stored state
     */
    @Override
    public void store() {
        storedMatrixIndex = currentMatrixIndex;
        storedBrLdIndex = currentBrLdIndex;
//        System.arraycopy(currentMatrixIndex, 0, storedMatrixIndex, 0, getNrOfNodes());
//        System.arraycopy(currentBrLdIndex, 0, storedBrLdIndex, 0, getNrOfNodes());
    }

    // ======= getters for unit tests =======

    public int getNrOfStates() {
        return nrOfStates;
    }

    public int getNrOfSites() {
        return nrOfSites;
    }

//    public int getNrOfTips() {
//        return tipStates.length;
//    }
//
//    public int getNrOfInterNodes() {
//        return internalNodeStates.getInternalNodeCount();
//    }
//
//    public int getNrOfNodes() {
//        return getNrOfTips() + getNrOfInterNodes();
//    }

    // nrOfCategories
    public int getNrOfCategories() {
        return nrOfCategories;
    }

    // = getNrOfSites();
    public int getBranchLdSize() {
        // branchLd[][root] == null
        return branchLd[0].length;
    }

    // transition probability matrix size = nrOfStates^2
    public int getMatrixSize() {
        return matrixSize;
    }


} // class
