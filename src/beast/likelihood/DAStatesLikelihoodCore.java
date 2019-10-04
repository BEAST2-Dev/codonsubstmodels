package beast.likelihood;


import beast.evolution.likelihood.LikelihoodCore;
import beast.tree.InternalNodeStates;

/**
 * standard data augmentation likelihood core, uses no caching *
 *
 *
 *
 */
public class DAStatesLikelihoodCore extends LikelihoodCore {
    protected int nrOfStates; // e.g. 64
    protected int nrOfNodes;
    protected int nrOfSites; // e.g. number of codons
    protected int partialsSize; // = nrOfSites * nrOfMatrices;
    protected int matrixSize; // nrOfStates^2
    protected int nrOfMatrices; // number of categories

    // partial likelihood matrix(ices),
    // 1st dimension is matrix index (current, stored),
    // 2nd is node index,
    // 3rd is improved to nrOfSites * nrOfMatrices
    protected double[][][] partials;

    // states in nodes: 0-63
    // 1st is node index, 2nd is site index
    protected int[][] states;

    // transition probability matrix(ices), P
    // 1st dimension is matrix index (current, stored),
    // 2nd is node index,
    // 3rd is nrOfMatrices * matrixSize
    protected double[][][] matrices;
    // store the matrix index, instead of different matrices
    protected int[] currentMatrixIndex;
    protected int[] storedMatrixIndex;
    protected int[] currentPartialsIndex;
    protected int[] storedPartialsIndex;

    protected boolean useScaling = false;

    protected double[][][] scalingFactors;

    private double scalingThreshold = 1.0E-100;
    double SCALE = 2;

    public DAStatesLikelihoodCore(int nrOfStates) {
        this.nrOfStates = nrOfStates;
    } // c'tor



    /**
     * Calculates pattern log likelihoods at a node.
     *
     * @param partials          the partials used to calculate the likelihoods
     * @param frequencies       an array of state frequencies
     * @param outLogLikelihoods an array into which the likelihoods will go
     */
    @Override
	public void calculateLogLikelihoods(double[] partials, double[] frequencies, double[] outLogLikelihoods) {
        int v = 0;
        for (int k = 0; k < nrOfSites; k++) {
//TODO validate index
            double sum = 0.0;
//            for (int i = 0; i < nrOfStates; i++) {
            // TODO check here
            int rootNr = nrOfNodes - 1;
            int i = states[rootNr][k];
                sum += frequencies[i] * partials[v];
                v++;
//            }
            outLogLikelihoods[k] = Math.log(sum) + getLogScalingFactor(k);
        }
    }


    /**
     * initializes likelihood arrays.
     * @param nodeCount           the number of nodes in the tree
     * @param siteCount        the number of patterns
     * @param categoryCount   the number of matrices (i.e., number of categories)
     * @param useAmbiguities  flag to indicate that sites containing ambiguous states should be handled instead of ignored
     */
    @Override
	public void initialize(int nodeCount, int siteCount, int categoryCount, boolean integrateCategories, boolean useAmbiguities) {

        this.nrOfNodes = nodeCount;
        this.nrOfSites = siteCount;
        this.nrOfMatrices = categoryCount;

//        this.integrateCategories = integrateCategories;
//        if (integrateCategories) { // always True
//            partialsSize = siteCount * nrOfStates * categoryCount;
//        } else {
//            partialsSize = siteCount * nrOfStates;
//        }
        // improve to site * category
        partialsSize = siteCount * categoryCount;
        partials = new double[2][nodeCount][];

        currentMatrixIndex = new int[nodeCount];
        storedMatrixIndex = new int[nodeCount];

        currentPartialsIndex = new int[nodeCount];
        storedPartialsIndex = new int[nodeCount];

        states = new int[nodeCount][];

        for (int i = 0; i < nodeCount; i++) {
            partials[0][i] = null;
            partials[1][i] = null;

            states[i] = null;
        }

        matrixSize = nrOfStates * nrOfStates;

        // 2 means current and stored
        matrices = new double[2][nodeCount][categoryCount * matrixSize];
    }

    /**
     * cleans up and deallocates arrays.
     */
    @Override
	public void finalize() throws Throwable {
        nrOfNodes = 0;
        nrOfSites = 0;
        nrOfMatrices = 0;

        partials = null;
        currentPartialsIndex = null;
        storedPartialsIndex = null;
        states = null;
        matrices = null;
        currentMatrixIndex = null;
        storedMatrixIndex = null;

        scalingFactors = null;
    }

    @Override
    public void setUseScaling(double scale) {
        useScaling = (scale != 1.0);

        if (useScaling) {
            scalingFactors = new double[2][nrOfNodes][nrOfSites];
        }
    }

    //============ states for a node ============

    /**
     * Allocates states for a node
     */
    public void createNodeStates(int nodeIndex) {

        this.states[nodeIndex] = new int[nrOfSites];
    }

    /**
     * Sets states for a node
     * TODO @param tips            sequences in tips, {@link final Alignment}
     * TODO @param ins             sequences in internal nodes, {@link InternalNodeStates}
     */
    @Override
	public void setNodeStates(int nodeIndex, int[] states) {

        if (this.states[nodeIndex] == null) {
            createNodeStates(nodeIndex);
        }
        System.arraycopy(states, 0, this.states[nodeIndex], 0, nrOfSites);
    }

    /**
     * Gets states for a node
     */
    @Override
	public void getNodeStates(int nodeIndex, int[] states) {
        System.arraycopy(this.states[nodeIndex], 0, states, 0, nrOfSites);
    }


    //============ transition probability matrix ============

    @Override
    public void setNodeMatrixForUpdate(int nodeIndex) {
        currentMatrixIndex[nodeIndex] = 1 - currentMatrixIndex[nodeIndex];

    }


    /**
     * Sets probability matrix for a node
     */
    @Override
	public void setNodeMatrix(int nodeIndex, int matrixIndex, double[] matrix) {
        System.arraycopy(matrix, 0, matrices[currentMatrixIndex[nodeIndex]][nodeIndex],
                matrixIndex * matrixSize, matrixSize);
    }

//    public void setPaddedNodeMatrices(int nodeIndex, double[] matrix) {
//        System.arraycopy(matrix, 0, matrices[currentMatrixIndex[nodeIndex]][nodeIndex],
//                0, nrOfMatrices * matrixSize);
//    }


    /**
     * Gets probability matrix for a node
     */
    @Override
	public void getNodeMatrix(int nodeIndex, int matrixIndex, double[] matrix) {
        System.arraycopy(matrices[currentMatrixIndex[nodeIndex]][nodeIndex],
                matrixIndex * matrixSize, matrix, 0, matrixSize);
    }

    //============ partial likelihood ============

    /**
     * Allocates partials for a node
     */
    @Override
    public void createNodePartials(int nodeIndex) {

        this.partials[0][nodeIndex] = new double[partialsSize];
        this.partials[1][nodeIndex] = new double[partialsSize];
    }

    /**
     * Sets partials for a node
     */
    @Override
    public void setNodePartials(int nodeIndex, double[] partials) {

        if (this.partials[0][nodeIndex] == null) {
            createNodePartials(nodeIndex);
        }
        if (partials.length < partialsSize) {
            int k = 0;
            for (int i = 0; i < nrOfMatrices; i++) {
                System.arraycopy(partials, 0, this.partials[0][nodeIndex], k, partials.length);
                k += partials.length;
            }
        } else {
            System.arraycopy(partials, 0, this.partials[0][nodeIndex], 0, partials.length);
        }
    }

    @Override
    public void getNodePartials(int nodeIndex, double[] partialsOut) {
        System.arraycopy(partials[currentPartialsIndex[nodeIndex]][nodeIndex], 0, partialsOut, 0, partialsOut.length);
    }

    @Override
    public void setNodePartialsForUpdate(int nodeIndex) {
        currentPartialsIndex[nodeIndex] = 1 - currentPartialsIndex[nodeIndex];
    }

    /**
     * Calculates partial likelihoods at a node.
     *
     * @param nodeIndex1 the 'child 1' node
     * @param nodeIndex2 the 'child 2' node
     * @param nodeIndex3 the 'parent' node
     */
    @Override
    public void calculatePartials(int nodeIndex1, int nodeIndex2, int nodeIndex3) {
        if (states[nodeIndex1] != null && states[nodeIndex2] != null && states[nodeIndex3] != null) {
            calculateStatesStates(
                    states[nodeIndex1], matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1],
                    states[nodeIndex2], matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2],
                    states[nodeIndex3], partials[currentPartialsIndex[nodeIndex3]][nodeIndex3]);
        } else {
            throw new IllegalArgumentException("Every node must has states ! \n" + "child 1 " +
                    (states[nodeIndex1] != null) + ", child 2 " + (states[nodeIndex2] != null) +
                    ", parent " + (states[nodeIndex3] != null));
        }

        if (useScaling) {
            throw new UnsupportedOperationException("in dev");
//            scalePartials(nodeIndex3);
        }

//
//        int k =0;
//        for (int i = 0; i < patternCount; i++) {
//            double f = 0.0;
//
//            for (int j = 0; j < stateCount; j++) {
//                f += partials[currentPartialsIndices[nodeIndex3]][nodeIndex3][k];
//                k++;
//            }
//            if (f == 0.0) {
//                Logger.getLogger("error").severe("A partial likelihood (node index = " + nodeIndex3 + ", pattern = "+ i +") is zero for all states.");
//            }
//        }
    }


    /**
     * Calculates partial likelihoods at a node given both children have states.
     * matrix P(t) is flattened to n = i * state + j where n is index of double[] matrices?.
     */
    protected void calculateStatesStates(int[] stateIndex1, double[] matrices1,
                                         int[] stateIndex2, double[] matrices2,
                                         int[] stateIndex3, double[] partials3) {
        int v = 0;

        for (int l = 0; l < nrOfMatrices; l++) {

            for (int k = 0; k < nrOfSites; k++) {

                int state1 = stateIndex1[k];
                int state2 = stateIndex2[k];
                int state3 = stateIndex3[k];

                int w = l * matrixSize;

                if (state1 < nrOfStates && state2 < nrOfStates) {

//                    for (int i = 0; i < nrOfStates; i++) {
//
//                        partials3[v] = matrices1[w + state1] * matrices2[w + state2];
//
//                        v++;
//                        w += nrOfStates;
//                    }


                    //TODO validate index
                    partials3[k] = matrices1[w + state1*nrOfStates + state3] * matrices2[w + state2*nrOfStates + state3];



//                } else if (state1 < nrOfStates) {
//                    // child 2 has a gap or unknown state so treat it as unknown
//
//                    for (int i = 0; i < nrOfStates; i++) {
//
//                        partials3[v] = matrices1[w + state1];
//
//                        v++;
//                        w += nrOfStates;
//                    }
//                } else if (state2 < nrOfStates) {
//                    // child 2 has a gap or unknown state so treat it as unknown
//
//                    for (int i = 0; i < nrOfStates; i++) {
//
//                        partials3[v] = matrices2[w + state2];
//
//                        v++;
//                        w += nrOfStates;
//                    }
                } else {
                    throw new UnsupportedOperationException("in dev");
                    // both children have a gap or unknown state so set partials to 1

//                    for (int j = 0; j < nrOfStates; j++) {
//                        partials3[v] = 1.0;
//                        v++;
//                    }
                }
            } // end k  nrOfSites
        } // end l nrOfMatrices
    }


    @Override
    public void integratePartials(int nodeIndex, double[] proportions, double[] outPartials) {
        calculateIntegratePartials(partials[currentPartialsIndex[nodeIndex]][nodeIndex], proportions, outPartials);
    }

    /**
     * Integrates partials across categories.
     *
     * @param inPartials  the array of partials to be integrated
     * @param proportions the proportions of sites in each category
     * @param outPartials an array into which the partials will go
     */
    @Override
    protected void calculateIntegratePartials(double[] inPartials, double[] proportions, double[] outPartials) {

        int u = 0;
        int v = 0;
        for (int k = 0; k < nrOfSites; k++) {
//TODO validate index
//            for (int i = 0; i < nrOfStates; i++) {

                outPartials[u] = inPartials[v] * proportions[0];
                u++;
                v++;
//            }
        }


        for (int l = 1; l < nrOfMatrices; l++) {
            u = 0;

//            for (int k = 0; k < nrOfSites; k++) {

                for (int i = 0; i < nrOfStates; i++) {

                    outPartials[u] += inPartials[v] * proportions[l];
                    u++;
                    v++;
                }
//            }
        }
    }


    /**
     * Scale the partials at a given node. This uses a scaling suggested by Ziheng Yang in
     * Yang (2000) J. Mol. Evol. 51: 423-432
     * <p/>
     * This function looks over the partial likelihoods for each state at each pattern
     * and finds the largest. If this is less than the scalingThreshold (currently set
     * to 1E-40) then it rescales the partials for that pattern by dividing by this number
     * (i.e., normalizing to between 0, 1). It then stores the log of this scaling.
     * This is called for every internal node after the partials are calculated so provides
     * most of the performance hit. Ziheng suggests only doing this on a proportion of nodes
     * but this sounded like a headache to organize (and he doesn't use the threshold idea
     * which improves the performance quite a bit).
     *
     * @param nodeIndex
     */
    protected void scalePartials(int nodeIndex) {

        int u = 0;

        for (int i = 0; i < nrOfSites; i++) {
//TODO validate index
            double scaleFactor = 0.0;
            int v = u;
            for (int k = 0; k < nrOfMatrices; k++) {
//                for (int j = 0; j < nrOfStates; j++) {
                    if (partials[currentPartialsIndex[nodeIndex]][nodeIndex][v] > scaleFactor) {
                        scaleFactor = partials[currentPartialsIndex[nodeIndex]][nodeIndex][v];
                    }
                    v++;
//                }
//                v += (nrOfSites - 1) * nrOfStates;
                v += (nrOfSites - 1);
            }

            if (scaleFactor < scalingThreshold) {

                v = u;
                for (int k = 0; k < nrOfMatrices; k++) {
//                    for (int j = 0; j < nrOfStates; j++) {
                        partials[currentPartialsIndex[nodeIndex]][nodeIndex][v] /= scaleFactor;
                        v++;
//                    }
//                    v += (nrOfSites - 1) * nrOfStates;
                    v += (nrOfSites - 1);
                }
                scalingFactors[currentPartialsIndex[nodeIndex]][nodeIndex][i] = Math.log(scaleFactor);

            } else {
                scalingFactors[currentPartialsIndex[nodeIndex]][nodeIndex][i] = 0.0;
            }
//            u += nrOfStates;


        }
    }

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
//                logScalingFactor += scalingFactors[currentPartialsIndex[i]][i][patternIndex_];
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
        int[] tmp1 = currentMatrixIndex;
        currentMatrixIndex = storedMatrixIndex;
        storedMatrixIndex = tmp1;

        int[] tmp2 = currentPartialsIndex;
        currentPartialsIndex = storedPartialsIndex;
        storedPartialsIndex = tmp2;
    }

    @Override
	public void unstore() {
        System.arraycopy(storedMatrixIndex, 0, currentMatrixIndex, 0, nrOfNodes);
        System.arraycopy(storedPartialsIndex, 0, currentPartialsIndex, 0, nrOfNodes);
    }

    /**
     * Restore the stored state
     */
    @Override
    public void store() {
        System.arraycopy(currentMatrixIndex, 0, storedMatrixIndex, 0, nrOfNodes);
        System.arraycopy(currentPartialsIndex, 0, storedPartialsIndex, 0, nrOfNodes);
    }



} // class
