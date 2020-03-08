package beast.evolution.operators;

import beast.core.Description;
import beast.evolution.likelihood.DABranchLikelihoodCore;
import beast.evolution.likelihood.DataAugTreeLikelihood;
import beast.evolution.tree.Node;
import beast.evolution.tree.NodeStatesArray;
import beast.util.RandomUtils;
import beast.util.Randomizer;

/**
 * Gibbs sampler to sample the internal node states.
 *
 * <p><code>w_i</code> is the ith state of total 60/61 states:
 * <p><code>w_i ~ P_z{w_i}(t) * P_{w_i}x(t) * P_{w_i}y(t)</code>
 * <p>At the root:
 * <p><code>w_i ~ P_{w_i}(t) * P_{w_i}x(t) * P_{w_i}y(t)</code>,
 * where P_{w_i}(t) is the equilibrium frequency.
 * <p>Renormalise all <code>w_i</code> and choose <code>w</code> from the distribution.
 *
 * @author Walter Xie
 */
@Description("Gibbs sampler to sample the internal node states.")
public class GibbsSampler {

    /*** they need to be updated each proposal ***/
    // DATreeLikelihood is calculation engine
    private DataAugTreeLikelihood daTreeLd;
    private Node node; // the node Gibbs operates on
    private NodeStatesArray nodesStates; // to get sequences, never change it here

    // control which node will be operated, if only working on a node
//    private int opNodeNr = -1;

    // avoid to create new T[]
    private int[] newStates; // proposed states
    private double[] cpd_w;  // cumulative probability distribution of states at w node

    private final int startInclusive, endExclusive;

    /**
     * Init caching arrays.
     * @param stateCount
     * @param startInclusive
     * @param endExclusive
     */
    public GibbsSampler(int stateCount, int startInclusive, int endExclusive) {
        assert endExclusive > startInclusive;
        this.startInclusive = startInclusive;
        this.endExclusive = endExclusive;

        newStates = new int[endExclusive - startInclusive];
        cpd_w = new double[stateCount];
    }

    /**
     * node, nodesStates need to be updated each proposal.
     * @param node
     */
    public void update(Node node, NodeStatesArray nodesStates, DataAugTreeLikelihood daTreeLd) {
        this.node = node;
        this.nodesStates = nodesStates;
        this.daTreeLd = daTreeLd;
    }

    /**
     * This method is used for multithreading.
     * Need to {@link #update(Node, NodeStatesArray, DataAugTreeLikelihood)} before this.
     * @return
     */
    public void gibbsSampling() {
        gibbsSampling(newStates, node, nodesStates, daTreeLd, startInclusive, endExclusive, cpd_w);
//        return newStates;
    }

    /**
     * node, nodesStates need to be updated each proposal.
     * @param node
     * @param nodesStates
     * @return
     */
    public int[] gibbsSampling(Node node, NodeStatesArray nodesStates, DataAugTreeLikelihood daTreeLd) {
        gibbsSampling(newStates, node, nodesStates, daTreeLd, startInclusive, endExclusive, cpd_w);
        return newStates;
    }

    public int[] getNewStates() {
        return newStates;
    }

    /**
     * Gibbs sampling all sites for a given node, and set new states.
     * @param newStates   proposed states from Gibbs sampling
     * @param node        {@link Node}
     * @param nodesStates {@link NodeStatesArray}
     * @param daTreeLd    {@link DataAugTreeLikelihood}
     * @param startInclusive   the site (codon) index to start in loop
     * @param endExclusive     the site (codon) index to end in loop
     * @param cpd_w        avoid to create new double[]
     * @see #gibbsSampling(int, int, DABranchLikelihoodCore, DABranchLikelihoodCore,
     *                     double[], double[], double[])
     * @see #gibbsSampling(int, int, int, DABranchLikelihoodCore, DABranchLikelihoodCore,
     *                     DABranchLikelihoodCore, double[], double[])
     */
    public void gibbsSampling(int[] newStates, final Node node, final NodeStatesArray nodesStates,
                              final DataAugTreeLikelihood daTreeLd, int startInclusive,
                              int endExclusive, double[] cpd_w) {
        // get parameters
        final double[] frequencies = daTreeLd.getSubstitutionModel().getFrequencies();
        final double[] proportions = daTreeLd.getSiteModel().getCategoryProportions(node);

        final int nodeNr = node.getNr();
        final int ch1Nr = node.getChild(0).getNr();
        final int ch2Nr = node.getChild(1).getNr();

        // w-x branch
        final DABranchLikelihoodCore wxBranchLd = daTreeLd.getDaBranchLdCores(ch1Nr);
        // w-y branch
        final DABranchLikelihoodCore wyBranchLd = daTreeLd.getDaBranchLdCores(ch2Nr);

        final boolean isRoot = node.isRoot();
        // z-w branch
        int parentNr = -1;
        DABranchLikelihoodCore zwBranchLd = null;
        if (!isRoot) {
            // z-w branch
            zwBranchLd = daTreeLd.getDaBranchLdCores(nodeNr);
            parentNr = node.getParent().getNr();
        }

        // sampling all sites for a given node
//        int oldState;
//        int changes = 0;
        int x,y,z,i=0;
        for (int k = startInclusive; k < endExclusive; k++) { // the site (codon) index
            // states at child nodes
            x = nodesStates.getState(ch1Nr, k);
            y = nodesStates.getState(ch2Nr, k);
            // match index i in newStates[] to the site index k in the whole sequence
            i = k - startInclusive;
            if (isRoot) {
                newStates[i] = gibbsSampling(x, y, wxBranchLd, wyBranchLd,
                        cpd_w, frequencies, proportions);
            } else {
                z = nodesStates.getState(parentNr, k);

                newStates[i] = gibbsSampling(x, y, z, wxBranchLd, wyBranchLd, zwBranchLd,
                        cpd_w, proportions);
            }

//            oldState = nodesStates.getState(nodeNr, k);
//            if (newStates[k] != oldState) {
//                System.out.println("Node " + nodeNr + " site " + k + " : state changes from " +
//                        oldState + " to " + newStates[k]);
//                changes++;
//            }
        }
//        System.out.println("Node " + nodeNr + " changed " + changes + " sites.");
    }

    /**
     * Gibbs sampling state at 1 site and the root.
     * As no parent, the equilibrium frequencies are used.
     * @param x            codon state at 1 site in child node x
     * @param y            codon state at 1 site in child node y
     * @param wxBranchLd   branch likelihood calculation core at branch w->x
     * @param wyBranchLd   branch likelihood calculation core at branch w->y
     * @param cpd_w        cumulative probability distribution of states at w node,
     *                     to avoid using <code>new double[]</code>.
     * @param frequencies  equilibrium frequencies
     * @param proportions  the array of proportions of sites in different categories.
     * @return the proposed state sampled prob from 3 branch likelihoods,
     *         or negative if something is wrong.
     */
    protected int gibbsSampling(final int x, final int y, final DABranchLikelihoodCore wxBranchLd,
                                final DABranchLikelihoodCore wyBranchLd, double[] cpd_w,
                                final double[] frequencies, final double[] proportions) {

        double pwx,pwy,pr_w = 0;
        // no z
        for (int w=0; w < cpd_w.length; w++) {
            // n = w + i * state + j
            pwx = wxBranchLd.calculateBranchLdAtSite(w, x, proportions);
            pwy = wyBranchLd.calculateBranchLdAtSite(w, y, proportions);
            // w_i ~ P_{w_i}(t) * P_{w_i}x(t) * P_{w_i}y(t)
            pr_w = frequencies[w] * pwx * pwy;
            // cumulate pr_w
            cumulatePr(cpd_w, pr_w, w);
        } // end w loop

        // choose final state w from the distribution
        double random = Randomizer.nextDouble() * cpd_w[cpd_w.length-1];

//        return RandomUtils.binarySearchSampling(cpd_w, random);
        // linear faster than binary search implementation nrOfState < 100
        return RandomUtils.randomChoice(cpd_w, random);
    }

    /**
     * Gibbs sampling state at 1 site 1 node. The node is not the root.
     * @param x            codon state at 1 site in child node x
     * @param y            codon state at 1 site in child node y
     * @param z            codon state at 1 site in parent node x
     * @param wxBranchLd   branch likelihood calculation core at branch w->x
     * @param wyBranchLd   branch likelihood calculation core at branch w->y
     * @param zwBranchLd   branch likelihood calculation core at branch z->w
     * @param cpd_w        cumulative probability distribution of states at w,
     *                     to avoid using <code>new double[]</code>.
     * @param proportions  the array of proportions of sites in different categories.
     * @return the proposed state sampled prob from 3 branch likelihoods,
     *         or negative if something is wrong.
     */
    protected int gibbsSampling( final int x, final int y, final int z,
                                 final DABranchLikelihoodCore wxBranchLd,
                                 final DABranchLikelihoodCore wyBranchLd,
                                 final DABranchLikelihoodCore zwBranchLd,
                                 double[] cpd_w, final double[] proportions) {
        double pzw,pwx,pwy,pr_w = 0;
        for (int w=0; w < cpd_w.length; w++) {
            pzw = zwBranchLd.calculateBranchLdAtSite(z, w, proportions);
            pwx = wxBranchLd.calculateBranchLdAtSite(w, x, proportions);
            pwy = wyBranchLd.calculateBranchLdAtSite(w, y, proportions);
            // w_i ~ P_z{w_i}(t) * P_{w_i}x(t) * P_{w_i}y(t)
            pr_w = pzw * pwx * pwy;
            // cumulate pr_w
            cumulatePr(cpd_w, pr_w, w);
        } // end w loop

        // choose final state w from the distribution
        double random = Randomizer.nextDouble() * cpd_w[cpd_w.length-1];

        return RandomUtils.binarySearchSampling(cpd_w, random);
    }

    // cumulate pr in cpd[], w is the index of cpd[]
    private void cumulatePr(double[] cpd, double pr, int w) {
        if (w == 0)
            cpd[0] = pr;
        else
            cpd[w] = cpd[w - 1] + pr;
    }

    /*** getter ***/
    public int getStartInclusive() {
        return startInclusive;
    }

    public int getEndExclusive() {
        return endExclusive;
    }
}
