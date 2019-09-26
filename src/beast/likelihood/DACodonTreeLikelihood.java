//package beast.likelihood;
//
//import beast.core.Input;
//import beast.core.State;
//import beast.core.util.Log;
//import beast.evolution.alignment.Alignment;
//import beast.evolution.alignment.CodonAlignment;
//import beast.evolution.branchratemodel.BranchRateModel;
//import beast.evolution.branchratemodel.StrictClockModel;
//import beast.evolution.likelihood.*;
//import beast.evolution.sitemodel.SiteModel;
//import beast.evolution.substitutionmodel.SubstitutionModel;
//import beast.evolution.tree.Node;
//import beast.evolution.tree.Tree;
//import beast.evolution.tree.TreeInterface;
//import beast.tree.InternalNodeSeqs;
//
//import java.io.PrintStream;
//import java.util.ArrayList;
//import java.util.Arrays;
//import java.util.List;
//import java.util.Random;
//
///**
// * Data augmentation to fast codon tree likelihood calculation.
// *
// *
// *
// * @author Walter Xie, Fabio Mendes
// */
//public class DACodonTreeLikelihood extends GenericTreeLikelihood {
//
//
//    final public Input<InternalNodeSeqs> internalNodeSeqsInput = new Input<>("internalNodeSeqs",
//            "The large 2-d matrix to store internal node sequences.", Input.Validate.REQUIRED);
//
//    InternalNodeSeqs internalNodeSeqs;
//
//    /**
//     * BEASTObject associated with inputs. Since none of the inputs are StateNodes, it
//     * is safe to link to them only once, during initAndValidate.
//     */
//    protected SubstitutionModel substitutionModel;
//    protected SiteModel.Base m_siteModel;
//    protected BranchRateModel.Base branchRateModel;
//
//    /**
//     * flag to indicate the
//     * // when CLEAN=0, nothing needs to be recalculated for the node
//     * // when DIRTY=1 indicates a node partial needs to be recalculated
//     * // when FILTHY=2 indicates the indices for the node need to be recalculated
//     * // (often not necessary while node partial recalculation is required)
//     */
//    protected int hasDirt;
//
//    /**
//     * Lengths of the branches in the tree associated with each of the nodes
//     * in the tree through their node  numbers. By comparing whether the
//     * current branch length differs from stored branch lengths, it is tested
//     * whether a node is dirty and needs to be recomputed (there may be other
//     * reasons as well...).
//     * These lengths take branch rate models in account.
//     */
//    protected double[] m_branchLengths;
//    protected double[] storedBranchLengths;
//
//
//    /**
//     * memory allocation for probability tables obtained from the SiteModel *
//     */
//    protected double[] probabilities;
//
//    protected int matrixSize;
//
//
//
//
//    @Override
//    public void initAndValidate() {
//
//        internalNodeSeqs = internalNodeSeqsInput.get();
//
//        CodonAlignment codonAlignment = CodonAlignment.toCodonAlignment(dataInput.get());
//
//        // sanity check: alignment should have same #taxa as tree
//        if (codonAlignment.getTaxonCount() != treeInput.get().getLeafNodeCount()) {
//            throw new IllegalArgumentException("The number of nodes in the tree does not match the number of sequences");
//        }
////        beagle = null;
////        beagle = new BeagleTreeLikelihood();
////        try {
////            beagle.initByName(
////                    "data", codonAlignment, "tree", treeInput.get(), "siteModel", siteModelInput.get(),
////                    "branchRateModel", branchRateModelInput.get(), "useAmbiguities", m_useAmbiguities.get(),
////                    "useTipLikelihoods", m_useTipLikelihoods.get(),"scaling", scaling.get().toString());
////            if (beagle.beagle != null) {
////                //a Beagle instance was found, so we use it
////                return;
////            }
////        } catch (Exception e) {
////            // ignore
////        }
////        // No Beagle instance was found, so we use the good old java likelihood core
////        beagle = null;
//
//        int nodeCount = treeInput.get().getNodeCount();
//        if (!(siteModelInput.get() instanceof SiteModel.Base)) {
//            throw new IllegalArgumentException("siteModel input should be of type SiteModel.Base");
//        }
//        m_siteModel = (SiteModel.Base) siteModelInput.get();
//        m_siteModel.setDataType(codonAlignment.getDataType());
//        substitutionModel = m_siteModel.substModelInput.get();
//
//        if (branchRateModelInput.get() != null) {
//            branchRateModel = branchRateModelInput.get();
//        } else {
//            branchRateModel = new StrictClockModel();
//        }
//        m_branchLengths = new double[nodeCount];
//        storedBranchLengths = new double[nodeCount];
//
//        int stateCount = codonAlignment.getMaxStateCount();
//
//
//        Log.info.println("  " + codonAlignment.toString(true));
//        // print startup messages via Log.print*
//
//
//        hasDirt = Tree.IS_FILTHY;
//
//
//        matrixSize = (stateCount + 1) * (stateCount + 1);
//        probabilities = new double[(stateCount + 1) * (stateCount + 1)];
//        Arrays.fill(probabilities, 1.0);
//
//    }
//
//
//    @Override
//    public double calculateLogP() {
//        //TODO beagle
//        logP = 0;
//
//        final TreeInterface tree = treeInput.get();
//
//        try {
//            if (traverse(tree.getRoot()) != Tree.IS_CLEAN)
//                for (int i = 0; i < codonAlignment.getPatternCount(); i++) {
//                    logP += patternLogLikelihoods[i] * codonAlignment.getPatternWeight(i);
//                }
//        }
//        catch (ArithmeticException e) {
//            return Double.NEGATIVE_INFINITY;
//        }
//
//        m_nScale++;
//        if (logP > 0 || (likelihoodCore.getUseScaling() && m_nScale > X)) {
////            System.err.println("Switch off scaling");
////            m_likelihoodCore.setUseScaling(1.0);
////            m_likelihoodCore.unstore();
////            m_nHasDirt = Tree.IS_FILTHY;
////            X *= 2;
////            traverse(tree.getRoot());
////            calcLogP();
////            return logP;
//        } else if (logP == Double.NEGATIVE_INFINITY && m_fScale < 10 && !scaling.get().equals(TreeLikelihood.Scaling.none)) { // && !m_likelihoodCore.getUseScaling()) {
//            m_nScale = 0;
//            m_fScale *= 1.01;
//            Log.warning.println("Turning on scaling to prevent numeric instability " + m_fScale);
//            likelihoodCore.setUseScaling(m_fScale);
//            likelihoodCore.unstore();
//            hasDirt = Tree.IS_FILTHY;
//            traverse(tree.getRoot());
//            calcLogP();
//            return logP;
//        }
//
//        return logP;
//    }
//
//
//    /* Assumes there IS a branch rate model as opposed to traverse() */
//    int traverse(final Node node) {
//
//        int update = (node.isDirty() | hasDirt);
//
//        final int nodeIndex = node.getNr();
//
//        final double branchRate = branchRateModel.getRateForBranch(node);
//        final double branchTime = node.getLength() * branchRate;
//
//        // First update the transition probability matrix(ices) for this branch
//        //if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_StoredBranchLengths[nodeIndex])) {
//        if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_branchLengths[nodeIndex])) {
//            m_branchLengths[nodeIndex] = branchTime;
//            final Node parent = node.getParent();
//            likelihoodCore.setNodeMatrixForUpdate(nodeIndex);
//            for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
//                final double jointBranchRate = m_siteModel.getRateForCategory(i, node) * branchRate;
//                substitutionModel.getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), jointBranchRate, probabilities);
//                //System.out.println(node.getNr() + " " + Arrays.toString(m_fProbabilities));
//                likelihoodCore.setNodeMatrix(nodeIndex, i, probabilities);
//            }
//            update |= Tree.IS_DIRTY;
//        }
//
//        // If the node is internal, update the partial likelihoods.
//        if (!node.isLeaf()) {
//
//            // Traverse down the two child nodes
//            final Node child1 = node.getLeft(); //Two children
//            final int update1 = traverse(child1);
//
//            final Node child2 = node.getRight();
//            final int update2 = traverse(child2);
//
//            // If either child node was updated then update this node too
//            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {
//
//                final int childNum1 = child1.getNr();
//                final int childNum2 = child2.getNr();
//
//                likelihoodCore.setNodePartialsForUpdate(nodeIndex);
//                update |= (update1 | update2);
//                if (update >= Tree.IS_FILTHY) {
//                    likelihoodCore.setNodeStatesForUpdate(nodeIndex);
//                }
//
//                if (m_siteModel.integrateAcrossCategories()) {
//                    likelihoodCore.calculatePartials(childNum1, childNum2, nodeIndex);
//                } else {
//                    throw new RuntimeException("Error TreeLikelihood 201: Site categories not supported");
//                    //m_pLikelihoodCore->calculatePartials(childNum1, childNum2, nodeNum, siteCategories);
//                }
//
//                if (node.isRoot()) {
//                    // No parent this is the root of the beast.tree -
//                    // calculate the pattern likelihoods
//                    final double[] frequencies = //m_pFreqs.get().
//                            substitutionModel.getFrequencies();
//
//                    final double[] proportions = m_siteModel.getCategoryProportions(node);
//                    likelihoodCore.integratePartials(node.getNr(), proportions, m_fRootPartials);
//
//                    if (constantPattern != null) { // && !SiteModel.g_bUseOriginal) {
//                        proportionInvariant = m_siteModel.getProportionInvariant();
//                        // some portion of sites is invariant, so adjust root partials for this
//                        for (final int i : constantPattern) {
//                            m_fRootPartials[i] += proportionInvariant;
//                        }
//                    }
//
//                    likelihoodCore.calculateLogLikelihoods(m_fRootPartials, frequencies, patternLogLikelihoods);
//                }
//
//            }
//        }
//        return update;
//    } // traverseWithBRM
//
//
//
//    @Override
//    public void sample(State state, Random random) {
//        super.sample(state, random);
//    }
//
//    @Override
//    protected boolean requiresRecalculation() {
//        return super.requiresRecalculation();
//    }
//
//    @Override
//    protected void accept() {
//        super.accept();
//    }
//
//    @Override
//    public void store() {
//        super.store();
//    }
//
//    @Override
//    public void restore() {
//        super.restore();
//    }
//
//    @Override
//    public void log(long sample, PrintStream out) {
//        super.log(sample, out);
//    }
//}
