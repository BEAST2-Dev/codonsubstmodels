package beast.evolution.tree;

import beast.util.RandomUtils;

import java.util.*;

/**
 * Reconstructing ancestral states using parsimony Fitch (1971) algorithm.
 * Equally to choose a state from the ambiguous set.
 *
 * @author Walter Xie
 */
public class RASParsimony1Site {

    // int[nodeNr]
    final protected int[] tipStates;
    final protected TreeInterface tree;
    // <nodeNr, Set<states>>
    private Map<Integer, Set<Integer>> ancestralStatesMap;

    /**
     * Reconstructing ancestral states using parsimony
     * Fitch (1971) algorithm, given tip states and tree.
     * Use {@link Node#getNr()} as the index of rows in the array.
     * @param tipStates   int[] tip states, [] is nodeNr
     * @param tree        {@link TreeInterface}
     */
    public RASParsimony1Site(int[] tipStates, TreeInterface tree) {
        this.tipStates = tipStates;
        this.tree = tree;
        assert tree.getLeafNodeCount() == getTipsCount();

        ancestralStatesMap = new HashMap<>();
    }

    /**
     * used to check if the tip state has been assigned to the correct tip.
     */
    public void printTipStates() {
        System.out.println("====== " + tipStates.length + " Tips States (ordered by node index) ======");
        for (int i=0; i < tipStates.length; i++) {
            String taxon = tree.getNode(i).getID();
            System.out.println(taxon + " (" + i + ")\t" + tipStates[i]);
        }
        System.out.println();
    }

    /**
     * the method to run Fitch (1971) algorithm.
     * @return   array of ancestral states, index = nodeNr - tipsCount.
     */
    public int[] reconstructAncestralStates() {
        //traverse 1: down pass optimization towards the root
        RASParsimony1Site downpass = new RASParsimony1Site(tipStates, tree);
        int score = 0;
        downpass.traverseTowardsRoot(tree.getRoot(), score);

        //traverse 2: up pass optimization away from the root
        RASParsimony1Site uppass = new RASParsimony1Site(tipStates, tree);
        uppass.traverseAwayRoot(tree.getRoot(), downpass.getAncestralStates());

        //3. final optimization, choose the state having greatest number
        int[] aS = new int[getInternalNodeCount()];
        traverseMPR(tree.getRoot(), downpass.getAncestralStates(), uppass.getAncestralStates(), aS);

        // 4. validate AS change == parsimony score
        int change = 0;
        countASChanges(tree.getRoot(), aS, change);
        if (change != score)
            System.err.println("Parsimony changes" + change + "in ancestral states should = " + score + " ! ");
        return aS;
    }

    /**
     * Recursively count the state change in the tree.
     * It should equal to the parsimony score.
     * @param node             start from root.
     * @param ancestralStates  array of ancestral states, index = nodeNr - tipsCount.
     * @param change           plus 1 if the parent node state != its child node state.
     */
    public void countASChanges(Node node, final int[] ancestralStates, int change) {

        int nr = node.getNr();
        assert nr >= getTipsCount();
        int idx = nr - getTipsCount();
        final int state = ancestralStates[idx];

        // Traverse down the two child nodes
        final Node child1 = node.getChild(0);
        if (child1.isLeaf()) {
            nr = child1.getNr();
            assert nr < getTipsCount();
            int childS = tipStates[nr];
            if (state != childS)
                change++;
        } else {
            countASChanges(child1, ancestralStates, change);
        }

        final Node child2 = node.getChild(1);
        if (child2.isLeaf()) {
            nr = child2.getNr();
            assert nr < getTipsCount();
            int childS = tipStates[nr];
            if (state != childS)
                change++;
        } else {
            countASChanges(child2, ancestralStates, change);
        }

    }

    /**
     * 1. Downpass optimization towards the root.
     * Add parsimony score in assignASAwayRoot().
     * @param node     start from root.
     * @param score    parsimony score
     */
    public void traverseTowardsRoot(Node node, int score) {
        // Traverse down the two child nodes
        final Node child1 = node.getChild(0);
        if (!child1.isLeaf()) traverseTowardsRoot(child1, score);

        final Node child2 = node.getChild(1);
        if (!child2.isLeaf()) traverseTowardsRoot(child2, score);

        assignASTowardsRoot(node, child1, child2, score);
    }

    /**
     * 2. uppass optimization away from the root.
     * @param node          start from root.
     * @param downpassMap   ancestral states from 1st step downpass
     */
    public void traverseAwayRoot(Node node, final Map<Integer, Set<Integer>> downpassMap) {
        final Node child1 = node.getChild(0);
        final Node child2 = node.getChild(1);

        if (node.isRoot()) {
            int nr = node.getNr();
            // use the root state assigned in down pass to start
            ancestralStatesMap.put(nr, downpassMap.get(nr));
        } else {
            final Node parent = node.getParent();
            assignASAwayRoot(node, parent, child1, child2, downpassMap);
        }

        if (!child1.isLeaf()) traverseAwayRoot(child1, downpassMap);
        if (!child2.isLeaf()) traverseAwayRoot(child2, downpassMap);
    }

    /**
     * 3. select final state using MPR: Most Parsimonious Reconstruction.
     * Randomly select 1 state if there are multiple transitions in MPR.
     * Fill in ancestralStates[], index is node.getNr() - getTipsCount().
     * @param node             start from root.
     * @param downpassMap      ancestral states from 1st step downpass
     * @param uppassMap        ancestral states from 2nd step uppass
     * @param ancestralStates  final ancestral states, index = nodeNr - tipsCount.
     */
    public void traverseMPR(Node node, final Map<Integer, Set<Integer>> downpassMap,
                            final Map<Integer, Set<Integer>> uppassMap, int[] ancestralStates) {
        // start from root
        int nr = node.getNr();
        assert nr >= getTipsCount();
        int idx = nr - getTipsCount();

        if (node.isRoot()) {
            assert downpassMap.size() == getInternalNodeCount() && uppassMap.size() == getInternalNodeCount();
            assert nr == getNodeCount() - 1;
            // downpass and uppass is same at root
            Set<Integer> states = downpassMap.get(nr);
            ancestralStates[idx] = RandomUtils.randomSelect(states);

        } else {
            final Node parent = node.getParent();
            int idxP = parent.getNr() - getTipsCount();
            int parentState = ancestralStates[idxP];
            if (parentState < 0)
                throw new IllegalArgumentException("parent state < 0 : " + parentState);

            Set<Integer> downpass = downpassMap.get(nr);
            Set<Integer> uppass = uppassMap.get(nr);
            ancestralStates[idx] = finalState(parentState, downpass, uppass);
        }

        final Node child1 = node.getChild(0);
        final Node child2 = node.getChild(1);
        if (!child1.isLeaf()) traverseMPR(child1, downpassMap, uppassMap, ancestralStates);
        if (!child2.isLeaf()) traverseMPR(child2, downpassMap, uppassMap, ancestralStates);
    }


    // retrieve node state per site,
    // if tip, take tipStates[nr].
    // If internal node, downpassMap.get(nr), which is Map<nodeNr, Set<states>>
    // to store ancestral states of all internal nodes at one site.
    // downpassMap is used by uppass, and is null during downpass.
    protected Set<Integer> getStatesFrom(Node node, final Map<Integer, Set<Integer>> downpassMap) {
        Set<Integer> stateSet = null;

        int nr = node.getNr();
        if (node.isLeaf()) { // tips
            assert nr < getTipsCount();
            final int state = tipStates[nr];
            stateSet = new HashSet<>();
            stateSet.add(state);
        } else {
            assert nr >= getTipsCount() && nr < getNodeCount();
            if (downpassMap != null) { // internal nodes
                stateSet = downpassMap.get(nr);
            } else {
                stateSet = ancestralStatesMap.get(nr);
            }
            // validate 1. get downpass states during uppass;
            // 2. no access internal node states during downpass
            if (stateSet == null)
                throw new IllegalArgumentException("Cannot find optimized states of internal node " + nr);
        }
        return stateSet;
    }

    // Input child1, child2, assign ancestral states to node in Map<Integer, Set<Integer>>.
    // Ancestral state is the intersection state of child1 and child2, otherwise take the union.
    // anceStates1SiteMap for node should be empty before this process.
    // if there is intersection, assign intersection, otherwise use union.
    protected void assignASTowardsRoot(Node node, Node child1, Node child2, int score) {
        Set<Integer> d1 = getStatesFrom(child1, null);
        final Set<Integer> d2 = getStatesFrom(child2, null);

        final int nr = node.getNr();
        if (ancestralStatesMap.containsKey(nr))
            throw new IllegalArgumentException("States of internal node " + nr + " cannot exist !");

        Set<Integer> intxn = new HashSet<>(d1);
        if (intxn.retainAll(d2)) // intersection
            ancestralStatesMap.put(nr, intxn);
        else { // union
            d1.addAll(d2);
            score++;
        }
        ancestralStatesMap.put(nr, d1);
    }

    // Input child1, child2, assign ancestral states to node in Map<Integer, Set<Integer>>.
    // Ancestral state is the intersection state of child1 and child2, otherwise take the union.
    // anceStates1SiteMap for node should be empty before this process.
    protected void assignASAwayRoot(Node node, Node parent, Node child1, Node child2,
                                    final Map<Integer, Set<Integer>> downpassMap) {
        final Set<Integer> upP = getStatesFrom(parent, null);
        final Set<Integer> downN = getStatesFrom(node, downpassMap);
        final Set<Integer> downD1 = getStatesFrom(child1, downpassMap);
        final Set<Integer> downD2 = getStatesFrom(child2, downpassMap);

        final int nr = node.getNr();
        if (ancestralStatesMap.containsKey(nr))
            throw new IllegalArgumentException("States of internal node " + nr + " cannot exist !");

        Set<Integer> tmpSet = new HashSet<>(downN);
        tmpSet.retainAll(upP);
        // 1. if d_n ∩ u_p = u_p
        if (tmpSet.containsAll(upP) && upP.containsAll(tmpSet)) {
            // then u_n = d_n ∩ u_p
            ancestralStatesMap.put(nr, tmpSet);
        } else {
            tmpSet = new HashSet<>(downD1);
            tmpSet.retainAll(downD2);
            // 2. if d_L ∩ d_R = ∅
            if (tmpSet.isEmpty()) {
                tmpSet = new HashSet<>(downN);
                tmpSet.addAll(upP); // then u_n = d_n ∪ u_p
                ancestralStatesMap.put(nr, tmpSet);
            } else {
                // 3. else u_n = d_n ∪ (u_p ∩ (d_L ∪ d_R))
                tmpSet = new HashSet<>(downD1);
                tmpSet.addAll(downD2); // d_L ∪ d_R
                tmpSet.retainAll(upP); // u_p ∩
                tmpSet.addAll(downN); // d_n ∪
                ancestralStatesMap.put(nr, tmpSet);
            }
        }
    }

    // select final state in MPR (Most Parsimonious Reconstruction) per branch given a state in parent node.
    // Set<Integer> downpass, Set<Integer> uppass for the child node of the branch.
    protected int finalState(int parentState, final Set<Integer> downpass, final Set<Integer> uppass) {
        // if parentState in the up-pass set of parent is in the down-pass set of child
        if (downpass.contains(parentState)) {
            return parentState;
        } else {
            // add all states in the down-pass set of child
            Set<Integer> tmpSet = new HashSet<>(downpass);
            // add parentState if parentState is in the up-pass set of child
            if (uppass.contains(parentState))
                tmpSet.add(parentState);

            int rand = RandomUtils.randomSelect(tmpSet);
            return rand;
        } // end if else
    }

    public String toString() {
        ArrayList<String>[] lineages = (ArrayList<String>[]) new ArrayList[getTipsCount()];
        nodeToState(tree.getRoot(), lineages, 0);

        StringBuilder str = new StringBuilder();
        for (int i = 0; i < lineages.length; i++) {
            ArrayList<String> lin = lineages[i];
            str.append(String.join("\t", lin));
            str.append("\n");
        }
        return str.toString();
    }

    protected void nodeToState(Node node, ArrayList<String>[] lineages, int row) {
        final int nr = node.getNr();
        if (node.isLeaf()) {
            lineages[row].add(nr + ":{" + tipStates[nr] + "}");
            row++;
        } else {
            String set = nr + ":{";
            for(Integer state : ancestralStatesMap.get(nr))
                set += state.toString() + ",";
            set += "}";
            lineages[row].add(set);
            int nrOfParents = lineages[row].size();

            // full lineage always in child 1
            final Node child1 = node.getChild(0);
            nodeToState(child1, lineages, row);

            // row will ++ when reach to leaf from child 1 traverse
            // add spaces for child 2
            for (int i = 0; i < nrOfParents; i++)
                lineages[row].add(set);
            final Node child2 = node.getChild(0);
            nodeToState(child2, lineages, row);
        }
    }


//+++++++ getters ++++++++

    public Map<Integer, Set<Integer>> getAncestralStates() {
        return ancestralStatesMap;
    }

    public int getTipsCount() {
        return tipStates.length;
    }

    /**
     * {@link #getTipsCount()} - 1
     */
    public int getInternalNodeCount() {
        return getTipsCount() - 1;
    }


    public int getNodeCount() {
        return 2*getTipsCount() - 1;
    }


//    private int[] sampleASUniform(Map<Integer, Set<Integer>> downpassMap,
//                                  Map<Integer, Set<Integer>> uppassMap, long seed) {
//        int[] anceStates1Site = new int[getInternalNodeCount()];
//
//        // internal nodes
//        for (int i = getTipsCount(); i < getNodeCount(); i++) {
//            Set<Integer> stateSetDown = downpassMap.get(i);
//            Set<Integer> stateSetUp = uppassMap.get(i);
//
//            if (stateSetDown == null)
//                throw new IllegalArgumentException("Cannot find optimized states of internal node " + i + " during downpass !");
//            // uppass does not assign ancestral states to the root
//            if (stateSetUp == null && i < getNodeCount()-1)
//                throw new IllegalArgumentException("Cannot find optimized states of internal node " + i + " during uppass !");
//
//            // use intersection and union to replace frequency count and maximization
//            // if intersection exists, it has that states having the greatest frequency,
//            // otherwise every states will be same frequent.
//            Set<Integer> candidateStates = new HashSet<>(stateSetDown);
//            if (stateSetUp != null) {
//                if (!candidateStates.retainAll(stateSetUp)) { // intersection
//                    candidateStates = new HashSet<>();
//                    // Union
//                    candidateStates.addAll(stateSetDown);
//                    candidateStates.addAll(stateSetUp);
//                }
//            }
//
//            // randomly choose the states having greatest number with equal probability
//            int finalState = chooseStateRandom(candidateStates, seed);
//            if (finalState < 0)
//                throw new IllegalArgumentException("Cannot select final states of internal node " + i +
//                        " from " + Arrays.toString(candidateStates.toArray()));
//
//            // final
//            anceStates1Site[i-getTipsCount()] = finalState;
//        }
//        return anceStates1Site;
//    }
//

}
