package beast.evolution.tree.parsimony;

import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;

import java.util.*;

/**
 * Reconstructing ancestral states using Wagner parsimony.
 * Cunningham et. al. 1998, Swofford & Maddison 1987, Kluge & Farris 1969.
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
     * Reconstructing ancestral states using parsimony given tip states and tree.
     * Use {@link Node#getNr()} as the index of rows in the array.
     * @param tipStates   int[] tip states, [] is nodeNr
     * @param tree        {@link TreeInterface}
     */
    public RASParsimony1Site(int[] tipStates, TreeInterface tree) {
        this.tipStates = tipStates;
        this.tree = tree;
        assert tree.getLeafNodeCount() == getTipsCount();
    }

    public int[] reconstructAncestralStates(int[] tipStates, TreeInterface tree) {
        int[] ancestralStates = new int[getInternalNodeCount()];

            //traverse 1: down pass optimization towards the root
        RASParsimony1Site downpass = new RASParsimony1Site(tipStates, tree);
        int score = 0;
        downpass.traverseTowardsRoot(tree.getRoot(), score);

            //traverse 2: up pass optimization away from the root
        RASParsimony1Site uppass = new RASParsimony1Site(tipStates, tree);
        uppass.traverseAwayRoot(tree.getRoot(), downpass.getAncestralStates());

            //final optimization, choose the state having greatest number


        return ancestralStates;
    }

    // downpass optimization towards the root.
    // input nodes are children nodes of the ancestral node which the optimized state is assigned to.
    // add parsimony score in assignASAwayRoot().
    public void traverseTowardsRoot(Node node, int score) {
        // Traverse down the two child nodes
        final Node child1 = node.getChild(0);
        if (!child1.isLeaf()) traverseTowardsRoot(child1, score);

        final Node child2 = node.getChild(1);
        if (!child2.isLeaf()) traverseTowardsRoot(child2, score);

        assignASTowardsRoot(node, child1, child2, score);
    }

    // uppass optimization away from the root
    // input nodes are the ancestor and its sister node.
    public void traverseAwayRoot(Node node, final Map<Integer, Set<Integer>> downpassMap) {
        final Node child1 = node.getChild(0);
        final Node child2 = node.getChild(1);

        if (node.isRoot()) {
            // use the root state assigned in down pass to start
            int nr = node.getNr();
            ancestralStatesMap.put(nr, downpassMap.get(nr));
        } else {
            final Node parent = node.getParent();
            assignASAwayRoot(node, parent, child1, child2, downpassMap);
        }

        if (!child1.isLeaf()) traverseAwayRoot(child1, downpassMap);
        if (!child2.isLeaf()) traverseAwayRoot(child2, downpassMap);
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
//    /**
//     * Randomly select the final state from a set of states.
//     * @param candidateStates   a set of states to choose
//     * @param seed              for {@link Random}
//     * @return
//     */
//    public int chooseStateRandom(Set<Integer> candidateStates, long seed) {
//        int size = candidateStates.size();
//        int item = new Random(seed).nextInt(size);
//        int s = 0;
//        for(Integer state : candidateStates) {
//            if (s == item)
//                return state;
//            s++;
//        }
//        return -1;
//    }

}
