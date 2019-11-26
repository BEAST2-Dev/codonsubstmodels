package beast.evolution.tree;

import beast.util.Randomizer;

import java.util.*;

/**
 * Reconstructing ancestral states using parsimony.
 * Snell & Childress 1987.
 * Cunningham et. al. 1998.
 * Equally to choose a state from the ambiguous set.
 *
 * @author Walter Xie
 */
public class RASParsimony {

    // int[nodeNr][site]
    final protected int[][] tipStates;
    protected TreeInterface tree;



    /**
     * Reconstructing ancestral states using parsimony given tip states and tree.
     * Use {@link Node#getNr()} as the index of rows in 2d array.
     * @param tipStates   int[][] tip states, 1st[] is nodeNr, 2nd[] is site.
     * @param tree        {@link TreeInterface}
     */
    public RASParsimony(int[][] tipStates, TreeInterface tree) {
        this.tipStates = tipStates;
        this.tree = tree;
        assert tree.getLeafNodeCount() == getTipsCount();
    }

    /**
     * Parsimony to init states Snell & Childress 1987.
     * Randomly choose a state from the ambiguous set with equal prob.
     * No ambiguous ancestral states.
     * @return int[][], rows are the internal node array index,
     *         which starts 0, and equals to nodeNr - tipsCount.
     *         Columns are sites.
     */
    public int[][] getAncestralStatesNoAmbiguous() {
        int[][] ancestralStates = new int[getInternalNodeCount()][getSiteCount()];
        for (int k = 0; k < getSiteCount(); k++) {
            //traverse 1: down pass optimization towards the root
            Map<Integer, Set<Integer>> downpass1SiteMap = new HashMap<>();
            traverseToRoot(tree.getRoot(), k, downpass1SiteMap);


            //traverse 2: up pass optimization away from the root
            Map<Integer, Set<Integer>> uppass1SiteMap = new HashMap<>();
            traverseToRoot(tree.getRoot(), k, uppass1SiteMap);

            //final optimization, choose the state having greatest number
            long seed = Randomizer.getSeed();
            int[] ancestralStates1Site = sampleASUniform(downpass1SiteMap, uppass1SiteMap, seed);
            for (int i = 0; i < getInternalNodeCount(); i++) {
                ancestralStates[i][k] = ancestralStates1Site[i];
            }
        }
        return ancestralStates;
    }

    // retrieve node state per site,
    // if tip, take tipStates[nr][siteNr], if internal node, anceStates1SiteMap.get(nr),
    // which is a Map<Integer, Set<Integer>> to store ancestral states of all internal nodes at one site.
    // Key is internal node Nr, value is the set of ancestral states at this node.
    private Set<Integer> getStatesFrom(Node node, final int siteNr,
                                       Map<Integer, Set<Integer>> anceStates1SiteMap) {
        Set<Integer> stateSet;

        int nr = node.getNr();
        if (node.isLeaf()) {
            assert nr < getTipsCount();
            final int state = tipStates[nr][siteNr];
            stateSet = new HashSet<>();
            stateSet.add(state);
        } else {
            assert nr >= getTipsCount() && nr < getNodeCount();
            stateSet = anceStates1SiteMap.get(nr);
            if (stateSet == null)
                throw new IllegalArgumentException("Cannot find optimized states of internal node " + nr);
        }
        return stateSet;
    }


    // Input node1, node2, assign ancestral states to node3 in Map<Integer, Set<Integer>>.
    // Ancestral state is the intersection state of node1 and node2, otherwise take the union.
    // anceStates1SiteMap for node3 should be empty before this process.
    private void assignAncestralStates(Node node1, Node node2, Node node3, final int siteNr,
                                       Map<Integer, Set<Integer>> anceStates1SiteMap) {
        Set<Integer> as1 = getStatesFrom(node1, siteNr, anceStates1SiteMap);
        Set<Integer> as2 = getStatesFrom(node2, siteNr, anceStates1SiteMap);
        final int nr = node3.getNr();
        if (anceStates1SiteMap.containsKey(nr))
            throw new IllegalArgumentException("States of internal node " + nr + " cannot exist !");

        Set<Integer> intersection = new HashSet<>(as1);
        if (intersection.retainAll(as2)) // intersection
            anceStates1SiteMap.put(nr, intersection);
        else // union
            as1.addAll(as2);
        anceStates1SiteMap.put(nr, as1);
    }

    // downpass optimization towards the root.
    // input nodes are children nodes of the ancestral node which the optimized state is assigned to.
    private void traverseToRoot(Node node, final int siteNr, Map<Integer, Set<Integer>> anceStates1SiteMap) {
        // Traverse down the two child nodes
        final Node child1 = node.getChild(0);
        if (!child1.isLeaf()) traverseToRoot(child1, siteNr, anceStates1SiteMap);

        final Node child2 = node.getChild(1);
        if (!child2.isLeaf()) traverseToRoot(child2, siteNr, anceStates1SiteMap);

        assignAncestralStates(child1, child2, node, siteNr, anceStates1SiteMap);
    }

    // uppass optimization away from the root
    // input nodes are the ancestor and its sister node.
    private void traverseAwayRoot(Node node, final int siteNr, Map<Integer, Set<Integer>> anceStates1SiteMap) {
        assert !node.isLeaf();

        // Traverse down the two child nodes
        final Node child1 = node.getChild(0);
        final Node child2 = node.getChild(1);

        assignAncestralStates(child1, child2, node, siteNr, anceStates1SiteMap);

        if (!child1.isLeaf()) traverseAwayRoot(child1, siteNr, anceStates1SiteMap);
        if (!child2.isLeaf()) traverseAwayRoot(child2, siteNr, anceStates1SiteMap);
    }


    private int[] sampleASUniform(Map<Integer, Set<Integer>> downpass1SiteMap,
                                    Map<Integer, Set<Integer>> uppass1SiteMap, long seed) {
        int[] anceStates1Site = new int[getInternalNodeCount()];

        // internal nodes
        for (int i = getTipsCount(); i < getNodeCount(); i++) {
            Set<Integer> stateSetDown = downpass1SiteMap.get(i);
            Set<Integer> stateSetUp = uppass1SiteMap.get(i);

            if (stateSetDown == null)
                throw new IllegalArgumentException("Cannot find optimized states of internal node " + i + " during downpass !");
            // uppass does not assign ancestral states to the root
            if (stateSetUp == null && i < getNodeCount()-1)
                throw new IllegalArgumentException("Cannot find optimized states of internal node " + i + " during uppass !");

            // use intersection and union to replace frequency count and maximization
            // if intersection exists, it has that states having the greatest frequency,
            // otherwise every states will be same frequent.
            Set<Integer> candidateStates = new HashSet<>(stateSetDown);
            if (stateSetUp != null) {
                if (!candidateStates.retainAll(stateSetUp)) { // intersection
                    candidateStates = new HashSet<>();
                    // Union
                    candidateStates.addAll(stateSetDown);
                    candidateStates.addAll(stateSetUp);
                }
            }

            // randomly choose the states having greatest number with equal probability
            int finalState = getRandomFinalState(candidateStates, seed);
            if (finalState < 0)
                throw new IllegalArgumentException("Cannot select final states of internal node " + i +
                        " from " + Arrays.toString(candidateStates.toArray()));

            // final
            anceStates1Site[i-getTipsCount()] = finalState;
        }
        return anceStates1Site;
    }

    public int getRandomFinalState (Set<Integer> candidateStates, long seed) {
        int size = candidateStates.size();
        int item = new Random(seed).nextInt(size);
        int s = 0;
        for(Integer state : candidateStates) {
            if (s == item)
                return state;
            s++;
        }
        return -1;
    }


    protected int getSiteCount() {
        return tipStates[0].length;
    }

    protected int getTipsCount() {
        return tipStates.length;
    }

    /**
     * {@link #getTipsCount()} - 1
     */
    protected int getInternalNodeCount() {
        return getTipsCount() - 1;
    }


    protected int getNodeCount() {
        return 2*getTipsCount() - 1;
    }

}
