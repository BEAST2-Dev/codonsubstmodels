package codonmodels.evolution.tree;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.util.Randomizer;
import codonmodels.util.RandomUtils;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * The test case is copied from Mark Holder's slides.
 * @author Walter Xie
 */
public class RASParsimonyTest {
    final TreeInterface tree = new Tree("(((A:1,D:1):1,(G:1,C:1):1):1,(B:1,(E:1,((H:1,I:1):1,F:1):1):1):1);");
    // A->0, C->1, G->2, ordered by tips: A, B, ..., I
    final int[] tipStates = new int[]{ 0, 1, 1, 0, 0, 1, 2, 0, 0 };

    RASParsimony1Site downpass;
    RASParsimony1Site uppass;

    @BeforeEach
    public void setUp() throws Exception {
        downpass = new RASParsimony1Site(tipStates, tree);
        uppass = new RASParsimony1Site(tipStates, tree);
    }

    @Test
    public void getters() {
        downpass.printTipStates();
        System.out.println(tree);

        assertEquals(9, downpass.getTipsCount(), "Tips");
        assertEquals(8, downpass.getInternalNodeCount(), "Internal nodes");
        assertEquals(17, downpass.getNodeCount(), "All nodes");
    }

    @Test
    public void traverseTowardsRoot() {
        // clear map before test
        downpass.clearAncestralStates();
        //traverse 1: down pass optimization towards the root
        int score = downpass.traverseTowardsRoot(tree.getRoot());
        System.out.println("Down-pass : \n" + downpass.toString());

        Map<Integer, Set<Integer>> asMap = downpass.getAncestralStates();
        assertEquals(new HashSet<>(Arrays.asList(0, 1)), asMap.get(16), "Root (16)");
        assertEquals(new HashSet<>(Arrays.asList(0, 1, 2)), asMap.get(11), "Node (11)");
        assertEquals(new HashSet<>(Arrays.asList(1, 2)), asMap.get(10), "Node (10)");
        assertEquals(new HashSet<>(Arrays.asList(0)), asMap.get(9), "Node (9)");

        System.out.println("Parsimony score = " + score);
        assertEquals(4, score, "Parsimony score");
    }


    @Test
    public void traverseAwayRoot() {
        downpass.clearAncestralStates();
        //traverse 1: down pass optimization towards the root
        int score = downpass.traverseTowardsRoot(tree.getRoot());

        uppass.clearAncestralStates();
        //traverse 2: up pass optimization away from the root
        uppass.traverseAwayRoot(tree.getRoot(), downpass.getAncestralStates());

        System.out.println("Up-pass : \n" + uppass.toString());
        Map<Integer, Set<Integer>> asMap = uppass.getAncestralStates();
        assertEquals(new HashSet<>(Arrays.asList(0, 1)), asMap.get(16), "Root (16)");
        assertEquals(new HashSet<>(Arrays.asList(0, 1)), asMap.get(11), "Node (11)");
        assertEquals(new HashSet<>(Arrays.asList(0, 1, 2)), asMap.get(10), "Node (10)");
        assertEquals(new HashSet<>(Arrays.asList(0)), asMap.get(9), "Node (9)");
    }

    @Test
    public void traverseMPR() {
        downpass.clearAncestralStates();
        //traverse 1: down pass optimization towards the root
        int score = downpass.traverseTowardsRoot(tree.getRoot());

        uppass.clearAncestralStates();
        //traverse 2: up pass optimization away from the root
        uppass.traverseAwayRoot(tree.getRoot(), downpass.getAncestralStates());

        //3. final optimization, choose the state having greatest number
        RASParsimony1Site mpr = new RASParsimony1Site(tipStates, tree);
        int[] aS = new int[mpr.getInternalNodeCount()];
        Arrays.fill(aS, -1);
        Randomizer.setSeed(1);
        mpr.traverseMPR(tree.getRoot(), downpass.getAncestralStates(), uppass.getAncestralStates(), aS);
        // choose 1 at root
        System.out.println("Ancestral states : " + Arrays.toString(aS) );
        assertArrayEquals(new int[]{0, 1, 1, 0, 1, 1, 1, 1}, aS, "Ancestral states");

        // 4. validate AS change == parsimony score
        int change = mpr.countASChanges(tree.getRoot(), aS);
        System.out.println("Parsimony score = " + score + ", state change = " + change );
        assertEquals(score, change, "State change");

        mpr.clearAncestralStates();
        aS = new int[mpr.getInternalNodeCount()];
        Arrays.fill(aS, -1);
        Randomizer.setSeed(2);
        mpr.traverseMPR(tree.getRoot(), downpass.getAncestralStates(), uppass.getAncestralStates(), aS);
        // choose 0 at root
        System.out.println("Ancestral states : " + Arrays.toString(aS) );
        assertArrayEquals(new int[]{0, 0, 0, 0, 0, 0, 0, 0}, aS, "Ancestral states");

        change = mpr.countASChanges(tree.getRoot(), aS);
        System.out.println("Parsimony score = " + score + ", state change = " + change );
        assertEquals(score, change, "State change");
    }

    @Test
    public void testIntxnUnion() {
        Set<Integer> d1 = new HashSet<>(Arrays.asList(0, 1));
        Set<Integer> d2 = new HashSet<>(Arrays.asList(1, 2));

        System.out.println("Set 1 : " + Arrays.toString(d1.toArray()) );
        System.out.println("Set 2 : " + Arrays.toString(d2.toArray()) );
        Set<Integer> intxn = new HashSet<>(d1);
        boolean changed = intxn.retainAll(d2);
        System.out.println("Intersection : " + " changed = " + changed + " " + Arrays.toString(intxn.toArray()) );

        Set<Integer> union = new HashSet<>(d1);
        changed = union.addAll(d2);
        System.out.println("Union : " + " changed = " + changed + " " + Arrays.toString(union.toArray()) );

        d1 = new HashSet<>(Arrays.asList(0, 3));
        System.out.println("Set 1 : " + Arrays.toString(d1.toArray()) );
        System.out.println("Set 2 : " + Arrays.toString(d2.toArray()) );
        intxn = new HashSet<>(d1);
        changed = intxn.retainAll(d2);
        System.out.println("Intersection : " + " changed = " + changed + " " + Arrays.toString(intxn.toArray()) );

        union = new HashSet<>(d1);
        changed = union.addAll(d2);
        System.out.println("Union : " + " changed = " + changed + " " + Arrays.toString(union.toArray()) );
    }

    @Test
    public void testRandomUtils() {
        Set<Integer> tmpSet = new HashSet<>(Arrays.asList(0));
        System.out.println("Set 1 : " + Arrays.toString(tmpSet.toArray()) );
        int rand = RandomUtils.randomSelect(tmpSet);
        System.out.println("Random state = " + rand);

        assertEquals(0, rand, "Random state");

        tmpSet = new HashSet<>(Arrays.asList(0, 1));
        System.out.println("Set 2 : " + Arrays.toString(tmpSet.toArray()) );
        Randomizer.setSeed(1);
        rand = RandomUtils.randomSelect(tmpSet);
        System.out.println("Random state = " + rand);
        assertEquals(1, rand, "Random state");

        Randomizer.setSeed(2);
        rand = RandomUtils.randomSelect(tmpSet);
        System.out.println("Random state = " + rand);
        assertEquals(0, rand, "Random state");
    }

}