package beast.evolution.tree;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

import java.io.PrintStream;

@Description("Logs tree annotated with metadata and/or rates")
public class NodeStatesLogger extends BEASTObject implements Loggable {
    final public Input<NodeStatesArray> nodesStatesInput = new Input<>("nodesStates",
            "Internal node states to be logged", Validate.REQUIRED);
    final public Input<Tree> treeInput = new Input<>("tree",
            "tree to provide node ID while logging.", Validate.REQUIRED);

//    final public Input<Boolean> compressInput = new Input<>("compress",
//            "if require to compress log", false, Validate.OPTIONAL);

    private NodeStatesArray nsa;
    private int tipCount;
    private int nodeCount;

    private boolean isCompress;

    @Override
    public void initAndValidate() {
        nsa = nodesStatesInput.get();
        tipCount = nsa.getTipsCount();
        nodeCount = nsa.getNodeCount();

//        isCompress = compressInput.get();
    }

    @Override
    public void init(PrintStream out) {
        out.println("#Internal Node Sequences\t" + nsa.getInternalNodeCount());
        out.println("#Sites\t" + nsa.getSiteCount());
        out.println("#States\t" + nsa.getStateCount());

        out.println("Sample\tNode\tStates");
        // state from 0 to 59/60

//        final int ncol = nsa.getSiteCount();
//        if (ncol == 1) {
//            out.print(getID() + "\t");
//        } else {
//            for (int i = 0; i < ncol; i++) {
//                out.print("s" + (i + 1) + "\t");
//            }
//        }
    }

    //TODO multithreading https://stackoverflow.com/questions/21632585/thread-safety-of-printstream-in-java
    @Override
    public void log(long sampleNr, PrintStream out) {
        NodeStatesArray nsaCurrent = (NodeStatesArray) nsa.getCurrent();
        Tree treeCurrent = (Tree) treeInput.get().getCurrent();
        out.print(sampleNr + "\t" + 0 + "\t");
        // log branches, parent node .. child, e.g. 3..1,5..3,3..2,5..4
        out.print(getNodesMap(treeCurrent));
        out.print("\n");
        // only internal nodes
        for (int i=tipCount; i < nodeCount; i++) {
            out.print("\t");
            nsaCurrent.getNodeStates(i).log(out); // state from 0 to 59/60
            out.print("\n");
        }
    }

    // parent node .. child, e.g. 3..1,5..3,3..2,5..4
    // Note: BEAST tree log makes nodeNr+1 in Newick
    private String getNodesMap(Tree tree) {
        StringBuilder nodesMap = new StringBuilder();
        for (int i = 0; i < tree.getNodeCount(); i++) {
            Node node = tree.getNode(i);
            if (!node.isRoot()) {
                Node parent = node.getParent();
                int pa = parent.getNr();
                // parent node .. child
                nodesMap.append(pa+1).append("..").append(i+1).append(",");
            }
        }
        // rm last ,
        if( nodesMap.length() > 0 )
            nodesMap.deleteCharAt( nodesMap.length() - 1 );
        return nodesMap.toString();
    }

    @Override
    public void close(PrintStream out) {
// nothing to do
    }

}

