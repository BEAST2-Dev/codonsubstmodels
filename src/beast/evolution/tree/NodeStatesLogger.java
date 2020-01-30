package beast.evolution.tree;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;

import java.io.PrintStream;

@Description("Logs tree annotated with metadata and/or rates")
public class NodeStatesLogger extends BEASTObject implements Loggable {
//    final public Input<Tree> treeInput = new Input<>("tree",
//            "tree to provide node ID while logging.", Validate.REQUIRED);
    final public Input<NodeStatesArray> nodesStatesInput = new Input<>("nodesStates",
            "Internal node states to be logged", Validate.REQUIRED);

    final public Input<Boolean> compressInput = new Input<>("compress",
            "if require to compress log", false);

    private NodeStatesArray nsa;
    private int tipCount;
    private int nodeCount;

    private boolean isCompress;

    @Override
    public void initAndValidate() {
        nsa = nodesStatesInput.get();
        tipCount = nsa.getTipsCount();
        nodeCount = nsa.getNodeCount();

        isCompress = compressInput.get();
    }

    @Override
    public void init(PrintStream out) {
        out.println("#Sequences\t" + nsa.getInternalNodeCount());
        out.println("#Sites\t" + nsa.getSiteCount());
        out.println("#States\t" + nsa.getStateCount());

        out.println("Sample\tNode\tStates");
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
        out.print(sampleNr);
        // only internal nodes
        for (int i=tipCount; i < nodeCount; i++) {
            out.print("\t");
            nsaCurrent.getNodeStates(i).log(out);
            out.print("\n");
        }
    }


    @Override
    public void close(PrintStream out) {
// nothing to do
    }

}

