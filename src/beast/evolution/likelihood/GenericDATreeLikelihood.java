package beast.evolution.likelihood;


import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.sitemodel.SiteModelInterface;
import beast.evolution.tree.NodesStates;
import beast.evolution.tree.TreeInterface;

import java.util.Collections;
import java.util.List;
import java.util.Random;


@Description("Generic data augmentation tree likelihood for data in all nodes " +
		"given a beast tree, a generic SiteModel and a branch rate model")
// Replace Alignment and TreeInterface into NodesStatesAndTree.
// Override Distribution.calculatLogP() to make this class functional.
public class GenericDATreeLikelihood extends Distribution {
    
//    final public Input<Alignment> dataInput = new Input<>("data", "sequence data for the beast.tree", Validate.REQUIRED);

	final public Input<NodesStates> nodesStatesInput = new Input<>("nodesStates",
			"States in all nodes with the beast.tree", Validate.REQUIRED);

	final public Input<TreeInterface> treeInput = new Input<>("tree",
			"phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);

	final public Input<SiteModelInterface> siteModelInput = new Input<>("siteModel",
			"site model for leafs in the beast.tree", Validate.REQUIRED);
    
    final public Input<BranchRateModel.Base> branchRateModelInput = new Input<>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.");

    /** data, tree and models **/
	protected NodesStates nodesStates;
	protected TreeInterface tree;
	protected SiteModel.Base siteModel;
	protected BranchRateModel.Base branchRateModel;


	/**
	 * @return a list of unique ids for the state nodes that form the argument
	 */
	@Override
	public List<String> getArguments() {
		return Collections.singletonList(nodesStates.getID());
	}

	/**
	 * @return a list of unique ids for the state nodes that make up the conditions
	 */
	@Override
	public List<String> getConditions() {
		return siteModel.getConditions();
	}

	@Override
	public void sample(State state, Random random) {}

}
