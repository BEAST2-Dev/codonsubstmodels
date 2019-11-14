package beast.evolution.likelihood;


import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.sitemodel.SiteModelInterface;
import beast.evolution.tree.NodesStatesAndTree;

import java.util.List;
import java.util.Random;


@Description("Generic data augmentation tree likelihood for data in all nodes " +
		"given a beast tree, a generic SiteModel and a branch rate model")
// Replace Alignment and TreeInterface into NodesStatesAndTree.
// Override Distribution.calculatLogP() to make this class functional.
public class GenericDATreeLikelihood extends Distribution {
    
//    final public Input<Alignment> dataInput = new Input<>("data", "sequence data for the beast.tree", Validate.REQUIRED);
//
//    final public Input<TreeInterface> treeInput = new Input<>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);

	final public Input<NodesStatesAndTree> dataAndTreeInput = new Input<>("dataAndTree",
			"States in all nodes with the beast.tree", Validate.REQUIRED);

    final public Input<SiteModelInterface> siteModelInput = new Input<>("siteModel",
			"site model for leafs in the beast.tree", Validate.REQUIRED);
    
    final public Input<BranchRateModel.Base> branchRateModelInput = new Input<>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.");

    
    
	@Override
	public List<String> getArguments() {return null;}

	@Override
	public List<String> getConditions() {return null;}

	@Override
	public void sample(State state, Random random) {}

}
