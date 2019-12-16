package beast.evolution.likelihood;


import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.sitemodel.SiteModelInterface;
import beast.evolution.substitutionmodel.SubstitutionModel;
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

	final public Input<NodesStates> nodesStatesInput = new Input<>("nodesStates",
			"States in all nodes with the beast.tree", Validate.REQUIRED);

	final public Input<TreeInterface> treeInput = new Input<>("tree",
			"phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);

	final public Input<SiteModelInterface> siteModelInput = new Input<>("siteModel",
			"site model for leafs in the beast.tree", Validate.REQUIRED);
    
    final public Input<BranchRateModel.Base> branchRateModelInput = new Input<>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.");

    /** data, tree and models **/
	/**
	 * Data : states in all nodes {@link NodesStates}
	 */
    protected NodesStates nodesStates;
	/**
	 * {@link TreeInterface}
	 */
	protected TreeInterface tree;
	/**
	 * {@link SiteModel}
	 */
	protected SiteModel.Base siteModel;
	/**
	 * Substitution model, which is from {@link SiteModel}.
	 */
	protected SubstitutionModel substitutionModel;
	/**
	 * molecular clock model, default to {@link StrictClockModel}
	 */
	protected BranchRateModel.Base branchRateModel;

	@Override
	public void initAndValidate() {
		// data
		nodesStates = nodesStatesInput.get();
		// tree
		tree = treeInput.get();

		// models
		if (!(siteModelInput.get() instanceof SiteModel.Base)) {
			throw new IllegalArgumentException("siteModel input should be of type SiteModel.Base");
		}
		siteModel = (SiteModel.Base) siteModelInput.get();
		// set data type from NodesStates
		siteModel.setDataType(nodesStates.getCodonDataType());
		substitutionModel = siteModel.getSubstitutionModel();

		if (branchRateModelInput.get() != null) {
			branchRateModel = branchRateModelInput.get();
		} else {
			branchRateModel = new StrictClockModel();
		}
	}



	public NodesStates getNodesStates() {
		return nodesStates;
	}

	public TreeInterface getTree() {
		return tree;
	}

	public SiteModel.Base getSiteModel() {
		return siteModel;
	}

	public SubstitutionModel getSubstitutionModel() {
		return substitutionModel;
	}

	public BranchRateModel.Base getBranchRateModel() {
		return branchRateModel;
	}

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
