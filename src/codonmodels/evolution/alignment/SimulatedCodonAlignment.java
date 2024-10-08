package codonmodels.evolution.alignment;

import java.io.FileNotFoundException;
import java.io.PrintStream;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.parser.XMLProducer;
import beast.base.util.Randomizer;
import codonmodels.evolution.datatype.Codon;
import codonmodels.evolution.datatype.GeneticCode;



/**
 * @author remco@cs.waikato.ac.nz
 */
@Description("An alignment containing sequences randomly generated using a"
        + "given site model down a given tree.")
public class SimulatedCodonAlignment extends CodonAlignment {
   // final public Input<Alignment> m_data = new Input<>("data", "alignment data which specifies datatype and taxa of the beast.tree", Validate.REQUIRED);
    final public Input<Tree> m_treeInput = new Input<>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);
    final public Input<SiteModel.Base> m_pSiteModelInput = new Input<>("siteModel", "site model for leafs in the beast.tree", Validate.REQUIRED);
    final public Input<BranchRateModel.Base> m_pBranchRateModelInput = new Input<>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.");
    final public Input<Integer> m_sequenceLengthInput = new Input<>("sequencelength", "nr of samples to generate (default 1000).", 1000);
    final public Input<String> m_outputFileNameInput = new Input<>(
            "outputFileName",
            "If provided, simulated alignment is additionally written to this file.");    

    /**
     * nr of samples to generate *
     */
    protected int m_sequenceLength;
    /**
     * tree used for generating samples *
     */
    protected Tree m_tree;
    /**
     * site model used for generating samples *
     */
    protected SiteModel.Base m_siteModel;
    /**
     * branch rate model used for generating samples *
     */
    protected BranchRateModel m_branchRateModel;
    /**
     * nr of categories in site model *
     */
    int m_categoryCount;
    /**
     * nr of states in site model *
     */
    int m_stateCount;

    /**
     * name of output file *
     */
    String m_outputFileName;    

    /**
     * an array used to transfer transition probabilities
     */
    protected double[][] m_probabilities;
    
    public SimulatedCodonAlignment() {
        
        // Override the sequence input requirement.
        sequenceInput.setRule(Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {
        m_tree = m_treeInput.get();
        m_siteModel = m_pSiteModelInput.get();
        m_branchRateModel = m_pBranchRateModelInput.get();
        m_sequenceLength = m_sequenceLengthInput.get();

        GeneticCode geneticCode = GeneticCode.findByName(geneticCodeInput.get());
        setGeneticCode(geneticCode);

        m_stateCount = geneticCode.getStateCount();
        m_categoryCount = m_siteModel.getCategoryCount();
        m_probabilities = new double[m_categoryCount][m_stateCount * m_stateCount];
        m_outputFileName = m_outputFileNameInput.get();
        
        sequenceInput.get().clear();
        
        simulate();        
        
        // Write simulated alignment to disk if requested:
        if (m_outputFileName != null) {
            PrintStream pstream;
			try {
				pstream = new PrintStream(m_outputFileName);
	            pstream.println(new XMLProducer().toRawXML(this));
	            pstream.close();
			} catch (FileNotFoundException e) {
				throw new IllegalArgumentException(e.getMessage());
			}
        }
        
        alignment = alignmentInput.get();
        alignment.sequenceInput.get().clear();
        alignment.sequenceInput.get().addAll(sequenceInput.get());
        
        super.initAndValidate();
    }

    /**
     * Convert integer representation of sequence into a Sequence
     *
     * @param seq  integer representation of the sequence
     * @param node used to determine taxon for sequence
     * @return Sequence
     */
    Sequence intArray2Sequence(int[] seq, Node node) {
        String seqString = m_dataType.encodingToString(seq);
        
        // Find taxon with same name if tree is labelled
        int taxonNum = node.getNr();
        if (node.getID() != null && !node.getID().isEmpty()) {
	        for (int i = 0; i < alignmentInput.get().getTaxaNames().size(); i ++) {
	        	if (alignmentInput.get().getTaxaNames().get(i).equals(node.getID())) {
	        		taxonNum=i;
	        		break;
	        	}
	        }
        }
        
        String taxon = alignmentInput.get().getTaxaNames().get(taxonNum);
        
        
        return new Sequence(taxon, seqString);
    } // intArray2Sequence

    /**
     * perform the actual sequence generation
     *
     * @return alignment containing randomly generated sequences for the nodes in the
     *         leaves of the tree
     */
    public void simulate() {
        Node root = m_tree.getRoot();


        double[] categoryProbs = m_siteModel.getCategoryProportions(root);
        int[] category = new int[m_sequenceLength];
        for (int i = 0; i < m_sequenceLength; i++) {
            category[i] = Randomizer.randomChoicePDF(categoryProbs);
        }

        double[] frequencies = m_siteModel.getSubstitutionModel().getFrequencies();
        int[] seq = new int[m_sequenceLength];
        for (int i = 0; i < m_sequenceLength; i++) {
            seq[i] = Randomizer.randomChoicePDF(frequencies);
        }

        traverse(root, seq, category);

    } // simulate

    /**
     * recursively walk through the tree top down, and add sequence to alignment whenever
     * a leave node is reached.
     *
     * @param node           reference to the current node, for which we visit all children
     * @param parentSequence randomly generated sequence of the parent node
     * @param category       array of categories for each of the sites
     * @param alignment
     */
    void traverse(Node node, int[] parentSequence, int[] category) {
        for (int childIndex = 0; childIndex < 2; childIndex++) {
            Node child = (childIndex == 0 ? node.getLeft() : node.getRight());
            for (int i = 0; i < m_categoryCount; i++) {
                getTransitionProbabilities(m_tree, child, i, m_probabilities[i]);
            }

            int[] seq = new int[m_sequenceLength];
            double[] cProb = new double[m_stateCount];
            for (int i = 0; i < m_sequenceLength; i++) {
                System.arraycopy(m_probabilities[category[i]], parentSequence[i] * m_stateCount, cProb, 0, m_stateCount);
                seq[i] = Randomizer.randomChoicePDF(cProb);
            }

            if (child.isLeaf()) {
                sequenceInput.setValue(intArray2Sequence(seq, child), this);
            } else {
                traverse(child, seq, category);
            }
        }
    } // traverse

    /**
     * get transition probability matrix for particular rate category *
     */
    void getTransitionProbabilities(Tree tree, Node node, int rateCategory, double[] probs) {

        Node parent = node.getParent();
        double branchRate = (m_branchRateModel == null ? 1.0 : m_branchRateModel.getRateForBranch(node));
        branchRate *= m_siteModel.getRateForCategory(rateCategory, node);

        m_siteModel.getSubstitutionModel().getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), branchRate, probs);

    } // getTransitionProbabilities


    @Override
    	public Codon getDataType() {
    		if (m_dataType == null) {
    			initDataType();
    		}
    		return super.getDataType();
    	}
    
} // class SequenceAlignment

