package beast.app.beauti;


import beast.app.draw.*;
import beast.core.*;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.core.util.CompoundDistribution;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.CodonAlignment;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.operators.DeltaExchangeOperator;
import beast.evolution.sitemodel.SiteModel;

import javax.swing.*;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import java.awt.*;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.List;

public class ConvertAlignmentInputEditor extends BEASTObjectInputEditor {
    private static final long serialVersionUID = 1L;

    IntegerInputEditor categoryCountEditor;
    JTextField categoryCountEntry;
    InputEditor gammaShapeEditor;
    ParameterInputEditor inVarEditor;

    // vars for dealing with mean-rate delta exchange operator
    JCheckBox fixMeanRatesCheckBox;
    DeltaExchangeOperator operator;
    protected SmallLabel fixMeanRatesValidateLabel;

	public ConvertAlignmentInputEditor(BeautiDoc doc) {
		super(doc);
	}

    @Override
    public Class<?> type() {
        return CodonAlignment.class;
    }

    @Override
    public void init(Input<?> input, BEASTInterface beastObject, int itemNr,
    		ExpandOption isExpandOption, boolean addButtons) {
    	super.init(input, beastObject, itemNr, isExpandOption, addButtons);

		List<Operator> operators = ((MCMC) doc.mcmc.get()).operatorsInput.get();

//    	operator = (DeltaExchangeOperator) doc.pluginmap.get("FixMeanMutationRatesOperator");
//    	if (operator == null) {
//    		operator = new DeltaExchangeOperator();
//    		try {
//    			operator.setID("FixMeanMutationRatesOperator");
//				operator.initByName("weight", 2.0, "delta", 0.75);
//			} catch (Throwable e1) {
//				// ignore initAndValidate exception
//			}
//    		doc.addPlugin(operator);
//    	}
		fixMeanRatesCheckBox.setSelected(operators.contains(operator));
		Box box = Box.createHorizontalBox();
		box.add(fixMeanRatesCheckBox);
		box.add(Box.createHorizontalGlue());
		fixMeanRatesValidateLabel = new SmallLabel("x", Color.GREEN);
		fixMeanRatesValidateLabel.setVisible(false);
		box.add(fixMeanRatesValidateLabel);

    	if (doc.alignments.size() >= 1) {
        	JComponent component = (JComponent) getComponents()[0];
    		component.add(box);
    	}
		setUpOperator();
    }

    public InputEditor createMutationRateEditor() {
    	SiteModel sitemodel = ((SiteModel) m_input.get());
        final Input<?> input = sitemodel.muParameterInput;
        ParameterInputEditor mutationRateEditor = new ParameterInputEditor(doc);
        mutationRateEditor.init(input, sitemodel, -1, ExpandOption.FALSE, true);
        mutationRateEditor.getEntry().setEnabled(!doc.autoUpdateFixMeanSubstRate);
        return mutationRateEditor;
    }

	public InputEditor createGammaCategoryCountEditor() {
    	SiteModel sitemodel = ((SiteModel) m_input.get());
        final Input<?> input = sitemodel.gammaCategoryCount;
        categoryCountEditor = new IntegerInputEditor(doc) {
			private static final long serialVersionUID = 1L;

			@Override
			public void validateInput() {
        		super.validateInput();
            	SiteModel sitemodel = (SiteModel) m_beastObject;
                if (sitemodel.gammaCategoryCount.get() < 2 && sitemodel.shapeParameterInput.get().isEstimatedInput.get()) {
                	m_validateLabel.m_circleColor = Color.orange;
                	m_validateLabel.setToolTipText("shape parameter is estimated, but not used");
                	m_validateLabel.setVisible(true);
                }
        	};
        };

        categoryCountEditor.init(input, sitemodel, -1, ExpandOption.FALSE, true);
        categoryCountEntry = categoryCountEditor.getEntry();
        categoryCountEntry.getDocument().addDocumentListener(new DocumentListener() {
            @Override
            public void removeUpdate(DocumentEvent e) {
                processEntry2();
            }

            @Override
            public void insertUpdate(DocumentEvent e) {
                processEntry2();
            }

            @Override
            public void changedUpdate(DocumentEvent e) {
                processEntry2();
            }
        });

       	categoryCountEditor.validateInput();
        return categoryCountEditor;
    }

    void processEntry2() {
        String categories = categoryCountEntry.getText();
        try {
            int categoryCount = Integer.parseInt(categories);
        	RealParameter shapeParameter = ((SiteModel) m_input.get()).shapeParameterInput.get();
            if (!gammaShapeEditor.getComponent().isVisible() && categoryCount >= 2) {
            	// we are flipping from no gamma to gamma heterogeneity accross sites
            	// so set the estimate flag on the shape parameter
            	shapeParameter.isEstimatedInput.setValue(true, shapeParameter);
            } else if (gammaShapeEditor.getComponent().isVisible() && categoryCount < 2) {
            	// we are flipping from with gamma to no gamma heterogeneity accross sites
            	// so unset the estimate flag on the shape parameter
            	shapeParameter.isEstimatedInput.setValue(false, shapeParameter);
            }
            Object o = ((ParameterInputEditor)gammaShapeEditor).getComponent();
            if (o instanceof ParameterInputEditor) {
	            ParameterInputEditor e = (ParameterInputEditor) o;
	            e.m_isEstimatedBox.setSelected(shapeParameter.isEstimatedInput.get());
            }
            gammaShapeEditor.getComponent().setVisible(categoryCount >= 2);
            repaint();
        } catch (NumberFormatException e) {
            // ignore.
        }
    }

    public InputEditor createShapeEditor() throws NoSuchMethodException, SecurityException, ClassNotFoundException, InstantiationException, IllegalAccessException, IllegalArgumentException, InvocationTargetException {
        final Input<?> input = ((SiteModel) m_input.get()).shapeParameterInput;
        gammaShapeEditor = doc.getInputEditorFactory().createInputEditor(input, (BEASTInterface) m_input.get(), doc);
        gammaShapeEditor.getComponent().setVisible(((SiteModel) m_input.get()).gammaCategoryCount.get() >= 2);
        return gammaShapeEditor;
    }


    public static boolean customConnector(BeautiDoc doc) {
 		try {
 	        DeltaExchangeOperator operator = (DeltaExchangeOperator) doc.pluginmap.get("FixMeanMutationRatesOperator");
 	        if (operator == null) {
 	        	return false;
 	        }

 	       	List<RealParameter> parameters = operator.parameterInput.get();
 	    	parameters.clear();
		   	//String weights = "";
		    CompoundDistribution likelihood = (CompoundDistribution) doc.pluginmap.get("likelihood");
		    boolean hasOneEstimatedRate = false;
		    List<String> rateIDs = new ArrayList<>();
		    List<Integer> weights = new ArrayList<>();
			for (Distribution d : likelihood.pDistributions.get()) {
				GenericTreeLikelihood treelikelihood = (GenericTreeLikelihood) d;
	    		Alignment data = treelikelihood.dataInput.get(); 
	    		int weight = data.getSiteCount();
	    		if (data.isAscertained) {
	    			weight -= data.getExcludedPatternCount();
	    		}
	    		if (treelikelihood.siteModelInput.get() instanceof SiteModel) {
		    		SiteModel siteModel = (SiteModel) treelikelihood.siteModelInput.get();
		    		RealParameter mutationRate = siteModel.muParameterInput.get();
		    		//clockRate.m_bIsEstimated.setValue(true, clockRate);
		    		if (mutationRate.isEstimatedInput.get()) {
		    			hasOneEstimatedRate = true;
		    			if (rateIDs.indexOf(mutationRate.getID()) == -1) {
			    			parameters.add(mutationRate);
			    			weights.add(weight);
			    			rateIDs.add(mutationRate.getID());
		    			} else {
		    				int k = rateIDs.indexOf(mutationRate.getID());
			    			weights.set(k,  weights.get(k) + weight);
		    			}
		    		}
	    		}
	    	}
			
			
		    IntegerParameter weightParameter;
			if (weights.size() == 0) {
		    	weightParameter = new IntegerParameter();
			} else {
				String weightString = "";
				for (int k : weights) {
					weightString += k + " ";
				}
		    	weightParameter = new IntegerParameter(weightString);
				weightParameter.setID("weightparameter");
				
			}
			weightParameter.isEstimatedInput.setValue(false, weightParameter);
	    	operator.parameterWeightsInput.setValue(weightParameter, operator);
	    	return hasOneEstimatedRate;
		} catch (Exception e) {
			
		}
		return false;
    }
    
    /** set up relative weights and parameter input **/
    public void setUpOperator() {

	}

}
