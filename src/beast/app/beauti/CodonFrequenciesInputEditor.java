package beast.app.beauti;


import beastfx.app.inputeditor.BEASTObjectInputEditor;
import beastfx.app.inputeditor.BeautiDoc;
import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import beast.evolution.alignment.CodonAlignment;
import codonmodels.CodonFrequencies;
import javafx.scene.control.ComboBox;
import javafx.scene.layout.Pane;

import java.util.Arrays;


public class CodonFrequenciesInputEditor extends BEASTObjectInputEditor {
    RealParameter freqsParameter;
    CodonAlignment alignment;


	public CodonFrequenciesInputEditor(BeautiDoc doc) {
		super(doc);
	}

    @Override
    public Class<?> type() {
        return CodonFrequencies.class;
        //return Frequencies.class;
    }

    @Override
    public void init(Input<?> input, BEASTInterface beastObject, int itemNr, ExpandOption isExpandOption, boolean addButtons) {
        super.init(input, beastObject, itemNr, isExpandOption, addButtons);
    } // init


    @Override
    /** suppress combobox **/
    protected void addComboBox(Pane box, Input<?> input, BEASTInterface beastObject0) {
        CodonFrequencies freqs = (CodonFrequencies) input.get();

        String[] freqsString = new String[]{"equal", "F1X4", "F3X4", "F60/F61"};
        ComboBox<String> comboBox = new ComboBox<>();
        comboBox.getItems().addAll(freqsString);
        
        String pi = freqs.piInput.get();
        int idx = Arrays.asList(freqsString).indexOf(pi);
        if (idx > 0) {
            comboBox.getSelectionModel().select(idx);
        } else {
            comboBox.getSelectionModel().select(0);
        }
        freqsParameter = freqs.frequenciesInput.get();
        alignment = (CodonAlignment) getCandidate(freqs.dataInput, freqs);

        comboBox.setOnAction(e -> {
                //@SuppressWarnings("unchecked")
                int selected = comboBox.getSelectionModel().getSelectedIndex();
                try {
                    freqs.piInput.setValue(freqsString[selected], freqs);
                    freqs.frequenciesInput.setValue(freqsParameter, freqs);
                } catch (Exception e2) {
                    e2.printStackTrace();
                }
                //System.err.println(freqs.frequencies.get() + " " + freqs.m_data.get() + " " + freqs.m_bEstimate.get());
            });
        box.getChildren().add(comboBox);
    }

    private BEASTInterface getCandidate(Input<?> input, CodonFrequencies freqs) {
        return getDoc().getPartition(freqs);
//		List<String> candidates = PluginPanel.getAvailablePlugins(input, freqs, null);
//		String id = candidates.get(0);
//		BEASTObject beastObject = PluginPanel.g_plugins.get(id);
//		return beastObject;
    }


    @Override
    /** suppress input label**/
    protected void addInputLabel() {
        super.addInputLabel();
    }

}
