package beast.app.beauti;


import beast.app.draw.BEASTObjectInputEditor;
import beast.core.BEASTInterface;
import beast.core.Input;
import codonmodels.CodonFrequencies;

import javax.swing.*;
import java.util.Arrays;


public class CodonFrequenciesInputEditor extends BEASTObjectInputEditor {
//    RealParameter freqsParameter;
//    CodonAlignment alignment;

    private static final long serialVersionUID = 1L;
//    boolean useDefaultBehavior;

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
    protected void addComboBox(JComponent box, Input<?> input, BEASTInterface beastObject) {
        CodonFrequencies freqs = (CodonFrequencies) input.get();

        String[] freqsString = new String[]{"equal", "F1X4", "F3X4", "F6n"};
        JComboBox<String> comboBox = new JComboBox<>(freqsString);
        String pi = freqs.piInput.get();
        int idx = Arrays.asList(freqsString).indexOf(pi);
        if (idx >= 0) {
            comboBox.setSelectedIndex(idx);
        } else {
            comboBox.setSelectedIndex(0);
        }
//        alignment = (CodonAlignment) getCandidate(null, freqs);
        comboBox.addActionListener(e -> {
                //@SuppressWarnings("unchecked")
                int selected = comboBox.getSelectedIndex();
                try {
                    freqs.piInput.setValue(freqsString[selected], freqs);
                } catch (Exception e2) {
                    e2.printStackTrace();
                }
                //System.err.println(freqs.frequencies.get() + " " + freqs.m_data.get() + " " + freqs.m_bEstimate.get());
            });
        box.add(comboBox);
    }

//    private BEASTInterface getCandidate(Input<?> input, CodonFrequencies freqs) {
//        return getDoc().getPartition(freqs);
////		List<String> candidates = PluginPanel.getAvailablePlugins(input, freqs, null);
////		String id = candidates.get(0);
////		BEASTObject beastObject = PluginPanel.g_plugins.get(id);
////		return beastObject;
//    }


    @Override
    /** suppress input label**/
    protected void addInputLabel() {
        super.addInputLabel();
    }

}
