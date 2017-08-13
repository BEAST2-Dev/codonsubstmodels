package beast.app.beauti;

import beast.app.draw.BEASTObjectInputEditor;
import beast.core.BEASTInterface;
import beast.core.Input;

import javax.swing.*;

public class CodonModelsInputEditor extends BEASTObjectInputEditor {
	private static final long serialVersionUID = 1L;

	public CodonModelsInputEditor(BeautiDoc doc) {
		super(doc);
	}

	@Override
	public Class<?> type() {
		return codonmodels.AbstractCodonModel.class;
	}

	@Override
	public void init(Input<?> input, BEASTInterface beastObject, int itemNr, ExpandOption isExpandOption,
			boolean addButtons) {
		super.init(input, beastObject, itemNr, isExpandOption, addButtons);
		Object o = getComponent(1);
		if (o instanceof Box) {
			o = ((Box)o).getComponent(0);
			if (o instanceof Box) {
				String label = ((BEASTInterface)input.get()).getDescription();
				label = label.replaceAll("\\n", "<br/>");
				((Box)o).add(new JLabel("<html>" + label + "</html>"));
			}
		}
	}

}
