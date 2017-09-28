package beast.app.beauti;

import beast.core.BEASTInterface;
import beast.core.Input;


public class CodonAlignmentListInputEditor extends AlignmentListInputEditor {
	private static final long serialVersionUID = 1L;

	public CodonAlignmentListInputEditor(BeautiDoc doc) {
		super(doc);
	}

	@Override
	@SuppressWarnings("unchecked")
	public void init(Input<?> input, BEASTInterface beastObject, int itemNr, ExpandOption isExpandOption, boolean addButtons) {
		super.init(input, beastObject, itemNr, isExpandOption, addButtons);

        splitButton.setVisible(false);
        replaceButton.setVisible(false);
	}


} // class AlignmentListInputEditor

