package beast.app.beauti;

import beast.core.BEASTInterface;
import beast.core.Input;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.CodonAlignment;


public class CodonAlignmentListInputEditor extends AlignmentListInputEditor {
	private static final long serialVersionUID = 1L;

	public CodonAlignmentListInputEditor(BeautiDoc doc) {
		super(doc);
	}

	@Override
	@SuppressWarnings("unchecked")
	public void init(Input<?> input, BEASTInterface beastObject, int itemNr, ExpandOption isExpandOption, boolean addButtons) {
		super.init(input, beastObject, itemNr, isExpandOption, addButtons);

		// TODO delItem() not working
		for (Alignment alignment : alignments) {
            if (alignment instanceof CodonAlignment) {
                splitButton.setVisible(false);
                replaceButton.setVisible(false);
                break;
            }
        }

	}


} // class AlignmentListInputEditor

