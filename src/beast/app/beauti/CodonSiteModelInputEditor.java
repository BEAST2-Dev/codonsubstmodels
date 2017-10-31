package beast.app.beauti;


import beast.core.BEASTInterface;
import beast.core.Input;
import codonmodels.CodonSiteModel;

public class CodonSiteModelInputEditor extends SiteModelInputEditor {
    private static final long serialVersionUID = 1L;

	public CodonSiteModelInputEditor(BeautiDoc doc) {
		super(doc);
	}

    @Override
    public void init(Input<?> input, BEASTInterface beastObject, int itemNr,
    		ExpandOption isExpandOption, boolean addButtons) {
    	super.init(input, beastObject, itemNr, isExpandOption, addButtons);
        fixMeanRatesCheckBox.setVisible(false);
    }

	@Override
    public Class<?> [] types() {
		Class<?>[] types = {CodonSiteModel.class};
		return types;
    }

}
