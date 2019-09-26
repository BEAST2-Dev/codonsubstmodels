package beast.tree;

import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.evolution.alignment.CodonAlignment;

import java.io.PrintStream;

/**
 * The large 2-d matrix to store internal node sequences.
 *
 * @author Walter Xie, Fabio Mendes
 */
public class InternalNodeSeqs extends IntegerParameter {

    final public Input<CodonAlignment> dataInput = new Input<>("data",
            "codon alignment to initialise internal node sequences", Input.Validate.REQUIRED);


    @Override
    public void initAndValidate() {
        CodonAlignment codonAlignment = dataInput.get();
        int stateCount = codonAlignment.getDataType().getStateCount(); // 64
        assert stateCount == 64;

        // 0 - 63, ignore lowerValueInput upperValueInput
        m_fLower = 0;
        m_fUpper = stateCount-1; // not deal unknown codon

        // internal nodes = n - 1
        minorDimension = codonAlignment.getTaxonCount()-1;
        int L = codonAlignment.getSiteCount();

        // length = (n-1)*L
        this.values = new Integer[minorDimension * L];
        this.storedValues = values.clone();

        m_bIsDirty = new boolean[values.length];

    }

    @Override
    public void log(long sampleNr, PrintStream out) {
        super.log(sampleNr, out);
    }
}
