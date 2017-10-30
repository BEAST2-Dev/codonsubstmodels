package beast.app.beauti;


import beast.core.BEASTInterface;
import beast.core.Description;
import beast.core.Input;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.CodonAlignment;
import beast.evolution.datatype.GeneticCode;

import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;


@Description("Class to convert nucleotide alignments to codon alignments,\n" +
        "using the same data and template but output codon data type.")
public class CodonAlignmentProvider extends BeautiAlignmentProvider {

    protected void addAlignments(BeautiDoc doc, List<BEASTInterface> selectedBEASTObjects) {
        for (BEASTInterface beastObject : selectedBEASTObjects) {
            if (beastObject instanceof Alignment &&
                    ((Alignment) beastObject).getDataType().getTypeDescription()=="nucleotide") {
                // ensure ID of alignment is unique
                int k = 0;
                String id = beastObject.getID();
                boolean found = true;
                while (doc.pluginmap.containsKey(id) && found) {
                    found = false;
                    for (Alignment data : doc.alignments) {
                        if (data.getID().equals(beastObject.getID())) {
                            found = true;
                            break;
                        }
                    }
                    if (found) {
                        k++;
                        id = beastObject.getID() + k;
                    } else {
                        BEASTInterface oldData = doc.pluginmap.get(beastObject.getID());
                        replaceItem(doc, oldData, beastObject);
                    }
                }
                beastObject.setID(id);
                // CodonAlignment wraps alignment
                CodonAlignment codonAlignment = new CodonAlignment((Alignment) beastObject, GeneticCode.UNIVERSAL);
                sortByTaxonName(((Alignment) beastObject).sequenceInput.get());
                if (getStartTemplate() != null) {
                    doc.addAlignmentWithSubnet(codonAlignment, getStartTemplate());
                }
            }
        }
    }

    @SuppressWarnings({ "rawtypes", "unchecked" })
    private void replaceItem(BeautiDoc doc, BEASTInterface oldData, BEASTInterface newData) {
        doc.pluginmap.remove(newData.getID());
        Set<BEASTInterface> outputs = new LinkedHashSet<>();
        outputs.addAll(oldData.getOutputs());
        for (BEASTInterface o : outputs) {
            for ( Input i : o.listInputs()) {
                if (i.get() == oldData) {
                    i.setValue(newData, o);
                } else if (i.get() != null && i.get() instanceof List) {
                    List list = (List) i.get();
                    int index = list.indexOf(oldData);
                    if (index >= 0) {
                        list.set(index, newData);
                        newData.getOutputs().add(o);
                    }
                }
            }
        }
    }


}
