package codonmodels.app.beauti;


import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.alignment.Alignment;
import beastfx.app.inputeditor.BeautiAlignmentProvider;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.util.Alert;
import codonmodels.evolution.alignment.CodonAlignment;
import codonmodels.evolution.datatype.GeneticCode;
import javafx.scene.control.ChoiceDialog;

import java.util.LinkedHashSet;
import java.util.List;
import java.util.Optional;
import java.util.Set;


@Description("Class to convert nucleotide alignments to codon alignments,\n" +
        "using the same data and template but output codon data type.")
public class CodonAlignmentProvider extends BeautiAlignmentProvider {

    protected void addAlignments(BeautiDoc doc, List<BEASTInterface> selectedBEASTObjects) {
        ChoiceDialog<String> dialog = new ChoiceDialog<>(
                GeneticCode.GENETIC_CODE_DESCRIPTIONS[GeneticCode.UNIVERSAL_ID],
                GeneticCode.GENETIC_CODE_DESCRIPTIONS);
        dialog.setHeaderText("Importing codon alignment");
        dialog.setTitle("Genetic Code");
        dialog.setContentText("Choose Genetic Code : ");

        Optional<String> choice = dialog.showAndWait();
        // GENETIC_CODE_DESCRIPTIONS
        loadCodonAlignment(choice.get(), doc, selectedBEASTObjects);

    }

    private void loadCodonAlignment(String geneticCodeDesc, BeautiDoc doc, List<BEASTInterface> selectedBEASTObjects) {
        if (geneticCodeDesc == null)
            Alert.showMessageDialog(null, "Fail to select a genetic code (null) ! ",
                    "Null Exception", Alert.ERROR_MESSAGE);
        else {
            GeneticCode geneticCode = GeneticCode.findByDescription(geneticCodeDesc);
            System.out.println("Choose genetic code " + geneticCode.getName());

            for (BEASTInterface beastObject : selectedBEASTObjects) {
                if (beastObject instanceof Alignment &&
                        ((Alignment) beastObject).getDataType().getTypeDescription() == "nucleotide") {
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
                    try {
                        // CodonAlignment wraps alignment
                        CodonAlignment codonAlignment = new CodonAlignment((Alignment) beastObject, geneticCode);
                        sortByTaxonName(((Alignment) beastObject).sequenceInput.get());
                        if (getStartTemplate() != null)
                            doc.addAlignmentWithSubnet(codonAlignment, getStartTemplate());

                    } catch (RuntimeException e) {
                        Alert.showMessageDialog(null,
                                "Cannot init codon alignment given genetic code " + geneticCode.getName() +
                                        " !\nPlease check if any stop codon exists in your alignment.\n\n" +
                                        e.getMessage());
                    }
                }
            }// end for
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

    @Override
    public
    void editAlignment(Alignment alignment, BeautiDoc doc) {
        if (alignment instanceof CodonAlignment codonAlignment) {
            try {
                CodonAlignmentViewer viewer = new CodonAlignmentViewer(codonAlignment);
                viewer.showInDialog();
            } catch (Exception e) {
                Alert.showMessageDialog(null,
                        "Something went wrong viewing the codon alignment: " + e.getMessage());
                e.printStackTrace();
            }
        } else {
            super.editAlignment(alignment, doc);
        }
    }
}
