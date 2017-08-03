/*
 * ConvertAlignment.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package beast.evolution.alignment;


import beast.core.Input;
import beast.evolution.datatype.*;

/**
 * An alignment class that takes another alignment and converts it on the fly
 * to a different dataType.
 * Currently only working on nucleotide => codon
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 *
 * Modified from BEAST 1 WrappedAlignment.
 */
public class WrappedAlignment extends Alignment {

    final public Input<Alignment> alignmentInput = new Input<>("data",
            "alignment to convert to new type specified by userDataType", Input.Validate.REQUIRED);

    final public Input<String> geneticCodeInput = new Input<>("geneticCode",
            "The rule to define how sequences of nucleotide triplets, " +
                    "called codons, specify which amino acid will be added next during protein synthesis.",
                    GeneticCode.GENETIC_CODE_NAMES[GeneticCode.UNIVERSAL_ID], Input.Validate.REQUIRED);

    protected Alignment alignment;
    protected GeneticCode geneticCode;

    public WrappedAlignment() {
        sequenceInput.setRule(Input.Validate.OPTIONAL);
        dataTypeInput.setRule(Input.Validate.OPTIONAL);
        siteWeightsInput.setRule(Input.Validate.FORBIDDEN); // ?
    }


    @Override
    public void initAndValidate() {

        alignment = alignmentInput.get();
        DataType originalType = alignment.getDataType();

        DataType newType = userDataTypeInput.get();

        //TODO generalise
        if (newType != null && newType instanceof Codon && originalType instanceof Nucleotide) {
            m_dataType = newType;
        } else {
            throw new UnsupportedOperationException("Currently only working on nucleotide => codon !");
        }

        this.geneticCode = GeneticCode.findByName(geneticCodeInput.get());

        calcPatterns();
        setupAscertainment();
    }


    /**
     * @return number of sites
     */
    public int getSiteCount() {
        if (alignment == null) throw new RuntimeException("WrappedAlignment has no alignment");

        DataType originalType = alignment.getDataType();
        int count = alignment.getSiteCount();

        if (originalType instanceof Nucleotide) {
            count /= 3;
        }

        return count;
    }

    /**
     * @return the sequence state at (taxon, site)
     */
    public int getState(int taxonIndex, int siteIndex) {
        if (alignment == null) throw new RuntimeException("WrappedAlignment has no alignment");

        DataType originalType = alignment.getDataType();
        DataType newType = getDataType();

        int state = 0;

        if (originalType instanceof Nucleotide) {
            int siteIndex3 = siteIndex * 3;
            int state1 = alignment.getState(taxonIndex, siteIndex3);
            int state2 = alignment.getState(taxonIndex, siteIndex3 + 1);
            int state3 = alignment.getState(taxonIndex, siteIndex3 + 2);

            if (newType instanceof Codon) {
                state = ((Codon)newType).getState(state1, state2, state3);
            } else { // newType instanceof Aminoacid
                state = geneticCode.getAminoAcidState(((Codon)newType).getCanonicalState(((Codon)newType).getState(state1, state2, state3)));
            }

        }
        else if (originalType instanceof Codon) {
            if (newType instanceof Aminoacid) {
                state = geneticCode.getAminoAcidState(alignment.getState(taxonIndex, siteIndex));
            } else { // newType instanceof Codon
                String string = Alignment.getSequence(alignment, taxonIndex);
                state = geneticCode.getNucleotideState(string.charAt(siteIndex));
            }
        }

        return state;
    }



}
