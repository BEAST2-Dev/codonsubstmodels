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
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.Nucleotide;
import beast.util.AddOnManager;
import beast.evolution.datatype.Codon;
import beast.evolution.datatype.GeneticCode;

import java.util.Arrays;
import java.util.List;

/**
 * An alignment class that takes another alignment and converts it on the fly
 * to a different dataType.
 * Currently only working on nucleotide => codon
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 *
 * Modified from BEAST 1 ConvertAlignment.
 */
public class ConvertAlignment extends Alignment { //TODO should have WrappedAlignment between ConvertAlignment and Alignment

    final public Input<Alignment> alignmentInput = new Input<>("data",
            "alignment to convert to new type specified by userDataType", Input.Validate.REQUIRED);

    final public Input<String> geneticCodeInput = new Input<>("geneticCode",
            "The rule to define how sequences of nucleotide triplets, " +
                    "called codons, specify which amino acid will be added next during protein synthesis.",
                    GeneticCode.GENETIC_CODE_NAMES[GeneticCode.UNIVERSAL_ID]);

    protected Alignment alignment;
//    protected GeneticCode geneticCode;

    public ConvertAlignment() {
        sequenceInput.setRule(Input.Validate.OPTIONAL);
        userDataTypeInput.setRule(Input.Validate.FORBIDDEN); // avoid confusion
        siteWeightsInput.setRule(Input.Validate.FORBIDDEN); // ?
    }


    @Override
    public void initAndValidate() {

        alignment = alignmentInput.get();
        DataType originalType = alignment.getDataType();

        initDataType();
        DataType newType = this.getDataType();
        //TODO generalise
        if (newType != null && newType instanceof Codon && originalType instanceof Nucleotide) {
            m_dataType = newType;
        } else {
            throw new UnsupportedOperationException("Currently only working on nucleotide => codon !");
        }

        // set geneticCode
        GeneticCode geneticCode = GeneticCode.findByName(geneticCodeInput.get());
        GeneticCode geneticCode2 = ((Codon) newType).getGeneticCode();
        if (! geneticCode.getName().equals( geneticCode2.getName() ) ) { //newType instanceof Codon &&
            ((Codon) newType).setGeneticCode(geneticCode);
        }

        convertCodonToState(true);

        if (alignmentInput.get().siteWeightsInput.get() != null) {
            String str = alignmentInput.get().siteWeightsInput.get().trim();
            String [] strs = str.split(",");
            siteWeights = new int[strs.length];
            for (int i = 0; i< strs.length; i++) {
                siteWeights[i] = Integer.parseInt(strs[i].trim());
            }
        }


        calcPatterns();
        setupAscertainment();
    }

    //TODO this should move to Alignment
    // Note: Codon uses default genetic code, need to change it if user specifies in xml
    protected void initDataType() {
        if (types.indexOf(dataTypeInput.get()) < 0) {
            throw new IllegalArgumentException("data type + '" + dataTypeInput.get() + "' cannot be found. " +
                    "Choose one of " + Arrays.toString(types.toArray(new String[0])));
        }
        // seems to spend forever in there??
        List<String> dataTypes = AddOnManager.find(beast.evolution.datatype.DataType.class, IMPLEMENTATION_DIR);
        for (String dataTypeName : dataTypes) {
            DataType dataType;
            try {
                dataType = (DataType) Class.forName(dataTypeName).newInstance();
                if (dataTypeInput.get().equals(dataType.getTypeDescription())) {
                    m_dataType = dataType;
                    break;
                }
            } catch (InstantiationException | IllegalAccessException | ClassNotFoundException e) {
                throw new IllegalArgumentException(e.getMessage());
            }
        }
    }

    protected void convertCodonToState(boolean log) {
        taxaNames.clear();
        stateCounts.clear();
        counts.clear();
        try {
            for (Sequence seq : alignment.sequences) {
                List<Integer> codonStates = seq.getSequence(getDataType()); // return mapCodeToStateSet indices i
                int tripletIndex = findStopCodon(codonStates);
                if (tripletIndex > -1)
                    throw new RuntimeException(seq.getTaxon() + " sequence contains stop codon at " +
                            (tripletIndex+1) + "th triplets ! \n" +
                            "Please either use codon alignment or provide a correct genetic code.");

                counts.add(codonStates);
                if (taxaNames.contains(seq.getTaxon())) {
                    throw new RuntimeException("Duplicate taxon found in alignment: " + seq.getTaxon());
                }
                taxaNames.add(seq.getTaxon());
                //TODO tipLikelihoods ?
//                tipLikelihoods.add(seq.getLikelihoods());
//                // if seq.isUncertain() == false then the above line adds 'null'
//                // to the list, indicating that this particular sequence has no tip likelihood information
//                usingTipLikelihoods |= (seq.getLikelihoods() != null);

                stateCounts.add(getDataType().getStateCount());

            }
            if (counts.size() == 0) {
                // no sequence data
                throw new RuntimeException("Sequence data expected, but none found");
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }
    }

    protected int findStopCodon(List<Integer> seqStates) {
        for (int i = 0; i < seqStates.size(); i++) {
            int codonState = seqStates.get(i);
            if (getGeneticCode().isStopCodon(codonState))
                return i;
        }
        return -1;
    }

    public GeneticCode getGeneticCode() {
        if (m_dataType instanceof Codon)
            return ((Codon) m_dataType).getGeneticCode();
        return null;
    }

    /**
     * @return number of sites
     */
    public int getSiteCount() {
        if (alignment == null) throw new RuntimeException("ConvertAlignment has no alignment");

        DataType originalType = alignment.getDataType();
        int count = alignment.getSiteCount();

        if (originalType instanceof Nucleotide) {
            count /= 3;
        }

        return count;
    }



}
