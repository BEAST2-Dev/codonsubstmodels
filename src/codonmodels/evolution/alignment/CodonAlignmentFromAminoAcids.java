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

package codonmodels.evolution.alignment;



import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.datatype.Aminoacid;
import beast.base.evolution.datatype.DataType;
import codonmodels.evolution.datatype.Codon;
import codonmodels.evolution.datatype.GeneticCode;

import java.util.ArrayList;
import java.util.List;

/**
 * An alignment class that takes nucleotide alignment and converts it
 * to a codon dataType.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 *
 * Modified from BEAST 1 ConvertAlignment.
 */
@Description("An alignment class that takes aminoacid alignment and converts it to an ambiguous codon dataType.")
public class CodonAlignmentFromAminoAcids extends CodonAlignment {

//    final public Input<Alignment> alignmentInput = new Input<>("data",
//            "Aminoacid alignment to convert into codon data type specified by dataType", Input.Validate.REQUIRED);
//
//    final public Input<String> geneticCodeInput = new Input<>("geneticCode",
//            "The rule to define how sequences of nucleotide triplets, " +
//                    "called codons, specify which amino acid will be added next during protein synthesis.",
//                    GeneticCode.GENETIC_CODE_NAMES[GeneticCode.UNIVERSAL_ID]);
//
//    final public Input<Boolean> verboseInput = new Input<>("verbose",
//            "Print the codon usage, etc.",
//            Boolean.TRUE);



    public CodonAlignmentFromAminoAcids() {
        sequenceInput.setRule(Input.Validate.OPTIONAL);
        userDataTypeInput.setRule(Input.Validate.FORBIDDEN); // avoid confusion
        siteWeightsInput.setRule(Input.Validate.FORBIDDEN); // ?
        // set default to Codon
        dataTypeInput.setValue(Codon.CODON, this);
    }

    public CodonAlignmentFromAminoAcids(Alignment alignment, GeneticCode geneticCode) {
        this();
        if (alignment.getID() != null)
            this.setID("ca." + alignment.getID()); // make CodonAlignment ID not null
        DataType alignmentType = alignment.getDataType();
        if (! (alignmentType instanceof Aminoacid) ) {
            throw new IllegalArgumentException("CodonAlignment currently only support to wrap the aminoacid alignment !");
        }

        alignmentInput.setValue(alignment, this);
        geneticCodeInput.setValue(geneticCode.getName(), this);

        if (this.m_dataType == null || !this.m_dataType.getTypeDescription().equals(Codon.CODON)) {
            this.m_dataType = new Codon(geneticCode);
            Log.info.println("Set CodonAlignment (" + this.getID() + ") to data type : " +
                    this.m_dataType + " - " + geneticCode.getDescription() + " !");
        }
        initAndValidate();
    }

    @Override
    public void initAndValidate() {
        alignment = alignmentInput.get();//Nucleotide
        DataType alignmentType = alignment.getDataType();

        // Note: Codon uses default genetic code, need to change it if user specifies in xml
        initDataType();

        DataType thisType = this.m_dataType;
        // only working for nucleotide => codon
        if (thisType != null && thisType instanceof Codon && alignmentType instanceof Aminoacid) {
            m_dataType = thisType;
        } else {
            throw new IllegalArgumentException("CodonAlignment only wraps the nucleotide alignment into codon alignment !");
        }

        // set geneticCode to geneticCodeInput.get() if different to dataType.getGeneticCode()
        GeneticCode geneticCode = GeneticCode.findByName(geneticCodeInput.get());
        setGeneticCode(geneticCode);

//        convertCodonToState(unknownCodeExceptionInput.get()); // default to true
        if (counts.size() < 1) {
            convertCodonToState();
        }

        // after convertCodonToState
        if (alignmentInput.get().siteWeightsInput.get() != null) {
            String str = alignmentInput.get().siteWeightsInput.get().trim();
            String [] strs = str.split(",");
            siteWeights = new int[strs.length];
            for (int i = 0; i< strs.length; i++) {
                siteWeights[i] = Integer.parseInt(strs[i].trim());
            }
        }

        // after checking siteWeightsInput
        // include calcPatterns and setup Ascertainment
        sanityCheckCalcPatternsSetUpAscertainment(true);

        Log.info.println("\nState count = " + getDataType().getStateCount() +
                ", using " + getGeneticCode().getDescription() + " genetic code, ");
        Log.info.println("ambiguous states start from " + getDataType().getStateCount() +
                " end at " + (getDataType().getStateCountAmbiguous()-1) );

        if (isAscertained)
            throw new UnsupportedOperationException("Ascertainment correction is not available !");
        
        
        
        
        String aaCodes = new Aminoacid().getCodeMap();
        String codonCodes = geneticCode.getCodeTableNoStopCodons();
        stateSet = new boolean[21][codonCodes.length()];
    	for (int i = 0; i < codonCodes.length(); i++) {
    		char codonChar = codonCodes.charAt(i);
    		int j = aaCodes.indexOf(codonChar);
    		stateSet[j][i] = true;
    	}
        for (int i = 0; i < codonCodes.length(); i++) {
        	stateSet[20][i] = true;
        }
    }
    
    private boolean [][] stateSet; 

    /**
     * modified from Alignment private initializeWithSequenceList(List<Sequence>, boolean)
     */
    protected void convertCodonToState() {
        taxaNames.clear();
        stateCounts.clear();
        counts.clear();
        try {
            for (Sequence seq : alignment.sequenceInput.get()) {
                // TODO need to distinguish ? and -, cannot just replace all ? to -
                String data = seq.getData().replaceAll("\\?", "-");
                seq.dataInput.setValue(data, seq);

                // no stop codon
                List<Integer> codonStates = null;
                try {
                // return mapCodeToStateSet indices i, also indices in Codon.CODON_TRIPLETS
                codonStates = seq.getSequence(alignment.getDataType());
                // unknownCodeException false to treat codons with partial ambiguities (-TA) as missing data
//                List<Integer> codonStates = getCodonStates(seq, getDataType(), unknownCodeException);
                } catch (IllegalArgumentException e) {
                    Log.err.println("Error: sequence " + seq.getTaxon() + " contains a stop codon : " +
                            e.getMessage() + " ! \n" +
                            "Please use the correct genetic code, current genetic code = " +
                            getGeneticCode().getDescription());
                    throw new RuntimeException(e);
                }

//                int stopCodon = findStopCodon(codonStates);
//                if (stopCodon > -1)
//                    Log.warning.println("Warning: " + seq.getTaxon() + " sequence contains a stop codon at " +
//                            (stopCodon+1) + "th triplets ! \n" +
//                            "Please either use a codon alignment or the correct genetic code.");

                counts.add(codonStates);
                if (taxaNames.contains(seq.getTaxon())) {
                    throw new RuntimeException("Duplicate taxon found in alignment: " + seq.getTaxon());
                }
                taxaNames.add(seq.getTaxon());

                if (seq.uncertainInput.get() != null && seq.uncertainInput.get())
                    throw new UnsupportedOperationException("Uncertain feature in sequence is not supported in this version !");

                    //TODO how to deal with uncertain? seq here is Nucleotide sequences not codons.
//                tipLikelihoods.add(seq.getLikelihoods());
//                // if seq.isUncertain() == false then the above line adds 'null'
//                // to the list, indicating that this particular sequence has no tip likelihood information
//                usingTipLikelihoods |= (seq.getLikelihoods() != null);

                stateCounts.add(getDataType().getStateCount());

            }
            Log.info.println(); // empty line in screen
            if (counts.size() == 0) {
                // no sequence data
                throw new RuntimeException("Sequence data expected, but none found");
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }
    }

     
    @Override
    public Codon getDataType() {
        if (!(m_dataType instanceof Codon))
            throw new UnsupportedOperationException("CodonAlignment only supports Codon data type !");
        return (Codon) m_dataType;
    }

    public GeneticCode getGeneticCode() {
        return getDataType().getGeneticCode();
    }

    public void setGeneticCode(GeneticCode geneticCode) {
        GeneticCode currentCode = getGeneticCode();
        if (! geneticCode.getName().equals( currentCode.getName() ) ) {
            Log.info.println("Change genetic code from " + currentCode.getName() + " to " + geneticCode.getName());
            ((Codon) m_dataType).setGeneticCode(geneticCode);
            geneticCodeInput.setValue(geneticCode.getName(), this);

            // refresh states after changing GeneticCode
            convertCodonToState();
        }
    }

    /**
     * @return number of codons
     */
    public int getSiteCount() {
        if (alignment == null) throw new RuntimeException("CodonAlignment has no alignment");
        int count = alignment.getSiteCount();
        return count;
    }

    /**
     * Retrieve the "weight" of a particular site.
     *
     * @param siteIndex_ Index into site array.
     * @return site weight
     */
    public int getSiteWeight(int siteIndex_) {
        return siteWeights[siteIndex_];
    }

    /**
     * Codon position * base (3x4) table, plus "overall" in last row.
     * The base order is "A", "C", "G", "T".
     * If any ambiguous, then add 0.25 to each Nucleotide count in 3 positions.
     * @return 2d frequency array [4][4]. Position is row, base is col. (3x4).
     *         The last row is sum over nucleotides, which can be used as 1x4.
     */
    public double[][] getObservedBaseFrequencies() {
        GeneticCode geneticCode = getGeneticCode();
        // position x base (3x4) table + average
        final int nrOfNuc = 4; // 4 nucleotides
        double[][] freqs = new double[4][nrOfNuc];
        for (int i = 0; i < counts.size(); i++) {
            List<Integer> codonStates = counts.get(i);
            for (int j = 0; j < codonStates.size(); j++) {
                int state = codonStates.get(j);
                // may have ambiguous nucleotide (i.e. 'A-T')
                String triplet = getDataType().encodingToString(new int[]{state});

                for (int pos = 0; pos < 3; pos++) {
                    // col index = nucState: A,C,G,T,-
                    int nucState = geneticCode.getNucleotideState(triplet.charAt(pos));
                    // if -, add 0.25 to each nucleotide
                    if (nucState >= nrOfNuc) {
                        for (int n = 0; n < nrOfNuc; n++)
                            freqs[pos][n] += 0.25;
                    } else {
                        freqs[pos][nucState] += 1;
                    }
                }
            }
        }
        // overall at last row
        for (int col = 0; col < nrOfNuc; col++) {
            freqs[3][col] = freqs[0][col] + freqs[1][col] + freqs[2][col];
        }
        // re-normalise
        for (int row = 0; row < 4; row++) {
            double rowSum = freqs[row][0];
            for (int col = 1; col < nrOfNuc; col++)
                rowSum += freqs[row][col];
            for (int col = 0; col < nrOfNuc; col++)
                freqs[row][col] = freqs[row][col] / rowSum;
        }
        return freqs;
    }

    @Override
    public boolean[] getStateSet(int state) {
    	if (state > 20) {
    		return stateSet[20];
    	}
    	if (state == 20) {
    		return stateSet[20];
    	}
    	return stateSet[state];
    }

}
