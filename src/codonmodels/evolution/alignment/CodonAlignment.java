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


import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.Nucleotide;
import codonmodels.evolution.datatype.Codon;
import codonmodels.evolution.datatype.GeneticCode;

import java.text.DecimalFormat;
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
public class CodonAlignment extends Alignment {

    final public Input<Alignment> alignmentInput = new Input<>("data",
            "Nucleotide alignment to convert into codon data type specified by dataType", Input.Validate.REQUIRED);

    final public Input<String> geneticCodeInput = new Input<>("geneticCode",
            "The rule to define how sequences of nucleotide triplets, " +
                    "called codons, specify which amino acid will be added next during protein synthesis.",
                    GeneticCode.GENETIC_CODE_NAMES[GeneticCode.UNIVERSAL_ID]);

    final public Input<Boolean> verboseInput = new Input<>("verbose",
            "Print the codon usage, etc.",
            Boolean.TRUE);

//    final public Input<Boolean> unknownCodeExceptionInput = new Input<>("unknownCodeException",
//            "Flag to throw exception for unknown code, but can be changed to 'false' to treat " +
//                    "partial ambiguities (-TA) as missing data.",
//            true, Input.Validate.OPTIONAL);

    protected Alignment alignment;

    public CodonAlignment() {
        sequenceInput.setRule(Input.Validate.OPTIONAL);
        userDataTypeInput.setRule(Input.Validate.FORBIDDEN); // avoid confusion
        siteWeightsInput.setRule(Input.Validate.FORBIDDEN); // ?
        // set default to Codon
        dataTypeInput.setValue(Codon.CODON, this);
    }

    public CodonAlignment(Alignment alignment, GeneticCode geneticCode) {
        this();
        if (alignment.getID() != null)
            this.setID("ca." + alignment.getID()); // make CodonAlignment ID not null
        DataType alignmentType = alignment.getDataType();
        if (! (alignmentType instanceof Nucleotide) ) {
            throw new IllegalArgumentException("CodonAlignment currently only support to wrap the nucleotide alignment !");
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
        if (thisType != null && thisType instanceof Codon && alignmentType instanceof Nucleotide) {
            m_dataType = thisType;
        } else {
            throw new IllegalArgumentException("CodonAlignment only wraps the nucleotide alignment into codon alignment !");
        }

        // set geneticCode to geneticCodeInput.get() if different to dataType.getGeneticCode()
        GeneticCode geneticCode = GeneticCode.findByName(geneticCodeInput.get());
        setGeneticCode(geneticCode);

//        convertCodonToState(unknownCodeExceptionInput.get()); // default to true
        if (counts.size() < 1)
            convertCodonToState();

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


        if (verboseInput.get()) {
            int[][] usage = getCodonUsage();
            printCodonUsage(usage);

            printCodonPositionBaseFrequencies();
            Log.info.println();
        }

        if (isAscertained)
            throw new UnsupportedOperationException("Ascertainment correction is not available !");
    }

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
                codonStates = seq.getSequence(getDataType());
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

    // if unknownCodeException is true, treat codons with partial ambiguities (-TA) as missing data;
    // if false, then use original code sequence.getSequence(dataType).
//    protected List<Integer> getCodonStates(Sequence sequence, DataType.Base dataType, boolean unknownCodeException) {
//        List<Integer> codonStates = new ArrayList<>();
//        if (unknownCodeException) {
//            // throw IllegalArgumentException, if codon has partial ambiguities (-TA)
//            codonStates = sequence.getSequence(dataType);
//        } else {
//            // allow partial ambiguities (-TA)
//            String data = sequence.getData();
//            // remove spaces
//            data = data.replaceAll("\\s", "");
//            data = data.toUpperCase();
//
//            // overwrite DataType.Base.stringToEncoding(data)
//            if (dataType.getCodeLength() == 3 || dataType.getCodeMap()==null) {
//                // use code map to resolve state codes
//                Map<String, Integer> map = new HashMap<>();
//                // fixed length code
//                for (int i = 0; i < dataType.getCodeMap().length(); i += dataType.getCodeLength()) {
//                    String code = dataType.getCodeMap().substring(i, i + dataType.getCodeLength());
//                    map.put(code, i / dataType.getCodeLength());
//                }
//
//                int unkn = 0;
//                for (int i = 0; i < data.length(); i += dataType.getCodeLength()) {
//                    String code = data.substring(i, i + dataType.getCodeLength()).toUpperCase();
//                    if (map.containsKey(code)) {
//                        codonStates.add(map.get(code));
//                    } else {
//                        // replace any unmapped triplets to ---
//                        codonStates.add(map.get("---"));
//                        unkn++;
//                    }
//                }
//                if (unkn > 0)
//                    Log.info.println("Replace " + unkn + " unmapped triplets to missing data in sequence " +
//                        sequence.getTaxon() + ".");
//            } else {
//                throw new IllegalArgumentException("Invalid data type !\n" +
//                        "Codon data type is required : " + dataType.getTypeDescription());
//            }
//        }
//        return codonStates;
//    }

//    @Deprecated
//    protected int findStopCodon(List<Integer> seqStates) {
//        for (int i = 0; i < seqStates.size(); i++) {
//            int codonState = seqStates.get(i);
//            if (getGeneticCode().isStopCodonIndex(codonState))
//                return i;
//        }
//        return -1;
//    }

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

        DataType originalType = alignment.getDataType();
        int count = alignment.getSiteCount();

        if (originalType instanceof Nucleotide) {
            count /= 3;
        }
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

    //============ Util ============

    /**
     * Util to cast Alignment to CodonAlignment
     * @param alignment
     * @return
     */
    public static CodonAlignment toCodonAlignment(Alignment alignment) {
        if (! (alignment instanceof CodonAlignment) )
            throw new IllegalArgumentException("Cannot parse to CodonAlignment class !");
        return (CodonAlignment) alignment;
    }

    //============ stats ============

    /**
     * The usage (counts) table of triplets
     * Use {@link Alignment#taxaNames taxaNames} as row indices, and
     * {@link Codon} states as column indices.
     * The current columns include ambiguous states.
     * @return matrix int[taxaNames.size()][getStateCountAmbiguous()]
     */
    public int[][] getCodonUsage() {
        if (taxaNames.size() != counts.size())
            throw new IllegalArgumentException("taxaNames.size() " + taxaNames.size() + " != counts.size() " + counts.size());

        // include --- ???
        int[][] usage = new int[taxaNames.size()][getDataType().getStateCountAmbiguous()];
        for (int i = 0; i < counts.size(); i++) {
            List<Integer> codonStates = counts.get(i);
            for (int j = 0; j < codonStates.size(); j++) {
                // states as column indices
                int state = codonStates.get(j);
                usage[i][state] += 1;
            }
        }
        return usage;
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

    //============ print ============

    /**
     * Codon usage in sequences
     */
    protected void printCodonUsage(int[][] usage) {

        Log.info.println("\n============ Codon usage in sequences ============");

        List<String> taxaNames = getTaxaNames();
        // header 1st cell to fill in spaces
        String firstTN = taxaNames.get(0);
        String spaceN = new String(new char[firstTN.length()+1]).replace('\0', ' ');

        final int stateCount = getDataType().getStateCount(); // 60/61

        // ambiguous check
        int ambiguous = 0;
        for (int i = 0; i < taxaNames.size(); i++) {
            for (int state = 0; state < usage[0].length; state++) {
                if (state >= stateCount)
                    ambiguous += usage[i][state];
            }
        }

        int colMax = getGeneticCode().getCodeTableLength();
        // not print ambiguous if no ambiguous
        if (ambiguous > 0)
            colMax += usage[0].length - stateCount;

        // header triplets
        Log.info.print(spaceN);
        // j is index not state
        for (int j = 0; j < colMax; j++)
            Log.info.print("\t" + getDataType().getTriplet(j));
        Log.info.println();

        // header AminoAcid
        Log.info.print(spaceN);
        // j is index not state
        for (int j = 0; j < colMax; j++)
            Log.info.print("\t" + getGeneticCode().getAminoAcidFromCodeTable(j));
        Log.info.println();

        // Codon Usage
        int[] colSums = new int[colMax];
        for (int i = 0; i < taxaNames.size(); i++) {
            Log.info.print(taxaNames.get(i));

            int state = 0;
            // j is index not state
            for (int j = 0; j < colMax; j++) {
                if (getGeneticCode().isStopCodonIndex(j)) {
                    Log.info.print("\t0");
                } else {
                    // usage[][state] not index
                    Log.info.print("\t" + usage[i][state]);
                    colSums[j] += usage[i][state];
                    state++;
                }
            }
            Log.info.println();
        }
//        Log.info.println();

        // overall
        Log.info.print("overall");
        for (int j = 0; j < colMax; j++)
            Log.info.print("\t" + colSums[j]);
        Log.info.println("\n");

        // print warning for ambiguous
        if (ambiguous > 0) {
            Log.info.println("Warning : find " + ambiguous + " ambiguous triplets in this alignment.");
            Log.info.println("Please be aware that the code should be non-ambiguous, " +
                    "a gap may be introduced by the sequencing error or bad alignment.");
        }

    }


    protected void printCodonPositionBaseFrequencies() {
        String[] rowNames = new String[]{"position 1 : ", "position 2 : ", "position 3 : ", "average : "};
        String[] colNames = new String[]{"A", "C", "G", "T"};
        double[][] freqs = getObservedBaseFrequencies();
        DecimalFormat df = new DecimalFormat("#");
        df.setMaximumFractionDigits(5);// 5 decimal places

        Log.info.println("\n============ Codon position * base (3x4) table + average ============");
        // header 1st cell to fill in spaces
        String firstTN = rowNames[0];
        String spaceN = new String(new char[firstTN.length()+1]).replace('\0', ' ');

        // header
        Log.info.print(spaceN);
        for (int j = 0; j < colNames.length; j++)
            Log.info.print("\t" + colNames[j]);
        Log.info.println();

        // freqs
        for (int i = 0; i < rowNames.length; i++) {
            Log.info.print(rowNames[i]);

            for (int j = 0; j < colNames.length; j++) {
                Log.info.print("\t" + df.format(freqs[i][j]));
            }
            Log.info.println();
        }
        Log.info.println();
    }


}
