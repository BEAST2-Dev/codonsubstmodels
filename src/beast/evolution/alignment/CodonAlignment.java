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
import beast.core.util.Log;
import beast.evolution.datatype.Codon;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.GeneticCode;
import beast.evolution.datatype.Nucleotide;
import beast.util.AddOnManager;

import java.text.DecimalFormat;
import java.util.Arrays;
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
public class CodonAlignment extends Alignment { //TODO should have WrappedAlignment between CodonAlignment and Alignment

    final public Input<Alignment> alignmentInput = new Input<>("data",
            "Nucleotide alignment to convert into codon data type specified by dataType", Input.Validate.REQUIRED);

    final public Input<String> geneticCodeInput = new Input<>("geneticCode",
            "The rule to define how sequences of nucleotide triplets, " +
                    "called codons, specify which amino acid will be added next during protein synthesis.",
                    GeneticCode.GENETIC_CODE_NAMES[GeneticCode.UNIVERSAL_ID]);

    protected Alignment alignment;
//    protected GeneticCode geneticCode;

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
            this.setID("CA." + alignment.getID());
        DataType alignmentType = alignment.getDataType();
        if (! (alignmentType instanceof Nucleotide) ) {
            throw new IllegalArgumentException("CodonAlignment currently only support to wrap the nucleotide alignment !");
        }

        alignmentInput.setValue(alignment, this);
        geneticCodeInput.setValue(geneticCode.getName(), this);

        if (this.m_dataType == null || !this.m_dataType.getTypeDescription().equals(Codon.CODON)) {
            DataType oldType = this.m_dataType;
            this.m_dataType = new Codon(geneticCode);
            Log.warning.println("Warning: CodonAlignment (" + this.getID() + ") original data type was " +
                    oldType + ", is corrected to " + this.m_dataType + " - " + geneticCode.getDescription() + " !");
        }
//        initAndValidate();
    }

    @Override
    public void initAndValidate() {
        alignment = alignmentInput.get();//Nucleotide
        DataType alignmentType = alignment.getDataType();
        GeneticCode geneticCode = GeneticCode.findByName(geneticCodeInput.get());

        initDataType(); //TODO need improve

        DataType thisType = this.m_dataType;
        // only working for nucleotide => codon
        if (thisType != null && thisType instanceof Codon && alignmentType instanceof Nucleotide) {
            m_dataType = thisType;
        } else {
            throw new IllegalArgumentException("CodonAlignment only wraps the nucleotide alignment into codon alignment !");
        }

        // set geneticCode to geneticCodeInput.get() if different
        GeneticCode geneticCode2 = ((Codon) thisType).getGeneticCode();
        if (! geneticCode.getName().equals( geneticCode2.getName() ) ) { //thisType instanceof Codon &&
            ((Codon) m_dataType).setGeneticCode(geneticCode);
        }

        convertCodonToState();

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

        Log.info.println("\nGenetic code is " + getGeneticCode().getDescription());

        int[][] usage = getCodonUsage();
        printCodonUsage(usage);

        printCodonPositionBaseFrequencies();

        double[] freqs = getCodonFrequenciesFromUsage(usage);
        printCodonFrequencies(freqs);
        Log.info.println();
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

    protected void convertCodonToState() {
        taxaNames.clear();
        stateCounts.clear();
        counts.clear();
        try {
            for (Sequence seq : alignment.sequences) {
                List<Integer> codonStates = seq.getSequence(getDataType()); // return mapCodeToStateSet indices i
                int tripletIndex = findStopCodon(codonStates);
                if (tripletIndex > -1)
                    Log.warning.println("Warning: " + seq.getTaxon() + " sequence contains a stop codon at " +
                            (tripletIndex+1) + "th triplets ! \n" +
                            "Please either use a codon alignment or the correct genetic code.");

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

    @Override
    public Codon getDataType() {
        if (!(m_dataType instanceof Codon))
            throw new UnsupportedOperationException("CodonAlignment only supports Codon data type !");
        return (Codon) m_dataType;
    }

    public GeneticCode getGeneticCode() {
        return getDataType().getGeneticCode();
    }

    /**
     * @return number of sites
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
     * Use {@link Alignment#taxaNames taxaNames} as row indices, and
     * {@link GeneticCode#GENETIC_CODE_TABLES GENETIC_CODE_TABLES} as column indices.
     * @return
     */
    protected int[][] getCodonUsage() {
        if (taxaNames.size() != counts.size())
            throw new IllegalArgumentException("taxaNames.size() " + taxaNames.size() + " != counts.size() " + counts.size());

        GeneticCode geneticCode = getGeneticCode();
        String code = geneticCode.getCodeTable();
        int[][] usage = new int[taxaNames.size()][code.length()];
        for (int i = 0; i < counts.size(); i++) {
            List<Integer> codonStates = counts.get(i);
            for (int j = 0; j < codonStates.size(); j++) {
                int state = codonStates.get(j);
                usage[i][state] += 1;
            }
        }
        return usage;
    }

    /**
     * Codon position * base (3x4) table, plus "overall" in last row.
     * The base order is "A", "C", "G", "T".
     * @return
     */
    public double[][] getCodonPositionBaseFrequencies() {
        GeneticCode geneticCode = getGeneticCode();
        // position x base (3x4) table + overall
        double[][] freqs = new double[4][4];
        for (int i = 0; i < counts.size(); i++) {
            List<Integer> codonStates = counts.get(i);
            for (int j = 0; j < codonStates.size(); j++) {
                int state = codonStates.get(j);
                String triplet = getDataType().state2string(new int[]{state});
                // position
                for (int pos = 0; pos < 3; pos++) {
                    // col index = nucState: A,C,G,T
                    int nucState = geneticCode.getNucleotideState(triplet.charAt(pos));
                    freqs[pos][nucState] += 1;
                }
            }
        }
        // overall
        for (int col = 0; col < 4; col++) {
            freqs[3][col] = freqs[0][col] + freqs[1][col] + freqs[2][col];
        }
        // compute frequencies
        for (int row = 0; row < 4; row++) {
            double rowSum = freqs[row][0];
            for (int col = 1; col < 4; col++)
                rowSum += freqs[row][col];
            for (int col = 0; col < 4; col++)
                freqs[row][col] = freqs[row][col] / rowSum;
        }
        return freqs;
    }

    protected double[] getCodonFrequenciesFromUsage(int[][] usage) {
        // usage last row is total
        double[] freqs = new double[usage[usage.length-1].length];
        double sum = 0;
        for(int i=0; i<freqs.length; i++) {
            freqs[i] = usage[usage.length-1][i];
            sum += freqs[i];
        }
        if (sum == 0)
            throw new IllegalArgumentException("Invalid codon usage, the total is 0 !");
        for(int i=0; i<freqs.length; i++)
            freqs[i] = freqs[i] / sum;
        return freqs;
    }

    /**
     * Codon frequencies from codon usage (AAA AAC AAG AAT ... TTT)
     * @return
     */
    public double[] getCodonFrequencies() {
        int[][] usage = getCodonUsage();
        return getCodonFrequenciesFromUsage(usage);
    }

    //============ print ============

    /**
     * Codon usage in sequences
     */
    protected void printCodonUsage(int[][] usage) {
        String codeTable = getGeneticCode().getCodeTable();
        List<String> taxaNames = getTaxaNames();

        Log.info.println("\n============ Codon Usage ============");
        // header 1st cell to fill in spaces
        String firstTN = taxaNames.get(0);
        String spaceN = new String(new char[firstTN.length()+1]).replace('\0', ' ');

        // header triplets
        Log.info.print(spaceN);
        for (int j = 0; j < codeTable.length(); j++)
            Log.info.print("\t" + getDataType().state2string(new int[]{j}));
        Log.info.println();

        // header AminoAcid
        Log.info.print(spaceN);
        for (int j = 0; j < codeTable.length(); j++)
            Log.info.print("\t" + codeTable.charAt(j));
        Log.info.println();

        // Codon Usage
        int[] colSums = new int[codeTable.length()];
        for (int i = 0; i < taxaNames.size(); i++) {
            Log.info.print(taxaNames.get(i));

            for (int j = 0; j < codeTable.length(); j++) {
                colSums[j] += usage[i][j];
                Log.info.print("\t" + usage[i][j]);
            }
            Log.info.println();
        }
//        Log.info.println();

        // overall
        Log.info.print("overall");
        for (int j = 0; j < codeTable.length(); j++)
            Log.info.print("\t" + colSums[j]);
        Log.info.println();
    }


    protected void printCodonPositionBaseFrequencies() {
        String[] rowNames = new String[]{"position 1 : ", "position 2 : ", "position 3 : ", "overall : "};
        String[] colNames = new String[]{"A", "C", "G", "T"};
        double[][] freqs = getCodonPositionBaseFrequencies();
        DecimalFormat df = new DecimalFormat("#");
        df.setMaximumFractionDigits(5);// 5 decimal places

        Log.info.println("\n============ Codon position * base (3x4) table + overall ============");
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
//        Log.info.println();
    }

    protected void printCodonFrequencies(double[] frequencies) {
        Log.info.println("\n============ Codon frequencies from usage (AAA AAC AAG AAT ... TTT) ============");
        DecimalFormat df = new DecimalFormat("#");
        df.setMaximumFractionDigits(8);
        for (int i = 0; i < frequencies.length; i++) {
            int state = getDataType().getStatesForCode(i)[0];
            if (i % 8 == 0) {
                Log.info.print("\n" + df.format(frequencies[state]));
            } else {
                Log.info.print("\t" + df.format(frequencies[state]));
            }
        }
        Log.info.println();
    }

}
