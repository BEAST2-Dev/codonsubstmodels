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
import beast.evolution.datatype.*;

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
            DataType oldType = this.m_dataType;
            this.m_dataType = new Codon(geneticCode);
            Log.warning.println("Warning: CodonAlignment (" + this.getID() + ") original data type was " +
                    oldType + ", is corrected to " + this.m_dataType + " - " + geneticCode.getDescription() + " !");
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
        sanityCheckCalcPatternsSetUpAscertainment(true);

        Log.info.println("\nGenetic code is " + getGeneticCode().getDescription());

        int[][] usage = getCodonUsage();
        printCodonUsage(usage);

        printCodonPositionBaseFrequencies();

        double[] freqs = getCodonFrequenciesByUsage(usage);
        printCodonFrequencies(freqs, "Codon frequencies by usage");
        Log.info.println();
    }

    //TODO need improve
    /**
     * modified from Alignment private initializeWithSequenceList(List<Sequence>, boolean)
     */
    protected void convertCodonToState() {
        taxaNames.clear();
        stateCounts.clear();
        counts.clear();
        try { //TODO treat codons with partial ambiguities (-TA) as missing data
            for (Sequence seq : alignment.sequences) {
                // return mapCodeToStateSet indices i, also indices in Codon.CODON_TRIPLETS
                List<Integer> tripletStates = seq.getSequence(getDataType());
                int stopCodon = findStopCodon(tripletStates);
                if (stopCodon > -1)
                    Log.warning.println("Warning: " + seq.getTaxon() + " sequence contains a stop codon at " +
                            (stopCodon+1) + "th triplets ! \n" +
                            "Please either use a codon alignment or the correct genetic code.");

                counts.add(tripletStates);
                if (taxaNames.contains(seq.getTaxon())) {
                    throw new RuntimeException("Duplicate taxon found in alignment: " + seq.getTaxon());
                }
                taxaNames.add(seq.getTaxon());

                if (seq.uncertain)
                    throw new UnsupportedOperationException("Uncertain sequence is not available in this version !");

                    //TODO how to deal with uncertain? seq here is Nucleotide sequences not codons.
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

    public void setGeneticCode(GeneticCode geneticCode) {
        GeneticCode geneticCode2 = getGeneticCode();
        if (! geneticCode.getName().equals( geneticCode2.getName() ) ) {
            ((Codon) m_dataType).setGeneticCode(geneticCode);
            geneticCodeInput.setValue(geneticCode.getName(), this);
            Log.info.println("Change genetic code from " + geneticCode2.getName() + " to " + geneticCode.getName());
        }
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
     * {@link Codon#CODON_TRIPLETS CODON_TRIPLETS} as column indices.
     * @return matrix int[taxaNames.size()][getDataType().ambiguousStateCount]
     */
    protected int[][] getCodonUsage() {
        if (taxaNames.size() != counts.size())
            throw new IllegalArgumentException("taxaNames.size() " + taxaNames.size() + " != counts.size() " + counts.size());

//        GeneticCode geneticCode = getGeneticCode();
//        String codeTable = geneticCode.getCodeTable();
//        int[][] usage = new int[taxaNames.size()][codeTable.length()];
        // include --- ???
        int[][] usage = new int[taxaNames.size()][getDataType().ambiguousStateCount];
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
     * If any ambiguous, then add 0.25 to each Nucleotide count.
     * @return
     */
    public double[][] getCodonPositionBaseFrequencies() {
        GeneticCode geneticCode = getGeneticCode();
        // position x base (3x4) table + overall
        int colMax = 4;
        double[][] freqs = new double[4][colMax];
        for (int i = 0; i < counts.size(); i++) {
            List<Integer> codonStates = counts.get(i);
            for (int j = 0; j < codonStates.size(); j++) {
                int state = codonStates.get(j);
                String triplet = getDataType().encodingToString(new int[]{state});
                // position
                for (int pos = 0; pos < 3; pos++) {
                    // col index = nucState: A,C,G,T,-
                    int nucState = geneticCode.getNucleotideState(triplet.charAt(pos));
                    // treat any ambiguous as a gap, add 1/4 to each freqs
                    if (nucState >= 4) {
                        for (int n = 0; n < 3; n++)
                            freqs[pos][n] += 0.25;
                    } else {
                        freqs[pos][nucState] += 1;
                    }
                }
            }
        }
        // overall
        for (int col = 0; col < colMax; col++) {
            freqs[3][col] = freqs[0][col] + freqs[1][col] + freqs[2][col];
        }
        // compute frequencies
        for (int row = 0; row < 4; row++) {
            double rowSum = freqs[row][0];
            for (int col = 1; col < colMax; col++)
                rowSum += freqs[row][col];
            for (int col = 0; col < colMax; col++)
                freqs[row][col] = freqs[row][col] / rowSum;
        }
        return freqs;
    }

    // int[][] usage has no total, and usage col has 2 ambiguous states count
    protected double[] getCodonFrequenciesByUsage(int[][] usage) {
        int stateMax = getDataType().getStateCount(); // 64

        double ambiguous = 0;
        for (int j = stateMax; j < usage[0].length; j++) {
            for (int i = 0; i < usage.length; i++) {
                ambiguous += usage[i][j];
            }
        }
//        if (ambiguous > 0)
//            Log.info.println("Find " + ambiguous + " ambiguous states in this alignment.");

        // sum up counts
        double[] freqs = new double[stateMax];
        double sum = 0;
        int actualStateCount = 0;
        for (int j = 0; j < stateMax; j++) {
            for (int i = 0; i < usage.length; i++) {
                freqs[j] += usage[i][j];
                sum += usage[i][j];
            }
            if (freqs[j] > 0)
                actualStateCount += 1;
        }
        if (sum == 0)
            throw new IllegalArgumentException("Invalid codon usage, the total is 0 !");

        // equally distribute ambiguous into the count of each of actual state
        if (ambiguous > 0) {
            for (int j = 0; j < stateMax; j++) {
                // freqs[j] == 0 is stop codon
                if (freqs[j] > 0)
                    freqs[j] += ambiguous / actualStateCount;
            }
            sum += ambiguous;
        }

        // calculate freqs
        for (int j = 0; j < stateMax; j++) {
            freqs[j] = freqs[j] / sum;
        }

        return freqs;
    }

    /**
     * Codon frequencies from codon usage (AAA AAC AAG AAT ... TTT)
     * @return
     */
    public double[] getCodonFrequencies() {
        int[][] usage = getCodonUsage();
        return getCodonFrequenciesByUsage(usage);
    }

    //============ print ============

    /**
     * Codon usage in sequences
     */
    protected void printCodonUsage(int[][] usage) {
//        String codeTable = getGeneticCode().getCodeTable();
//        int colMax = codeTable.length(); // no - ?
        int colMax = getDataType().ambiguousStateCount;
        List<String> taxaNames = getTaxaNames();

        Log.info.println("\n============ Codon Usage ============");
        // header 1st cell to fill in spaces
        String firstTN = taxaNames.get(0);
        String spaceN = new String(new char[firstTN.length()+1]).replace('\0', ' ');

        // header triplets
        Log.info.print(spaceN);
        for (int j = 0; j < colMax; j++)
            Log.info.print("\t" + getDataType().encodingToString(new int[]{j}));
        Log.info.println();

        // header AminoAcid
        Log.info.print(spaceN);
        for (int j = 0; j < colMax; j++)
            Log.info.print("\t" + getGeneticCode().getAminoAcid(j));
        Log.info.println();

        // Codon Usage
        int[] colSums = new int[colMax];
        for (int i = 0; i < taxaNames.size(); i++) {
            Log.info.print(taxaNames.get(i));

            for (int j = 0; j < colMax; j++) {
                colSums[j] += usage[i][j];
                Log.info.print("\t" + usage[i][j]);
            }
            Log.info.println();
        }
//        Log.info.println();

        // overall
        Log.info.print("overall");
        for (int j = 0; j < colMax; j++)
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

    public void printCodonFrequencies(double[] frequencies, String title) {
        Log.info.println("\n============ " + title + " (AAA AAC AAG AAT ... TTA TTC TTG TTT) ============");
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
