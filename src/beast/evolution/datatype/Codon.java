/*
 * Codons.java
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

package beast.evolution.datatype;

import beast.core.Input;
import beast.util.StringUtils;

/**
 * Implements DataType for codons
 * <p/>
 * Codon have tree different representations:
 * State numbers - 0-63 + 64, 65 as unknown and gap
 * State chars - the above converted into chars, starting at 'A'
 * and '?' + '-' for unknown and gap
 * Strings or chars of three nucleotide characters
 *
 * Modified from BEAST 1 Codon.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class Codon extends DataType.Base {

    final public Input<String> geneticCodeInput = new Input<>("geneticCode",
            "The rule to define how sequences of nucleotide triplets, " +
                    "called codons, specify which amino acid will be added next during protein synthesis.",
            GeneticCode.GENETIC_CODE_NAMES[GeneticCode.UNIVERSAL_ID], Input.Validate.REQUIRED);


//    public static final Codon UNIVERSAL = new Codon(GeneticCode.UNIVERSAL);
//    public static final Codon VERTEBRATE_MT = new Codon(GeneticCode.VERTEBRATE_MT);
//    public static final Codon YEAST = new Codon(GeneticCode.YEAST);
//    public static final Codon MOLD_PROTOZOAN_MT = new Codon(GeneticCode.MOLD_PROTOZOAN_MT);
//    public static final Codon MYCOPLASMA = new Codon(GeneticCode.MYCOPLASMA);
//    public static final Codon INVERTEBRATE_MT = new Codon(GeneticCode.INVERTEBRATE_MT);
//    public static final Codon CILIATE = new Codon(GeneticCode.CILIATE);
//    public static final Codon ECHINODERM_MT = new Codon(GeneticCode.ECHINODERM_MT);
//    public static final Codon EUPLOTID_NUC = new Codon(GeneticCode.EUPLOTID_NUC);
//    public static final Codon BACTERIAL = new Codon(GeneticCode.BACTERIAL);
//    public static final Codon ALT_YEAST = new Codon(GeneticCode.ALT_YEAST);
//    public static final Codon ASCIDIAN_MT = new Codon(GeneticCode.ASCIDIAN_MT);
//    public static final Codon FLATWORM_MT = new Codon(GeneticCode.FLATWORM_MT);
//    public static final Codon BLEPHARISMA_NUC = new Codon(GeneticCode.BLEPHARISMA_NUC);
//    public static final Codon NO_STOPS = new Codon(GeneticCode.NO_STOPS);

    public static final int UNKNOWN_STATE = 64;
    public static final int GAP_STATE = 65;

    public static final String[] CODON_TRIPLETS = {
            "AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT",
            "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT",
            "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT",
            "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT",
            "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT",
            "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT",
            "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT",
            "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT",
            "???", "---"
    };

    /**
     * This character represents the amino acid equivalent of a stop codon to cater for
     * situations arising from converting coding DNA to an amino acid sequence.
     */
    public static final char STOP_CHARACTER = '*';

    public static final int STOP_STATE = 23;

    public final static String CODON = "codon";

    @Override
    public void initAndValidate() {
        GeneticCode geneticCode = GeneticCode.findByName(geneticCodeInput.get());
        setGeneticCode(geneticCode);
    }

    protected GeneticCode geneticCode;
    /**
     * stateMap saves the states 0-63, but keeps the states of stop codon in the last.
     * It is different to states from {@link DataType#string2state(String) string2state}
     * using codeMap
     */
//    protected int[] stateMap; //TODO merge stateMap to codeMap
//    protected int[] reverseMap;

    //TODO move to DataType
    protected int ambiguousStateCount;

    // this constructor used by xml not java
    public Codon() {
        this(GeneticCode.UNIVERSAL);
    }

    public Codon(GeneticCode geneticCode) {
        setGeneticCode(geneticCode);
    }

    // offer to change GeneticCode later, especially after default constructor called by CodonAlignment.initDataType()
    public void setGeneticCode(GeneticCode geneticCode) {
        this.geneticCode = geneticCode;

//        stateCount = 64 - geneticCode.getStopCodonCount();
        stateCount = 64;
        codeLength = 3;
        codeMap = StringUtils.concatenateToString(CODON_TRIPLETS);

        ambiguousStateCount = 66;
        mapCodeToStateSet = new int[ambiguousStateCount][];
        for (int i = 0; i < stateCount; i++) {
            mapCodeToStateSet[i] = new int[1];
            mapCodeToStateSet[i][0] = i;
        }
        int[] all = new int[stateCount];
        for (int i = 0; i < stateCount; i++) {
            all[i] = i;
        }
        mapCodeToStateSet[64] = all;
        mapCodeToStateSet[65] = all;


//        int j = 0;
//        int k = 64 - geneticCode.getStopCodonCount();
//        // i is CODON_TRIPLETS also codonTable index, and [stateCount, 64] are stop codons state numbers
//        for (int i = 0; i < 64; i++) {
//            if (!geneticCode.isStopCodon(i)) {
//                mapCodeToStateSet[i] = new int[]{j};
//                j++;
//            } else {
//                mapCodeToStateSet[i] = new int[]{k};
//                k++;
//            }
//
//        }
//        int[] all = new int[64];
//        for (int i = 0; i < 64; i++) {
//            all[i] = i;
//        }
//        for (int i = 64; i < ambiguousStateCount; i++) {
//            mapCodeToStateSet[i] = all;
//        }
//
//        // TODO stateMap, reverseMap => mapCodeToStateSet, where reverseMap mapCodeToStateSet are same
//
//        stateMap = new int[ambiguousStateCount];
//        reverseMap = new int[ambiguousStateCount];
//
//        j = 0;
//        k = 64 - geneticCode.getStopCodonCount();
//        for (int i = 0; i < 64; i++) {
//            if (!geneticCode.isStopCodon(i)) {
//                stateMap[j] = i;
//                reverseMap[i] = j;
//                j++;
//            } else {
//                stateMap[k] = i;
//                reverseMap[i] = k;
//                k++;
//            }
//        }
//        for (int i = 64; i < ambiguousStateCount; i++) {
//            stateMap[i] = i;
//            reverseMap[i] = i;
//        }
    }

    @Override
    public String getTypeDescription() {
        return CODON;
    }

    //=========== for codons ==========
    /**
     * index in {@link Nucleotide#codeMap codeMap}
     */
    public static final int NUCLEOTIDE_UNKNOWN_STATE = 18; // ?
    public static final int NUCLEOTIDE_GAP_STATE = 17; // -

    /**
     * Get state corresponding to a nucleotide triplet
     *
     * @param nuc1 the codon triplet as chars
     * @param nuc2 the codon triplet as chars
     * @param nuc3 the codon triplet as chars
     * @return state
     */
    public final int getCodonState(char nuc1, char nuc2, char nuc3) {
        char ns1 = geneticCode.getNucleotideChar(nuc1);
        char ns2 = geneticCode.getNucleotideChar(nuc2);
        char ns3 = geneticCode.getNucleotideChar(nuc3);

        if (ns1 == NUCLEOTIDE_GAP_STATE || ns2 == NUCLEOTIDE_GAP_STATE ||
                ns3 == NUCLEOTIDE_GAP_STATE )
            return GAP_STATE;
        if (isAmbiguousState(ns1) || isAmbiguousState(ns2) || isAmbiguousState(ns3))
            return UNKNOWN_STATE;

        int codonState = string2state("" + nuc1 + nuc2 + nuc3).get(0);

        return codonState;
    }

    /**
     * @return the canonical state (in standard combinatorial order)
     *         of a funny codon state.
     */
//    public final int getCanonicalState(int funnyState) {
//        if (funnyState >= stateMap.length)
//            throw new IllegalArgumentException("Invalid state >= length of stateMap ! ");
//        return stateMap[funnyState];
//    }

    /**
     * Get codon state corresponding to a nucleotide triplet
     *
     * @param ns1 the codon triplet as states
     * @param ns2 the codon triplet as states
     * @param ns3 the codon triplet as states
     * @return state
     */
    public final int getCodonState(int ns1, int ns2, int ns3) {
        if (ns1 == NUCLEOTIDE_GAP_STATE || ns2 == NUCLEOTIDE_GAP_STATE ||
                ns3 == NUCLEOTIDE_GAP_STATE )
            return GAP_STATE;
        if (isAmbiguousState(ns1) || isAmbiguousState(ns2) || isAmbiguousState(ns3))
            return UNKNOWN_STATE;

//        int canonicalState = (ns1 * 16) + (ns2 * 4) + ns3; // cannot use BEAST1 nice design
        char nuc1 = geneticCode.getNucleotideChar(ns1);
        char nuc2 = geneticCode.getNucleotideChar(ns2);
        char nuc3 = geneticCode.getNucleotideChar(ns3);

        int codonState = string2state("" + nuc1 + nuc2 + nuc3).get(0);

        return codonState;
    }

    /**
     * Get character corresponding to a given state
     *
     * @param state state
     *              <p/>
     *              return corresponding character
     */
    public final char getChar(int state) {
        throw new IllegalArgumentException("Codon datatype cannot be expressed as char");
    }

    /**
     * Get triplet string corresponding to a given state
     * indexed by {@link Codon#codeMap codeMap}
     *
     * @param state state
     *              <p/>
     *              return corresponding triplet string
     */
    public final String getTriplet(int state) {
//        return CODON_TRIPLETS[stateMap[state]]; // states from stateMap
        return state2string(new int[]{state}); // states from codeMap
    }

    /**
     * Get an array of three nucleotide states making this codon state,
     * where nucleotide states are indexed by {@link Nucleotide#codeMap codeMap}.
     *
     * @see GeneticCode#getNucleotideState
     * @param state state
     *              <p/>
     *              return corresponding triplet string
     */
    public final int[] getTripletStates(int state) {
        int[] triplet = new int[3];

        triplet[0] = geneticCode.getNucleotideState(getTriplet(state).charAt(0));
        triplet[1] = geneticCode.getNucleotideState(getTriplet(state).charAt(1));
        triplet[2] = geneticCode.getNucleotideState(getTriplet(state).charAt(2));

        return triplet;
    }


    /**
     * @return the genetic code
     */
    public GeneticCode getGeneticCode() {
        return geneticCode;
    }


    /**
     * Parse a text string to return a genetic code
     */
//    public static Codon findByName(String codeStr) {
//        Codon codon = null;
//        if (codeStr.equals(GeneticCode.UNIVERSAL.getName())) {
//            codon = Codon.UNIVERSAL;
//        } else if (codeStr.equals(GeneticCode.VERTEBRATE_MT.getName())) {
//            codon = Codon.VERTEBRATE_MT;
//        } else if (codeStr.equals(GeneticCode.YEAST.getName())) {
//            codon = Codon.YEAST;
//        } else if (codeStr.equals(GeneticCode.MOLD_PROTOZOAN_MT.getName())) {
//            codon = Codon.MOLD_PROTOZOAN_MT;
//        } else if (codeStr.equals(GeneticCode.INVERTEBRATE_MT.getName())) {
//            codon = Codon.INVERTEBRATE_MT;
//        } else if (codeStr.equals(GeneticCode.CILIATE.getName())) {
//            codon = Codon.CILIATE;
//        } else if (codeStr.equals(GeneticCode.ECHINODERM_MT.getName())) {
//            codon = Codon.ECHINODERM_MT;
//        } else if (codeStr.equals(GeneticCode.EUPLOTID_NUC.getName())) {
//            codon = Codon.EUPLOTID_NUC;
//        } else if (codeStr.equals(GeneticCode.BACTERIAL.getName())) {
//            codon = Codon.BACTERIAL;
//        } else if (codeStr.equals(GeneticCode.ALT_YEAST.getName())) {
//            codon = Codon.ALT_YEAST;
//        } else if (codeStr.equals(GeneticCode.ASCIDIAN_MT.getName())) {
//            codon = Codon.ASCIDIAN_MT;
//        } else if (codeStr.equals(GeneticCode.FLATWORM_MT.getName())) {
//            codon = Codon.FLATWORM_MT;
//        } else if (codeStr.equals(GeneticCode.BLEPHARISMA_NUC.getName())) {
//            codon = Codon.BLEPHARISMA_NUC;
//        } else if (codeStr.equals(GeneticCode.NO_STOPS.getName())) {
//            codon = Codon.NO_STOPS;
//        } else {
//            throw new RuntimeException("Unknown genetics code");
//        }
//        return codon;
//    }

}
