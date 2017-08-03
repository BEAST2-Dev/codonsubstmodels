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
import util.StringUtils;

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


    @Override
    public void initAndValidate() {
        GeneticCode geneticCode = GeneticCode.findByName(geneticCodeInput.get());
        initStats(geneticCode);
    }

    //TODO move DataType
    protected int ambiguousStateCount;

    // need this constructor to pass Input validation
    public Codon() {
        this(GeneticCode.UNIVERSAL);
    }

    public Codon(GeneticCode geneticCode) {
        initStats(geneticCode);
    }

    private void initStats(GeneticCode geneticCode) {
        this.geneticCode = geneticCode;

        stateCount = 64 - geneticCode.getStopCodonCount();
        codeLength = 3;
        codeMap = StringUtils.concatenateToString(CODON_TRIPLETS);

        ambiguousStateCount = 66;
        mapCodeToStateSet = new int[ambiguousStateCount][];

        int j = 0;
        int k = stateCount;
        for (int i = 0; i < 64; i++) {
            if (!geneticCode.isStopCodon(i)) {
                mapCodeToStateSet[i] = new int[]{j};
                j++;
            } else {
                mapCodeToStateSet[i] = new int[]{k};
                k++;
            }

        }
        int[] all = new int[64];
        for (int i = 0; i < 64; i++) {
            all[i] = i;
        }
        for (int i = 64; i < ambiguousStateCount; i++) {
            mapCodeToStateSet[i] = all;
        }

        // TODO stateMap, reverseMap => mapCodeToStateSet

        stateMap = new int[ambiguousStateCount];
        reverseMap = new int[ambiguousStateCount];

        j = 0;
        k = stateCount;
        for (int i = 0; i < 64; i++) {
            if (!geneticCode.isStopCodon(i)) {
                stateMap[j] = i;
                reverseMap[i] = j;
                j++;
            } else {
                stateMap[k] = i;
                reverseMap[i] = k;
                k++;
            }
        }
        for (int i = 64; i < ambiguousStateCount; i++) {
            stateMap[i] = i;
            reverseMap[i] = i;
        }
    }

    @Override
    public String getTypeDescription() {
        return "codon";
    }

    //=========== for codons ==========

    /**
     * Get state corresponding to a nucleotide triplet
     *
     * @param nuc1 the codon triplet as chars
     * @param nuc2 the codon triplet as chars
     * @param nuc3 the codon triplet as chars
     * @return state
     */
    public final int getState(char nuc1, char nuc2, char nuc3) {
        return getState(geneticCode.getNucleotideState(nuc1),
                geneticCode.getNucleotideState(nuc2),
                geneticCode.getNucleotideState(nuc3));
    }

    /**
     * @return the canonical state (in standard combinatorial order)
     *         of a funny codon state.
     */
    public final int getCanonicalState(int funnyState) {
        if (funnyState >= stateMap.length)
            throw new IllegalArgumentException("Invalid state >= length of stateMap ! ");
        return stateMap[funnyState];
    }

//    public static final int NUCLEOTIDE_UNKNOWN_STATE = 16;
    public static final int NUCLEOTIDE_GAP_STATE = 17;

    /**
     * Get state corresponding to a nucleotide triplet
     *
     * @param nuc1 the codon triplet as states
     * @param nuc2 the codon triplet as states
     * @param nuc3 the codon triplet as states
     * @return state
     */
    public final int getState(int nuc1, int nuc2, int nuc3) {
        if (nuc1 == NUCLEOTIDE_GAP_STATE || nuc2 == NUCLEOTIDE_GAP_STATE ||
                nuc3 == NUCLEOTIDE_GAP_STATE ) {
            return GAP_STATE;
        }

        if (isAmbiguousState(nuc1) || isAmbiguousState(nuc2) || isAmbiguousState(nuc3)) {
            return UNKNOWN_STATE;
        }

        int canonicalState = (nuc1 * 16) + (nuc2 * 4) + nuc3;

        return reverseMap[canonicalState];
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
     *
     * @param state state
     *              <p/>
     *              return corresponding triplet string
     */
    public final String getTriplet(int state) {
        return CODON_TRIPLETS[stateMap[state]];
    }

    /**
     * Get an array of three nucleotide states making this codon state
     *
     * @param state state
     *              <p/>
     *              return corresponding triplet string
     */
    public final int[] getTripletStates(int state) {
        int[] triplet = new int[3];

        triplet[0] = geneticCode.getNucleotideState(CODON_TRIPLETS[stateMap[state]].charAt(0));
        triplet[1] = geneticCode.getNucleotideState(CODON_TRIPLETS[stateMap[state]].charAt(1));
        triplet[2] = geneticCode.getNucleotideState(CODON_TRIPLETS[stateMap[state]].charAt(2));

        return triplet;
    }


    /**
     * @return the genetic code
     */
    public GeneticCode getGeneticCode() {
        return geneticCode;
    }

    /**
     * Takes non-canonical state and returns true if it represents stop codon
     *
     * @param state
     * @return true if the given state is a stop codon
     */
    public final boolean isStopCodon(int state) {
        return geneticCode.isStopCodon(stateMap[state]);
    }

    public static byte[] constructRateMap(Codon codonDataType) {
        final int stateCount = codonDataType.getStateCount();
        final int rateCount = (stateCount * (stateCount - 1)) / 2;
        return constructRateMap(rateCount, stateCount, codonDataType, codonDataType.getGeneticCode());
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

    /**
	 * Construct a map of the rate classes in the rate matrix using the current
	 * genetic code. Classes:
	 *		0: codon changes in more than one codon position (or stop codons)
	 *		1: synonymous transition
	 *		2: synonymous transversion
	 *		3: non-synonymous transition
	 *		4: non-synonymous transversion
	 */
	public static byte[] constructRateMap(int rateCount, int stateCount, Codon codonDataType, GeneticCode geneticCode)	{
		int u, v, i1, j1, k1, i2, j2, k2;
		byte rateClass;
		int[] codon;
		int cs1, cs2, aa1, aa2;

		int i = 0;

		byte[] rateMap = new byte[rateCount];

		for (u = 0; u < stateCount; u++) {

			codon = codonDataType.getTripletStates(u);
			i1 = codon[0];
			j1 = codon[1];
			k1 = codon[2];

			cs1 = codonDataType.getState(i1, j1, k1);
			aa1 = geneticCode.getAminoAcidState(codonDataType.getCanonicalState(cs1));

			for (v = u + 1; v < stateCount; v++) {

				codon = codonDataType.getTripletStates(v);
				i2 = codon[0];
				j2 = codon[1];
				k2 = codon[2];

				cs2 = codonDataType.getState(i2, j2, k2);
				aa2 = geneticCode.getAminoAcidState(codonDataType.getCanonicalState(cs2));

				rateClass = -1;
				if (i1 != i2) {
					if ( (i1 == 0 && i2 == 2) || (i1 == 2 && i2 == 0) || // A <-> G
						 (i1 == 1 && i2 == 3) || (i1 == 3 && i2 == 1) ) { // C <-> T
						rateClass = 1; // Transition at position 1
					} else {
						rateClass = 2; // Transversion at position 1
					}
				}
				if (j1 != j2) {
					if (rateClass == -1) {
						if ( (j1 == 0 && j2 == 2) || (j1 == 2 && j2 == 0) || // A <-> G
							 (j1 == 1 && j2 == 3) || (j1 == 3 && j2 == 1) ) { // C <-> T
							rateClass = 1; // Transition
						} else {
							rateClass = 2; // Transversion
						}
					} else
						rateClass = 0; // Codon changes at more than one position
				}
				if (k1 != k2) {
					if (rateClass == -1) {
						if ( (k1 == 0 && k2 == 2) || (k1 == 2 && k2 == 0) || // A <-> G
							 (k1 == 1 && k2 == 3) || (k1 == 3 && k2 == 1) ) { // C <-> T
							rateClass = 1; // Transition
						} else {
							rateClass = 2; // Transversion
						}
					} else
						rateClass = 0; // Codon changes at more than one position
				}

	 			if (rateClass != 0) {
					if (aa1 != aa2) {
						rateClass += 2; // Is a non-synonymous change
					}
				}

				rateMap[i] = rateClass;
				i++;
			}

		}
        return rateMap;
	}

    // Private members

    protected GeneticCode geneticCode;
    protected int[] stateMap;
    protected int[] reverseMap;
}
