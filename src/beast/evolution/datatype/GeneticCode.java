/*
 * GeneticCode.java
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


import beast.core.Description;

/**
 * A set of standard genetic codes.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 *
 * Merged and modified from BEAST 1 GeneticCode, AminoAcids, Nucleotides.
 */
@Description("A set of standard genetic codes.")
public final class GeneticCode {

//    public static final String GENETIC_CODE = "geneticCode";
    
    /**
     * Constants used to refer to the built in code tables
     */
    public static final int UNIVERSAL_ID = 0;
    public static final int VERTEBRATE_MT_ID = 1;
    public static final int YEAST_ID = 2;
    public static final int MOLD_PROTOZOAN_MT_ID = 3;
    public static final int MYCOPLASMA_ID = 4;
    public static final int INVERTEBRATE_MT_ID = 5;
    public static final int CILIATE_ID = 6;
    public static final int ECHINODERM_MT_ID = 7;
    public static final int EUPLOTID_NUC_ID = 8;
    public static final int BACTERIAL_ID = 9;
    public static final int ALT_YEAST_ID = 10;
    public static final int ASCIDIAN_MT_ID = 11;
    public static final int FLATWORM_MT_ID = 12;
    public static final int BLEPHARISMA_NUC_ID = 13;
    public static final int NO_STOPS_ID = 14;

    /**
     * Standard genetic code tables from GENBANK
     * Nucleotides go A, C, G, T - Note: this is not the order used by the Genbank web site
     * With the first codon position most significant (i.e. AAA, AAC, AAG, AAT, ACA, etc.).
     */
    public static final String[] GENETIC_CODE_TABLES = {
        // Universal
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF",
        // Vertebrate Mitochondrial
        "KNKNTTTT*S*SMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
        // Yeast
        "KNKNTTTTRSRSMIMIQHQHPPPPRRRRTTTTEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
        // Mold Protozoan Mitochondrial
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
        // Mycoplasma
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
        // Invertebrate Mitochondrial
        "KNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
        // Ciliate
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSS*CWCLFLF",
        // Echinoderm Mitochondrial
        "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
        // Euplotid Nuclear
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSCCWCLFLF",
        // Bacterial
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF",
        // Alternative Yeast
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLSLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF",
        // Ascidian Mitochondrial
        "KNKNTTTTGSGSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
        // Flatworm Mitochondrial
        "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYY*YSSSSWCWCLFLF",
        // Blepharisma Nuclear
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YQYSSSS*CWCLFLF",
        // No stops
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYYQYSSSSWCWCLFLF"
    };

    /**
     * Names of the standard genetic code tables from GENBANK
     */
    public static final String[] GENETIC_CODE_NAMES = {
        "universal", "vertebrateMitochondrial", "yeast", "moldProtozoanMitochondrial",
        "mycoplasma", "invertebrateMitochondrial", "ciliate", "echinodermMitochondrial",
        "euplotidNuclear", "bacterial", "alternativeYeast", "ascidianMitochondrial",
        "flatwormMitochondrial", "blepharismaNuclear", "noStops"
    };

    /**
     * Descriptions of the standard genetic code tables from GENBANK
     */
    public static final String[] GENETIC_CODE_DESCRIPTIONS = {
        "Universal", "Vertebrate Mitochondrial", "Yeast", "Mold Protozoan Mitochondrial",
        "Mycoplasma", "Invertebrate Mitochondrial", "Ciliate", "Echinoderm Mitochondrial",
        "Euplotid Nuclear", "Bacterial", "Alternative Yeast", "Ascidian Mitochondrial",
        "Flatworm Mitochondrial", "Blepharisma Nuclear", "Test case with no stop codons"
    };

    public static final GeneticCode UNIVERSAL = new GeneticCode(UNIVERSAL_ID);
    public static final GeneticCode VERTEBRATE_MT = new GeneticCode(VERTEBRATE_MT_ID);
    public static final GeneticCode YEAST = new GeneticCode(YEAST_ID);
    public static final GeneticCode MOLD_PROTOZOAN_MT = new GeneticCode(MOLD_PROTOZOAN_MT_ID);
    public static final GeneticCode MYCOPLASMA = new GeneticCode(MYCOPLASMA_ID);
    public static final GeneticCode INVERTEBRATE_MT = new GeneticCode(INVERTEBRATE_MT_ID);
    public static final GeneticCode CILIATE = new GeneticCode(CILIATE_ID);
    public static final GeneticCode ECHINODERM_MT = new GeneticCode(ECHINODERM_MT_ID);
    public static final GeneticCode EUPLOTID_NUC = new GeneticCode(EUPLOTID_NUC_ID);
    public static final GeneticCode BACTERIAL = new GeneticCode(BACTERIAL_ID);
    public static final GeneticCode ALT_YEAST = new GeneticCode(ALT_YEAST_ID);
    public static final GeneticCode ASCIDIAN_MT = new GeneticCode(ASCIDIAN_MT_ID);
    public static final GeneticCode FLATWORM_MT = new GeneticCode(FLATWORM_MT_ID);
    public static final GeneticCode BLEPHARISMA_NUC = new GeneticCode(BLEPHARISMA_NUC_ID);
    public static final GeneticCode NO_STOPS = new GeneticCode(NO_STOPS_ID);

//    public static final GeneticCode[] GENETIC_CODES = {
//        UNIVERSAL, VERTEBRATE_MT, YEAST, MOLD_PROTOZOAN_MT, MYCOPLASMA, INVERTEBRATE_MT,
//        CILIATE, ECHINODERM_MT, EUPLOTID_NUC, BACTERIAL, ALT_YEAST, ASCIDIAN_MT,
//        FLATWORM_MT, BLEPHARISMA_NUC, NO_STOPS
//    };

    protected Nucleotide nucleotide = new Nucleotide(); // char <=> state

    public GeneticCode(int geneticCodeId) {
        this.geneticCodeId = geneticCodeId;
        codeTable = GENETIC_CODE_TABLES[geneticCodeId];
    }
    
    /**
     * Returns the name of the genetic code
     */
    public String getName() {
        return GENETIC_CODE_NAMES[geneticCodeId];
    }
    
    /**
     * Returns the description of the genetic code
     */
    public String getDescription() {
        return GENETIC_CODE_DESCRIPTIONS[geneticCodeId];
    }

    public static GeneticCode findByName(String codeStr) {
        GeneticCode geneticCode = null;
        if (codeStr.equals(GeneticCode.UNIVERSAL.getName())) {
            geneticCode = GeneticCode.UNIVERSAL;
        } else if (codeStr.equals(GeneticCode.VERTEBRATE_MT.getName())) {
            geneticCode = GeneticCode.VERTEBRATE_MT;
        } else if (codeStr.equals(GeneticCode.YEAST.getName())) {
            geneticCode = GeneticCode.YEAST;
        } else if (codeStr.equals(GeneticCode.MOLD_PROTOZOAN_MT.getName())) {
            geneticCode = GeneticCode.MOLD_PROTOZOAN_MT;
        } else if (codeStr.equals(GeneticCode.INVERTEBRATE_MT.getName())) {
            geneticCode = GeneticCode.INVERTEBRATE_MT;
        } else if (codeStr.equals(GeneticCode.CILIATE.getName())) {
            geneticCode = GeneticCode.CILIATE;
        } else if (codeStr.equals(GeneticCode.ECHINODERM_MT.getName())) {
            geneticCode = GeneticCode.ECHINODERM_MT;
        } else if (codeStr.equals(GeneticCode.EUPLOTID_NUC.getName())) {
            geneticCode = GeneticCode.EUPLOTID_NUC;
        } else if (codeStr.equals(GeneticCode.BACTERIAL.getName())) {
            geneticCode = GeneticCode.BACTERIAL;
        } else if (codeStr.equals(GeneticCode.ALT_YEAST.getName())) {
            geneticCode = GeneticCode.ALT_YEAST;
        } else if (codeStr.equals(GeneticCode.ASCIDIAN_MT.getName())) {
            geneticCode = GeneticCode.ASCIDIAN_MT;
        } else if (codeStr.equals(GeneticCode.FLATWORM_MT.getName())) {
            geneticCode = GeneticCode.FLATWORM_MT;
        } else if (codeStr.equals(GeneticCode.BLEPHARISMA_NUC.getName())) {
            geneticCode = GeneticCode.BLEPHARISMA_NUC;
        } else if (codeStr.equals(GeneticCode.NO_STOPS.getName())) {
            geneticCode = GeneticCode.NO_STOPS;
        } else {
            throw new RuntimeException("Unknown genetics code");
        }
        return geneticCode;
    }


    /**
     * Get Nucleotide states from {@link Nucleotide#codeMap codeMap}.
     * @see DataType#string2state(String)
     * @param c
     * @return
     */
    public int getNucleotideState(char c) {
        return nucleotide.string2state(String.valueOf(c)).get(0);
    }

    /**
     * Infer char from {@link Nucleotide#codeMap codeMap}.
     * @see DataType#state2string(int[])
     * @param state
     * @return
     */
    public char getNucleotideChar(int state) {
        return nucleotide.state2string(new int[]{state}).charAt(0);
    }

    /**
     * This table maps amino acid characters into state codes (0-25).
     * Amino Acids go ACDEFGHIKLMNPQRSTVWYBZX*?-,
     * Other letters; j, o, and u are mapped to ?
     * *, ? and - are mapped to themselves
     * All other chars are mapped to -
     */
    public static final int[] AMINOACID_STATES = {
            25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,    // 0-15
            25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,    // 16-31
    //                                     *        -
            25,25,25,25,25,25,25,25,25,25,23,25,25,25,25,25,    // 32-47
    //                                                    ?
            25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,24,    // 48-63
    //          A  B  C  D  E  F  G  H  I  j  K  L  M  N  o
            25, 0,20, 1, 2, 3, 4, 5, 6, 7,24, 8, 9,10,11,24,    // 64-79
    //       P  Q  R  S  T  u  V  W  X  Y  Z
            12,13,14,15,16,24,17,18,22,19,21,25,25,25,25,25,    // 80-95
    //          A  B  C  D  E  F  G  H  I  j  K  L  M  N  o
            25, 0,20, 1, 2, 3, 4, 5, 6, 7,24, 8, 9,10,11,24,    // 96-111
    //       P  Q  R  S  T  u  V  W  X  Y  Z
            12,13,14,15,16,24,17,18,22,19,21,25,25,25,25,25     // 112-127
    };

    /**
     * A table to translate state numbers (0-25) into one letter codes.
     */
    public static final char[] AMINOACID_CHARS= {
            'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R',
            'S','T','V','W','Y','B','Z','X',
            Codon.STOP_CHARACTER,DataType.MISSING_CHAR,DataType.GAP_CHAR
    };

    /**
     * Get character corresponding to a given state mapping to
     * {@link #AMINOACID_STATES AMINOACID_STATES}.
     * @return
     */
    public char getAminoAcidChar(int state) {
        return AMINOACID_CHARS[state];
    }

    /**
     * Returns the state associated with AminoAcid ({@link #AMINOACID_STATES AMINOACID_STATES})
     * represented by codonState.
     * Note: the state is the canonical state (generated combinatoriall),
     * which is different to BEAST 2 Aminoacid states, and which can be
     * from {@link Codon#getCanonicalState(int)} getCanonicalState}.
     *
     * @see #AMINOACID_STATES
     * @return '?' if codon unknown
     */
    public int getAminoAcidCodonState(int codonState) {
        if (codonState == Codon.UNKNOWN_STATE)
            return Aminoacid.MISSING_CHAR;
        else if (codonState == Codon.GAP_STATE)
            return Aminoacid.GAP_CHAR;

        return AMINOACID_STATES[codeTable.charAt(codonState)];
    }

//    /**
//     * Note that the state is the canonical state (generated combinatorially)
//     * @return whether the codonState is a stop codon
//     */
//    public boolean isStopCodon(int codonState) {
//        return (getAminoAcidCodonState(codonState) == Codon.STOP_STATE);
//    }


    public boolean isStopCodon(int codonState) {
        return (codeTable.charAt(codonState) == Codon.STOP_CHARACTER);
    }

    /**
     * @return the codon states of stop amino acids.
     */
    public int[] getStopCodonIndices() {
    
        int i, j, n = getStopCodonCount();
        int[] indices = new int[n];
        
        j = 0;
        for (i = 0; i < 64; i++) {
            if (codeTable.charAt(i) == Codon.STOP_CHARACTER) {
                indices[j] = i;
                j++;
            }
        }
        
        return indices;
    }

    /**
     * Returns the number of terminator amino acids.
     */
    public int getStopCodonCount() {
        int i, count = 0;
        
        for (i = 0; i < 64; i++) {
            if (codeTable.charAt(i) == Codon.STOP_CHARACTER)
                count++;
        }
        
        return count;
    }
    
    private int geneticCodeId;
    private String codeTable;

}
