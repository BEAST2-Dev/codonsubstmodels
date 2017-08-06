/*
 * AbstractCodonModel.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package beast.evolution.substitutionmodel.codons;


import beast.core.Description;
import beast.core.Input;
import beast.evolution.alignment.ConvertAlignment;
import beast.evolution.datatype.Codon;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.GeneticCode;
import beast.evolution.substitutionmodel.GeneralSubstitutionModel;

import java.lang.reflect.InvocationTargetException;

/**
 * Imported from BEAST 1.
 * @author Marc A. Suchard
 */
@Description("Abstract codon model")
public abstract class AbstractCodonModel extends GeneralSubstitutionModel {

    //TODO I only need DataType, how to pass it simply?
    final public Input<ConvertAlignment> convertAlignmentInput = new Input<>("data",
            "Converted alignment to provide codon data type", Input.Validate.REQUIRED);

    protected byte[] rateMap;

    protected Codon codonDataType;
    protected GeneticCode geneticCode;
    protected int rateCount;

    public AbstractCodonModel() {
        ratesInput.setRule(Input.Validate.FORBIDDEN); // only use internally
//        try {
//            ratesInput.setValue(null, this);
//        } catch (Exception e) {
//            e.printStackTrace();
//            // TODO: handle exception
//        }
    }

    @Override
    public void initAndValidate() {
        this.frequencies = frequenciesInput.get();

        ConvertAlignment alignment = convertAlignmentInput.get();
        DataType dataType = alignment.getDataType();

        if (! (dataType instanceof Codon) )
            throw new IllegalArgumentException("Codon data type is required !");
        this.codonDataType = (Codon) dataType;

        this.geneticCode = codonDataType.getGeneticCode();

        //====== init states and rates ======
        updateMatrix = true;
        nrOfStates = frequencies.getFreqs().length;

        try {
            eigenSystem = createEigenSystem();
        } catch (SecurityException | ClassNotFoundException | InstantiationException | IllegalAccessException | IllegalArgumentException
                | InvocationTargetException e) {
            throw new IllegalArgumentException(e.getMessage());
        }
        //eigenSystem = new DefaultEigenSystem(m_nStates);

        rateMatrix = new double[nrOfStates][nrOfStates];

        rateCount = getRateCount(nrOfStates);
        if (rateCount < 1)
            throw new RuntimeException("Incorrect rate count = " + rateCount);
        relativeRates = new double[rateCount];
        storedRelativeRates = new double[rateCount];

        constructRateMap();
    }

    //TODO move to SubstitutionModel.Base
    protected int getRateCount(int stateCount) {
        return ((stateCount - 1) * stateCount);
    }


    /**
     * Construct a map of the rate classes in the rate matrix using the current
     * genetic code. Classes:
     * 0: codon changes in more than one codon position (or stop codons)
     * 1: synonymous transition
     * 2: synonymous transversion
     * 3: non-synonymous transition
     * 4: non-synonymous transversion
     */
    protected void constructRateMap() {
        // Refactored into static function, since CodonProductChains need this functionality
        rateMap = constructRateMap(rateCount, nrOfStates, codonDataType, geneticCode);
    }

    /**
     * Construct a map of the rate classes in the rate matrix using the current
     * genetic code. Classes:
     *		0: codon changes in more than one codon position (or stop codons)
     *		1: synonymous transition
     *		2: synonymous transversion
     *		3: non-synonymous transition
     *		4: non-synonymous transversion
     */
    public byte[] constructRateMap(int rateCount, int stateCount, Codon codonDataType, GeneticCode geneticCode)	{
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

    @Override
    public boolean canHandleDataType(DataType dataType) {
        return dataType instanceof Codon;
    }

//    public void printRateMap() {
//        int u, v, i1, j1, k1, i2, j2, k2;
//        byte rateClass;
//        int[] codon;
//        int cs1, cs2, aa1, aa2;
//
//        System.out.print("\t");
//        for (v = 0; v < stateCount; v++) {
//            codon = codonDataType.getTripletStates(v);
//            i2 = codon[0];
//            j2 = codon[1];
//            k2 = codon[2];
//
//            System.out.print("\t" + Nucleotides.INSTANCE.getChar(i2));
//            System.out.print(Nucleotides.INSTANCE.getChar(j2));
//            System.out.print(Nucleotides.INSTANCE.getChar(k2));
//        }
//        System.out.println();
//
//        System.out.print("\t");
//        for (v = 0; v < stateCount; v++) {
//            codon = codonDataType.getTripletStates(v);
//            i2 = codon[0];
//            j2 = codon[1];
//            k2 = codon[2];
//
//            cs2 = codonDataType.getState(i2, j2, k2);
//            aa2 = geneticCode.getAminoAcidState(codonDataType.getCanonicalState(cs2));
//            System.out.print("\t" + AminoAcids.INSTANCE.getChar(aa2));
//        }
//        System.out.println();
//
//        for (u = 0; u < stateCount; u++) {
//
//            codon = codonDataType.getTripletStates(u);
//            i1 = codon[0];
//            j1 = codon[1];
//            k1 = codon[2];
//
//            System.out.print(Nucleotides.INSTANCE.getChar(i1));
//            System.out.print(Nucleotides.INSTANCE.getChar(j1));
//            System.out.print(Nucleotides.INSTANCE.getChar(k1));
//
//            cs1 = codonDataType.getState(i1, j1, k1);
//            aa1 = geneticCode.getAminoAcidState(codonDataType.getCanonicalState(cs1));
//
//            System.out.print("\t" + AminoAcids.INSTANCE.getChar(aa1));
//
//            for (v = 0; v < stateCount; v++) {
//
//                codon = codonDataType.getTripletStates(v);
//                i2 = codon[0];
//                j2 = codon[1];
//                k2 = codon[2];
//
//                cs2 = codonDataType.getState(i2, j2, k2);
//                aa2 = geneticCode.getAminoAcidState(codonDataType.getCanonicalState(cs2));
//
//                rateClass = -1;
//                if (i1 != i2) {
//                    if ((i1 == 0 && i2 == 2) || (i1 == 2 && i2 == 0) || // A <-> G
//                            (i1 == 1 && i2 == 3) || (i1 == 3 && i2 == 1)) { // C <-> T
//                        rateClass = 1; // Transition at position 1
//                    } else {
//                        rateClass = 2; // Transversion at position 1
//                    }
//                }
//                if (j1 != j2) {
//                    if (rateClass == -1) {
//                        if ((j1 == 0 && j2 == 2) || (j1 == 2 && j2 == 0) || // A <-> G
//                                (j1 == 1 && j2 == 3) || (j1 == 3 && j2 == 1)) { // C <-> T
//                            rateClass = 1; // Transition
//                        } else {
//                            rateClass = 2; // Transversion
//                        }
//                    } else
//                        rateClass = 0; // Codon changes at more than one position
//                }
//                if (k1 != k2) {
//                    if (rateClass == -1) {
//                        if ((k1 == 0 && k2 == 2) || (k1 == 2 && k2 == 0) || // A <-> G
//                                (k1 == 1 && k2 == 3) || (k1 == 3 && k2 == 1)) { // C <-> T
//                            rateClass = 1; // Transition
//                        } else {
//                            rateClass = 2; // Transversion
//                        }
//                    } else
//                        rateClass = 0; // Codon changes at more than one position
//                }
//
//                if (rateClass != 0) {
//                    if (aa1 != aa2) {
//                        rateClass += 2; // Is a non-synonymous change
//                    }
//                }
//
//                System.out.print("\t" + rateClass);
//
//            }
//            System.out.println();
//
//        }
//    }
}