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
import beast.core.Param;
import beast.evolution.alignment.CodonAlignment;
import beast.evolution.datatype.Codon;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.GeneticCode;
import beast.evolution.substitutionmodel.GeneralSubstitutionModel;

import java.lang.reflect.InvocationTargetException;

/**
 * Modified from BEAST 1 AbstractCodonModel.
 *
 * @author Marc A. Suchard
 * @author Walter Xie
 */
@Description("Abstract codon model to construct an array of rate classes " +
        "using the current genetic code.")
public class AbstractCodonModel extends GeneralSubstitutionModel {

    final public Input<CodonAlignment> convertAlignmentInput = new Input<>("data",
            "Converted alignment to provide codon data type", Input.Validate.REQUIRED);
    final public Input<Boolean> printRateMapInput = new Input<>("printRateMap",
            "Print the rate classes in the rate matrix using the current genetic code",
            Boolean.FALSE);

    protected byte[] rateMap;

    protected Codon codonDataType;
    protected GeneticCode geneticCode;
    protected int rateCount;

    public AbstractCodonModel() {
        ratesInput.setRule(Input.Validate.FORBIDDEN); // only use internally
    }

    public AbstractCodonModel(@Param(name="mode", description="model identifier is a 6 digit number describing which parameters are linked. "
            + "For instance 010010 = HKY, 012345 = GTR. "
            + "Order of digits are for rates AC, AG, AT, CG, CT, GT respectively ",
            defaultValue="000000") String modelID) {
        this();
    }

    @Override
    public void initAndValidate() {
        this.frequencies = frequenciesInput.get();

        CodonAlignment alignment = convertAlignmentInput.get();
        DataType dataType = alignment.getDataType();

        if (! (dataType instanceof Codon) )
            throw new IllegalArgumentException("Codon data type is required !");
        this.codonDataType = (Codon) dataType;

        this.geneticCode = codonDataType.getGeneticCode();
        System.out.println("Genetic code is " + geneticCode.getDescription());

        //====== init states and rates ======
        updateMatrix = true;
        nrOfStates = frequencies.getFreqs().length;

        try {
            eigenSystem = createEigenSystem();
        } catch (SecurityException | ClassNotFoundException | InstantiationException |
                IllegalAccessException | IllegalArgumentException | InvocationTargetException e) {
            throw new IllegalArgumentException(e.getMessage());
        }
        //eigenSystem = new DefaultEigenSystem(m_nStates);

        rateMatrix = new double[nrOfStates][nrOfStates];

        rateCount = getRateCount(nrOfStates);
        if (rateCount < 1)
            throw new RuntimeException("Incorrect rate count = " + rateCount);
        relativeRates = new double[rateCount];
        storedRelativeRates = new double[rateCount];

        rateMap = constructRateMap(rateCount, nrOfStates, codonDataType, geneticCode);

        if (printRateMapInput.get())
            printRateMap(); // debug
    }

    //TODO move to GeneralSubstitutionModel ?
    // relativeRates.length = nrOfStates * (nrOfStates - 1) in GeneralSubstitutionModel
    protected int getRateCount(int stateCount) {
        return ((stateCount - 1) * stateCount);
    }

    /**
     *
     * Construct a map of the rate classes in the rate matrix
     * using the current genetic code. Classes:
     *		0: codon changes in more than one codon position (or stop codons)
     *		1: synonymous transition
     *		2: synonymous transversion
     *		3: non-synonymous transition
     *		4: non-synonymous transversion
     *
     * The map is in an array whose length = nrOfStates * (nrOfStates - 1).
     * @param rateCount  get from {@link #getRateCount(int)} getRateCount}.
     * @param stateCount  nrOfStates
     * @param codonDataType  {@link Codon Codon}
     * @param geneticCode  {@link GeneticCode#GENETIC_CODE_NAMES GENETIC_CODE_NAMES}
     * @return
     */
    protected byte[] constructRateMap(int rateCount, int stateCount, Codon codonDataType, GeneticCode geneticCode)	{
        byte rateClass;
        byte[] rateMap = new byte[rateCount];

        // this needs to match rateMatrix[i][j] <= relativeRates[] in setupRateMatrix()
        for (int i = 0; i < stateCount; i++) {

            // i1, j1, k1, aa1
            int[] ids1 = getCodonStatesForRateClass(i, codonDataType, geneticCode);

            for (int j = 0; j < i; j++) {
                // i2, j2, k2, aa2
                int[] ids2 = getCodonStatesForRateClass(j, codonDataType, geneticCode);

                rateClass = getRateClass(ids1[0], ids1[1], ids1[2], ids2[0], ids2[1], ids2[2], ids1[3], ids2[3]);
                rateMap[i * (stateCount - 1) + j] = rateClass;

            }
            for (int j = i + 1; j < stateCount; j++) {
                // i2, j2, k2, aa2
                int[] ids2 = getCodonStatesForRateClass(j, codonDataType, geneticCode);

                rateClass = getRateClass(ids1[0], ids1[1], ids1[2], ids2[0], ids2[1], ids2[2], ids1[3], ids2[3]);
                rateMap[i * (stateCount - 1) + j - 1] = rateClass;
            }
        }

        return rateMap;
    }

    // return i1, j1, k1, aa1 given i
    private int[] getCodonStatesForRateClass(int i, Codon codonDataType, GeneticCode geneticCode) {
        int i1, j1, k1, cs1, aa1;
        int[] codon  = codonDataType.getTripletStates(i);
        i1 = codon[0];
        j1 = codon[1];
        k1 = codon[2];

        cs1 = codonDataType.getState(i1, j1, k1);
        aa1 = geneticCode.getAminoAcidCodonState(codonDataType.getCanonicalState(cs1));
        return new int[]{i1, j1, k1, aa1};
    }

    // get rateClass for constructRateMap
    private byte getRateClass(int i1, int j1, int k1, int i2, int j2, int k2, int aa1, int aa2) {
        byte rateClass = -1;
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
        return rateClass;
    }

    @Override
    protected void setupRelativeRates() {
        throw new UnsupportedOperationException("setupRelativeRates needs to override in codon models !");
    }


    @Override
    public boolean canHandleDataType(DataType dataType) {
        return dataType instanceof Codon;
    }

    /** rateClass:
     *		0: codon changes in more than one codon position (or stop codons)
     *		1: synonymous transition
     *		2: synonymous transversion
     *		3: non-synonymous transition
     *		4: non-synonymous transversion
     */
    public void printRateMap() {
        byte rateClass;
        int stateCount = nrOfStates;

        System.out.println("\n============ rate matrix ============");
        System.out.println("  0: codon changes in more than one codon position (or stop codons)");
        System.out.println("  1: synonymous transition");
        System.out.println("  2: synonymous transversion");
        System.out.println("  3: non-synonymous transition");
        System.out.println("  4: non-synonymous transversion");
        System.out.print("\t");
        for (int j = 0; j < stateCount; j++) {
            // i2, j2, k2, aa2
            int[] ids2 = getCodonStatesForRateClass(j, codonDataType, geneticCode);

            System.out.print("\t" + geneticCode.getNucleotideChar(ids2[0]));
            System.out.print(geneticCode.getNucleotideChar(ids2[1]));
            System.out.print(geneticCode.getNucleotideChar(ids2[2]));
        }
        System.out.println();

        System.out.print("\t");
        for (int j = 0; j < stateCount; j++) {
            int[] ids2 = getCodonStatesForRateClass(j, codonDataType, geneticCode);
            System.out.print("\t" + geneticCode.getAminoAcidChar(ids2[3]));
        }
        System.out.println();

        for (int i = 0; i < stateCount; i++) {

            // i1, j1, k1, aa1
            int[] ids1 = getCodonStatesForRateClass(i, codonDataType, geneticCode);

            System.out.print(geneticCode.getNucleotideChar(ids1[0]));
            System.out.print(geneticCode.getNucleotideChar(ids1[1]));
            System.out.print(geneticCode.getNucleotideChar(ids1[2]));
            System.out.print("\t" + geneticCode.getAminoAcidChar(ids1[3]));

            for (int j = 0; j < i; j++) {
                // i2, j2, k2, aa2
                int[] ids2 = getCodonStatesForRateClass(j, codonDataType, geneticCode);

                rateClass = getRateClass(ids1[0], ids1[1], ids1[2], ids2[0], ids2[1], ids2[2], ids1[3], ids2[3]);
                System.out.print("\t" + rateClass);

            }
            System.out.print("\t.");// i=j
            for (int j = i + 1; j < stateCount; j++) {
                // i2, j2, k2, aa2
                int[] ids2 = getCodonStatesForRateClass(j, codonDataType, geneticCode);

                rateClass = getRateClass(ids1[0], ids1[1], ids1[2], ids2[0], ids2[1], ids2[2], ids1[3], ids2[3]);
                System.out.print("\t" + rateClass);
            }
            System.out.println();
        }

        System.out.println("\n============ rates in array ============");
        int col = rateMap.length / (stateCount - 1);
        for (int i = 0; i < col; i++) {
            System.out.print("rateMap[" + i * (stateCount - 1) + " - " + (i+1) * (stateCount - 1) + "] = ");
            for (int j = 0; j < (stateCount - 1); j++) {
                System.out.print("\t" + rateMap[i * (stateCount - 1) + j]);
            }
            System.out.println();
        }
    }

}