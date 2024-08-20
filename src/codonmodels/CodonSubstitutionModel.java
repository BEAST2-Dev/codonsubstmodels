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

package codonmodels;


import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import codonmodels.evolution.alignment.CodonAlignment;
import codonmodels.evolution.datatype.Codon;
import codonmodels.evolution.datatype.GeneticCode;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.base.evolution.tree.Node;

import java.lang.reflect.InvocationTargetException;

/**
 * Modified from BEAST 1 AbstractCodonModel.
 *
 * This also uses new code to compute transition-prob matrix,
 * which is faster than
 * {@link GeneralSubstitutionModel#getTransitionProbabilities(Node, double, double, double, double[], boolean)}.
 *
 * @author Alexei Drummond
 * @author Marc A. Suchard
 * @author Walter Xie
 */
@Description("Codon substitution model to construct an array of rate classes " +
        "using the current genetic code.")
public class CodonSubstitutionModel extends GeneralSubstitutionModel {

    final public Input<Boolean> verboseInput = new Input<>("verbose",
            "Print the rate classes in the rate matrix, etc.",
            Boolean.TRUE);

    protected byte[] rateMap;

    protected Codon codonDataType;
//    protected GeneticCode geneticCode; // get from codon
    protected int rateCount;

    public CodonSubstitutionModel() {
        ratesInput.setRule(Input.Validate.FORBIDDEN); // only use internally
    }

    @Override
    public void initAndValidate() {
        this.frequencies = frequenciesInput.get();
//        if (! (frequencies instanceof CodonFrequencies) )
//            throw new IllegalArgumentException("Codon frequencies is required by CodonSubstitutionModel !");

        Alignment data = null;
        if (frequencies instanceof CodonFrequencies) {
        	data = frequencies.dataInput.get();
        } else {
        	for (BEASTInterface o : getOutputs()) {
        		if (o instanceof SiteModel) {
        			SiteModel sitemodel = (SiteModel) o;
                	for (BEASTInterface o2 : sitemodel.getOutputs()) {
                		if (o2 instanceof GenericTreeLikelihood) {
                			GenericTreeLikelihood tl = (GenericTreeLikelihood) o2;
                			data = tl.dataInput.get();
                		}
                	}        			
        		}
        	}
        }
        
        CodonAlignment alignment = getCodonAlignment(data);
        this.codonDataType = alignment.getDataType();

        //====== init states and rates ======
        updateMatrix = true;
        double[] freqs = frequencies.getFreqs();
        nrOfStates = freqs.length;

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

        rateMap = constructRateMap(rateCount, nrOfStates, codonDataType);

        if (verboseInput.get())
            printRateMap(codonDataType); // debug

    }

    /**
     * Faster code to replace
     * {@link GeneralSubstitutionModel#getTransitionProbabilities(Node, double, double, double, double[], boolean)}.
     *
     * @param startTime parent.getHeight()
     * @param endTime   node.getHeight()
     * @param rate      joint rate = rate for a site category * mean branch rate.
     * @param iexp      iexp, without creating a new matrix each call.
     * @param matrix    P(t), without creating a new matrix each call.
     */
    public void getTransiProbs(double startTime, double endTime, double rate,
                               double[] iexp, double[] matrix) {//, boolean normalized) {
        double distance = (startTime - endTime) * rate;

        getTransiProbs(distance, iexp, matrix);
    }

    /**
     * Faster code to replace
     * {@link GeneralSubstitutionModel#getTransitionProbabilities(Node, double, double, double, double[], boolean)}.
     * @param distance  distance = (startTime - endTime) * mean branch rate * rate for a site category.
     * @param iexp      iexp, without creating a new matrix each call.
     * @param matrix    P(t), without creating a new matrix each call.
     */
    public void getTransiProbs(double distance, double[] iexp, double[] matrix) {
        int i, j, k;
        double temp;
        // this must be synchronized to avoid being called simultaneously by
        // two different likelihood threads - AJD
        synchronized (this) {
            if (updateMatrix) {
                setupRelativeRates();
//                if (normalized) {
                setupRateMatrix();
//                } else {
//                    setupRateMatrixUnnormalized();
//                }
                eigenDecomposition = eigenSystem.decomposeMatrix(rateMatrix);
                updateMatrix = false;
            }
        }

        // is the following really necessary?
        // implemented a pool of iexp matrices to support multiple threads
        // a quick timing experiment shows no difference - RRB
//        double[] iexp = new double[nrOfStates * nrOfStates];
        double[] Evec = eigenDecomposition.getEigenVectors();
        // inverse Eigen vectors
        double[] Ievc = eigenDecomposition.getInverseEigenVectors();
        // Eigen values
        double[] Eval = eigenDecomposition.getEigenValues();

        // faster computation, reviewed by AJD
        int x = 0;
        for (i = 0; i < nrOfStates; i++) {
            temp = Math.exp(distance * Eval[i]);
            for (j = 0; j < nrOfStates; j++) {
                // iexp[i * nrOfStates + j] = Ievc[i * nrOfStates + j] * temp;
                iexp[x] = Ievc[x] * temp;
                x++; // save time, plus once
            }
        }

        x = 0;
        int y;
        int u = 0;
        for (i = 0; i < nrOfStates; i++) {
            for (j = 0; j < nrOfStates; j++) {
                y = j;
                temp = 0.0;
                for (k = 0; k < nrOfStates; k++) {
//                    temp += Evec[i * nrOfStates + k] * iexp[k * nrOfStates + j];
                    temp += Evec[x + k] * iexp[y];
                    y += nrOfStates; // y = k + nrOfStates
                }

                matrix[u] = Math.abs(temp);
                u++; // u = u + j
            }
            x += nrOfStates; // x = i + nrOfStates
        }
    }

    protected CodonAlignment getCodonAlignment(Alignment data) {
        return CodonAlignment.toCodonAlignment(data);
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
     * @param stateCount  nrOfStates 64
     * @param codon  {@link Codon Codon}
     * @return
     */
    protected byte[] constructRateMap(int rateCount, int stateCount, Codon codon)	{
        byte rateClass;
        byte[] rateMap = new byte[rateCount];

        // this needs to match rateMatrix[i][j] <= relativeRates[] in setupRateMatrix()
        // i j is codonState
        for (int i = 0; i < stateCount; i++) {

            // i1, j1, k1, aa1
            int[] ids1 = getCodonStatesForRateClass(i, codon);

            for (int j = 0; j < i; j++) {
                // i2, j2, k2, aa2
                int[] ids2 = getCodonStatesForRateClass(j, codon);

                rateClass = getRateClass(ids1[0], ids1[1], ids1[2], ids2[0], ids2[1], ids2[2], ids1[3], ids2[3]);
                rateMap[i * (stateCount - 1) + j] = rateClass;

            }
            for (int j = i + 1; j < stateCount; j++) {
                // i2, j2, k2, aa2
                int[] ids2 = getCodonStatesForRateClass(j, codon);

                rateClass = getRateClass(ids1[0], ids1[1], ids1[2], ids2[0], ids2[1], ids2[2], ids1[3], ids2[3]);
                rateMap[i * (stateCount - 1) + j - 1] = rateClass;
            }
        }

        return rateMap;
    }

    // return i1, j1, k1, aa1, given codonState
    private int[] getCodonStatesForRateClass(int codonState, Codon codon) {
        int i1, j1, k1, cs1, aa1;
//        int codonState = codon.getCodonState(stateMapIndex);
        int[] nucStates  = codon.getTripletNucStates(codonState);
        i1 = nucStates[0];
        j1 = nucStates[1];
        k1 = nucStates[2];

        cs1 = codon.getCodonState(i1, j1, k1);

        GeneticCode geneticCode = codon.getGeneticCode();
        aa1 = geneticCode.getAminoAcidState(cs1);
        return new int[]{i1, j1, k1, aa1};
    }

    // get rateClass for constructRateMap, i,j = {A,C,G,T}
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

//    @Override
//    protected void setupRelativeRates() {
//        throw new UnsupportedOperationException("setupRelativeRates needs to be overridden " +
//                "in the child class of CodonSubstitutionModel !");
//    }


    @Override
    public boolean canHandleDataType(DataType dataType) {
        return dataType instanceof Codon;
    }

    //============ print ============

    /** rateClass:
     *		0: codon changes in more than one codon position (or stop codons)
     *		1: synonymous transition
     *		2: synonymous transversion
     *		3: non-synonymous transition
     *		4: non-synonymous transversion
     * @param codon
     */
    protected void printRateMap(Codon codon) {
        byte rateClass;
        int stateCount = nrOfStates;
        GeneticCode geneticCode = codon.getGeneticCode();

        Log.info.println("\n============ Rate Matrix ============");
        Log.info.println("  0: codon changes in more than one codon position (or stop codons)");
        Log.info.println("  1: synonymous transition");
        Log.info.println("  2: synonymous transversion");
        Log.info.println("  3: non-synonymous transition");
        Log.info.println("  4: non-synonymous transversion");
        Log.info.print("\t");
        for (int j = 0; j < stateCount; j++) {
            // i2, j2, k2, aa2
            int[] ids2 = getCodonStatesForRateClass(j, codon);

            Log.info.print("\t" + geneticCode.getNucleotideChar(ids2[0]));
            Log.info.print(geneticCode.getNucleotideChar(ids2[1]));
            Log.info.print(geneticCode.getNucleotideChar(ids2[2]));
        }
        Log.info.println();

        Log.info.print("\t");
        for (int j = 0; j < stateCount; j++) {
            int[] ids2 = getCodonStatesForRateClass(j, codon);
            Log.info.print("\t" + geneticCode.getAminoAcidChar(ids2[3]));
        }
        Log.info.println();

        for (int i = 0; i < stateCount; i++) {

            // i1, j1, k1, aa1
            int[] ids1 = getCodonStatesForRateClass(i, codon);

            Log.info.print(geneticCode.getNucleotideChar(ids1[0]));
            Log.info.print(geneticCode.getNucleotideChar(ids1[1]));
            Log.info.print(geneticCode.getNucleotideChar(ids1[2]));
            Log.info.print("\t" + geneticCode.getAminoAcidChar(ids1[3]));
            // lower triangle
            for (int j = 0; j < i; j++) {
                // i2, j2, k2, aa2
                int[] ids2 = getCodonStatesForRateClass(j, codon);

                rateClass = getRateClass(ids1[0], ids1[1], ids1[2], ids2[0], ids2[1], ids2[2], ids1[3], ids2[3]);
                Log.info.print("\t" + rateClass);

            }
            Log.info.print("\t.");// i=j
            // upper triangle
            for (int j = i + 1; j < stateCount; j++) {
                // i2, j2, k2, aa2
                int[] ids2 = getCodonStatesForRateClass(j, codon);

                rateClass = getRateClass(ids1[0], ids1[1], ids1[2], ids2[0], ids2[1], ids2[2], ids1[3], ids2[3]);
                Log.info.print("\t" + rateClass);
            }
            Log.info.println();
        }

//        Log.info.println("\n============ rates in array ============");
//        int col = rateMap.length / (stateCount - 1);
//        for (int i = 0; i < col; i++) {
//            Log.info.print("rateMap[" + i * (stateCount - 1) + " - " + (i+1) * (stateCount - 1) + "] = ");
//            for (int j = 0; j < (stateCount - 1); j++) {
//                Log.info.print("\t" + rateMap[i * (stateCount - 1) + j]);
//            }
//            Log.info.println();
//        }
    }

}