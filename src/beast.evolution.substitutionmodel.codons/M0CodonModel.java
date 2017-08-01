///*
// * GY94CodonModel.java
// *
// * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
// *
// * This file is part of BEAST.
// * See the NOTICE file distributed with this work for additional
// * information regarding copyright ownership and licensing.
// *
// * BEAST is free software; you can redistribute it and/or modify
// * it under the terms of the GNU Lesser General Public License as
// * published by the Free Software Foundation; either version 2
// * of the License, or (at your option) any later version.
// *
// *  BEAST is distributed in the hope that it will be useful,
// *  but WITHOUT ANY WARRANTY; without even the implied warranty of
// *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// *  GNU Lesser General Public License for more details.
// *
// * You should have received a copy of the GNU Lesser General Public
// * License along with BEAST; if not, write to the
// * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
// * Boston, MA  02110-1301  USA
// */
//
//package beast.evolution.substitutionmodel.codons;
//
//import beast.core.Citation;
//import beast.core.Description;
//
//
//import java.util.Collections;
//import java.util.List;
//
///**
// * Yang model of codon evolution
// *
// * @author Andrew Rambaut
// * @author Alexei Drummond
// * @author Marc A. Suchard
// * @version $Id: YangCodonModel.java,v 1.21 2005/05/24 20:25:58 rambaut Exp $
// */
//@Citation("Nick Goldman and Ziheng Yang. A codon-based model of nucleotide substitution for protein-coding DNA sequences. " +
//        "Molecular biology and evolution 11.5 (1994): 725-736.")
//@Description("M0 codon model, also called as GY94, published by Goldman and Yang 1994")
//public class M0CodonModel extends AbstractCodonModel {
//    /**
//     * kappa
//     */
//    protected Parameter kappaParameter;
//
//    /**
//     * omega
//     */
//    protected Parameter omegaParameter;
//
//    public M0CodonModel(Codons codonDataType, Parameter omegaParameter, Parameter kappaParameter,
//                          FrequencyModel freqModel) {
//        this(codonDataType, omegaParameter, kappaParameter, freqModel,
//                new DefaultEigenSystem(codonDataType.getStateCount()));
//    }
//
//    public M0CodonModel(Codons codonDataType,
//                          Parameter omegaParameter,
//                          Parameter kappaParameter,
//                          FrequencyModel freqModel, EigenSystem eigenSystem) {
//
//        super("GY94", codonDataType, freqModel, eigenSystem);
//
//        this.omegaParameter = omegaParameter;
//
//        int dim = omegaParameter.getDimension();
//        double value = omegaParameter.getParameterValue(dim - 1);
//        if(value < 0) {
//        	throw new RuntimeException("Negative Omega parameter value " + value);
//        }//END: negative check
//
//        addVariable(omegaParameter);
//        omegaParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0,
//                omegaParameter.getDimension()));
//
//        this.kappaParameter = kappaParameter;
//
//        dim = kappaParameter.getDimension();
//        value = kappaParameter.getParameterValue(dim - 1);
//        if(value < 0) {
//        	throw new RuntimeException("Negative kappa parameter value value " + value);
//        }//END: negative check
//
//        addVariable(kappaParameter);
//        kappaParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0,
//                kappaParameter.getDimension()));
//
//        // Assuming it's always the same dim
//
//
//        addStatistic(synonymousRateStatistic);
//    }
//
//    /**
//     * set kappa
//     *
//     * @param kappa kappa
//     */
//    public void setKappa(double kappa) {
//        kappaParameter.setParameterValue(0, kappa);
//        updateMatrix = true;
//    }
//
//    /**
//     * @return kappa
//     */
//    public double getKappa() {
//        return kappaParameter.getParameterValue(0);
//    }
//
//    /**
//     * set dN/dS
//     *
//     * @param omega omega
//     */
//    public void setOmega(double omega) {
//        omegaParameter.setParameterValue(0, omega);
//        updateMatrix = true;
//    }
//
//    /**
//     * @return dN/dS
//     */
//    public double getOmega() {
//        return omegaParameter.getParameterValue(0);
//    }
//
//    public double getSynonymousRate() {
//        double k = getKappa();
//        double o = getOmega();
//        return ((31.0 * k) + 36.0) / ((31.0 * k) + 36.0 + (138.0 * o) + (58.0 * o * k));
//    }
//
//    public double getNonSynonymousRate() {
//        return 0;
//    }
//
//    protected void setupRelativeRates(double[] rates) {
//
//        double kappa = getKappa();
//        double omega = getOmega();
//        for (int i = 0; i < rateCount; i++) {
//            switch (rateMap[i]) {
//                case 0:
//                    rates[i] = 0.0;
//                    break;            // codon changes in more than one codon position
//                case 1:
//                    rates[i] = kappa;
//                    break;        // synonymous transition
//                case 2:
//                    rates[i] = 1.0;
//                    break;            // synonymous transversion
//                case 3:
//                    rates[i] = kappa * omega;
//                    break;// non-synonymous transition
//                case 4:
//                    rates[i] = omega;
//                    break;        // non-synonymous transversion
//            }
//        }
//    }
//
//
//}