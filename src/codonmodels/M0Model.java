/*
 * GY94CodonModel.java
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

import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;

/**
 * Yang model of codon evolution
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Marc A. Suchard
 * @author Walter Xie
 */
@Citation(value= "Goldman N, & Yang Z (1994). A codon-based model of nucleotide substitution for\n" +
        "protein-coding DNA sequences. Molecular biology and evolution 11(5), 725-736.", year = 2014,
         DOI = "10.1093/oxfordjournals.molbev.a040153")
@Description("M0 codon substitution model, also called as GY94, published by Goldman and Yang 1994")
public class M0Model extends CodonSubstitutionModel {//implements Loggable {
    final public Input<RealParameter> kappaInput = new Input<>("kappa",
            "kappa parameter for transition-transversion rate ratio", Input.Validate.REQUIRED);

    final public Input<RealParameter> omegaInput = new Input<>("omega",
            "omega parameter to represent the nonsynonymous-synonymous rate ratio", Input.Validate.REQUIRED);


    public M0Model() {
        super();
    }

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        double value = omegaInput.get().getValue();
        if(value < 0) {
            throw new RuntimeException("Negative Omega parameter value " + value);
        }//END: negative check

        omegaInput.get().setBounds(Math.max(0.0, omegaInput.get().getLower()), omegaInput.get().getUpper());

        value = kappaInput.get().getValue();
        if(value < 0) {
            throw new RuntimeException("Negative kappa parameter value value " + value);
        }//END: negative check

        kappaInput.get().setBounds(Math.max(0.0, kappaInput.get().getLower()), kappaInput.get().getUpper());

    }

    @Override
    protected void setupRelativeRates() {
        double kappa = kappaInput.get().getValue();
        double omega = omegaInput.get().getValue();

//        this.synonymousRate = getSynonymousRate(kappa, omega);//not here

        // multiply pi in GeneralSubstitutionModel#setupRateMatrix()
        for (int i = 0; i < rateCount; i++) {
            switch (rateMap[i]) {
                case 0:
                    relativeRates[i] = 0.0;   // q_ij = 0
                    break;            // codon changes in more than one codon position
                case 1:
                    relativeRates[i] = kappa; // q_ij = pi_j * kappa
                    break;            // synonymous transition
                case 2:
                    relativeRates[i] = 1.0;   // q_ij = pi_j
                    break;            // synonymous transversion
                case 3:
                    relativeRates[i] = kappa * omega; // q_ij = pi_j * kappa * omega
                    break;            // non-synonymous transition
                case 4:
                    relativeRates[i] = omega; // q_ij = pi_j * omega
                    break;            // non-synonymous transversion
            }
        }
    }


//    @Override
//    public void init(PrintStream out) {
//        out.print("synonymousRate" + "\t");
//    }
//
//    @Override
//    public void log(int sample, PrintStream out) {
//        out.print(synonymousRate + "\t");
//    }
//
//    @Override
//    public void close(PrintStream out) {
// nothing to do
//    }
}