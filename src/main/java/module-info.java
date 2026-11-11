open module codonmodels {
    requires beast.base;
    requires beast.pkgmgmt;
    requires org.apache.commons.math4.legacy;

    // BEAUti support is optional at runtime
    requires static beast.fx;
    requires static javafx.controls;
    requires static javafx.swing;
    requires static java.desktop;

    exports codonmodels;
    exports codonmodels.app.beauti;
    exports codonmodels.evolution.alignment;
    exports codonmodels.evolution.datatype;
    exports codonmodels.evolution.likelihood;
    exports codonmodels.evolution.operators;
    exports codonmodels.evolution.tree;
    exports codonmodels.util;

    provides beast.base.core.BEASTInterface with
        codonmodels.app.beauti.CodonAlignmentProvider,
        codonmodels.ApproximateSubstModel,
        codonmodels.ApproxP_dist_CubicSpline,
        codonmodels.ApproxP_dist_Piecewise,
        codonmodels.ApproxP_dist_Polynomial,
        codonmodels.CodonFrequencies,
        codonmodels.CodonSubstitutionModel,
        codonmodels.evolution.alignment.CodonAlignment,
        codonmodels.evolution.alignment.CodonAlignmentFromAminoAcids,
        codonmodels.evolution.alignment.SimulatedCodonAlignment,
        codonmodels.evolution.datatype.Codon,
        codonmodels.evolution.likelihood.DataAugTreeLikelihood,
        codonmodels.evolution.likelihood.GenericDATreeLikelihood,
        codonmodels.evolution.operators.GibbsSamplingOperator,
        codonmodels.evolution.operators.ScaleOperatorGibbs,
        codonmodels.evolution.operators.UniformGibbs,
        codonmodels.evolution.tree.NodeStatesArray,
        codonmodels.evolution.tree.NodeStatesLogger,
        codonmodels.util.Report,
        codonmodels.M0Model,
        codonmodels.GeneralCodonFrequencies,
        codonmodels.GeneralCodonSubstitutionModel,
        codonmodels.EmpiricalCodonModel,
        codonmodels.EmpiricalAACodonModel;

    provides beast.base.evolution.datatype.DataType with
        codonmodels.evolution.datatype.Codon;
}
