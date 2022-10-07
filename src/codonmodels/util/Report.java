package codonmodels.util;


import beast.base.inference.Distribution;
import beast.base.core.Input;
import beast.base.inference.Logger;
import beast.base.core.Log;

import java.util.ArrayList;
import java.util.List;

/**
 * Create a report.
 *
 * examples/testCodonLikelihood.xml
 *
 * @author Walter Xie
 */
public class Report extends beast.base.inference.Runnable {

    public Input<Distribution> distributionInput = new Input<>("distribution",
            "Distribution to calculate log likelihood.",
            Input.Validate.REQUIRED);

    public Input<List<Logger>> loggersInput = new Input<>("logger",
            "report parameters and trees.",
            new ArrayList<>());

    Distribution distribution;
    List<Logger> loggers;

//    TreeLikelihood likelihood;

    public void initAndValidate() {
        distribution = distributionInput.get();
        loggers = loggersInput.get();

//        likelihood = new TreeLikelihood();
    }


    public void run() throws Exception {
        final double logLikelihood = distribution.calculateLogP();

        Log.info.println("Distribution (" + distribution.getID() + ") log likelihood = " + logLikelihood);

//        if () {
//
//        }

        // Initialize loggers, default to System.out
        for (Logger logger : loggers)
            logger.init();

        for (Logger logger : loggers) {
            logger.log(0);
        }

        // Finalize loggers
        for (Logger logger: loggers)
            logger.close();

    }


}