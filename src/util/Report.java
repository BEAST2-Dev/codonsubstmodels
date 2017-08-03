package util;


import beast.core.Input;
import beast.core.Logger;

import java.util.ArrayList;
import java.util.List;

/**
 * Create a report.
 *
 * examples/testCodonLikelihood.xml
 *
 * @author Walter Xie
 */
public class Report extends beast.core.Runnable {

    public Input<List<Logger>> loggersInput = new Input<>("logger",
            "report parameters and trees.",
            new ArrayList<>());

    List<Logger> loggers;

    public void initAndValidate() {
        loggers = loggersInput.get();
    }


    public void run() throws Exception {
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