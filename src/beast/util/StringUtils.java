package beast.util;


import beast.evolution.datatype.Codon;

import java.text.DecimalFormat;

public class StringUtils {

    // not use set because hard to get element given index
    public static <E> String concatenateToString(E[] array) {
        String concat = "";
        for (E a : array) {
            concat+= a.toString();
        }
        return concat;
    }

    public static Double[] roundDoubleArrays(Double[] numbers, int decimal) {
        Double[] numbersNew = new Double[numbers.length];
        DecimalFormat df = new DecimalFormat("#");
        df.setMaximumFractionDigits(decimal);
        for(int i = 0; i < numbers.length; i++ ){
            numbersNew[i] = Double.valueOf(df.format(numbers[i]));
        }
        return numbersNew;
    }

    public static double[] roundDoubleArrays(double[] numbers, int decimal) {
        double[] numbersNew = new double[numbers.length];
        DecimalFormat df = new DecimalFormat("#");
        df.setMaximumFractionDigits(decimal);
        for(int i = 0; i < numbers.length; i++ ){
            numbersNew[i] = Double.valueOf(df.format(numbers[i]));
        }
        return numbersNew;
    }

    /**
     * for print and debug view, matrixSize = nrOfStates^2
     * @param matrix transition probability matrix
     * @return
     */
    public static String get2DMatrixString(double[] matrix, Codon codon) {
        StringBuilder output = new StringBuilder();
        int dim = (int) Math.sqrt(matrix.length);

        if (matrix.length != dim * dim)
            throw new IllegalArgumentException("The array length " + matrix.length + " != " + dim + " * " + dim);

        if (codon != null) {
            // add header
            for (int j = 0; j < dim; j++)
                output.append("\t").append(codon.encodingToString(new int[]{j}));
            output.append("\n");
        }
        for (int i = 0; i < dim; i++) {
            if (codon != null) {
                // add row names
                output.append(codon.encodingToString(new int[]{i})).append("\t");
            }
            for (int j = 0; j < dim; j++) {
                output.append(matrix[i * dim + j]).append("\t");
            }
            output.append("\n");
        }
        return output.toString();
    }

    public static String get2DMatrixString(double[][] matrix, Codon codon) {
        StringBuilder output = new StringBuilder();
        int dim = matrix.length;

        if (matrix[0].length != dim)
            throw new IllegalArgumentException("The array length " + matrix.length + " != " + dim + " * " + dim);

        if (codon != null) {
            // add header
            for (int j = 0; j < dim; j++)
                output.append("\t").append(codon.encodingToString(new int[]{j}));
            output.append("\n");
        }
        for (int i = 0; i < dim; i++) {
            if (codon != null) {
                // add row names
                output.append(codon.encodingToString(new int[]{i})).append("\t");
            }
            for (int j = 0; j < dim; j++) {
                output.append(matrix[i][j]).append("\t");
            }
            output.append("\n");
        }
        return output.toString();
    }

}
