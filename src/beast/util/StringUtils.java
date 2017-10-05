package beast.util;


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

}
