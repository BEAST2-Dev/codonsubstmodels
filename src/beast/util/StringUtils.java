package beast.util;


public class StringUtils {

    // not use set because hard to get element given index
    public static <E> String concatenateToString(E[] array) {
        String concat = "";
        for (E a : array) {
            concat+= a.toString();
        }
        return concat;
    }


}
