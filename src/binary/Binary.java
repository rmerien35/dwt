package binary;

import java.util.*;

/**
	Binary Class
*/

public class Binary {

	public static final int DUMP_LINE_SIZE = 16;

	// ----------------------------------------------------------------------------

	static public String toSingleBinaryString(byte b) {
    String str = "";

    if ((b == 0x00)) str = "0";
    else             str = "1";

    return str;
  }

  // ----------------------------------------------------------------------------
  
  static public String toBinaryString(byte b) {
    String str = "";
    int p = 0x00;
    /*
    [0] =  ((b&0x01) == 1);
    [1] =  ((b&0x02) == 2);
    [2] =  ((b&0x04) == 4);
    [3] =  ((b&0x08) == 8);
    [4] =  ((b&0x0F) == 16);
    [5] =  ((b&0x20) == 32);
    [6] =  ((b&0x40) == 64);
    [7] =  ((b&0x80) == 128);
    */
    for (int i=0; i<8; i++) {
     	p = ((int)b >> i) & 0x01;
    	str = String.valueOf(p) + str;
    }
    return str;
  }
  
  static public String toBinaryString(int b, int size) {
    String str = "";
    int p = 0x00;
   
    for (int i=0; i<size; i++) {
     	p = (b >> i) & 0x01;
    	str = String.valueOf(p) + str;
    }
    return str;
  }
  
  // ----------------------------------------------------------------------------

  static public String toDecimalString(byte b) {
  	int i = ((int)b) & 0xFF;
  	
    return String.valueOf(i);
  }

  // ----------------------------------------------------------------------------

  static public String toHexaString(byte b) {
  	int i = ((int)b) & 0xFF;
    
  	if (i < 0x10)	return("0" + (Integer.toHexString(i)).toUpperCase());
    else 					return((Integer.toHexString(i)).toUpperCase());
	}

	// ----------------------------------------------------------------------------

  static public String toString(byte b) {
    return String.valueOf((char) b);
  }
  
  // ----------------------------------------------------------------------------
  
  
}