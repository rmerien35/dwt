package binary;

import java.io.*;

/**
	BinaryInputStream
*/
      
public class BinaryInputStream extends FilterInputStream {

	protected int currentValue;
	protected byte reference;
    
	protected int bits = 0;
   
	protected boolean statistics = false;
	public int[] tabStatistics_1b = null;
	public int[] tabStatistics_8b = null;

	protected int ch = -1;
    
	public BinaryInputStream(InputStream in) {
		super(in);
      currentValue = 0x00;
      reference    = 0;
      bits         = 0;
      
      if (statistics) {
        tabStatistics_1b = new int[2];
        tabStatistics_1b[0] = 0;
        tabStatistics_1b[1] = 0;
               
        tabStatistics_8b = new int[256];
        for (int i=0; i<256; i++) tabStatistics_8b[i] =0;
      }
    }

    /**
	    Read byte as
      Single Binary Digit
    */
    public byte readBit() throws IOException {
      if (reference==0) {
        ch = in.read();
        if (ch < 0) throw new EOFException();
        currentValue = ch;
      }

      byte b = (byte) ( (currentValue >> reference) & 0x01 );
      			
      if (reference==7) {
      	/*
			if (statistics) {
            int p = -1;
            
            for (int i=0; i<8; i++) {
                p = ((int)currentValue >> i) & 0x01;
                tabStatistics_1b[p]++;
            }
                    
            p = ((int)currentValue) & 0xFF;
            tabStatistics_8b[p]++; 
        }
		  */

        // currentValue = 0x00;
        reference = 0;
      }
      else reference++;

      bits++;
      return b;
    }

	 public int readBit(int size) throws IOException {
        int p = 0x00;
        try {
            for (int i=0; i<size; i++) {
            	p = (((int) readBit()) << i) | p;
    	      }
    	      return p;
        }
        catch (IOException e) {
        		throw new IOException();
        }    
    }

    /**
	    Read byte as
      Unsigned 8 bits integer
    */
    public byte readByte() throws IOException {
      if (reference == 0) {
      	ch = in.read();
        if (ch < 0) throw new EOFException();
				byte b = (byte) ch;
				
        if (statistics) {
            int p = -1;
            
            for (int i=0; i<8; i++) {
                p = ((int)b >> i) & 0x01;
                tabStatistics_1b[p]++;
            }
                      
            p = ((int)b) & 0xFF;
            tabStatistics_8b[p]++; 
        }
		
        bits = bits + 8;
        return b;
      }
      else throw new IOException();
    }

    /**
     * Returns the number of bits read from this data input stream.
     */
    public int getBits() {
	    return bits;
    }

    /**
     * Returns the number of bytes read from this data input stream.
     */
    public int size() {
	    return (bits/8);
    }
    
    // --------------------------------------------------------------
    
    public void setStatistics(boolean s) {
        statistics = s;
    }
    
    public boolean isStatistics() {
        return statistics;
    }
    
    public void genStatistics() {
        if (statistics && (getBits() > 0)) {
            System.out.println("Generating statistics ...");
        
            System.out.println("Nb bits -> " + getBits());
            System.out.println("0 -> " + tabStatistics_1b[0]);
            System.out.println("1 -> " + tabStatistics_1b[1]);
                       
            System.out.println("Nb bytes -> " + size());
            for (int i=0; i<256; i++) {
                System.out.println(Binary.toHexaString((byte)i) + " -> " + tabStatistics_8b[i]);
            }
        }       
    }
    
}