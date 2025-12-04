import java.io.*;
import java.lang.Math;
import java.util.*;
import binary.*;

/**
*	Image Compression, DWT / IDWT

*	@author Ronan Merien <rmerien@hotmail.com>
*
*/
public class SimpleDWT {

	static public int COLS = -1;
	static public int ROWS = -1;
	static public int N = 8;
	static public double quality = 2;
	
	// Masques de Quantification pour 5 niveaux DWT (LL5, Details L5, Details L4, Details L3, Details L2, Details L1)
	// Les pas sont choisis pour être plus petits pour Y (luminance) et plus petits pour les basses fréquences (LL5).
	// L'index 0 est pour LL5; les index 1 à 5 sont pour les détails de L=5 (i=16) à L=1 (i=1).
	
	//static private float[] QUANT_MASK_Y    = {1.0f, 2.0f, 3.0f, 4.0f, 6.0f, 8.0f};
	//static private float[] QUANT_MASK_CRCB = {4.0f, 6.0f, 8.0f, 10.0f, 12.0f, 14.0f};
	
	//static private float[] QUANT_MASK_Y    = {2.0f, 4.0f, 6.0f, 8.0f, 10.0f, 12.0f};
	//static private float[] QUANT_MASK_CRCB = {8.0f, 10.0f, 12.0f, 14.0f, 18.0f, 22.0f};

	static private float[] QUANT_MASK_Y    = {4.0f, 6.0f, 8.0f, 12.0f, 16.0f, 20.0f};
	static private float[] QUANT_MASK_CRCB = {12.0f, 16.0f, 20.0f, 24.0f, 30.0f, 36.0f};
	
	/*
	LL : factor = 1.0
	LH, HL : factor = 0.7
	HH : factor = 0.6
	 */
	float bandFactor_LL = 1.0f;   
	float bandFactor_LH = 0.7f;  
	float bandFactor_HL = 0.7f;  
	float bandFactor_HH = 0.6f; // pour HH (moins quantifier)
	
	
	// x = i = row ; y = j = col

	static float[][] imageY;
	static float[][] imageCr;
	static float[][] imageCb;
	
	static double[][] dwtTemp;

	float[][] C = new float[N][N];
	float[][] Ct = new float[N][N];
	float[][] temp = new float[N][N];

	float[][] inputY = new float[N][N];
	float[][] inputCr = new float[N][N];
	float[][] inputCb = new float[N][N];

	float[][] outputY = new float[N][N];
	float[][] outputCr = new float[N][N];
	float[][] outputCb = new float[N][N];


	public SimpleDWT() {
	}

	protected int round_byte(double a) {
		if (a < 0) return 0;
		else if (a > 255) return 255;
		else return (int) Math.round(a);
	}


	// ----------------------------------------------------------------------------------------------------


	/**
	* Applique la quantification scalaire uniforme sur une région spécifiée de la matrice.
	* La matrice de coefficients est mise à jour "in place" (float -> int -> float).
	*/
	

    public void quantizeRegionInPlace(float[][] coefficients, int startRow, int endRow, int startCol, int endCol, float stepSize) {
	    
        for (int i = startRow; i < endRow; i++) {
            for (int j = startCol; j < endCol; j++) {
                float c = coefficients[i][j];

                int quantizedValue;
                quantizedValue = (int) Math.round(c / stepSize);
                
                // Mise à jour de la matrice float avec la valeur entière quantifiée
                coefficients[i][j] = (float) quantizedValue;
            }
        }
        
    }
   
    /**
     * Applique la déquantification scalaire uniforme sur une région spécifiée de la matrice.
     * Elle multiplie les coefficients entiers stockés (float) par le pas de quantification (stepSize).
     */
  
     public void dequantizeRegionInPlace(float[][] coefficients, int startRow, int endRow, int startCol, int endCol, float stepSize) {

         for (int i = startRow; i < endRow; i++) {
             for (int j = startCol; j < endCol; j++) {
                 // Le coefficient stocké est un entier (résultat de la quantification).

            	 int c_q = Math.round(coefficients[i][j]);
         
                 // Requantification: c_r = c_q * Delta_b
            	 float reconstructedValue = (float) c_q * stepSize;
                 
                 // Mise à jour de la matrice avec la valeur flottante reconstruite
                 coefficients[i][j] = reconstructedValue;
             }
         }
     }
    

	public void quantizeAllSubbands() {
	    
	    // LL5 (Approximation) - Facteur de résolution 32 (taille R/32 x C/32)
	    int LL_i = 32; 
	    int LL_R = ROWS / LL_i;
	    int LL_C = COLS / LL_i;
	    
	    // 1. LL5 Quantification (Top-Left-most band) - Index 0 du masque
	    quantizeRegionInPlace(imageY, 0, LL_R, 0, LL_C, QUANT_MASK_Y[0]*bandFactor_LL);
	    quantizeRegionInPlace(imageCr, 0, LL_R, 0, LL_C, QUANT_MASK_CRCB[0]*bandFactor_LL);
	    quantizeRegionInPlace(imageCb, 0, LL_R, 0, LL_C, QUANT_MASK_CRCB[0]*bandFactor_LL);
	    
	    // 2. Quantification des Bandes de Détail L=5 (i=16) à L=1 (i=1) - Index 1 à 5
	    
	    int[] resolution_i = {16, 8, 4, 2, 1}; // i pour les détails des niveaux L=5 à L=1
	    
	    for (int idx = 0; idx < 5; idx++) {
	        int i = resolution_i[idx]; 
	        int mask_idx = idx + 1; // Index 1 à 5
	        
	        // La région où se trouvent les détails LH_L, HL_L, HH_L
	        int R_block = ROWS / i;
	        int C_block = COLS / i;
	        
	        // Taille des sous-bandes de détail à ce niveau (R/2i x C/2i)
	        int R_sub = R_block / 2;
	        int C_sub = C_block / 2;
	        
			// --- Quantification des 3 bandes de détail (LH, HL, HH) ---
	        
	        // LH_L: Détail Horizontal (Top-Right)
	        quantizeRegionInPlace(imageY, 0, R_sub, C_sub, C_block, QUANT_MASK_Y[mask_idx]*bandFactor_LH);
	        quantizeRegionInPlace(imageCr, 0, R_sub, C_sub, C_block, QUANT_MASK_CRCB[mask_idx]*bandFactor_LH);
	        quantizeRegionInPlace(imageCb, 0, R_sub, C_sub, C_block, QUANT_MASK_CRCB[mask_idx]*bandFactor_LH);
	        
	        // HL_L: Détail Vertical (Bottom-Left)
	        quantizeRegionInPlace(imageY, R_sub, R_block, 0, C_sub, QUANT_MASK_Y[mask_idx]*bandFactor_HL);
	        quantizeRegionInPlace(imageCr, R_sub, R_block, 0, C_sub, QUANT_MASK_CRCB[mask_idx]*bandFactor_HL);
	        quantizeRegionInPlace(imageCb, R_sub, R_block, 0, C_sub, QUANT_MASK_CRCB[mask_idx]*bandFactor_HL);
	        
	        // HH_L: Détail Diagonal (Bottom-Right)
	        quantizeRegionInPlace(imageY, R_sub, R_block, C_sub, C_block, QUANT_MASK_Y[mask_idx]*bandFactor_HH);
	        quantizeRegionInPlace(imageCr, R_sub, R_block, C_sub, C_block, QUANT_MASK_CRCB[mask_idx]*bandFactor_HH);
	        quantizeRegionInPlace(imageCb, R_sub, R_block, C_sub, C_block, QUANT_MASK_CRCB[mask_idx]*bandFactor_HH);
	    }
	}
	

	/**
	    * Applique la déquantification avec le masque sur toutes les sous-bandes DWT (5 niveaux).
	    */
	    public void dequantizeAllSubbands() {
	        
	        // LL5 (Approximation) - Facteur de résolution 32
	        int LL_i = 32; 
	        int LL_R = ROWS / LL_i;
	        int LL_C = COLS / LL_i;
	        
	        // 1. LL5 Déquantification (Top-Left-most band) - Index 0 du masque
	        dequantizeRegionInPlace(imageY, 0, LL_R, 0, LL_C, QUANT_MASK_Y[0]*bandFactor_LL);
	        dequantizeRegionInPlace(imageCr, 0, LL_R, 0, LL_C, QUANT_MASK_CRCB[0]*bandFactor_LL);
	        dequantizeRegionInPlace(imageCb, 0, LL_R, 0, LL_C, QUANT_MASK_CRCB[0]*bandFactor_LL);

	        // 2. Déquantification des Bandes de Détail L=5 (i=16) à L=1 (i=1) - Index 1 à 5
	        
	        int[] resolution_i = {16, 8, 4, 2, 1}; // i pour les détails des niveaux L=5 à L=1
	        
	        for (int idx = 0; idx < 5; idx++) {
	            int i = resolution_i[idx]; 
	            int mask_idx = idx + 1; // Index 1 à 5
	            
	            // La région où se trouvent les détails LH_L, HL_L, HH_L
	            int R_block = ROWS / i;
	            int C_block = COLS / i;
	            
	            // Taille des sous-bandes de détail à ce niveau (R/2i x C/2i)
	            int R_sub = R_block / 2;
	            int C_sub = C_block / 2;
	            
	            // --- Déquantification des 3 bandes de détail (LH, HL, HH) ---
	            
	            // LH_L: Détail Horizontal (Top-Right)
	            dequantizeRegionInPlace(imageY, 0, R_sub, C_sub, C_block, QUANT_MASK_Y[mask_idx]*bandFactor_LH);
	            dequantizeRegionInPlace(imageCr, 0, R_sub, C_sub, C_block, QUANT_MASK_CRCB[mask_idx]*bandFactor_LH);
	            dequantizeRegionInPlace(imageCb, 0, R_sub, C_sub, C_block, QUANT_MASK_CRCB[mask_idx]*bandFactor_LH);
	            
	            // HL_L: Détail Vertical (Bottom-Left)
	            dequantizeRegionInPlace(imageY, R_sub, R_block, 0, C_sub, QUANT_MASK_Y[mask_idx]*bandFactor_HL);
	            dequantizeRegionInPlace(imageCr, R_sub, R_block, 0, C_sub, QUANT_MASK_CRCB[mask_idx]*bandFactor_HL);
	            dequantizeRegionInPlace(imageCb, R_sub, R_block, 0, C_sub, QUANT_MASK_CRCB[mask_idx]*bandFactor_HL);
	            
	            // HH_L: Détail Diagonal (Bottom-Right)
	            dequantizeRegionInPlace(imageY, R_sub, R_block, C_sub, C_block, QUANT_MASK_Y[mask_idx]*bandFactor_HH);
	            dequantizeRegionInPlace(imageCr, R_sub, R_block, C_sub, C_block, QUANT_MASK_CRCB[mask_idx]*bandFactor_HH);
	            dequantizeRegionInPlace(imageCb, R_sub, R_block, C_sub, C_block, QUANT_MASK_CRCB[mask_idx]*bandFactor_HH);
	        }
	    }
	    
	    
	// ----------------------------------------------------------------------------------------------------

	    /**
	    * l'image est suppose etre multiple de 64
	    * cad 16 (la resolution maximale pour la DWT) par 8 
	    * dans notre exemple : test.bmp 384 * 256 = (3*8*16) * (2*8*16)
	    * 
	    * LL : première zone
	    * LH : zone en haut-droite
	    * HL : zone en bas-gauche
	    * HH : zone en bas-droite
	    */
	    
    
	    public void forwardDWT_level(float[][] image, int rows, int cols) {
	        int row, col;
	        double sqrt2 = (float)Math.sqrt(2.0);

	        // --- Transformation Horizontale ---
	        for (row = 0; row < rows; row++) {
	            for (col = 0; col < cols/2; col++) {
	                // L = (A + B) / √2
	                dwtTemp[row][col] = (image[row][2*col] + image[row][2*col+1]) / sqrt2; 
	                // H = (A - B) / √2
	                dwtTemp[row][col + cols/2] = (image[row][2*col] - image[row][2*col+1]) / sqrt2;
	            }
	        }

	        for (row = 0; row < rows; row++) {
	            for (col = 0; col < cols; col++) {
	                image[row][col] = (float) dwtTemp[row][col];
	            }
	        }

	        // --- Transformation Verticale ---
	        for (row = 0; row < rows/2; row++) {
	            for (col = 0; col < cols; col++) {
	                // L = (A + B) / √2
	                dwtTemp[row][col] = (image[2*row][col] + image[2*row+1][col]) / sqrt2;
	                // H = (A - B) / √2
	                dwtTemp[row + rows/2][col] = (image[2*row][col] - image[2*row+1][col]) / sqrt2;
	            }
	        }

	        for (row = 0; row < rows; row++) {
	            for (col = 0; col < cols; col++) {
	                image[row][col] = (float) dwtTemp[row][col];
	            }
	        }
	    }

	   
	/**
	 * Applique l'IDWT de Haar sur la région de travail spécifiée (rows x cols).
	 * Reconstruit la région LL de niveau supérieur.
	 */

	    public void inverseDWT_level(float[][] image, int rows, int cols) {
	        int row, col;
	        double sqrt2 = (float)Math.sqrt(2.0);
	        
	        // --- IDWT Verticale ---
	        for (row = 0; row < rows/2; row++) {
	            for (col = 0; col < cols; col++) {
	                // A = (L + H) / √2
	                dwtTemp[2*row][col] = (image[row][col] + image[row + rows/2][col]) / sqrt2; 
	                // B = (L - H) / √2
	                dwtTemp[2*row+1][col] = (image[row][col] - image[row + rows/2][col]) / sqrt2;
	            }
	        }
	        
	        for (row = 0; row < rows; row++) {
	            for (col = 0; col < cols; col++) {
	                image[row][col] = (float) dwtTemp[row][col];
	            }
	        }
	        
	        // --- IDWT Horizontale ---
	        for (row = 0; row < rows; row++) {
	            for (col = 0; col < cols/2; col++) {
	                // A = (L + H) / √2
	                dwtTemp[row][2*col] = (image[row][col] + image[row][col + cols/2]) / sqrt2;
	                // B = (L - H) / √2
	                dwtTemp[row][2*col+1] = (image[row][col] - image[row][col + cols/2]) / sqrt2;
	            }
	        }
	        
	        for (row = 0; row < rows; row++) {
	            for (col = 0; col < cols; col++) {
	                image[row][col] = (float) dwtTemp[row][col];
	            }
	        }
	    }
	    	

	// ---------------------------------------------------------------------------------------

	public void compressFile(String inFile) throws Exception {
	//public void compressFile(String inFile, String outFile) throws Exception {
		try {
			int row, col;

			FileInputStream fis = new FileInputStream(inFile+".bmp");
			BinaryInputStream bis = new BinaryInputStream(new BufferedInputStream(fis));

			// Reading the input BMP imagefile
			int p;

			// BitmapFileHeader (14 bytes)

			if ((bis.readByte() != (byte) 'B') || (bis.readByte() != (byte) 'M')) throw new Exception("Not a Bitmap file"); // header = 'BM' (2 bytes)
			p = bis.readBit(32); // BMP file size (4 bytes)
			p = bis.readBit(64); // Reserved & Offset (8 bytes)

			// BitmapInfoHeader (40 bytes)

			p = bis.readBit(32); // info header size = 40 (4 bytes)
			COLS = bis.readBit(32); // width (4 bytes)
			ROWS = bis.readBit(32); // height (4 bytes)
			p = bis.readBit(16); // planes = 1 (2 bytes)
			p = bis.readBit(16); // bitcount = 24 (2 bytes)
			if (p != 24) throw new Exception("Not a 24bits Bitmap file");
			p = bis.readBit(32); // compression = 0 (4 bytes)
			p = bis.readBit(32); // image size (4 bytes)
			p = bis.readBit(32); // parameters (4 bytes)
			p = bis.readBit(32); // parameters (4 bytes)
			p = bis.readBit(32); // parameters (4 bytes)
			p = bis.readBit(32); // parameters (4 bytes)

			System.out.println("COLS="+COLS);
			System.out.println("ROWS="+ROWS);

			imageY = new float[ROWS][COLS];
			imageCr = new float[ROWS][COLS];
			imageCb = new float[ROWS][COLS];

			dwtTemp = new double[ROWS][COLS];

			int blue, green, red;
			float y, cb, cr;

			for (row = ROWS - 1; row >= 0; row--) {
				for (col = 0; col < COLS; col++) {

					blue = bis.readBit(8);
					green = bis.readBit(8);
					red = bis.readBit(8);

					y =  (float) (0.299 * red + 0.587 * green + 0.114 * blue);
					cb = (float) (-0.1687 * red - 0.3313 * green + 0.5 * blue); // cb = (float) (0.564*(blue - y));
					cr = (float) (0.5 * red - 0.4187 * green - 0.0813 * blue); // cr = (float) (0.713*(red - y));

					imageY[row][col] = y- 128.0f;
					imageCr[row][col] = cr  ;
					imageCb[row][col] = cb ;

				}
			}

			exportBMP_RGB("compress"+String.valueOf(0));
			
			int current_cols = COLS;
			int current_rows = ROWS;
		    
		    // DWT FORWARD (5 niveaux, de i=1 à i=16)
		    int i =1; // (i=1; i<=16; i = 2*i)
		    for (int k = 0; k < 5; k++) {
		    	
		    	 System.out.println("i="+i);
			     System.out.println("current_rows="+current_rows);
			     System.out.println("current_cols="+current_cols);
				
				// Appliquer la DWT uniquement à la région LL courante
			    System.out.println("forwardDWT_level");
		        forwardDWT_level(imageY, current_rows, current_cols);
		        forwardDWT_level(imageCr, current_rows, current_cols);
		        forwardDWT_level(imageCb, current_rows, current_cols);
		        
		        // La nouvelle bande LL (approximation) est divisée par 2
		        current_rows /= 2;
		        current_cols /= 2;
		        				
				exportBMP_RGB("compress"+String.valueOf(i));
				i = 2*i;
			}

		    
			quantizeAllSubbands();
			
			Huffman trans = new Huffman();
			trans.compressFile(inFile+".dwt");
/*			
			trans.expandFile(inFile+".dwt");
		
			dequantizeAllSubbands();
				
			//importBMP_YCrCb(inFile+".dwt"); // Charge les matrice imageY, imageCr, imageCb

			current_cols = COLS/32;
			current_rows = ROWS/32;
			
			// DWT INVERSE (5 niveaux, de i=16 à i=1)
			i=16; //i=16; i>=1; i = i/2
			while (i >= 1) {
		        // Réaugmenter la taille pour la reconstruction
		        current_rows *= 2;
		        current_cols *= 2;
		        
		    	System.out.println("i="+i);
			    System.out.println("current_rows="+current_rows);
			    System.out.println("current_cols="+current_cols);
		        
		        exportBMP_RGB("expand"+String.valueOf(i));
		        	
		        // Appliquer l'IDWT pour reconstruire la région LL précédente
		        System.out.println("inverseDWT_level");
		        inverseDWT_level(imageY, current_rows, current_cols);
		        inverseDWT_level(imageCr, current_rows, current_cols);
		        inverseDWT_level(imageCb, current_rows, current_cols);
		        
		        i=i/2;
		    }

		    exportBMP_RGB("expand0");
*/
		} catch (Exception e) {
			System.out.println(e);
		}


	}

	// ---------------------------------------------------------------------------------
	public void expandFile(String inFile) throws Exception {
	//public void expandFile(String inFile, String outFile) throws Exception {
		try {
			
			Huffman trans = new Huffman();
			trans.expandFile(inFile+".dwt");
			
			dequantizeAllSubbands();
			
			int current_cols = COLS/32;
			int current_rows = ROWS/32;
		    
			// DWT INVERSE (5 niveaux, de i=16 à i=1)
			int i=16; //i=16; i>=1; i = i/2
			while (i >= 1) {
		        // Réaugmenter la taille pour la reconstruction
		        current_rows *= 2;
		        current_cols *= 2;
		        
		        System.out.println("i="+i);
		        System.out.println("current_rows="+current_rows);
		        System.out.println("current_cols="+current_cols);
		        
		        exportBMP_RGB("expand"+String.valueOf(i));
		        
		        // Appliquer l'IDWT pour reconstruire la région LL précédente
		        System.out.println("inverseDWT_level");
		        inverseDWT_level(imageY, current_rows, current_cols);
		        inverseDWT_level(imageCr, current_rows, current_cols);
		        inverseDWT_level(imageCb, current_rows, current_cols);
		        
		        i=i/2;
		    }

		    exportBMP_RGB(inFile+".idwt.bmp");
			
		}
		catch (EOFException e) {}
		catch (IOException e) {System.out.println("expandFile :" + e.getStackTrace());}
	}

	/* Export a BMP image with YCrCb values
	 * 
	 */
	public void exportBMP_YCrCb(String name) throws Exception
	{
		int row, col;
		//int blue, green, red;
		int y, cb, cr;

		FileOutputStream fos = new FileOutputStream(name + ".bmp");
		BinaryOutputStream bos = new BinaryOutputStream(new BufferedOutputStream(fos));

		// 14 bytes
		bos.writeByte((byte) 'B');
		bos.writeByte((byte) 'M');
		bos.writeBit(COLS * ROWS * 3 + 54, 32); // BMP file length

		bos.writeBit(0, 32); // Reserved
		bos.writeByte((byte) 54); // Offset
		bos.writeByte((byte) 0);
		bos.writeByte((byte) 0);
		bos.writeByte((byte) 0);

		// 40 bytes
		bos.writeBit(40, 32); // 40 bytes
		bos.writeBit(COLS, 32); // largeur
		bos.writeBit(ROWS, 32); // hauteur
		bos.writeBit(1, 16);


		bos.writeBit(24, 16); // bits by pixel

		bos.writeBit(0, 32);

		bos.writeBit(COLS * ROWS * 3, 32); // image size

		bos.writeBit(0, 32);
		bos.writeBit(0, 32);
		bos.writeBit(0, 32);
		bos.writeBit(0, 32);

		for (row = ROWS - 1; row >= 0; row--) {
			for (col = 0; col < COLS; col++) {

				y = round_byte(imageY[row][col]+ 128.0f);
				cr = round_byte(imageCr[row][col]+ 128.0f);
				cb = round_byte(imageCb[row][col]+ 128.0f);

				bos.writeBit(y, 8);
				bos.writeBit(cr, 8);
				bos.writeBit(cb, 8);
			}
		}

		bos.flush();
		fos.close();

	}

	/* Export a BMP image with RGB values
	 * 
	 */
	public void exportBMP_RGB(String name) throws Exception
	{
		int row, col;
		int blue, green, red;
		double y, cb, cr;

		FileOutputStream fos = new FileOutputStream(name + ".bmp");
		BinaryOutputStream bos = new BinaryOutputStream(new BufferedOutputStream(fos));

		// 14 bytes
		bos.writeByte((byte) 'B');
		bos.writeByte((byte) 'M');
		bos.writeBit(COLS * ROWS * 3 + 54, 32); // BMP file length

		bos.writeBit(0, 32); // Reserved
		bos.writeByte((byte) 54); // Offset
		bos.writeByte((byte) 0);
		bos.writeByte((byte) 0);
		bos.writeByte((byte) 0);

		// 40 bytes
		bos.writeBit(40, 32); // 40 bytes
		bos.writeBit(COLS, 32); // largeur
		bos.writeBit(ROWS, 32); // hauteur
		bos.writeBit(1, 16);


		bos.writeBit(24, 16); // bits by pixel

		bos.writeBit(0, 32);

		bos.writeBit(COLS * ROWS * 3, 32); // image size

		bos.writeBit(0, 32);
		bos.writeBit(0, 32);
		bos.writeBit(0, 32);
		bos.writeBit(0, 32);

		for (row = ROWS - 1; row >= 0; row--) {
			for (col = 0; col < COLS; col++) {

				y = (float) imageY[row][col] + 128.0f;
				cr = (float) imageCr[row][col] ;
				cb = (float) imageCb[row][col] ;

				blue = round_byte(y + 1.772 * cb);
				green = round_byte(y - 0.34414 * cb - 0.71414 * cr);
				red = round_byte(y + 1.402 * cr);

				bos.writeBit(blue, 8);
				bos.writeBit(green, 8);
				bos.writeBit(red, 8);
			}
		}

		bos.flush();
		fos.close();

	}

	/* Load Matrix imageY, imageCr, imageCb
	 * from a BMP image
	 */
	 
	public void importBMP_YCrCb(String name) throws Exception
	{
		int row, col;

		FileInputStream fis = new FileInputStream(name + ".bmp");
		BinaryInputStream bis = new BinaryInputStream(new BufferedInputStream(fis));

		// Reading the input BMP imagefile
		int p;

		// BitmapFileHeader (14 bytes)

		if ((bis.readByte() != (byte) 'B') || (bis.readByte() != (byte) 'M')) throw new Exception("Not a Bitmap file"); // header = 'BM' (2 bytes)
		p = bis.readBit(32); // BMP file size (4 bytes)
		p = bis.readBit(64); // Reserved & Offset (8 bytes)

		// BitmapInfoHeader (40 bytes)

		p = bis.readBit(32); // info header size = 40 (4 bytes)
		COLS = bis.readBit(32); // width (4 bytes)
		ROWS = bis.readBit(32); // height (4 bytes)
		p = bis.readBit(16); // planes = 1 (2 bytes)
		p = bis.readBit(16); // bitcount = 24 (2 bytes)
		if (p != 24) throw new Exception("Not a 24bits Bitmap file");
		p = bis.readBit(32); // compression = 0 (4 bytes)
		p = bis.readBit(32); // image size (4 bytes)
		p = bis.readBit(32); // parameters (4 bytes)
		p = bis.readBit(32); // parameters (4 bytes)
		p = bis.readBit(32); // parameters (4 bytes)
		p = bis.readBit(32); // parameters (4 bytes)

		int blue, green, red; 
		float y, cb, cr;


		for (row = ROWS - 1; row >= 0; row--) {
			for (col = 0; col < COLS; col++) {
				
				blue = bis.readBit(8);
				green = bis.readBit(8);
				red = bis.readBit(8);

				y =  (float) (0.299 * red + 0.587 * green + 0.114 * blue);
				cb = (float) (-0.1687 * red - 0.3313 * green + 0.5 * blue) +128; // cb = (float) (0.564*(blue - y));
				cr = (float) (0.5 * red - 0.4187 * green - 0.0813 * blue) +128; // cr = (float) (0.713*(red - y));

				imageY[row][col] = (float) y-128.0f;
				imageCr[row][col] = (float) cr-128.0f;
				imageCb[row][col] = (float) cb-128.0f;
			}
		}

		bis.close();
		fis.close();	
	}

	// ---------------------------------------------------------------------------------------------

	static public void help()
	{
		System.out.println("SimpleDWT 1.1 powered by Ronan Merien <rmerien@hotmail.com>");
		System.out.println("[YCbCr 4:2:2, DWT 1/2/4/8/16 DCT 8*8, quantization]");
		System.out.println("Usage: java SimpleDWT [options] imagefile");
		System.out.println("-c		compress bmp_to_dwtfile");
		System.out.println("-e		expand dwt_to_bmpfile");
		System.out.println();
		//System.out.println("-quality n	quality factor between 1 and 100 (default 20)");
		System.out.println("Examples:  java SimpleDWT -c photo.bmp");
		System.out.println("           java SimpleDWT -e photo.dwt");
		System.out.println("");
	}

	static public void main(String args[])
	{
		String option = "";
		String imagefile = "";

		int argc = args.length;
		int argIndex = 0;
		String arg;

		while (argc > 0) {
			arg = args[argIndex];

			if (arg.startsWith("-")) {

				if (arg.equals("-quality")) {
					--argc;
					++argIndex;
					arg = args[argIndex];
					int p = Integer.parseInt(arg);
					if ((p >= 1) && (p <= 100)) {
						SimpleDWT.quality = p / 10;
					}
				}
				else option = arg;

			}
			else if (imagefile.equals("")) {
				imagefile = arg;
			}

			--argc;
			++argIndex;
		}

		if (imagefile.indexOf(".") != -1) {
			StringTokenizer f = new StringTokenizer(imagefile, ".");
			imagefile = f.nextToken();
		}

		SimpleDWT trans = new SimpleDWT();

		// compress bmp_to_dwtfile
		if (option.equals("-c"))
		{
			System.out.println("compress " + imagefile + ".bmp to " + imagefile + ".dwt");

			try {
				trans.compressFile(imagefile);
				//trans.compressFile(imagefile + ".bmp", imagefile + ".dwt");
			}
			catch (Exception e) { System.out.println(e); }

		}

		// expand dwt_to_bmpfile
		else if (option.equals("-e"))
		{
			System.out.println("expand " + imagefile + ".dwt to " + imagefile + ".idwt.bmp");

			try {
				trans.expandFile(imagefile);
				//trans.expandFile(imagefile + ".dwt", imagefile + ".idwt.bmp");
			}
			catch (Exception e) { System.out.println(e); }

		}

		// help instructions
		else {
			help();
		}
	}

}