import java.io.*;
import java.lang.Math;
import java.util.*;
import binary.*;

/**
*	Image Compression, DCT / IDCT
*	[YCbCr 4:2:2, DCT 8*8, quantization, DC differential, zigzag sequence, RLE & Huffman]

*	@author Ronan Merien <rmerien@hotmail.com>
*
*/
public class Huffman {


	static public int COLS = -1;
	static public int ROWS = -1;
	static public int N = 8;
	static public double quality = 2;
	
	double[][] C = new double[N][N];
	double[][] Ct = new double[N][N];
	double[][] temp = new double[N][N];

	int[][] inputY = new int[N][N];
	int[][] inputCr = new int[N][N];
	int[][] inputCb = new int[N][N];

	int[][] outputY = new int[N][N];
	int[][] outputCr = new int[N][N];
	int[][] outputCb = new int[N][N];

	protected int inputRLE = 0;
	protected int outputRLE = 0;

	protected int previousDC = 0;
	protected int nextDC = 0;
	protected boolean encodingDC = false;

	int code_ascii;

	static int maxBits = 20;
	int[] count = new int[maxBits];
	int[] next_code = new int[maxBits];

	BinaryTree[] tabStatistics_rle;
	BinaryTree[] tabStatistics_len;

	static int rle_max = N * N;
	static int len_max = 258; // 256 rle_mode, 257 out of range

	BinaryTree tree_len = new BinaryTree(-1, 0);
	BinaryTree tree_rle = new BinaryTree(-1, 0);

	@SuppressWarnings("rawtypes")
	Vector gen;

	public Huffman() {
	}

	protected int round_byte(double a) {
		if (a < 0) return 0;
		else if (a > 255) return 255;
		else return (int) Math.round(a);
	}


	// ----------------------------------------------------------------------------------------------------

	public void outputCode(BinaryOutputStream bos, int code) throws Exception {
		try {
			if ((code == 0) && !encodingDC) {
				outputRLE++;
				return;
			}

			while (outputRLE > 0) {
				tabStatistics_len[256].frequence++;
				gen.addElement(new MatchLength(256));
				// bos.writeBit((byte) 0x00); // RLE mode
				//System.out.print("0, RLE ");

				if (outputRLE <= rle_max) {
					tabStatistics_rle[outputRLE - 1].frequence++;
					gen.addElement(new MatchLength(outputRLE - 1));
					// bos.writeBit(outputRLE-1,6);
					//System.out.println(outputRLE);

					outputRLE = 0;
				} else {
					tabStatistics_rle[rle_max - 1].frequence++;
					gen.addElement(new MatchLength(rle_max - 1));
					// bos.writeBit(64-1,8);
					//System.out.println(64);

					outputRLE = outputRLE - 64;
				}
			}

			// bos.writeBit((byte) 0x01); // Huffman mode

			if (encodingDC) {
				// differential DC encoding
				encodingDC = false;
				nextDC = code - previousDC;
				previousDC = code;
				code = nextDC;
			}

			//System.out.print("1, Huffman ");

			if ((code >= -128) && (code <= 127)) {
				tabStatistics_len[code + 128].frequence++;
				gen.addElement(new MatchLength(code + 128));
				//System.out.println(code);
			} else if ((code < -128) || (code > 127)) {
				//System.out.println("out of range " + code);
				tabStatistics_len[257].frequence++;
				gen.addElement(
					new MatchLength(257, code + 128 * N, (int) (10 + N / 8)));
			} else {
				System.out.print("??? ");
			}

			// bos.writeBit(code + 1024, 11);

		} catch (Exception e) {
			System.out.println(e);
			throw e;
		}
	}

	// ---------------------------------------------------------------------------------
	public void outputCodeFlush(BinaryOutputStream bos) throws Exception {
		try {
			while (outputRLE > 0) {
				tabStatistics_len[256].frequence++;
				gen.addElement(new MatchLength(256));
				// bos.writeBit((byte) 0x00); // RLE mode
				//System.out.print("0, RLE ");

				if (outputRLE <= rle_max) {
					tabStatistics_rle[outputRLE - 1].frequence++;
					gen.addElement(new MatchLength(outputRLE - 1));
					// bos.writeBit(outputRLE-1,8);
					//System.out.println(outputRLE);

					outputRLE = 0;
				} else {
					tabStatistics_rle[rle_max - 1].frequence++;
					gen.addElement(new MatchLength(rle_max - 1));
					// bos.writeBit(rle_max-1,8);
					//System.out.println(rle_max);

					outputRLE = outputRLE - 64;
				}
			}

		} catch (Exception e) {
			System.out.println(e);
			throw e;
		}
	}
	
	// ---------------------------------------------------------------------------------
	
	/**
	 * Linéarise les coefficients DWT 2D en une séquence 1D 
	 * dans l'ordre de la plus basse fréquence (LL5) vers la plus haute fréquence (LH1, HL1, HH1).
	 *
	 * @param image Matrice de coefficients quantifiés (float, mais contenant des entiers)
	 * @return Une liste 1D de coefficients entiers.
	 */
	
	private List<Integer> linearizeDWT(float[][] image) {
	    List<Integer> coeffs = new ArrayList<>();
	    int rows = image.length;
	    int cols = image[0].length;
	    int max_level = 5; // Assumant 5 niveaux de décomposition (LL5 est le plus petit)

	    // --- A. Traitement du Bloc d'Approximation Final (LL5) - Inclus le DC ---
	    int ll_rows = rows / (1 << max_level); // Ex: R/32
	    int ll_cols = cols / (1 << max_level); // Ex: C/32

	    // 1. LL5 : Lecture en Raster Scan (pour les coefficients DC/LL5)
	    for (int r = 0; r < ll_rows; r++) {
	        for (int c = 0; c < ll_cols; c++) {
	            coeffs.add((int) Math.rint(image[r][c]));
	        }
	    }

	    // --- B. Traitement des Bandes de Détail (AC) : du niveau L=5 au niveau L=1 ---
	    // L'ordre peut être ajusté, ici L=5 (coarsest) à L=1 (finest)
	    
	    for (int L = max_level; L >= 1; L--) {
	        // Taille de la sous-bande au niveau L (e.g., L=5 -> R/32 x C/32)
	        int half_rows = rows / (1 << L);
	        int half_cols = cols / (1 << L);

	        // Taille de la région contenant LL(L-1) : R/2^(L-1) x C/2^(L-1)
	        int size_rows = half_rows * 2;
	        int size_cols = half_cols * 2;
	        
	        // 2. HL(L) band (Horizontal Low-pass / Vertical High-pass) - Détail Horizontal
	        // Région: [half_rows..size_rows-1] x [0..half_cols-1]
	        for (int r = half_rows; r < size_rows; r++) {
	            for (int c = 0; c < half_cols; c++) {
	                coeffs.add((int) Math.rint(image[r][c]));
	            }
	        }
	        
	        // 3. LH(L) band (Horizontal High-pass / Vertical Low-pass) - Détail Vertical
	        // Région: [0..half_rows-1] x [half_cols..size_cols-1]
	        for (int r = 0; r < half_rows; r++) {
	            for (int c = half_cols; c < size_cols; c++) {
	                coeffs.add((int) Math.rint(image[r][c]));
	            }
	        }
	        
	        // 4. HH(L) band (Diagonal Detail) - Détail Diagonal
	        // Région: [half_rows..size_rows-1] x [half_cols..size_cols-1]
	        for (int r = half_rows; r < size_rows; r++) {
	            for (int c = half_cols; c < size_cols; c++) {
	                coeffs.add((int) Math.rint(image[r][c]));
	            }
	        }
	    }
	    
	    return coeffs;
	}
	

	/**
	 * Reconstruit la matrice de coefficients 2D (float[][]) à partir 
	 * de la séquence 1D (List<Integer>), en suivant l'ordre inverse de linearizeDWT.
	 * * @param coefficients Liste 1D des coefficients entiers décodés.
	 * @param image Matrice de destination float[][] (imageY, imageCr, ou imageCb).
	 */
	private void unlinearizeDWT(List<Integer> coefficients, float[][] image) {
	    if (image.length == 0 || coefficients.isEmpty()) return;
	    
	    int rows = image.length;
	    int cols = image[0].length;
	    int max_level = 5; // Supposant 5 niveaux de décomposition
	    int index = 0;     // Index de lecture dans la liste 1D
	    
	    // --- A. Reconstruction du Bloc d'Approximation Final (LL5) ---
	    // Ordre : Raster Scan sur R/32 x C/32
	    int ll_rows = rows / (1 << max_level);
	    int ll_cols = cols / (1 << max_level);

	    for (int r = 0; r < ll_rows; r++) {
	        for (int c = 0; c < ll_cols; c++) {
	            // Lecture de l'entier et écriture dans le tableau float
	            image[r][c] = (float) coefficients.get(index++);
	        }
	    }

	    // --- B. Reconstruction des Bandes de Détail (AC) : du niveau L=5 au niveau L=1 ---
	    // L'ordre est L=5 à L=1, avec le même ordre d'itération que dans linearizeDWT.
	    
	    for (int L = max_level; L >= 1; L--) {
	        int half_rows = rows / (1 << L);
	        int half_cols = cols / (1 << L);

	        int size_rows = half_rows * 2;
	        int size_cols = half_cols * 2;
	        
	        // 2. HL(L) band - Détail Horizontal
	        // Région: [half_rows..size_rows-1] x [0..half_cols-1]
	        for (int r = half_rows; r < size_rows; r++) {
	            for (int c = 0; c < half_cols; c++) {
	                image[r][c] = (float) coefficients.get(index++);
	            }
	        }
	        
	        // 3. LH(L) band - Détail Vertical
	        // Région: [0..half_rows-1] x [half_cols..size_cols-1]
	        for (int r = 0; r < half_rows; r++) {
	            for (int c = half_cols; c < size_cols; c++) {
	                image[r][c] = (float) coefficients.get(index++);
	            }
	        }
	        
	        // 4. HH(L) band - Détail Diagonal
	        // Région: [half_rows..size_rows-1] x [half_cols..size_cols-1]
	        for (int r = half_rows; r < size_rows; r++) {
	            for (int c = half_cols; c < size_cols; c++) {
	                image[r][c] = (float) coefficients.get(index++);
	            }
	        }
	    }
	}
	
	// ---------------------------------------------------------------------------------
	public int inputCode(BinaryInputStream bis) throws Exception {
		int code = 0;

		BinaryTree node = null;
		boolean find = false;
		int value = -1;
		byte b = 0;

		try {
			if (inputRLE > 0) {
				inputRLE--;
				return 0;
			}

			find = false;
			node = tree_len;

			while (!find) {
				b = bis.readBit();
				// System.out.print(Binary.toSingleBinaryString(b));

				if (b == 0x00)
					node = node.node_0;
				else
					node = node.node_1;

				if (node == null)
					throw new Exception();

				value = node.codeAscii;
				if (value != -1)
					find = true;
			}
			code = value;
			//System.out.println("value = " + value);

			if (code == 256) {
				// RLE mode
				//System.out.print("0, RLE ");

				find = false;
				node = tree_rle;

				while (!find) {
					b = bis.readBit();
					// System.out.print(Binary.toSingleBinaryString(b));

					if (b == 0x00)
						node = node.node_0;
					else
						node = node.node_1;

					if (node == null)
						throw new Exception();

					value = node.codeAscii;
					if (value != -1)
						find = true;
				}
				inputRLE = value;
				//System.out.println(inputRLE);

				return 0;
			}

			// 1 : Huffman mode

			if (code == 257) {
				//System.out.println("out of range 257");
				code = bis.readBit((int) (10 + N / 8)) - 128 * N;
			} else
				code = code - 128;


			if (encodingDC) {
				// differential DC encoding
				encodingDC = false;
				nextDC = code + previousDC;
				code = nextDC;
				previousDC = nextDC;
			}

			return code;

		} catch (Exception e) {
			System.out.println("inputCode " + e);
			throw e;
		}
	}

	// ---------------------------------------------------------------------------------
	@SuppressWarnings("rawtypes")
	public void compressFile(String inFile, String outFile) throws Exception {
		try {
			int row, col, i, j, k;
			
			FileInputStream fis = new FileInputStream(inFile);
			BinaryInputStream bis = new BinaryInputStream(new BufferedInputStream(fis));

			FileOutputStream fos = new FileOutputStream(outFile);
			BinaryOutputStream bos = new BinaryOutputStream(new BufferedOutputStream(fos));

			tabStatistics_len = new BinaryTree[len_max];
			for (i = 0; i < len_max; i++) {
				tabStatistics_len[i] = new BinaryTree(i, 0);
			}

			tabStatistics_rle = new BinaryTree[rle_max];
			for (i = 0; i < rle_max; i++) {
				tabStatistics_rle[i] = new BinaryTree(i, 0);
			}

			gen = new Vector(fis.available(), 100000);


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


			bos.writeBit(COLS, 16);
			bos.writeBit(ROWS, 16);

			int blue, green, red, y, cb, cr;
		
			
			for (row = ROWS - 1; row >= 0; row--) {
				for (col = 0; col < COLS; col++) {

					blue = bis.readBit(8);
					green = bis.readBit(8);
					red = bis.readBit(8);

					y = (int) (0.299 * red + 0.587 * green + 0.114 * blue);
					cb = (int) (-0.1687 * red - 0.3313 * green + 0.5 * blue); // cb = (int) (0.564*(blue - y));
					cr = (int) (0.5 * red - 0.41874 * green - 0.08130 * blue); // cr = (int) (0.713*(red - y));
					
					SimpleDWT.imageY[row][col] = y - 128;
					SimpleDWT.imageCr[row][col] = cr;
					SimpleDWT.imageCb[row][col] = cb;

				}
			}
			
			List<Integer> coefficientsY = linearizeDWT(SimpleDWT.imageY); 
			System.out.println("Encodage coefficientsY");
			for (i = 0; i < coefficientsY.size(); i++) {
				if (i == 0) encodingDC = true;
				int val_lue = coefficientsY.get(i);
				outputCode(bos, val_lue);
			}

			System.out.println("Encodage coefficientsCr");
			List<Integer> coefficientsCr = linearizeDWT(SimpleDWT.imageCr); 
			for (i = 0; i < coefficientsCr.size(); i++) {
				if (i == 0) encodingDC = true;
				int val_lue = coefficientsCr.get(i);
				outputCode(bos, val_lue);
			}

			System.out.println("Encodage coefficientsCb");
			List<Integer> coefficientsCb = linearizeDWT(SimpleDWT.imageCb); 
			for (i = 0; i < coefficientsCb.size(); i++) {
				if (i == 0) encodingDC = true;
				int val_lue = coefficientsCb.get(i);
				outputCode(bos, val_lue);
			}
			
			outputCodeFlush(bos);

			bis.close();

			// build the vector
			Vector list = new Vector();
			for (i = 0; i < len_max; i++) {
				if (tabStatistics_len[i].frequence > 0) {
					list.addElement(tabStatistics_len[i]);
				}
			}
			sort(list, 0, list.size() - 1);

			Vector list_rle = new Vector();
			for (i = 0; i < rle_max; i++) {
				if (tabStatistics_rle[i].frequence > 0) {
					list_rle.addElement(tabStatistics_rle[i]);
				}
			}
			sort(list_rle, 0, list_rle.size() - 1);

			// build the huffman tree
			for (i = 0; i < maxBits; i++) {
				count[i] = 0;
				next_code[i] = 0;
			}

			int size = list.size();

			BinaryTree node_0 = null;
			BinaryTree node_1 = null;
			BinaryTree tree = null;

			for (k = 0; k < size - 1; k++) {

				node_0 = (BinaryTree) list.elementAt(0);
				node_1 = (BinaryTree) list.elementAt(1);

				int sum = node_0.frequence + node_1.frequence;

				tree = new BinaryTree(-1, sum, node_0, node_1);

				list.removeElementAt(0);
				list.removeElementAt(0);

				list.insertElementAt(tree, 0);
				sort(list, 0, list.size() - 1);
			}

			convert_tree_to_code(tree, 0);


			int compressCode = 0;
			count[0] = 0;
			for (int bits = 1; bits < maxBits; bits++) {
				compressCode = (compressCode + count[bits - 1]) << 1;
				next_code[bits] = compressCode;
			}

			BinaryTree node = null;
			for (i = 0; i < len_max; i++) {
				node = tabStatistics_len[i];

				if (node != null) {
					int len = node.nbBits;
					if (len != 0) {
						node.compressCode = reverse(next_code[len], len);
						next_code[len]++;
					}
				}
			}

			
			for (i = 0; i < maxBits; i++) {
				count[i] = 0;
				next_code[i] = 0;
			}

			node_0 = null;
			node_1 = null;
			tree = null;

			size = list_rle.size();

			for (k = 0; k < size - 1; k++) {

				node_0 = (BinaryTree) list_rle.elementAt(0);
				node_1 = (BinaryTree) list_rle.elementAt(1);

				int sum = node_0.frequence + node_1.frequence;

				tree = new BinaryTree(-1, sum, node_0, node_1);

				list_rle.removeElementAt(0);
				list_rle.removeElementAt(0);

				list_rle.insertElementAt(tree, 0);
				sort(list_rle, 0, list_rle.size() - 1);
			}

			if (tree != null)
				convert_tree_to_code(tree, 0);

			compressCode = 0;
			count[0] = 0;
			for (int bits = 1; bits < maxBits; bits++) {
				compressCode = (compressCode + count[bits - 1]) << 1;
				next_code[bits] = compressCode;
			}

			for (i = 0; i < rle_max; i++) {
				node = tabStatistics_rle[i];

				if (node != null) {
					int len = node.nbBits;
					if (len != 0) {
						node.compressCode = reverse(next_code[len], len);
						next_code[len]++;
					}
				}
			}

			
			boolean mode_rle = false;
			int compteur_rle = 0;

			// second pass : generating huffman codes

			// next, write huffman statistics for length into the output file
			for (i = 0; i < len_max; i++) {
				if (tabStatistics_len[i].frequence > 0) {
					if (mode_rle) {
						bos.writeBit(0x00, 4);
						bos.writeBit(compteur_rle - 1, 8);
						mode_rle = false;
						compteur_rle = 0;
					}

					BinaryTree element = tabStatistics_len[i];
					if (element.nbBits < 15) {
						bos.writeBit(element.nbBits, 4);
					} else if (
						(element.nbBits >= 15) && (element.nbBits < 18)) {
						bos.writeBit(15, 4);
						bos.writeBit(element.nbBits - 15, 2);
						// mode (0-> 15, 1-> 16, 2-> 17, 3-> 18+ bits)
					} else if (element.nbBits < maxBits) {
						bos.writeBit(15, 4);
						bos.writeBit(3, 2);
						bos.writeBit(element.nbBits - 18, 2);
						// mode (0-> 18, 1-> 19 bits)
					}

				} else {
					mode_rle = true;
					compteur_rle++;
				}
			}

			if (mode_rle) {
				bos.writeBit(0x00, 4);
				bos.writeBit(compteur_rle - 1, 8);

				mode_rle = false;
				compteur_rle = 0;
			}

			// next, write huffman statistics for rle encoding into the output file
			for (i = 0; i < rle_max; i++) {
				if (tabStatistics_rle[i].frequence > 0) {
					if (mode_rle) {
						bos.writeBit(0x00, 4);
						bos.writeBit(compteur_rle - 1, 5);

						mode_rle = false;
						compteur_rle = 0;
					}

					BinaryTree element = tabStatistics_rle[i];
					bos.writeBit(element.nbBits, 4);

				} else {
					mode_rle = true;
					compteur_rle++;
				}
			}

			if (mode_rle) {
				bos.writeBit(0x00, 4);
				bos.writeBit(compteur_rle - 1, 5);

				mode_rle = false;
				compteur_rle = 0;
			}
			
		
			MatchLength len = null;

			BinaryTree node_len = null;
			BinaryTree node_rle = null;

			// next, write contains of the "gen" vector into the output file
			for (Enumeration e = gen.elements(); e.hasMoreElements();) {
				len = (MatchLength) e.nextElement();
				
				if (len.value == 256) {
								
					node_len = tabStatistics_len[len.value];
					bos.writeBit(node_len.compressCode, node_len.nbBits);

					len = (MatchLength) e.nextElement();
					//System.out.println("rle = " + len.value);
					node_rle = tabStatistics_rle[len.value];
					bos.writeBit(node_rle.compressCode, node_rle.nbBits);
				} else {
					
					node_len = tabStatistics_len[len.value];
					bos.writeBit(node_len.compressCode, node_len.nbBits);

					// even, generating an extra parameter
					if (len.nbExtraBits > 0) {
						bos.writeBit(len.extraValue, len.nbExtraBits);

					}
				}
			}
			
			//System.out.println("compress tabStatistics_rle="+tabStatistics_rle.length);
			
			bos.writeEOF();

			bos.flush();
		} catch (EOFException e) {
			System.out.println(e);
		} catch (IOException e) {
			System.out.println("compressFile :" + e);
		}
	}

	
	public void encodeFile(String inFile, String outFile) throws Exception {
		try {
			int row, col, i, j, k;
			
			FileInputStream fis = new FileInputStream(inFile);
			BinaryInputStream bis = new BinaryInputStream(new BufferedInputStream(fis));

			FileOutputStream fos = new FileOutputStream(outFile);
			BinaryOutputStream bos = new BinaryOutputStream(new BufferedOutputStream(fos));

			tabStatistics_len = new BinaryTree[len_max];
			for (i = 0; i < len_max; i++) {
				tabStatistics_len[i] = new BinaryTree(i, 0);
			}

			tabStatistics_rle = new BinaryTree[rle_max];
			for (i = 0; i < rle_max; i++) {
				tabStatistics_rle[i] = new BinaryTree(i, 0);
			}

			gen = new Vector(fis.available(), 100000);


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


			bos.writeBit(COLS, 16);
			bos.writeBit(ROWS, 16);

			SimpleDWT.imageY = new float[ROWS][COLS];
			SimpleDWT.imageCr = new float[ROWS][COLS];
			SimpleDWT.imageCb = new float[ROWS][COLS];

			int blue, green, red, y, cb, cr;
		
			
			for (row = ROWS - 1; row >= 0; row--) {
				for (col = 0; col < COLS; col++) {

					blue = bis.readBit(8);
					green = bis.readBit(8);
					red = bis.readBit(8);

					y = (int) (0.299 * red + 0.587 * green + 0.114 * blue);
					cb = (int) (-0.1687 * red - 0.3313 * green + 0.5 * blue); // cb = (int) (0.564*(blue - y));
					cr = (int) (0.5 * red - 0.41874 * green - 0.08130 * blue); // cr = (int) (0.713*(red - y));
					
					SimpleDWT.imageY[row][col] = y - 128;
					SimpleDWT.imageCr[row][col] = cr;
					SimpleDWT.imageCb[row][col] = cb;

				}
			}
			
			List<Integer> coefficientsY = linearizeDWT(SimpleDWT.imageY); 
			System.out.println("Encodage coefficientsY");
			for (i = 0; i < coefficientsY.size(); i++) {
				if (i == 0) encodingDC = true;
				Integer val_lue = coefficientsY.get(i);
				//outputCode(bos, val_lue);
				bos.writeBit(coefficientsY.get(i).byteValue());
			}

			System.out.println("Encodage coefficientsCr");
			List<Integer> coefficientsCr = linearizeDWT(SimpleDWT.imageCr); 
			for (i = 0; i < coefficientsCr.size(); i++) {
				if (i == 0) encodingDC = true;
				Integer val_lue = coefficientsCr.get(i);
				//outputCode(bos, val_lue);
				bos.writeBit(val_lue.byteValue());
			}

			System.out.println("Encodage coefficientsCb");
			List<Integer> coefficientsCb = linearizeDWT(SimpleDWT.imageCb); 
			for (i = 0; i < coefficientsCb.size(); i++) {
				if (i == 0) encodingDC = true;
				Integer val_lue = coefficientsCb.get(i);
				//outputCode(bos, val_lue);
				bos.writeBit(val_lue.byteValue());
			}
			
			outputCodeFlush(bos);

			bis.close();		
			bos.writeEOF();

			bos.flush();
		} catch (EOFException e) {
			System.out.println(e);
		} catch (IOException e) {
			System.out.println("compressFile :" + e);
		}
	}
	
	// ---------------------------------------------------------------------------------
	public void expandFile(String inFile, String outFile) throws Exception {
		try {
			int row, col, i, j, k;

			FileInputStream fis = new FileInputStream(inFile);
			BinaryInputStream bis =	new BinaryInputStream(new BufferedInputStream(fis));

			FileOutputStream fos = new FileOutputStream(outFile);
			BinaryOutputStream bos = new BinaryOutputStream(new BufferedOutputStream(fos));

			byte b = 0;

			COLS = bis.readBit(16);	
			ROWS = bis.readBit(16);
			
			//System.out.println("COLS="+COLS);
			//System.out.println("ROWS="+ROWS);

			// get statistics from input stream

			tabStatistics_len = new BinaryTree[len_max];
			for (i = 0; i < len_max; i++) {
				tabStatistics_len[i] = new BinaryTree(i, 0);
			}

			tabStatistics_rle = new BinaryTree[rle_max];
			for (i = 0; i < rle_max; i++) {
				tabStatistics_rle[i] = new BinaryTree(i, 0);
			}

			for (i = 0; i < maxBits; i++) {
				count[i] = 0;
				next_code[i] = 0;
			}

			for (i = 0; i < len_max; i++) {
				int value = bis.readBit(4);
				
				if (value == 0x00) {
					i = i + bis.readBit(8);
				} else {
					if (value == 15)
						value = value + bis.readBit(2);
					if (value == 18)
						value = value + bis.readBit(2);
					tabStatistics_len[i].nbBits = value;
					count[value]++;
				}
			}
			
			int code = 0;
			count[0] = 0;
			for (int bits = 1; bits < maxBits; bits++) {
				code = (code + count[bits - 1]) << 1;
				next_code[bits] = code;
			}

			for (i = 0; i < len_max; i++) {
				int len = tabStatistics_len[i].nbBits;
				if (len != 0) {
					tabStatistics_len[i].compressCode = next_code[len];
					convert_code_to_tree(tree_len, tabStatistics_len[i], len);
					next_code[len]++;
				}
			}

			
			// --------------------------------------------------------

			for (i = 0; i < maxBits; i++) {
				count[i] = 0;
				next_code[i] = 0;
			}

			for (i = 0; i < rle_max; i++) {
				int value = bis.readBit(4);

				if (value == 0x00) {
					i = i + bis.readBit(5);
				} else {
					tabStatistics_rle[i].nbBits = value;
					count[value]++;
				}
			}
			
			code = 0;
			count[0] = 0;
			for (int bits = 1; bits < maxBits; bits++) {
				code = (code + count[bits - 1]) << 1;
				next_code[bits] = code;
			}
		
			for (i = 0; i < rle_max; i++) {
				int len = tabStatistics_rle[i].nbBits;
				if (len != 0) {
					tabStatistics_rle[i].compressCode = next_code[len];
					convert_code_to_tree(tree_rle, tabStatistics_rle[i], len);
					next_code[len]++;
				}
			}
			
			//System.out.println("expand tabStatistics_rle="+tabStatistics_rle.length);

			SimpleDWT.imageY = new float[ROWS][COLS];
			SimpleDWT.imageCr = new float[ROWS][COLS];
			SimpleDWT.imageCb = new float[ROWS][COLS];

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


			int blue, green, red;
			float y, cb, cr;
			
			List<Integer> coefficientsY = new ArrayList();
			for(k=0; k<COLS * ROWS; k++ ) {
				if (k == 0) encodingDC = true;
				
				int val_lue = inputCode(bis);
				coefficientsY.add(k,val_lue);
			}
			unlinearizeDWT(coefficientsY, SimpleDWT.imageY);
			
			List<Integer> coefficientsCr = new ArrayList();
			for(k=0; k<COLS * ROWS; k++ ) {
				if (k == 0) encodingDC = true;
				int val_lue = inputCode(bis);
				coefficientsCr.add(k,val_lue);
			}
			unlinearizeDWT(coefficientsCr, SimpleDWT.imageCr);
			
			List<Integer> coefficientsCb = new ArrayList();
			for(k=0; k<COLS * ROWS; k++ ) {
				if (k == 0) encodingDC = true;
				int val_lue = inputCode(bis);
				coefficientsCb.add(k,val_lue);
			}
			unlinearizeDWT(coefficientsCb, SimpleDWT.imageCb);
			
	
			for (row = ROWS - 1; row >= 0; row--) {
				for (col = 0; col < COLS; col++) {

					y = SimpleDWT.imageY[row][col] + 128;
					cr = SimpleDWT.imageCr[row][col];
					cb = SimpleDWT.imageCb[row][col];

					blue = round_byte(y + 1.773 * cb);
					green = round_byte(y - 0.34414 * cb - 0.71414 * cr);
					red = round_byte(y + 1.402 * cr);

					bos.writeBit(blue, 8);
					bos.writeBit(green, 8);
					bos.writeBit(red, 8);
				}
			}

			bos.flush();
			fos.close();
			fis.close();
		} 
		catch (EOFException e) {}
		catch (IOException e) {System.out.println("expandFile :" + e);}

		//exportBMP_RGB("RGB_"+outFile);
	}

	
	public void decodeFile(String inFile, String outFile) throws Exception {
		try {
			int row, col, i, j, k;

			FileInputStream fis = new FileInputStream(inFile);
			BinaryInputStream bis =	new BinaryInputStream(new BufferedInputStream(fis));

			FileOutputStream fos = new FileOutputStream(outFile);
			BinaryOutputStream bos = new BinaryOutputStream(new BufferedOutputStream(fos));

			byte b = 0;

			COLS = bis.readBit(16);	
			ROWS = bis.readBit(16);
			
			//System.out.println("COLS="+COLS);
			//System.out.println("ROWS="+ROWS);

			SimpleDWT.imageY = new float[ROWS][COLS];
			SimpleDWT.imageCr = new float[ROWS][COLS];
			SimpleDWT.imageCb = new float[ROWS][COLS];

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


			int blue, green, red;
			float y, cb, cr;
			
			List<Integer> coefficientsY = new ArrayList();
			for(k=0; k<COLS * ROWS; k++ ) {
				//if (k == 0) encodingDC = true;
				//coefficientsY.add(k,inputCode(bis));
				int val_lue = bis.readBit();
				coefficientsY.add(k,val_lue);
			}
			unlinearizeDWT(coefficientsY, SimpleDWT.imageY);
			
			List<Integer> coefficientsCr = new ArrayList();
			for(k=0; k<COLS * ROWS; k++ ) {
				//if (k == 0) encodingDC = true;
				//coefficientsCr.add(k,inputCode(bis));
				int val_lue = bis.readBit();
				coefficientsCr.add(k,val_lue);
			}
			unlinearizeDWT(coefficientsCr, SimpleDWT.imageCr);
			
			List<Integer> coefficientsCb = new ArrayList();
			for(k=0; k<COLS * ROWS; k++ ) {
				//if (k == 0) encodingDC = true;
				//coefficientsCb.add(k,inputCode(bis));
				int val_lue = bis.readBit();
				coefficientsCb.add(k,val_lue);
			}
			unlinearizeDWT(coefficientsCb, SimpleDWT.imageCb);
			
			
	
			for (row = ROWS - 1; row >= 0; row--) {
				for (col = 0; col < COLS; col++) {

					y = SimpleDWT.imageY[row][col] + 128;
					cr = SimpleDWT.imageCr[row][col];
					cb = SimpleDWT.imageCb[row][col];

					blue = round_byte(y + 1.773 * cb);
					green = round_byte(y - 0.34414 * cb - 0.71414 * cr);
					red = round_byte(y + 1.402 * cr);

					bos.writeBit(blue, 8);
					bos.writeBit(green, 8);
					bos.writeBit(red, 8);
				}
			}

			bos.flush();
			fos.close();
			fis.close();
		} 
		catch (EOFException e) {}
		catch (IOException e) {System.out.println("expandFile :" + e);}

		//exportBMP_RGB("RGB_"+outFile);
	}
	
	// ---------------------------------------------------------------------------------


// -------------------------------------------------------------------------------------------
	public int reverse(int value, int size) 
	{
		int result = 0x00;
		for (int i = 0; i < size; i++) {
			result = (result << 1) | ((value >> i) & 0x01);
		}
		return result;
	}


	public void convert_tree_to_code(BinaryTree tree, int nbBits) 
	{
		if (tree.isLeafNode()) {
			tree.nbBits = nbBits;
			count[nbBits]++;
			return;
		}

		nbBits++;

		convert_tree_to_code(tree.node_0, nbBits);
		convert_tree_to_code(tree.node_1, nbBits);
	}


	public BinaryTree convert_code_to_tree(BinaryTree tree, BinaryTree node, int nbBits) 
	{
		if (nbBits == 0) {
			return node;
		}

		if (tree == null) tree = new BinaryTree(-1, 0);

		int currentBit = ((node.compressCode >> (nbBits - 1)) & 0x01);
		nbBits--;

		if (currentBit == 0x00)	tree.node_0 = convert_code_to_tree(tree.node_0, node, nbBits);
		else tree.node_1 = convert_code_to_tree(tree.node_1, node, nbBits);

		return tree;
	}


	static public void sort(Vector list, int low, int high) 
	{
		if (!list.isEmpty()) {
			BinaryTree previous = null;
			BinaryTree current = null;

			int l = low;
			int h = high;
			int pivot =
				((BinaryTree) list.elementAt((low + high) / 2)).frequence;

			while (l <= h) {
				while (((BinaryTree) list.elementAt(l)).frequence < pivot) l++;
				while (pivot < ((BinaryTree) list.elementAt(h)).frequence) h--;

				if (l <= h) {

					previous = (BinaryTree) list.elementAt(l);
					current = (BinaryTree) list.elementAt(h);

					list.removeElementAt(l);
					list.insertElementAt(current, l);

					list.removeElementAt(h);
					list.insertElementAt(previous, h);

					l++;
					h--;
				}
			}

			// Sort the low part of the list :
			if (low < h)
				sort(list, low, h);

			// Sort the high part of the list :
			if (l < high)
				sort(list, l, high);
		}
	}


}