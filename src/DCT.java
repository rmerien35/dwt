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
public class DCT {


	static public int COLS = -1;
	static public int ROWS = -1;
	static public int N = 8;
	static public double quality = 2;

	int[][] imageY;
	int[][] imageCr;
	int[][] imageCb;

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

	Vector gen;

	public void DCT() {
	}

	protected int round_byte(double a) {
		if (a < 0) return 0;
		else if (a > 255) return 255;
		else return (int) Math.round(a);
	}

	// ----------------------------------------------------------------------------------------------------

	public class Pixel {
		int col, row;

		public Pixel(int i, int j) {
			col = i;
			row = j;
		}
	}

	Pixel[] zigzag8 =
	{
			new Pixel(0, 0),
			new Pixel(0, 1), new Pixel(1, 0),
			new Pixel(2, 0), new Pixel(1, 1), new Pixel(0, 2),
			new Pixel(0, 3), new Pixel(1, 2), new Pixel(2, 1), new Pixel(3, 0),
			new Pixel(4, 0), new Pixel(3, 1), new Pixel(2, 2), new Pixel(1, 3), new Pixel(0, 4), 
			new Pixel(0, 5), new Pixel(1, 4), new Pixel(2, 3), new Pixel(3, 2), new Pixel(4, 1), new Pixel(5, 0), 
			new Pixel(6, 0), new Pixel(5, 1), new Pixel(4, 2), new Pixel(3, 3), new Pixel(2, 4), new Pixel(1, 5), new Pixel(0, 6),
			new Pixel(0, 7), new Pixel(1, 6), new Pixel(2, 5), new Pixel(3, 4), new Pixel(4, 3), new Pixel(5, 2), new Pixel(6, 1), new Pixel(7, 0),
			new Pixel(7, 1), new Pixel(6, 2), new Pixel(5, 3), new Pixel(4, 4), new Pixel(3, 5), new Pixel(2, 6), new Pixel(1, 7), 
			new Pixel(2, 7), new Pixel(3, 6), new Pixel(4, 5), new Pixel(5, 4), new Pixel(6, 3), new Pixel(7, 2),
			new Pixel(7, 3), new Pixel(6, 4), new Pixel(5, 5), new Pixel(4, 6), new Pixel(3, 7),
			new Pixel(4, 7), new Pixel(5, 6), new Pixel(6, 5), new Pixel(7, 4),
			new Pixel(7, 5), new Pixel(6, 6), new Pixel(5, 7),
			new Pixel(6, 7), new Pixel(7, 6),
			new Pixel(7, 7)
	};

	/** 
	* zigzag sequence
	*
	*/

	protected Pixel zigzag(int k) {
		return zigzag8[k];
	}

	// ----------------------------------------------------------------------------------------------------

	double[][] quantumY = 
	{ 
		{ 15, 10, 10, 10, 10, 10, 10, 10 }, 
		{ 10, 10, 10, 10, 10, 10, 10, 10 }, 
		{ 10, 10, 10, 10, 10, 10, 10, 10 }, 
		{ 10, 10, 10, 10, 10, 10, 10, 10 }, 
		{ 10, 10, 10, 10, 10, 10, 10, 10 }, 
		{ 10, 10, 10, 10, 10, 10, 10, 10 }, 
		{ 10, 10, 10, 10, 10, 10, 10, 10 }, 
		{ 10, 10, 10, 10, 10, 10, 10, 10 }
	};

	double[][] quantumCrCb = 
	{ 
		{ 15, 10, 10, 10, 10, 10, 10, 10 }, 
		{ 10, 10, 10, 10, 10, 10, 10, 10 }, 
		{ 10, 10, 10, 10, 10, 10, 10, 10 }, 
		{ 10, 10, 10, 10, 10, 10, 10, 10 }, 
		{ 10, 10, 10, 10, 10, 10, 10, 10 }, 
		{ 10, 10, 10, 10, 10, 10, 10, 10 }, 
		{ 10, 10, 10, 10, 10, 10, 10, 10 }, 
		{ 10, 10, 10, 10, 10, 10, 10, 10 }
	};

	/** 
	* initialyze
	*
	* quantified_DCT[i][j] = DCT[i][j] / quantum[i][j]
	*
	* Cosinus Transform Matrix is C
	* C Transposed Matrix is Ct
	*/

	public void initialyze() 
	{
		int i, j;

		for (j = 0; j < N; j++) {

			C[0][j] = 1.0 / Math.sqrt(N);
			Ct[j][0] = C[0][j];

		}

		for (i = 1; i < N; i++) {
			for (j = 0; j < N; j++) {

				C[i][j] = Math.sqrt(2.0 / N) * Math.cos(((2 * j + 1) * i * Math.PI) / (2.0 * N));
				Ct[j][i] = C[i][j];

			}
		}
	}

	// ----------------------------------------------------------------------------------------------------

	/** 
	* forwardDCT
   *
	* DCT[i][j] = (1/sqrt(2N)) C(x) C(y) SUM(x,0,N-1) SUM(y,0,N-1) pixel[x][y] cos((i*pi*(2x+1))/2N) cos((j*pi*(2y+1))/2N)
	* C(0) = 1/sqrt(2) & C(x) = 1 if x>0
	*
   * DCT = C * image * Ct
	* where temp = image * Ct
	* and   DCT = C * temp
	*	
	* input byte from -128 to +127
	*/

	public void forwardDCT(int[][] image, int[][] DCT) 
	{
		int i, j, k;

		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++) {

				temp[i][j] = 0.0;
				for (k = 0; k < N; k++) temp[i][j] += image[i][k] * Ct[k][j];

			}
		}

		double vtemp;

		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++)	{

				vtemp = 0.0;
				for (k = 0; k < N; k++)	vtemp += C[i][k] * temp[k][j];

				DCT[i][j] = (int) Math.round(vtemp);

			}
		}
	}

	// ----------------------------------------------------------------------------------------------------

	/** 
	* inverseDCT
	*
	* IDCT[x][y] = (1/sqrt(2N)) SUM(i,0,N-1) SUM(j,0,N-1) C(i) C(j) DCT[i][j] cos((i*pi*(2x+1))/2N) cos((j*pi*(2y+1))/2N)
	* C(0) = 1/sqrt(2) & C(x) = 1 if x>0
	*
	* IDCT = Ct * DCT * C
	* where temp = image * C
	* and IDCT = Ct * temp
	*
	*/

	public void inverseDCT(int[][] DCT, int[][] IDCT)
	{
		int i, j, k;
		double vtemp;

		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++) {

				temp[i][j] = 0.0;
				for (k = 0; k < N; k++) temp[i][j] += DCT[i][k] * C[k][j];

			}
		}

		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++) {

				vtemp = 0.0;
				for (k = 0; k < N; k++)	vtemp += Ct[i][k] * temp[k][j];

				// output byte from -128 to +127
				if (vtemp < -128) 		IDCT[i][j] = -128;
				else if (vtemp > 127) 	IDCT[i][j] = 127;
				else							IDCT[i][j] = (int) Math.round(vtemp);

			}
		}
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
			} else if ((code < 128) || (code > 127)) {
				System.out.println("out of range " + code);
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
				System.out.print("out of range ");
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

			imageY = new int[ROWS][COLS];
			imageCr = new int[ROWS][COLS];
			imageCb = new int[ROWS][COLS];

			initialyze();

			int blue, green, red, y, cb, cr;

			for (row = ROWS - 1; row >= 0; row--) {
				for (col = 0; col < COLS; col++) {
/*
					blue = bis.readBit(8);
					green = bis.readBit(8);
					red = bis.readBit(8);

					y = (int) (0.299 * red + 0.587 * green + 0.114 * blue);
					cb = (int) (-0.1687 * red - 0.3313 * green + 0.5 * blue); // cb = (int) (0.564*(blue - y));
					cr = (int) (0.5 * red - 0.41874 * green - 0.08130 * blue); // cr = (int) (0.713*(red - y));
					
					imageY[row][col] = y - 128;
					imageCr[row][col] = cr;
					imageCb[row][col] = cb;
*/

					y = bis.readBit(8);
					cr = bis.readBit(8);
					cb = bis.readBit(8);

					imageY[row][col] = y-128;
					imageCr[row][col] = cr -128;
					imageCb[row][col] = cb -128;

				}
			}

			for (row = 0; row < ROWS; row += N) {
				for (col = 0; col < COLS; col += N) {

					for (i = 0; i < N; i++) {
						for (j = 0; j < N; j++) {
							inputY[i][j] = imageY[row + i][col + j];
							inputCr[i][j] = imageCr[row + i][col + j];
							inputCb[i][j] = imageCb[row + i][col + j];
						}
					}

					forwardDCT(inputY, outputY);

					forwardDCT(inputCr, outputCr);

					forwardDCT(inputCb, outputCb);

	
					for (k = 0; k < N * N; k++) {
						i = zigzag(k).col;
						j = zigzag(k).row;
						outputY[i][j] =
							(int) Math.round(outputY[i][j] / quantumY[i][j]);

						if (k == 0)
							encodingDC = true;
						outputCode(bos, outputY[i][j]);
					}

					for (k = 0; k < N * N; k++) {
						i = zigzag(k).col;
						j = zigzag(k).row;
						outputCr[i][j] = (int) Math.round(outputCr[i][j] / quantumCrCb[i][j]);

						if (k == 0) encodingDC = true;
						outputCode(bos, outputCr[i][j]);
					}

					for (k = 0; k < N * N; k++) {
						i = zigzag(k).col;
						j = zigzag(k).row;
						outputCb[i][j] = (int) Math.round(outputCb[i][j] / quantumCrCb[i][j]);

						if (k == 0)	encodingDC = true;
						outputCode(bos, outputCb[i][j]);
					}
				}
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

			COLS = bis.readBit(16);	// System.out.println("COLS = " + COLS);
			ROWS = bis.readBit(16); // System.out.println("ROWS = " + ROWS);

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


			imageY = new int[ROWS][COLS];
			imageCr = new int[ROWS][COLS];
			imageCb = new int[ROWS][COLS];

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


			initialyze();

			int blue, green, red, y, cb, cr;

			for (row = 0; row < ROWS; row += N) {
				for (col = 0; col < COLS; col += N) {

					for (k = 0; k < N * N; k++) {
						i = zigzag(k).col;
						j = zigzag(k).row;

						if (k == 0) encodingDC = true;
						inputY[i][j] = inputCode(bis);
						inputY[i][j] = (int) Math.round(inputY[i][j] * quantumY[i][j]);
					}

					for (k = 0; k < N * N; k++) {
						i = zigzag(k).col;
						j = zigzag(k).row;

						if (k == 0) encodingDC = true;
						inputCr[i][j] = inputCode(bis);
						inputCr[i][j] = (int) Math.round(inputCr[i][j] * quantumCrCb[i][j]);
					}

					for (k = 0; k < N * N; k++) {
						i = zigzag(k).col;
						j = zigzag(k).row;

						if (k == 0)
							encodingDC = true;
						inputCb[i][j] = inputCode(bis);
						inputCb[i][j] = (int) Math.round(inputCb[i][j] * quantumCrCb[i][j]);
					}

					inverseDCT(inputY, outputY);

					inverseDCT(inputCr, outputCr);

					inverseDCT(inputCb, outputCb);

					for (i = 0; i < N; i++) {
						for (j = 0; j < N; j++) {
							imageY[row + i][col + j] = outputY[i][j];
						}
					}
				
					for (i = 0; i < N; i++) {
						for (j = 0; j < N; j++) {
							imageCr[row + i][col + j] = outputCr[i][j];
						}
					}

					for (i = 0; i < N; i++) {
						for (j = 0; j < N; j++) {
							imageCb[row + i][col + j] = outputCb[i][j];
						}
					}
				}
			}

			for (row = ROWS - 1; row >= 0; row--) {
				for (col = 0; col < COLS; col++) {
/*
					y = imageY[row][col] + 128;
					cr = imageCr[row][col];
					cb = imageCb[row][col];

					blue = round_byte(y + 1.773 * cb);
					green = round_byte(y - 0.34414 * cb - 0.71414 * cr);
					red = round_byte(y + 1.402 * cr);

					bos.writeBit(blue, 8);
					bos.writeBit(green, 8);
					bos.writeBit(red, 8);
*/
				y = (int) imageY[row][col]+ 128;
				cr = (int) imageCr[row][col]+ 128;
				cb = (int) imageCb[row][col]+ 128;

				bos.writeBit(y, 8);
				bos.writeBit(cr, 8);
				bos.writeBit(cb, 8);
				}
			}

			bos.flush();
			fos.close();
			fis.close();
		} 
		catch (EOFException e) {}
		catch (IOException e) {System.out.println("expandFile :" + e);}

		exportBMP_RGB("RGB_"+outFile);
	}

	// ---------------------------------------------------------------------------------

	public void exportBMP_RGB(String name) throws Exception
	{
		int row, col;
		int blue, green, red, y, cb, cr;

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

				y = (int) imageY[row][col] + 128;
				cr = (int) imageCr[row][col];
				cb = (int) imageCb[row][col];

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

	}

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

	// ---------------------------------------------------------------------------------------------

	static public void help() 
	{
		System.out.println("DCT 1.1 powered by Ronan Merien <rmerien@hotmail.com>");
		System.out.println("under the LGPL open source licence");
		System.out.println("[YCbCr 4:2:2, DCT 8*8, quantization, zigzag sequence, RLE & Huffman]");
		System.out.println("Usage: java DCT [options] imagefile");
		System.out.println("-c		compress bmp_to_dctfile");
		System.out.println("-e		expand dct_to_bmpfile");
		System.out.println();
		//System.out.println("-quality n	quality factor between 1 and 100 (default 20)");
		System.out.println("Examples:  java DCT -c photo.bmp");
		System.out.println("           java DCT -e photo.dct");
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
						DCT.quality = p / 10;
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

		DCT trans = new DCT();

		// compress bmp_to_dctfile 
		if (option.equals("-c")) 
		{
			System.out.println("compress " + imagefile + ".bmp to " + imagefile + ".dct");

			try {	
				trans.compressFile(imagefile + ".bmp", imagefile + ".dct"); 
			} 
			catch (Exception e) { System.out.println(e); }

		}

		// expand dct_to_bmpfile
		else if (option.equals("-e")) 
		{
			System.out.println("expand " + imagefile + ".dct to " + imagefile + ".idct.bmp");
			
			try {
				trans.expandFile(imagefile + ".dct", imagefile + ".idct.bmp");
			}
			catch (Exception e) { System.out.println(e); }
			
		}

		// help instructions 
		else {
			help();
		}
	}

}