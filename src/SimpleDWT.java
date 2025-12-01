import java.io.*;
import java.lang.Math;
import java.util.*;
import binary.*;

/**
*	Image Compression, DWT / IDWT
*	[YCbCr 4:2:2, DCT 8*8, zigzag]

*	@author Ronan Merien <rmerien@hotmail.com>
*
*/
public class SimpleDWT {

	static public int COLS = -1;
	static public int ROWS = -1;
	static public int N = 8;
	static public double quality = 2;

	// x = i = row ; y = j = col

	float[][] imageY;
	float[][] imageCr;
	float[][] imageCb;

	float[][] dwtTemp;

	float[][] C = new float[N][N];
	float[][] Ct = new float[N][N];
	float[][] temp = new float[N][N];

	float[][] inputY = new float[N][N];
	float[][] inputCr = new float[N][N];
	float[][] inputCb = new float[N][N];

	float[][] outputY = new float[N][N];
	float[][] outputCr = new float[N][N];
	float[][] outputCb = new float[N][N];


	public void SimpleDWT() {
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

	float[][] quantumY =
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

	float[][] quantumCrCb =
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
/*
	double[][] quantumY =
	{
		{ 16, 11, 10, 16, 24, 40, 51, 61 },
		{ 12, 12, 14, 19, 26, 58, 60, 55 },
		{ 14, 13, 16, 24, 40, 57, 69, 56 },
		{ 14, 17, 22, 29, 51, 87, 80, 62 },
		{ 18, 22, 37, 56, 68, 109, 103, 77 },
		{ 24, 35, 59, 64, 81, 104, 113, 92 },
		{ 49, 64, 78, 87, 103, 121, 120, 101 },
		{ 72, 92, 95, 98, 112, 100, 103, 99 }
	};

	double[][] quantumCrCb =
	{
		{ 17, 18, 24, 47, 99, 99, 99, 99 },
		{ 18, 21, 26, 66, 99, 99, 99, 99 },
		{ 24, 26, 56, 99, 99, 99, 99, 99 },
		{ 47, 99, 99, 99, 99, 99, 99, 99 },
		{ 99, 99, 99, 99, 99, 99, 99, 99 },
		{ 99, 99, 99, 99, 99, 99, 99, 99 },
		{ 99, 99, 99, 99, 99, 99, 99, 99 },
		{ 99, 99, 99, 99, 99, 99, 99, 99 }
	};
*/
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

			C[0][j] = (float) (1.0f / Math.sqrt(N));
			Ct[j][0] = C[0][j];

		}

		for (i = 1; i < N; i++) {
			for (j = 0; j < N; j++) {

				C[i][j] = (float) (Math.sqrt(2.0f / N) * Math.cos(((2 * j + 1) * i * Math.PI) / (2.0f * N)));
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

	public void forwardDCT(float[][] image, float[][] DCT)
	{
		int i, j, k;

		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++) {

				temp[i][j] = 0.0f;
				for (k = 0; k < N; k++) temp[i][j] += image[i][k] * Ct[k][j];

			}
		}

		float vtemp;

		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++)	{

				vtemp = 0.0f;
				for (k = 0; k < N; k++)	vtemp += C[i][k] * temp[k][j];

				DCT[i][j] = (int) Math.round(vtemp);

			}
		}
	}

	// ----------------------------------------------------------------------------------------------------

/**
* l'image est suppos e etre multiple de 64
* cad 16 (la resolution maximale pour la DWT) par 8 (taille d'un bloc elementaire pour appliquer la DCT)
* dans notre exemple : test.bmp 384 * 256 = (3*8*16) * (2*8*16)
*/
	public void forwardDWT(float[][] image, int resolution)
	{
		int row, col;
		int rows = (int) ROWS/resolution;
		int cols = (int) COLS/resolution;

		System.out.println("rows=" + rows);
		System.out.println("cols=" + cols);

		// dwtTemp = f1 | fh

		// image en sortie :
		// f11 | fh1
		// f1h | fhh

		for (row = 0; row < rows; row ++) {
			for (col = 0; col < cols/2; col ++) {


	// Moyenner les pixels de l'image originale deux   deux suivant l'axe horizontal
	// f1(x,y) = (image(x,y) + image(x,y+1)) / 2

				dwtTemp[row][col] = ((image[row][2*col]
					+ image[row][2*col+1]) / 2);

	// Calculer l'erreur entre l'image originale et l'image sous- chantillonn es dans le sens horizontal
	// fh(x,y) = (image(x,y) - image(x,y+1)) / 2
				dwtTemp[row][col+cols/2] = ((image[row][2*col]
					- image[row][2*col+1]) / 2);

			}
		}

		for (row = 0; row < rows; row ++) {
			for (col = 0; col < cols; col ++) {
				image[row][col] = dwtTemp[row][col];
			}
		}

		for (row = 0; row < rows/2; row ++) {
			for (col = 0; col < cols/2; col ++) {

	//Pour chacune des deux images interm diaires, moyenner les pixels deux   deux suivant l'axe vertical
	//f11(x,y) = (f1(x,y) + f1(x+1,y)) / 2
	//fh1(x,y) = (fh(x,y) + fh(x+1,y)) / 2

	//Pour chacune des deux images interm diaires, calculer l'erreur suivant l'axe vertical
	//f1h(x,y) = (f1(x,y) - f1(x+1,y)) / 2
	//fhh(x,y) = (fh(x,y) - fh(x+1,y)) / 2

				image[row][col] =  ((dwtTemp[2*row][col]
					+ dwtTemp[2*row+1][col]) / 2);

				image[row][col+cols/2] =  ((dwtTemp[2*row][col+cols/2]
					+ dwtTemp[2*row+1][col+cols/2]) / 2);

				image[row+rows/2][col]  =  ((dwtTemp[2*row][col]
					- dwtTemp[2*row+1][col]) / 2);

				image[row+rows/2][col+cols/2] = ((dwtTemp[2*row][col+cols/2]
					- dwtTemp[2*row+1][col+cols/2]) / 2);

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

	public void inverseDCT(float[][] DCT, float[][] IDCT)
	{
		int i, j, k;
		float vtemp;

		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++) {

				temp[i][j] = 0.0f;
				for (k = 0; k < N; k++) temp[i][j] += DCT[i][k] * C[k][j];

			}
		}

		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++) {

				vtemp = 0.0f;
				for (k = 0; k < N; k++)	vtemp += Ct[i][k] * temp[k][j];

				// output byte from -128 to +127
				if (vtemp < -128) 		IDCT[i][j] = -128;
				else if (vtemp > 127) 	IDCT[i][j] = 127;
				else					IDCT[i][j] = (int) Math.round(vtemp);
				
				//IDCT[i][j] = vtemp; 

			}
		}
	}

	// ---------------------------------------------------------------------------------

	public void inverseDWT(float[][] image, int resolution)
	{
		int row, col;
		int rows = (int) ROWS/resolution;
		int cols = (int) COLS/resolution;

		System.out.println("rows=" + rows);
		System.out.println("cols=" + cols);

		// dwtTemp = f1 | fh

		// image en sortie :
		// f11 | fh1
		// f1h | fhh

		for (row = 0; row < rows/2; row ++) {
			for (col = 0; col < cols/2; col ++) {

				//System.out.println("row=" + row + " col=" + col);

				// f1(x+y) = f11 + f1h
				dwtTemp[2*row][col] = image[row][col] + image[row + rows/2][col];
				dwtTemp[2*row+1][col] = image[row][col] - image[row + rows/2][col];

				// fh(x+y) = fh1 + fhh
				dwtTemp[2*row][col+cols/2] = image[row][col+cols/2] + image[row + rows/2][col+cols/2];
				dwtTemp[2*row+1][col+cols/2] = image[row][col+cols/2] - image[row + rows/2][col+cols/2];
			}
		}


		for (row = 0; row < rows; row ++) {
			for (col = 0; col < cols/2; col ++) {

				//System.out.println("row=" + row + " col=" + col);

				//image(x+y) = f1 + fh

				image[row][2*col] = dwtTemp[row][col] + dwtTemp[row][col + cols/2];

				image[row][2*col+1] = dwtTemp[row][col] - dwtTemp[row][col + cols/2];

			}
		}


	}

	// ---------------------------------------------------------------------------------------

	public void compressFile(String inFile) throws Exception {
	//public void compressFile(String inFile, String outFile) throws Exception {
		try {
			int row, col, i, j, k;

			FileInputStream fis = new FileInputStream(inFile+".bmp");
			BinaryInputStream bis = new BinaryInputStream(new BufferedInputStream(fis));
			
			//FileOutputStream fos = new FileOutputStream(outFile);
			//BinaryOutputStream bos = new BinaryOutputStream(new BufferedOutputStream(fos));


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


			//bos.writeBit(COLS, 16);
			//bos.writeBit(ROWS, 16);

			imageY = new float[ROWS][COLS];
			imageCr = new float[ROWS][COLS];
			imageCb = new float[ROWS][COLS];

			dwtTemp = new float[ROWS][COLS];

			initialyze();

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

					imageY[row][col] = y- 128.0f;
					imageCr[row][col] = cr - 128.0f ;
					imageCb[row][col] = cb -128.0f ;

				}
			}

			exportBMP_YCrCb("YCrCb_before");
			exportBMP_RGB("compress"+String.valueOf(0));
			for (i=1; i<=16; i = 2*i) {
				forwardDWT(imageY, i);
				forwardDWT(imageCr, i);
				forwardDWT(imageCb, i);
				
				exportBMP_RGB("compress"+String.valueOf(i));
			}

			exportBMP_YCrCb("compress16YCrCb");

			DCT trans = new DCT();

			trans.compressFile("compress16YCrCb.bmp", inFile+".dwt");
			trans.expandFile(inFile+".dwt", "expand16YCrCb.bmp");
			importBMP_YCrCb("expand16YCrCb");
			importBMP_YCrCb("compress16YCrCb");



			for (i=16; i>=1; i = i/2) {
				exportBMP_RGB("expand"+String.valueOf(i));
				inverseDWT(imageY, i);
				inverseDWT(imageCr, i);
				inverseDWT(imageCb, i);
			}
			exportBMP_RGB("expand"+String.valueOf(i));



		} catch (Exception e) {
			System.out.println(e);
		}


	}

	// ---------------------------------------------------------------------------------
	public void expandFile(String inFile) throws Exception {
	//public void expandFile(String inFile, String outFile) throws Exception {
		try {
			int row, col, i, j, k;

			FileInputStream fis = new FileInputStream(inFile+".bmp");
			BinaryInputStream bis =	new BinaryInputStream(new BufferedInputStream(fis));

			String outFile = inFile+".idwt.bmp";
			FileOutputStream fos = new FileOutputStream(outFile);
			BinaryOutputStream bos = new BinaryOutputStream(new BufferedOutputStream(fos));

			byte b = 0;

			COLS = bis.readBit(16);
			ROWS = bis.readBit(16);

			imageY = new float[ROWS][COLS];
			imageCr = new float[ROWS][COLS];
			imageCb = new float[ROWS][COLS];

			dwtTemp = new float[ROWS][COLS];

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

			int blue, green, red;
			double y, cb, cr;

			for (row = 0; row < ROWS; row += N) {
				for (col = 0; col < COLS; col += N) {


					for (k = 0; k < N * N; k++) {
						i = zigzag(k).col;
						j = zigzag(k).row;

						inputY[i][j] = (float) (bis.readBit(11) - 1024);
						inputY[i][j] = (float) (inputY[i][j] * quantumY[i][j]);
					}

					inverseDCT(inputY, outputY);

					for (i = 0; i < N; i++) {
						for (j = 0; j < N; j++) {
							imageY[row + i][col + j] = outputY[i][j];
						}
					}

					for (k = 0; k < N * N; k++) {
						i = zigzag(k).col;
						j = zigzag(k).row;

						inputCr[i][j] = (float) (bis.readBit(11) - 1024);
						inputCr[i][j] = (float) (inputCr[i][j] * quantumCrCb[i][j]);
					}

					inverseDCT(inputCr, outputCr);

					for (i = 0; i < N; i++) {
						for (j = 0; j < N; j++) {
							imageCr[row + i][col + j] = outputCr[i][j];
						}
					}

					for (k = 0; k < N * N; k++) {
						i = zigzag(k).col;
						j = zigzag(k).row;

						inputCb[i][j] = (float) (bis.readBit(11) - 1024);
						inputCb[i][j] = (float) (inputCb[i][j] * quantumCrCb[i][j]);
					}

					inverseDCT(inputCb, outputCb);

					for (i = 0; i < N; i++) {
						for (j = 0; j < N; j++) {
							imageCb[row + i][col + j] = outputCb[i][j];
						}
					}
				}
			}

			for (row = 0; row < ROWS; row ++) {
				for (col = 0; col < COLS; col ++) {
					imageY[row][col] = imageY[row][col]/4;
					imageCr[row][col] = imageCr[row][col]/4;
					imageCb[row][col] = imageCb[row][col]/4;
				}
			}

			for (i=16; i>=1; i = i/2) {
				exportBMP_RGB("expand"+String.valueOf(i));
				inverseDWT(imageY, i);
				inverseDWT(imageCr, i);
				inverseDWT(imageCb, i);
			}
			exportBMP_RGB("expand"+String.valueOf(0));


			for (row = ROWS - 1; row >= 0; row--) {
				for (col = 0; col < COLS; col++) {

					y = (float) imageY[row][col] + 128.0f;
					cr = (float) imageCr[row][col];
					cb = (float) imageCb[row][col];

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
			fis.close();
		}
		catch (EOFException e) {}
		catch (IOException e) {System.out.println("expandFile :" + e);}
	}

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

		//int blue, green, red; 
		int y, cb, cr;


		for (row = ROWS - 1; row >= 0; row--) {
			for (col = 0; col < COLS; col++) {

				y = bis.readBit(8);
				cr = bis.readBit(8);
				cb = bis.readBit(8);

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