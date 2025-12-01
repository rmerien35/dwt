package binary;

/*
	Zip++ application
	Copyright 2000, 2003 Zip-Technology.com

	@Author Ronan Merien <ronan.merien@zip-technology.com>
	@Version 1.2, 2003/10/15
*/

import binary.*;

/**
	<p>A Huffman code can be represented as a binary tree.</p>
	<p>This class is used to construct Huffman codes storing words from the dictionnary,
	length and position (distance) for the LZSS algorithm.</p>
*/
public class BinaryTree {

	public String word = null;
	public BinaryTree[] tab = null;

	public int codeAscii = -1;
	public int frequence = 0;
	public int compressCode = 0x00;
	public int expandCode = 0x00;
	public int nbBits = 0;

	// public BinaryTree parent = null;
	public BinaryTree node_0 = null;
	public BinaryTree node_1 = null;

 	public BinaryTree ()
	{
		tab = new BinaryTree[26];
	}


	public BinaryTree (int b, String w, int f)
	{
		codeAscii = b;
		word = w;
		frequence = f;
	}

	public BinaryTree (int b, int f)
	{
		codeAscii = b;
		frequence = f;
	}

	public BinaryTree (int b, int f, BinaryTree n0, BinaryTree n1)
	{
		codeAscii = b;
		frequence = f;
		node_0 = n0;
		node_1 = n1;
	}


	public boolean isLeafNode() {
		if ((node_0 == null) || (node_1 == null)) return true;
	   else return false;
	}

	public void print() {
		if (codeAscii != -1) toString();

		if (node_0 != null) {
			node_0.print();
		}
		if (node_1 != null) {
			node_1.print();
		}
	}

	public java.lang.String toString() {

		System.out.println((codeAscii != -1 ? "codeAscii = " + (char) codeAscii : "")
						+ " , frequence = " + frequence
						+ " , compressCode = " + Binary.toBinaryString(compressCode,nbBits)
						+ " , expandCode = " + Binary.toBinaryString(expandCode,nbBits)
						+ " , nbBits = " + nbBits
						);
		return null;
	}

}