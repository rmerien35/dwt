package binary;

/*
	Zip++ application
	Copyright 2000, 2003 Zip-Technology.com

	@Author Ronan Merien <ronan.merien@zip-technology.com>
	@Version 1.2, 2003/10/15
*/

/**
	<p>A MatchDistance object contain the distance of a string matching.</p>
*/
public class MatchDistance {

	public int value = -1;
	public int nbBits = -1;

	public int extraValue = 0;
	public int nbExtraBits = 0;

	public MatchDistance ()
	{
	}

	public MatchDistance (int _value, int _nbBits)
	{
		value = _value;
		nbBits = _nbBits;
	}

	public MatchDistance (int _value, int _nbBits, int _extraValue, int _nbExtraBits)
	{
		value = _value;
		nbBits = _nbBits;

		extraValue = _extraValue;
		nbExtraBits = _nbExtraBits;
	}
}