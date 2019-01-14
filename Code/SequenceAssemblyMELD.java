import java.io.*;
import javax.imageio.*;
import java.awt.image.*;
import java.awt.*;
import java.lang.*;
import java.lang.Math.*;
import java.util.*;
import java.util.Arrays.*;
import java.util.regex.*;

//-------------------------
//Java code executed from the command line for sequence assembly of MELD protein digest peptide sequences (based on mass spectrometry).
//This code takes as input the MS confidences scores and digested peptide sequences and automatically identifies sutable seed sequences.
//The end result of this code is an elongated sequence tag for each seed sequence 
//Created during 2011-2018 by Tyler. A. Zimmerman, Ph.D. 
//-------------------------

public class SequenceAssemblyMELD
{
	public static void main(String args[ ]) throws IOException
		{

		String[] filename = new String[2];
		filename[0]="C:\\Documents\\";
		filename[1]="input.csv";
		ReadLines(filename);

		}

static void ReadLines(String[] filename) throws IOException
		{

		String WorkingPath = filename[0];

                if (WorkingPath.endsWith("\\") == false | WorkingPath.endsWith("/") == false)
		{
		WorkingPath = WorkingPath + "\\";
		}

//------------------------------------------------------
//Reading the raw data
//------------------------------------------------------
                
        BufferedReader inLines = new BufferedReader(new FileReader(WorkingPath + filename[1]));

		int lineCount = 0;
			
		do
		{
	        String line = inLines.readLine();
			if (line == null) break;
	        	lineCount++;	
		} while (true);

		String[] PeptideComplete = new String[lineCount];  //PeptideComplete includes all variable PTMs
		String[] Peptide = new String[lineCount];			//Peptide has eliminated the variable PTMs (such as N-glycosylation) and replaced them with the base residues
		String[] Confidence = new String[lineCount];
		int[] ALC = new int[lineCount];
		double[] PPM = new double[lineCount];
		String YesNo = "Yes";


        BufferedReader inLines1 = new BufferedReader(new FileReader(WorkingPath + filename[1]));

		int Rownum = 0;

		do
		{
			String Line1 = inLines1.readLine();			

			if (Line1 == null) break;			
			
			if (Line1.indexOf("local confidence") == -1)
			{
			
				System.out.println(Line1);
			
				PeptideComplete[Rownum] = Line1.split(",")[1].trim();
				Peptide[Rownum] = Line1.split(",")[1].replaceAll("N\\(\\+1606\\.59\\)", "N").replaceAll("N\\(\\+1444\\.53\\)", "N").replaceAll("N\\(\\+1768\\.64\\)", "N").replaceAll("C\\(\\+57\\.02\\)", "C").replaceAll("N\\(\\+\\.98\\)", "D").replaceAll("M\\(\\+15\\.99\\)", "M").replaceAll("Q\\(\\+\\.98\\)", "E").replaceAll("I", "I");
				Confidence[Rownum] = Line1.split(",")[9];
				ALC[Rownum] = Integer.valueOf(Line1.split(",")[3]);
				PPM[Rownum] = Double.parseDouble(Line1.split(",")[8]);
				
				if (Peptide[Rownum] == null) break;

				Rownum++;
			}	

		} while (true && Rownum < lineCount);


//------------------------------------------------------
//Defining the seed sequences of 5 aa, as taken from 
//the C-terminal constant region of the IgG
//------------------------------------------------------

		String SeedSeqsTemp[] = new String[Peptide.length*2];
		double SeedSeqsConfidenceTemp[][] = new double[Peptide.length*2][5];
		int SeedCount = 0;

		for (int i = 0; i <= Peptide.length-1; i++)
		{
			if (Peptide[i] != null && Confidence[i] != null)
			{
			if (Peptide[i].length() > 6)
			{
				String ConfidenceSeed[] = new String[100];

				ConfidenceSeed = Confidence[i].split(" ");
				
//-------------------------
//Here is defined the MinConfidence, Default = 70
//-------------------------

				int MinConfidence = 70;

				if (ConfidenceSeed.length > 7)
				{
				try
				{
				for (int j = 0; j <= ConfidenceSeed.length-6; j++)
				{
					if (Integer.valueOf(ConfidenceSeed[j]) >  MinConfidence && Integer.valueOf(ConfidenceSeed[j+1]) >  MinConfidence && Integer.valueOf(ConfidenceSeed[j+2]) >  MinConfidence && Integer.valueOf(ConfidenceSeed[j+3]) >  MinConfidence && Integer.valueOf(ConfidenceSeed[j+4]) >  MinConfidence)
					{
						SeedSeqsTemp[SeedCount] = Peptide[i].substring(j, j+5);
						SeedSeqsConfidenceTemp[SeedCount][0] = Double.parseDouble(ConfidenceSeed[j]);
						SeedSeqsConfidenceTemp[SeedCount][1] = Double.parseDouble(ConfidenceSeed[j+1]);
						SeedSeqsConfidenceTemp[SeedCount][2] = Double.parseDouble(ConfidenceSeed[j+2]);
						SeedSeqsConfidenceTemp[SeedCount][3] = Double.parseDouble(ConfidenceSeed[j+3]);
						SeedSeqsConfidenceTemp[SeedCount][4] = Double.parseDouble(ConfidenceSeed[j+4]);
						SeedCount++;
						break;
					}				
				}
				} catch (Exception e) {System.out.println(e);}
				}
			}
			}
		}

//------------------------------------------------------
//Removing redundant seed sequences
//------------------------------------------------------

		String SeedSeqsTemp2[] = new String[Peptide.length*2];
		int SeedSeqsCount = 0;

		for (int i = 0; i <= SeedSeqsTemp.length-1; i++)
		{
			if (SeedSeqsTemp[i] != null)
			{	

			int Already = 0;
			int AlreadyCount = 0;

			for (int j = 0; j <= SeedSeqsTemp2.length-1; j++)
			{
				if (SeedSeqsTemp2[j] != null)
				{
					if (SeedSeqsTemp[i].equals(SeedSeqsTemp2[j]) == true)
					{
						Already = 1;

					}
				}
			}

			if (Already == 0)
			{
				SeedSeqsTemp2[SeedSeqsCount] = SeedSeqsTemp[i];
				SeedSeqsCount++;
			}

			}
		}

		String SeedSeqs[] = new String[SeedSeqsCount];
		double SeedSeqsConfidence[][] = new double[SeedSeqsCount][5];

		for (int i = 0; i <= SeedSeqsTemp2.length-1; i++)
		{
			if (SeedSeqsTemp2[i] != null)
			{
				SeedSeqs[i] = SeedSeqsTemp2[i];
			}
		}

//------------------------------------------------------
//Calculating the seed sequence redundancy
//------------------------------------------------------

		int Redundancy[] = new int[SeedSeqs.length];

		for (int i = 0; i <= SeedSeqs.length-1; i++)
		{
			for (int j = 0; j <= SeedSeqsTemp.length-1; j++)
			{
				if (SeedSeqsTemp[j] != null)
				{
				if (SeedSeqsTemp[j].equals(SeedSeqs[i]) == true)
				{
					Redundancy[i]++;
				}
				}
			}
		}

//------------------------------------------------------
//Calculating the seed sequence confidence, 
//averaged over all instances of a seed sequence
//------------------------------------------------------

		int SeedSeqsConfidenceCount[] = new int[SeedSeqsConfidence.length];

		for (int i = 0; i <= SeedSeqs.length-1; i++)
		{
			for (int j = 0; j <= SeedSeqsTemp.length-1; j++)
			{
				if (SeedSeqs[i].equals(SeedSeqsTemp[j]) )
				{
					SeedSeqsConfidenceCount[i]++;
				}
			}
		}


		for (int i = 0; i <= SeedSeqs.length-1; i++)
		{
			for (int j = 0; j <= SeedSeqsTemp.length-1; j++)
			{
				if (SeedSeqs[i].equals(SeedSeqsTemp[j]) )
				{
					for (int k = 0; k <= SeedSeqsConfidenceTemp[0].length-1; k++)
					{

						SeedSeqsConfidence[i][k] += SeedSeqsConfidenceTemp[j][k];
					}
				}
			}
		}

		for (int i = 0; i <= SeedSeqsConfidence.length-1; i++)
		{
			for (int k = 0; k <= SeedSeqsConfidenceTemp[0].length-1; k++)
			{
				SeedSeqsConfidence[i][k] = SeedSeqsConfidence[i][k] / SeedSeqsConfidenceCount[i];
			}
		}

//------------------------------------------------------
//Initalizing some variables
//------------------------------------------------------

		String SeedSeq = null;
		int SeedRedundancy = 0;
		double SeedSeqsConfidenceOutput[] = new double[5];

		int Manual = 0;

		Manual = 1;

		if (Manual == 1)
		{
			//SeedSeq = "HEAL";
			//SeedSeq = "NVNH";

			//SeedSeq = "HATKHK";
			//SeedSeq = "FRND";	
			//SeedSeq = "PETL";
			//SeedSeq = "LSDG";
			//SeedSeq = filename[2];
		}

//------------------------------------------------------
//Output file name
//------------------------------------------------------

	FileOutputStream out;
	PrintStream p;
	out = new FileOutputStream(WorkingPath + "Results_Reduundancy30.txt");
	p = new PrintStream( out );


//------------------------------------------------------
//Beginning of the big for loop
//------------------------------------------------------

System.out.println(Arrays.toString(SeedSeqs));
System.out.println(Arrays.deepToString(SeedSeqsConfidence));

for (int I = 0; I <= SeedSeqs.length-1; I++)
{
		SeedSeq = SeedSeqs[I];
		SeedRedundancy = Redundancy[I];
		SeedSeqsConfidenceOutput = SeedSeqsConfidence[I];

//------------------------------------------------------
//Defining the output parameters
//------------------------------------------------------


		String CurrentSeqL;
		String CurrentSeqR;
		CurrentSeqL = SeedSeq;
		CurrentSeqR = SeedSeq;

		String TotalSeq;
		TotalSeq = SeedSeq;


		String ConfidenceOutputFreq = "X X X X X";
		String CountCorrectOutput = "X X X X X";
		String AvgConfidenceOutput = "X X X X X";
		String TotalConfidenceOutput = "X X X X X";

		int ExtensionCountL = 1000000;
		int ExtensionCountR = 1000000;

//------------------------------------------------------
//Extending the seed sequence, beginning of the big "do-while" loop
//------------------------------------------------------

int ALCthreshold = -1;

int StopLeft = 0;
int StopRight = 0;

if (SeedRedundancy > 30)
{


while (ExtensionCountL > 0 | ExtensionCountR > 0 && StopLeft == 0 | StopRight == 0 && TotalSeq.length() < 100)
{
		CurrentSeqL = TotalSeq.substring(0, 3);
		CurrentSeqR = TotalSeq.substring(TotalSeq.length()-3, TotalSeq.length());

		System.out.println(CurrentSeqL);
		System.out.println(CurrentSeqR);

		String[] GrowthSeqLeft = new String[Peptide.length];
		String[] GrowthSeqLeftComplete = new String[Peptide.length];
		String[] GrowthSeqRight = new String[Peptide.length];
		String[][] GrowthConfidenceLeft = new String[Peptide.length][100];
		String[][] GrowthConfidenceRight = new String[Peptide.length][100];

		int CountL = 0;	
		int CountR = 0;	

		for (int i = 0; i <= Peptide.length-1; i++)
		{
				if (Peptide[i] != null)
				{
				if (Peptide[i].contains(CurrentSeqL) == true)
				{
				if (Peptide[i].indexOf(CurrentSeqL) - 1 >= 0)
				{			
					GrowthSeqLeft[CountL] = Peptide[i].substring(Peptide[i].indexOf(CurrentSeqL) - 1, Peptide[i].indexOf(CurrentSeqL) + CurrentSeqL.length());	
					GrowthSeqLeftComplete[CountL] = PeptideComplete[i].substring(Peptide[i].indexOf(CurrentSeqL) - 1, Peptide[i].indexOf(CurrentSeqL) + CurrentSeqL.length());
					GrowthConfidenceLeft[CountL] = Arrays.copyOfRange(Confidence[i].split(" "), Peptide[i].indexOf(CurrentSeqL) - 1, Peptide[i].indexOf(CurrentSeqL) + CurrentSeqL.length());	
					CountL++;	
				}
				}
				}

				if (Peptide[i] != null)
				{
				if (Peptide[i].contains(CurrentSeqR) == true)
				{
				if (Peptide[i].indexOf(CurrentSeqR) + CurrentSeqR.length() + 1 <= Peptide[i].length())
				{			
					GrowthSeqRight[CountR] = Peptide[i].substring(Peptide[i].indexOf(CurrentSeqR), Peptide[i].indexOf(CurrentSeqR) + CurrentSeqR.length() + 1);	
					GrowthConfidenceRight[CountR] = Arrays.copyOfRange(Confidence[i].split(" "), Peptide[i].indexOf(CurrentSeqR), Peptide[i].indexOf(CurrentSeqR) + CurrentSeqR.length() + 1);
					CountR++;
				}
				}
				}
				
				//System.out.println(Arrays.toString(GrowthConfidenceLeft[CountL]));
				//System.out.println(CountL);
				//System.out.println(Arrays.toString(GrowthConfidenceRight[CountR]));
				//System.out.println(GrowthSeqRight[CountR]);				
		}

				System.out.println(CurrentSeqL + " is found " + CountL + " times");
				System.out.println(CurrentSeqR + " is found " + CountR + " times");

//------------------------------------------------------
//Counting the frequency of the one amino acid adjacent 
//to the seed sequence, for the left (N-terminal) side
//------------------------------------------------------

		ExtensionCountL = 0;

		for (int i = 0; i <= Peptide.length-1; i++)
		{
			if (Peptide[i] != null)
			{
			if (Peptide[i].contains(CurrentSeqL) == true && Peptide[i].indexOf(CurrentSeqL) - 1 >= 0)
			{
				ExtensionCountL++;
			}
			}
		}

		if (ExtensionCountL <= 5)
		{
			StopLeft = 1;
		}

if (ExtensionCountL > 5 && StopLeft == 0)
{
		int[] NumValuesLeft = new int[CountL];
		char[] CharValuesLeft = new char[CountL];


		for (int i = 0; i <= CountL-1; i++)
		{
			if (GrowthSeqLeft[i] != null)
			{
				char Current1Left = GrowthSeqLeft[i].charAt(0);  // the first amino acid that was concatenated above.
				NumValuesLeft[i] = Character.getNumericValue(Current1Left);  // Converting char to numeric values
				CharValuesLeft[i] = Current1Left;  // The corresponding char value
			}
		}

		//System.out.println(Arrays.toString(NumValuesLeft));
		//System.out.println(Arrays.toString(CharValuesLeft));
		
		// maxKey is the letter identity of the amino acid
		// maxValue is its frequency of occurrance

		Map<Integer, Integer> m = new HashMap<Integer, Integer>();

	        int maxKey = -1;  
	        int maxValue = -1;
	        int maxKey2nd = -1;  
	        int maxValue2nd = -1;
	        int maxKey3rd = -1;  
	        int maxValue3rd = -1;   
	        int maxKey4th = -1;  
	        int maxValue4th = -1; 			

		for (int i = 0; i <= NumValuesLeft.length-1; i++)
		{
			if (!m.containsKey(NumValuesLeft[i]))   // If m does not contain it yet, add it to m.
			{
				m.put(NumValuesLeft[i], 1);
			}
			else
			{
				m.put(NumValuesLeft[i], m.get(NumValuesLeft[i])+1);  // Otherwise, if it does contain it, add it with ...
			}
		}

	    for (Map.Entry<Integer, Integer> entry : m.entrySet())  
		{  
            if (entry.getValue() > maxValue) 
			{  
            	maxKey = entry.getKey();  // Getting the max frequency occuring value
	            maxValue = entry.getValue();  
		    }  
		}

	    for (Map.Entry<Integer, Integer> entry1 : m.entrySet())  
		{  
        		if (entry1.getKey() != maxKey && entry1.getKey() != 0 && entry1.getValue() > maxValue2nd) 
			{  
                		maxKey2nd = entry1.getKey();  
	        	        maxValue2nd = entry1.getValue();  
		        }  
		}

	    for (Map.Entry<Integer, Integer> entry2 : m.entrySet()) 
		{  
        		if (entry2.getKey() != maxKey && entry2.getKey() != maxKey2nd && entry2.getKey() != 0 && entry2.getValue() > maxValue3rd) 
			{  
                		maxKey3rd = entry2.getKey();
	        	        maxValue3rd = entry2.getValue();
		        }  
		}

	    for (Map.Entry<Integer, Integer> entry3 : m.entrySet())  
		{  
        	if (entry3.getKey() != maxKey && entry3.getKey() != maxKey2nd && entry3.getKey() != maxKey3rd && entry3.getKey() != 0 && entry3.getValue() > maxValue4th) 
			{  
          		maxKey4th = entry3.getKey();  
	            maxValue4th = entry3.getValue();  
		    }  
		}

		System.out.println("The winner is residue "+Character.forDigit(maxKey, Character.MAX_RADIX)+" its frequency of occurrence is "+maxValue);  
		System.out.println("The 2nd place winner is residue "+Character.forDigit(maxKey2nd, Character.MAX_RADIX)+" its frequency of occurrence is "+maxValue2nd);  
		System.out.println("The 3rd place winner is residue "+Character.forDigit(maxKey3rd, Character.MAX_RADIX)+" its frequency of occurrence is "+maxValue3rd);  
		System.out.println("The 4th place winner is residue "+Character.forDigit(maxKey4th, Character.MAX_RADIX)+" its frequency of occurrence is "+maxValue4th);  

		// Converting numeric values back to their char equivalents
		
		String maxKeyString = String.valueOf(Character.forDigit(maxKey, Character.MAX_RADIX));
		maxKeyString = maxKeyString.toUpperCase();
		String maxKey2ndString = String.valueOf(Character.forDigit(maxKey2nd, Character.MAX_RADIX));
		maxKey2ndString = maxKey2ndString.toUpperCase();
		String maxKey3rdString = String.valueOf(Character.forDigit(maxKey3rd, Character.MAX_RADIX));
		maxKey3rdString = maxKey3rdString.toUpperCase();
		String maxKey4thString = String.valueOf(Character.forDigit(maxKey4th, Character.MAX_RADIX));
		maxKey4thString = maxKey4thString.toUpperCase();

		// Doing the same thing for 4 potential concatenated residues, but this time using confidence score rather than total frequency count
		// Counting the total occurances (?) of each potential concatenated residue (4 total concatenated residues) 

		// The Correct tallies seem to be checking if the 3 amino acids (CurrentSeqL) of the seed sequence exist in the peptide & and 
		// frist 4 amino acids of the TotalSeq (Seed Sequence) are there. So, it's looking for slightly more amino acids, one extra.
		
		Double CountCorrect1st = 0.0;
		Double CountCorrect2nd = 0.0;
		Double CountCorrect3rd = 0.0;
		Double CountCorrect4th = 0.0;

		for (int i = 0; i <= Peptide.length-1; i++)
		{
			if (Peptide[i] != null)
			{
				if (Peptide[i].contains(TotalSeq.substring(0, 4)) == true && Peptide[i].contains(maxKeyString.concat(CurrentSeqL)) == true)
				{
				CountCorrect1st++;
				}

				if (Peptide[i].contains(TotalSeq.substring(0, 4)) == true && Peptide[i].contains(maxKey2ndString.concat(CurrentSeqL)) == true)
				{
				CountCorrect2nd++;
				}

				if (Peptide[i].contains(TotalSeq.substring(0, 4)) == true && Peptide[i].contains(maxKey3rdString.concat(CurrentSeqL)) == true)
				{
				CountCorrect3rd++;
				}

				if (Peptide[i].contains(TotalSeq.substring(0, 4)) == true && Peptide[i].contains(maxKey4thString.concat(CurrentSeqL)) == true)
				{
				CountCorrect4th++;
				}
			}
		}

		System.out.println("CountCorrect for " + maxKeyString + " = " + CountCorrect1st);
		System.out.println("CountCorrect for " + maxKey2ndString + " = " + CountCorrect2nd);
		System.out.println("CountCorrect for " + maxKey3rdString + " = " + CountCorrect3rd);
		System.out.println("CountCorrect for " + maxKey4thString + " = " + CountCorrect4th);
	

		Double TotalConfidenceLeft = 0.0;
		Double TotalConfidenceLeftUn = 0.0;
		Double TotalZerosLeft = 0.0;
		Double TotalHighLeft = 0.0;
		

		for (int i = 0; i <= CountL-1; i++)
		{
			if (GrowthSeqLeft[i] != null && GrowthSeqLeft[i].startsWith(maxKeyString))
			{
				TotalConfidenceLeftUn += Integer.valueOf(GrowthConfidenceLeft[i][0]);
			
				if (Integer.valueOf(GrowthConfidenceLeft[i][0]) == 0)
				{
				TotalZerosLeft++;  // Counting the total number of zero value confidences for this residue, to be used a another measure of confidence
				}

				if (Integer.valueOf(GrowthConfidenceLeft[i][0]) > 95)
				{
				TotalHighLeft++;
				}
			}
		}
		
		Double AvgConfidenceLeft = (double)TotalConfidenceLeftUn / maxValue;
		TotalConfidenceLeft = TotalConfidenceLeftUn / (TotalZerosLeft+1) * TotalHighLeft;

		System.out.println(TotalConfidenceLeftUn + " " + TotalConfidenceLeft);
		System.out.println(TotalZerosLeft);
		System.out.println(TotalHighLeft);

		Double TotalConfidenceLeft2nd = 0.0;
		Double TotalConfidenceLeft2ndUn = 0.0;
		Double TotalZerosLeft2nd = 0.0;
		Double TotalHighLeft2nd = 0.0;

		for (int i = 0; i <= CountL-1; i++)
		{
			if (GrowthSeqLeft[i] != null && GrowthSeqLeft[i].startsWith(maxKey2ndString))
			{
				TotalConfidenceLeft2ndUn += Integer.valueOf(GrowthConfidenceLeft[i][0]);

				if (Integer.valueOf(GrowthConfidenceLeft[i][0]) == 0)
				{
				TotalZerosLeft2nd++;
				}

				if (Integer.valueOf(GrowthConfidenceLeft[i][0]) > 95)
				{
				TotalHighLeft2nd++;
				}
			}
		}

		Double AvgConfidenceLeft2nd = (double)TotalConfidenceLeft2ndUn / maxValue2nd;
		TotalConfidenceLeft2nd = TotalConfidenceLeft2ndUn / (TotalZerosLeft2nd+1) * TotalHighLeft2nd;

		System.out.println(TotalConfidenceLeft2ndUn + " " + TotalConfidenceLeft2nd);
		System.out.println(TotalZerosLeft2nd);
		System.out.println(TotalHighLeft2nd);

		Double TotalConfidenceLeft3rd = 0.0;
		Double TotalConfidenceLeft3rdUn = 0.0;
		Double TotalZerosLeft3rd = 0.0;
		Double TotalHighLeft3rd = 0.0;

		for (int i = 0; i <= CountL-1; i++)
		{
			if (GrowthSeqLeft[i] != null && GrowthSeqLeft[i].startsWith(maxKey3rdString))
			{
				TotalConfidenceLeft3rdUn += Integer.valueOf(GrowthConfidenceLeft[i][0]);

				if (Integer.valueOf(GrowthConfidenceLeft[i][0]) == 0)
				{
				TotalZerosLeft3rd++;
				}

				if (Integer.valueOf(GrowthConfidenceLeft[i][0]) > 95)
				{
				TotalHighLeft3rd++;
				}
			}
		}

		Double AvgConfidenceLeft3rd = (double)TotalConfidenceLeft3rdUn / maxValue3rd;
		TotalConfidenceLeft3rd = TotalConfidenceLeft3rdUn / (TotalZerosLeft3rd+1) * TotalHighLeft3rd;

		System.out.println(TotalConfidenceLeft3rdUn + " " + TotalConfidenceLeft3rd);
		System.out.println(TotalZerosLeft3rd);
		System.out.println(TotalHighLeft3rd);

		Double TotalConfidenceLeft4th = 0.0;
		Double TotalConfidenceLeft4thUn = 0.0;
		Double TotalZerosLeft4th = 0.0;
		Double TotalHighLeft4th = 0.0;

		for (int i = 0; i <= CountL-1; i++)
		{
			if (GrowthSeqLeft[i] != null && GrowthSeqLeft[i].startsWith(maxKey4thString))
			{
				TotalConfidenceLeft4thUn += Integer.valueOf(GrowthConfidenceLeft[i][0]);

				if (Integer.valueOf(GrowthConfidenceLeft[i][0]) == 0)
				{
				TotalZerosLeft4th++;
				}

				if (Integer.valueOf(GrowthConfidenceLeft[i][0]) > 95)
				{
				TotalHighLeft4th++;
				}
			}
		}

		Double AvgConfidenceLeft4th = (double)TotalConfidenceLeft4thUn / maxValue4th;
		TotalConfidenceLeft4th = TotalConfidenceLeft4thUn / (TotalZerosLeft4th+1) * TotalHighLeft4th;

		System.out.println(TotalConfidenceLeft4thUn + " " + TotalConfidenceLeft4th);
		System.out.println(TotalZerosLeft4th);
		System.out.println(TotalHighLeft4th);

//------------------------------------------------------
//Concatenation of the left (N-terminal) residue
//ocurrs here
//------------------------------------------------------

	int YesNo0 = 0;
	int YesNo1 = 0;
	int YesNo2 = 0;
	int YesNo3 = 0;
	int YesNo4 = 0;
	int YesNo5 = 0;

	if (CountCorrect1st > 10 | CountCorrect2nd > 10 | CountCorrect3rd > 10 | CountCorrect4th > 10)
	{
		if (CountCorrect1st*maxValue > CountCorrect2nd*maxValue2nd && CountCorrect1st*maxValue > CountCorrect3rd*maxValue3rd && CountCorrect1st*maxValue > CountCorrect4th*maxValue4th && TotalConfidenceLeft > 1)
		{
			 TotalSeq = maxKeyString.concat(TotalSeq);
			 YesNo0 = 1;
		}

		if (CountCorrect2nd*maxValue2nd > CountCorrect1st*maxValue && CountCorrect2nd*maxValue2nd > CountCorrect3rd*maxValue3rd && CountCorrect2nd*maxValue2nd > CountCorrect4th*maxValue4th && TotalConfidenceLeft2nd > 1)
		{
			 TotalSeq = maxKey2ndString.concat(TotalSeq);
			 YesNo0 = 1;
		}

		if (CountCorrect3rd*maxValue3rd > CountCorrect1st*maxValue && CountCorrect3rd*maxValue3rd > CountCorrect2nd*maxValue2nd && CountCorrect3rd*maxValue3rd > CountCorrect4th*maxValue4th && TotalConfidenceLeft3rd > 1)
		{
			 TotalSeq = maxKey3rdString.concat(TotalSeq);
			 YesNo0 = 1;
		}


		if (CountCorrect4th*maxValue4th > CountCorrect1st*maxValue && CountCorrect4th*maxValue4th > CountCorrect2nd*maxValue2nd && CountCorrect4th*maxValue4th > CountCorrect3rd*maxValue3rd && TotalConfidenceLeft4th > 1)
		{
			 TotalSeq = maxKey4thString.concat(TotalSeq);
			 YesNo0 = 1;
		}
	}

	if (YesNo0 == 0)
	{
		if (CountCorrect1st - 2 > CountCorrect2nd && CountCorrect1st - 2 > CountCorrect3rd && CountCorrect1st - 2 > CountCorrect4th && maxValue - 5 > maxValue2nd && maxValue - 5 > maxValue3rd && maxValue - 5 > maxValue4th)
		{
			 TotalSeq = maxKeyString.concat(TotalSeq);
			 YesNo1 = 1;
		}

		if (CountCorrect2nd - 2 > CountCorrect1st && CountCorrect2nd - 2 > CountCorrect3rd && CountCorrect2nd - 2 > CountCorrect4th && maxValue2nd - 5 > maxValue && maxValue2nd - 5 > maxValue3rd && maxValue2nd - 5 > maxValue4th)
		{
			 TotalSeq = maxKey2ndString.concat(TotalSeq);
			 YesNo1 = 1;
		}

		if (CountCorrect3rd - 2 > CountCorrect1st && CountCorrect3rd - 2 > CountCorrect2nd && CountCorrect3rd - 2 > CountCorrect4th && maxValue3rd - 5 > maxValue && maxValue3rd - 5 > maxValue2nd && maxValue3rd - 5 > maxValue4th)
		{
			 TotalSeq = maxKey3rdString.concat(TotalSeq);
			 YesNo1 = 1;
		}

		if (CountCorrect4th - 2 > CountCorrect1st && CountCorrect4th - 2 > CountCorrect2nd && CountCorrect4th - 2 > CountCorrect3rd && maxValue4th - 5 > maxValue && maxValue4th - 5 > maxValue2nd && maxValue4th - 5 > maxValue3rd)
		{
			 TotalSeq = maxKey4thString.concat(TotalSeq);
			 YesNo1 = 1;
		}
	}


	if (YesNo0 == 0 && YesNo1 == 0 && TotalConfidenceLeft > 200.0 | TotalConfidenceLeft2nd > 200.0 | TotalConfidenceLeft3rd > 200.0 | TotalConfidenceLeft4th > 200.0)
	{
		if (CountCorrect1st*TotalConfidenceLeft > 0.0 && CountCorrect2nd*TotalConfidenceLeft2nd <= 0.0 && CountCorrect3rd*TotalConfidenceLeft3rd <= 0.0 && CountCorrect4th*TotalConfidenceLeft4th <= 0.0 && CountCorrect1st*TotalConfidenceLeft > 50.0)
		{
			 TotalSeq = maxKeyString.concat(TotalSeq);
			 YesNo2 = 1;
		}

		if (CountCorrect1st*TotalConfidenceLeft <= 0.0 && CountCorrect2nd*TotalConfidenceLeft2nd > 0.0 && CountCorrect3rd*TotalConfidenceLeft3rd <= 0.0 && CountCorrect4th*TotalConfidenceLeft4th <= 0.0 && CountCorrect2nd*TotalConfidenceLeft2nd > 50.0)
		{
			 TotalSeq = maxKey2ndString.concat(TotalSeq);
			 YesNo2 = 1;
		}

		if (CountCorrect1st*TotalConfidenceLeft <= 0.0 && CountCorrect2nd*TotalConfidenceLeft2nd <= 0.0 && CountCorrect3rd*TotalConfidenceLeft3rd > 0.0 && CountCorrect4th*TotalConfidenceLeft4th <= 0.0 && CountCorrect3rd*TotalConfidenceLeft3rd > 50.0)
		{
			 TotalSeq = maxKey3rdString.concat(TotalSeq);
			 YesNo2 = 1;
		}

		if (CountCorrect1st*TotalConfidenceLeft <= 0.0 && CountCorrect2nd*TotalConfidenceLeft2nd <= 0.0 && CountCorrect3rd*TotalConfidenceLeft3rd <= 0.0 && CountCorrect4th*TotalConfidenceLeft4th > 0.0 && CountCorrect4th*TotalConfidenceLeft4th > 50.0)
		{
			 TotalSeq = maxKey4thString.concat(TotalSeq);
			 YesNo2 = 1;
		}
	}

	if (YesNo0 == 0 && YesNo1 == 0 && YesNo2 == 0)
	{
		if (CountCorrect1st > CountCorrect2nd && CountCorrect1st > CountCorrect3rd && CountCorrect1st > CountCorrect4th && TotalConfidenceLeft > 1)
		{
			TotalSeq = maxKeyString.concat(TotalSeq);
			YesNo3 = 1;
		}

		if (CountCorrect2nd > CountCorrect1st && CountCorrect2nd > CountCorrect3rd && CountCorrect2nd > CountCorrect4th && TotalConfidenceLeft2nd > 1)
		{
			TotalSeq = maxKey2ndString.concat(TotalSeq);
			YesNo3 = 1;
		}

		if (CountCorrect3rd > CountCorrect1st && CountCorrect3rd > CountCorrect2nd && CountCorrect3rd > CountCorrect4th && TotalConfidenceLeft3rd > 1)
		{
			TotalSeq = maxKey3rdString.concat(TotalSeq);
			YesNo3 = 1;
		}

		if (CountCorrect4th > CountCorrect1st && CountCorrect4th > CountCorrect2nd && CountCorrect4th > CountCorrect3rd && TotalConfidenceLeft4th > 1)
		{
			TotalSeq = maxKey4thString.concat(TotalSeq);
			YesNo3 = 1;
		}
	}

	if (YesNo0 == 0 && YesNo1 == 0 && YesNo2 == 0 && YesNo3 == 0 && TotalConfidenceLeft > 200.0 | TotalConfidenceLeft2nd > 200.0 | TotalConfidenceLeft3rd > 200.0 | TotalConfidenceLeft4th > 200.0)
	{
		if (CountCorrect1st*TotalConfidenceLeft > CountCorrect2nd*TotalConfidenceLeft2nd && CountCorrect1st*TotalConfidenceLeft > CountCorrect3rd*TotalConfidenceLeft3rd && CountCorrect1st*TotalConfidenceLeft > CountCorrect4th*TotalConfidenceLeft4th && CountCorrect1st*TotalConfidenceLeft > 50)
		{
			TotalSeq = maxKeyString.concat(TotalSeq);
			YesNo4 = 1;
		}

		if (CountCorrect2nd*TotalConfidenceLeft2nd > CountCorrect1st*TotalConfidenceLeft && CountCorrect2nd*TotalConfidenceLeft2nd > CountCorrect3rd*TotalConfidenceLeft3rd && CountCorrect2nd*TotalConfidenceLeft2nd > CountCorrect4th*TotalConfidenceLeft4th && CountCorrect2nd*TotalConfidenceLeft2nd > 50)
		{
			TotalSeq = maxKey2ndString.concat(TotalSeq);
			YesNo4 = 1;
		}

		if (CountCorrect3rd*TotalConfidenceLeft3rd > CountCorrect1st*TotalConfidenceLeft && CountCorrect3rd*TotalConfidenceLeft3rd > CountCorrect2nd*TotalConfidenceLeft2nd && CountCorrect3rd*TotalConfidenceLeft3rd > CountCorrect4th*TotalConfidenceLeft4th && CountCorrect3rd*TotalConfidenceLeft3rd > 50)
		{
			TotalSeq = maxKey3rdString.concat(TotalSeq);
			YesNo4 = 1;
		}

		if (CountCorrect4th*TotalConfidenceLeft4th > CountCorrect1st*TotalConfidenceLeft && CountCorrect4th*TotalConfidenceLeft4th > CountCorrect2nd*TotalConfidenceLeft2nd && CountCorrect4th*TotalConfidenceLeft4th > CountCorrect3rd*TotalConfidenceLeft3rd && CountCorrect4th*TotalConfidenceLeft4th > 50)
		{
			TotalSeq = maxKey4thString.concat(TotalSeq);
			YesNo4 = 1;
		}
	}

	if (YesNo0 == 0 && YesNo1 == 0 && YesNo2 == 0 && YesNo3 == 0 && YesNo4 == 0 && TotalConfidenceLeft > 200.0 | TotalConfidenceLeft2nd > 200.0 | TotalConfidenceLeft3rd > 200.0 | TotalConfidenceLeft4th > 200.0)
	{
		if (CountCorrect1st*TotalConfidenceLeft > CountCorrect2nd*TotalConfidenceLeft2nd && CountCorrect1st*TotalConfidenceLeft > CountCorrect3rd*TotalConfidenceLeft3rd && CountCorrect1st*TotalConfidenceLeft > CountCorrect4th*TotalConfidenceLeft4th)
		{
			TotalSeq = maxKeyString.concat(TotalSeq);
			YesNo5 = 1;
		}

		if (CountCorrect2nd*TotalConfidenceLeft2nd > CountCorrect1st*TotalConfidenceLeft && CountCorrect2nd*TotalConfidenceLeft2nd > CountCorrect3rd*TotalConfidenceLeft3rd && CountCorrect2nd*TotalConfidenceLeft2nd > CountCorrect4th*TotalConfidenceLeft4th)
		{
			TotalSeq = maxKey2ndString.concat(TotalSeq);
			YesNo5 = 1;
		}

		if (CountCorrect3rd*TotalConfidenceLeft3rd > CountCorrect1st*TotalConfidenceLeft && CountCorrect3rd*TotalConfidenceLeft3rd > CountCorrect2nd*TotalConfidenceLeft2nd && CountCorrect3rd*TotalConfidenceLeft3rd > CountCorrect4th*TotalConfidenceLeft4th)
		{
			TotalSeq = maxKey3rdString.concat(TotalSeq);
			YesNo5 = 1;
		}

		if (CountCorrect4th*TotalConfidenceLeft4th > CountCorrect1st*TotalConfidenceLeft && CountCorrect4th*TotalConfidenceLeft4th > CountCorrect2nd*TotalConfidenceLeft2nd && CountCorrect4th*TotalConfidenceLeft4th > CountCorrect3rd*TotalConfidenceLeft3rd)
		{
			TotalSeq = maxKey4thString.concat(TotalSeq);
			YesNo5 = 1;
		}

		if (CountCorrect4th*TotalConfidenceLeft4th == 0 && CountCorrect3rd*TotalConfidenceLeft3rd == 0 && CountCorrect2nd*TotalConfidenceLeft2nd == 0 && CountCorrect1st*TotalConfidenceLeft == 0)
		{
			if (CountCorrect1st > CountCorrect2nd && CountCorrect1st > CountCorrect3rd && CountCorrect1st > CountCorrect4th)
			{
				TotalSeq = maxKeyString.concat(TotalSeq);
				YesNo5 = 1;
			}

			if (CountCorrect2nd > CountCorrect1st && CountCorrect2nd > CountCorrect3rd && CountCorrect2nd > CountCorrect4th)
			{
				TotalSeq = maxKey2ndString.concat(TotalSeq);
				YesNo5 = 1;
			}

			if (CountCorrect3rd > CountCorrect1st && CountCorrect3rd > CountCorrect2nd && CountCorrect3rd > CountCorrect4th)
			{
				TotalSeq = maxKey3rdString.concat(TotalSeq);
				YesNo5 = 1;
			}

			if (CountCorrect4th > CountCorrect1st && CountCorrect4th > CountCorrect2nd && CountCorrect4th > CountCorrect3rd)
			{
				TotalSeq = maxKey4thString.concat(TotalSeq);
				YesNo5 = 1;
			}
		}

	}

	System.out.println("YesNo1 = " + YesNo1);
	System.out.println("YesNo2 = " + YesNo2);
	System.out.println("YesNo3 = " + YesNo3);
	System.out.println("YesNo4 = " + YesNo4);
	System.out.println("YesNo5 = " + YesNo5);

	if (YesNo0 == 0 && YesNo1 == 0 && YesNo2 == 0 && YesNo3 == 0 && YesNo4 == 0 && YesNo5 == 0)
	{
		StopLeft = 1;
	}

	if (TotalSeq.startsWith("KKKK") | TotalSeq.startsWith("LLLL") )
	{
		StopLeft = 1;
	}

	if (TotalSeq.substring(0,1).equals(maxKeyString)) 
	{
		ConfidenceOutputFreq = (String.valueOf(maxValue) + " ").concat(ConfidenceOutputFreq);
		CountCorrectOutput = (String.valueOf(CountCorrect1st.intValue()) + " ").concat(CountCorrectOutput);
		AvgConfidenceOutput = (String.valueOf(AvgConfidenceLeft.intValue()) + " ").concat(AvgConfidenceOutput);
		TotalConfidenceOutput = (String.valueOf(TotalConfidenceLeftUn.intValue()) + " ").concat(TotalConfidenceOutput); 
	}

	if (TotalSeq.substring(0,1).equals(maxKey2ndString))
	{
		ConfidenceOutputFreq = (String.valueOf(maxValue2nd) + " ").concat(ConfidenceOutputFreq);
		CountCorrectOutput = (String.valueOf(CountCorrect2nd.intValue()) + " ").concat(CountCorrectOutput);
		AvgConfidenceOutput = (String.valueOf(AvgConfidenceLeft2nd.intValue()) + " ").concat(AvgConfidenceOutput);
		TotalConfidenceOutput = (String.valueOf(TotalConfidenceLeft2ndUn.intValue()) + " ").concat(TotalConfidenceOutput); 
	}

	if (TotalSeq.substring(0,1).equals(maxKey3rdString))
	{
		ConfidenceOutputFreq = (String.valueOf(maxValue3rd) + " ").concat(ConfidenceOutputFreq);
		CountCorrectOutput = (String.valueOf(CountCorrect3rd.intValue()) + " ").concat(CountCorrectOutput);
		AvgConfidenceOutput = (String.valueOf(AvgConfidenceLeft3rd.intValue()) + " ").concat(AvgConfidenceOutput);
		TotalConfidenceOutput = (String.valueOf(TotalConfidenceLeft3rdUn.intValue()) + " ").concat(TotalConfidenceOutput); 
	}

	if (TotalSeq.substring(0,1).equals(maxKey4thString))
	{
		ConfidenceOutputFreq = (String.valueOf(maxValue4th) + " ").concat(ConfidenceOutputFreq);
		CountCorrectOutput = (String.valueOf(CountCorrect4th.intValue()) + " ").concat(CountCorrectOutput);
		AvgConfidenceOutput = (String.valueOf(AvgConfidenceLeft4th.intValue()) + " ").concat(AvgConfidenceOutput);
		TotalConfidenceOutput = (String.valueOf(TotalConfidenceLeft4thUn.intValue()) + " ").concat(TotalConfidenceOutput); 
	}

System.out.println(TotalSeq);	

}


//------------------------------------------------------
//Counting the frequency of the amino acids adjacent 
//to the seed sequence, for the right (C-terminal) side
//------------------------------------------------------


		ExtensionCountR = 0;

		for (int i = 0; i <= Peptide.length-1; i++)
		{
			if (Peptide[i] != null)
			{
			if (Peptide[i].contains(CurrentSeqR) == true && Peptide[i].indexOf(CurrentSeqR) + 5 <= Peptide[i].length() && ALC[i] > ALCthreshold)
			{
				ExtensionCountR++;
			}
			}
		}

		if (ExtensionCountR <= 5)
		{
			StopRight = 1;
		}

if (ExtensionCountR > 5 && StopRight == 0 && TotalSeq.length() < 100)
{
		int[] NumValuesRight = new int[CountR];
		char[] CharValuesRight = new char[CountR];

		for (int i = 0; i <= CountR-1; i++)
		{
			if (GrowthSeqRight[i] != null)
			{
				char Current1Right = GrowthSeqRight[i].charAt(GrowthSeqRight[i].length()-1);
				NumValuesRight[i] = Character.getNumericValue(Current1Right);
				CharValuesRight[i] = Current1Right;
			}
		}

		//System.out.println(Arrays.toString(NumValuesRight));
		//System.out.println(Arrays.toString(CharValuesRight));

		Map<Integer, Integer> m1 = new HashMap<Integer, Integer>();

	        int maxKeyR = -1;  
	        int maxValueR = -1;
	        int maxKey2ndR = -1;  
	        int maxValue2ndR = -1;  
	        int maxKey3rdR = -1;  
	        int maxValue3rdR = -1;
	        int maxKey4thR = -1;  
	        int maxValue4thR = -1;  

		for (int i = 0; i <= NumValuesRight.length-1; i++)
		{
			if (!m1.containsKey(NumValuesRight[i]))
			{
				m1.put(NumValuesRight[i], 1);
			}
			else
			{
				m1.put(NumValuesRight[i], m1.get(NumValuesRight[i])+1);
			}
		}

	        for (Map.Entry<Integer, Integer> entry2R : m1.entrySet())  
		{  
        		if (entry2R.getValue() > maxValueR) 
			{  
                		maxKeyR = entry2R.getKey();  
	        	        maxValueR = entry2R.getValue();  
		        }  
		}

	        for (Map.Entry<Integer, Integer> entry3R : m1.entrySet())  
		{  
        		if (entry3R.getKey() != maxKeyR && entry3R.getKey() != 0 && entry3R.getValue() > maxValue2ndR) 
			{  
                		maxKey2ndR = entry3R.getKey();  
	        	        maxValue2ndR = entry3R.getValue();  
		        }  
		}

	        for (Map.Entry<Integer, Integer> entry4R : m1.entrySet())  
		{  
        		if (entry4R.getKey() != maxKeyR && entry4R.getKey() != maxKey2ndR && entry4R.getKey() != 0 && entry4R.getValue() > maxValue3rdR) 
			{  
                		maxKey3rdR = entry4R.getKey();  
	        	        maxValue3rdR = entry4R.getValue();  
		        }  
		}

	        for (Map.Entry<Integer, Integer> entry5R : m1.entrySet())  
		{  
        		if (entry5R.getKey() != maxKeyR && entry5R.getKey() != maxKey2ndR && entry5R.getKey() != maxKey3rdR  && entry5R.getKey() != 0 && entry5R.getValue() > maxValue4thR) 
			{  
                		maxKey4thR = entry5R.getKey();  
	        	        maxValue4thR = entry5R.getValue();  
		        }  
		}

		System.out.println("The winner (right) is residue "+Character.forDigit(maxKeyR, Character.MAX_RADIX)+" its frequency of occurrence is "+maxValueR);  
		System.out.println("The 2nd place winner (Right) is residue "+Character.forDigit(maxKey2ndR, Character.MAX_RADIX)+" its frequency of occurrence is "+maxValue2ndR);  
		System.out.println("The 3rd place winner (Right) is residue "+Character.forDigit(maxKey3rdR, Character.MAX_RADIX)+" its frequency of occurrence is "+maxValue3rdR);  
		System.out.println("The 4th place winner (Right) is residue "+Character.forDigit(maxKey4thR, Character.MAX_RADIX)+" its frequency of occurrence is "+maxValue4thR);  

		String maxKeyStringR = String.valueOf(Character.forDigit(maxKeyR, Character.MAX_RADIX));
		maxKeyStringR = maxKeyStringR.toUpperCase();
		String maxKey2ndStringR = String.valueOf(Character.forDigit(maxKey2ndR, Character.MAX_RADIX));
		maxKey2ndStringR = maxKey2ndStringR.toUpperCase();
		String maxKey3rdStringR = String.valueOf(Character.forDigit(maxKey3rdR, Character.MAX_RADIX));
		maxKey3rdStringR = maxKey3rdStringR.toUpperCase();
		String maxKey4thStringR = String.valueOf(Character.forDigit(maxKey4thR, Character.MAX_RADIX));
		maxKey4thStringR = maxKey4thStringR.toUpperCase();

		int CountCorrect1stR = 0;
		int CountCorrect2ndR = 0;
		int CountCorrect3rdR = 0;
		int CountCorrect4thR = 0;

		for (int i = 0; i <= Peptide.length-1; i++)
		{
			if (Peptide[i] != null)
			{
				if (Peptide[i].contains(TotalSeq.substring(TotalSeq.length()-4, TotalSeq.length())) == true && Peptide[i].contains(CurrentSeqR.concat(maxKeyStringR)) == true)
				{
				CountCorrect1stR++;
				}

				if (Peptide[i].contains(TotalSeq.substring(TotalSeq.length()-4, TotalSeq.length())) == true && Peptide[i].contains(CurrentSeqR.concat(maxKey2ndStringR)) == true)
				{
				CountCorrect2ndR++;
				}

				if (Peptide[i].contains(TotalSeq.substring(TotalSeq.length()-4, TotalSeq.length())) == true && Peptide[i].contains(CurrentSeqR.concat(maxKey3rdStringR)) == true)
				{
				CountCorrect3rdR++;
				}

				if (Peptide[i].contains(TotalSeq.substring(TotalSeq.length()-4, TotalSeq.length())) == true && Peptide[i].contains(CurrentSeqR.concat(maxKey4thStringR)) == true)
				{
				CountCorrect4thR++;
				}
			}
		}

		System.out.println("CountCorrectRight for " + maxKeyStringR + " = " + CountCorrect1stR);
		System.out.println("CountCorrectRight for " + maxKey2ndStringR + " = " + CountCorrect2ndR);
		System.out.println("CountCorrectRight for " + maxKey3rdStringR + " = " + CountCorrect3rdR);
		System.out.println("CountCorrectRight for " + maxKey4thStringR + " = " + CountCorrect4thR);

		Double TotalConfidenceR = 0.0;
		Double TotalZerosR = 0.0;
		Double TotalHighR = 0.0;		

		for (int i = 0; i <= CountR-1; i++)
		{
			if (GrowthSeqRight[i] != null && GrowthSeqRight[i].endsWith(maxKeyStringR))
			{			
				TotalConfidenceR += Integer.valueOf(GrowthConfidenceRight[i][GrowthConfidenceRight[i].length-1]);

				if (Integer.valueOf(GrowthConfidenceRight[i][GrowthConfidenceRight[i].length-1]) == 0)
				{
				TotalZerosR++;
				}

				if (Integer.valueOf(GrowthConfidenceRight[i][GrowthConfidenceRight[i].length-1]) > 95)
				{
				TotalHighR++;
				}
			}
		}

		Double AvgConfidenceR = (double)TotalConfidenceR / maxValueR;
		TotalConfidenceR = TotalConfidenceR / (TotalZerosR+1) * TotalHighR;

		System.out.println(TotalConfidenceR);
		System.out.println(TotalZerosR);
		System.out.println(TotalHighR);

		Double TotalConfidenceR2nd = 0.0;
		Double TotalZerosR2nd = 0.0;
		Double TotalHighR2nd = 0.0;

		for (int i = 0; i <= CountR-1; i++)
		{
			if (GrowthSeqRight[i] != null && GrowthSeqRight[i].endsWith(maxKey2ndStringR))
			{
				TotalConfidenceR2nd += Integer.valueOf(GrowthConfidenceRight[i][GrowthConfidenceRight[i].length-1]);

				if (Integer.valueOf(GrowthConfidenceRight[i][GrowthConfidenceRight[i].length-1]) == 0)
				{
				TotalZerosR2nd++;
				}

				if (Integer.valueOf(GrowthConfidenceRight[i][GrowthConfidenceRight[i].length-1]) > 95)
				{
				TotalHighR2nd++;
				}
			}
		}

		Double AvgConfidenceR2nd = (double)TotalConfidenceR2nd / maxValue2ndR;
		TotalConfidenceR2nd = TotalConfidenceR2nd / (TotalZerosR2nd+1) * TotalHighR2nd;

		System.out.println(TotalConfidenceR2nd);
		System.out.println(TotalZerosR2nd);
		System.out.println(TotalHighR2nd);

		Double TotalConfidenceR3rd = 0.0;
		Double TotalZerosR3rd = 0.0;
		Double TotalHighR3rd = 0.0;

		for (int i = 0; i <= CountR-1; i++)
		{
			if (GrowthSeqRight[i] != null && GrowthSeqRight[i].endsWith(maxKey3rdStringR))
			{
				TotalConfidenceR3rd += Integer.valueOf(GrowthConfidenceRight[i][GrowthConfidenceRight[i].length-1]);

				if (Integer.valueOf(GrowthConfidenceRight[i][GrowthConfidenceRight[i].length-1]) == 0)
				{
				TotalZerosR3rd++;
				}

				if (Integer.valueOf(GrowthConfidenceRight[i][GrowthConfidenceRight[i].length-1]) > 95)
				{
				TotalHighR3rd++;
				}
			}
		}

		Double AvgConfidenceR3rd = (double)TotalConfidenceR3rd / maxValue3rdR;
		TotalConfidenceR3rd = TotalConfidenceR3rd / (TotalZerosR3rd+1) * TotalHighR3rd;

		System.out.println(TotalConfidenceR3rd);
		System.out.println(TotalZerosR3rd);
		System.out.println(TotalHighR3rd);

		Double TotalConfidenceR4th = 0.0;
		Double TotalZerosR4th = 0.0;
		Double TotalHighR4th = 0.0;

		for (int i = 0; i <= CountR-1; i++)
		{
			if (GrowthSeqRight[i] != null && GrowthSeqRight[i].endsWith(maxKey4thStringR))
			{
				TotalConfidenceR4th += Integer.valueOf(GrowthConfidenceRight[i][GrowthConfidenceRight[i].length-1]);

				if (Integer.valueOf(GrowthConfidenceRight[i][GrowthConfidenceRight[i].length-1]) == 0)
				{
				TotalZerosR4th++;
				}

				if (Integer.valueOf(GrowthConfidenceRight[i][GrowthConfidenceRight[i].length-1]) > 95)
				{
				TotalHighR4th++;
				}
			}
		}

		Double AvgConfidenceR4th = (double)TotalConfidenceR4th / maxValue4thR;
		TotalConfidenceR4th = TotalConfidenceR4th / (TotalZerosR4th+1) * TotalHighR4th;

		System.out.println(TotalConfidenceR4th);
		System.out.println(TotalZerosR4th);
		System.out.println(TotalHighR4th);

//------------------------------------------------------
//Concatenation of the right (C-terminal) residue
//ocurrs here
//------------------------------------------------------

	int YesNoR0 = 0;
	int YesNoR1 = 0;
	int YesNoR2 = 0;
	int YesNoR3 = 0;
	int YesNoR4 = 0;
	int YesNoR5 = 0;


	if (CountCorrect1stR > 10 | CountCorrect2ndR > 10 | CountCorrect3rdR > 10 | CountCorrect4thR > 10)
	{
		if (CountCorrect1stR*maxValueR > CountCorrect2ndR*maxValue2ndR && CountCorrect1stR*maxValueR > CountCorrect3rdR*maxValue3rdR && CountCorrect1stR*maxValueR > CountCorrect4thR*maxValue4thR && TotalConfidenceR > 1)
		{
			TotalSeq = TotalSeq.concat(maxKeyStringR);
			YesNoR0 = 1;
		}

		if (CountCorrect2ndR*maxValue2ndR > CountCorrect1stR*maxValueR && CountCorrect2ndR*maxValue2ndR > CountCorrect3rdR*maxValue3rdR && CountCorrect2ndR*maxValue2ndR > CountCorrect4thR*maxValue4thR && TotalConfidenceR2nd > 1)
		{
			TotalSeq = TotalSeq.concat(maxKey2ndStringR);
			YesNoR0 = 1;
		}

		if (CountCorrect3rdR*maxValue3rdR > CountCorrect1stR*maxValueR && CountCorrect3rdR*maxValue3rdR > CountCorrect2ndR*maxValue2ndR && CountCorrect3rdR*maxValue3rdR > CountCorrect4thR*maxValue4thR  && TotalConfidenceR3rd > 1)
		{
			TotalSeq = TotalSeq.concat(maxKey3rdStringR);
			YesNoR0 = 1;
		}

		if (CountCorrect4thR*maxValue4thR > CountCorrect1stR*maxValueR && CountCorrect4thR*maxValue4thR > CountCorrect2ndR*maxValue2ndR && CountCorrect4thR*maxValue3rdR > CountCorrect4thR*maxValue3rdR && TotalConfidenceR4th > 1)
		{
			TotalSeq = TotalSeq.concat(maxKey4thStringR);
			YesNoR0 = 1;
		}
	}	

	if (YesNoR0 == 0)
	{
		if (CountCorrect1stR - 2 > CountCorrect2ndR && CountCorrect1stR - 2 > CountCorrect3rdR && CountCorrect1stR - 2 > CountCorrect4thR && maxValueR - 5 > maxValue2ndR && maxValueR - 5 > maxValue3rdR && maxValueR - 5 > maxValue4thR)
		{
			TotalSeq = TotalSeq.concat(maxKeyStringR);
			YesNoR1 = 1;
		}

		if (CountCorrect2ndR - 2 > CountCorrect1stR && CountCorrect2ndR - 2 > CountCorrect3rdR && CountCorrect2ndR - 2 > CountCorrect4thR && maxValue2ndR - 5 > maxValueR && maxValue2ndR - 5 > maxValue3rdR && maxValue2ndR - 5 > maxValue4thR)
		{
			TotalSeq = TotalSeq.concat(maxKey2ndStringR);
			YesNoR1 = 1;
		}

		if (CountCorrect3rdR - 2 > CountCorrect1stR && CountCorrect3rdR - 2 > CountCorrect2ndR && CountCorrect3rdR - 2 > CountCorrect4thR && maxValue3rdR - 5 > maxValueR && maxValue3rdR - 5 > maxValue2ndR && maxValue3rdR - 5 > maxValue4thR)
		{
			TotalSeq = TotalSeq.concat(maxKey3rdStringR);
			YesNoR1 = 1;
		}

		if (CountCorrect4thR - 2 > CountCorrect1stR && CountCorrect4thR - 2 > CountCorrect2ndR && CountCorrect4thR - 2 > CountCorrect3rdR && maxValue4thR - 5 > maxValueR && maxValue4thR - 5 > maxValue2ndR && maxValue4thR - 5 > maxValue3rdR)
		{
			TotalSeq = TotalSeq.concat(maxKey4thStringR);
			YesNoR1 = 1;
		}
	}


	if (YesNoR0 == 0 && YesNoR1 == 0 && TotalConfidenceR > 200.0 | TotalConfidenceR2nd > 200.0 | TotalConfidenceR3rd > 200.0 && TotalConfidenceR4th > 200.0)
	{
		if (CountCorrect1stR*TotalConfidenceR > 0.0 && CountCorrect2ndR*TotalConfidenceR2nd <= 0.0 && CountCorrect3rdR*TotalConfidenceR3rd <= 0.0 && CountCorrect4thR*TotalConfidenceR4th <= 0.0)
		{
			 TotalSeq = TotalSeq.concat(maxKeyStringR);
			 YesNoR2 = 1;
		}

		if (CountCorrect1stR*TotalConfidenceR <= 0.0 && CountCorrect2ndR*TotalConfidenceR2nd > 0.0 && CountCorrect3rdR*TotalConfidenceR3rd <= 0.0 && CountCorrect4thR*TotalConfidenceR4th <= 0.0)
		{
			 TotalSeq = TotalSeq.concat(maxKey2ndStringR);
			 YesNoR2 = 1;
		}

		if (CountCorrect1stR*TotalConfidenceR <= 0 && CountCorrect2ndR*TotalConfidenceR2nd <= 0.0 && CountCorrect3rdR*TotalConfidenceR3rd > 0.0 && CountCorrect4thR*TotalConfidenceR4th <= 0.0)
		{
			 TotalSeq = TotalSeq.concat(maxKey3rdStringR);
			 YesNoR2 = 1;
		}

		if (CountCorrect1stR*TotalConfidenceR <= 0 && CountCorrect2ndR*TotalConfidenceR2nd <= 0.0 && CountCorrect3rdR*TotalConfidenceR3rd <= 0.0 && CountCorrect4thR*TotalConfidenceR4th > 0.0)
		{
			 TotalSeq = TotalSeq.concat(maxKey4thStringR);
			 YesNoR2 = 1;
		}
	}

	if (YesNoR0 == 0 && YesNoR1 == 0 && YesNoR2 == 0)
	{
		if (CountCorrect1stR > CountCorrect2ndR && CountCorrect1stR > CountCorrect3rdR && CountCorrect1stR > CountCorrect4thR && TotalConfidenceR > 1)
		{
			TotalSeq = TotalSeq.concat(maxKeyStringR);
			YesNoR3 = 1;
		}

		if (CountCorrect2ndR > CountCorrect1stR && CountCorrect2ndR > CountCorrect3rdR && CountCorrect2ndR > CountCorrect4thR && TotalConfidenceR2nd > 1)
		{
			TotalSeq = TotalSeq.concat(maxKey2ndStringR);
			YesNoR3 = 1;
		}

		if (CountCorrect3rdR > CountCorrect1stR && CountCorrect3rdR > CountCorrect2ndR && CountCorrect3rdR > CountCorrect4thR && TotalConfidenceR3rd > 1)
		{
			TotalSeq = TotalSeq.concat(maxKey3rdStringR);
			YesNoR3 = 1;
		}

		if (CountCorrect4thR > CountCorrect1stR && CountCorrect4thR > CountCorrect2ndR && CountCorrect4thR > CountCorrect3rdR && TotalConfidenceR4th > 1)
		{
			TotalSeq = TotalSeq.concat(maxKey4thStringR);
			YesNoR3 = 1;
		}
	}

	if (YesNoR0 == 0 && YesNoR1 == 0 && YesNoR2 == 0 && YesNoR3 == 0 && TotalConfidenceR > 200.0 | TotalConfidenceR2nd > 200.0 | TotalConfidenceR3rd > 200.0 && TotalConfidenceR4th > 200.0)
	{
		if (CountCorrect1stR*TotalConfidenceR > CountCorrect2ndR*TotalConfidenceR2nd && CountCorrect1stR*TotalConfidenceR > CountCorrect3rdR*TotalConfidenceR3rd && CountCorrect1stR*TotalConfidenceR > CountCorrect4thR*TotalConfidenceR4th)
		{
			TotalSeq = TotalSeq.concat(maxKeyStringR);
			YesNoR4 = 1;
		}

		if (CountCorrect2ndR*TotalConfidenceR2nd > CountCorrect1stR*TotalConfidenceR && CountCorrect2ndR*TotalConfidenceR2nd > CountCorrect3rdR*TotalConfidenceR3rd && CountCorrect2ndR*TotalConfidenceR2nd > CountCorrect4thR*TotalConfidenceR4th)
		{
			TotalSeq = TotalSeq.concat(maxKey2ndStringR);
			YesNoR4 = 1;
		}

		if (CountCorrect3rdR*TotalConfidenceR3rd > CountCorrect1stR*TotalConfidenceR && CountCorrect3rdR*TotalConfidenceR3rd > CountCorrect2ndR*TotalConfidenceR2nd && CountCorrect3rdR*TotalConfidenceR3rd > CountCorrect4thR*TotalConfidenceR4th)
		{
			TotalSeq = TotalSeq.concat(maxKey3rdStringR);
			YesNoR4 = 1;
		}

		if (CountCorrect4thR*TotalConfidenceR4th > CountCorrect1stR*TotalConfidenceR && CountCorrect4thR*TotalConfidenceR4th > CountCorrect2ndR*TotalConfidenceR2nd && CountCorrect4thR*TotalConfidenceR4th > CountCorrect3rdR*TotalConfidenceR3rd)
		{
			TotalSeq = TotalSeq.concat(maxKey4thStringR);
			YesNoR4 = 1;
		}

	}

	if (YesNoR0 == 0 && YesNoR1 == 0 && YesNoR2 == 0 && YesNoR3 == 0 && YesNoR4 == 0)
	{
	if (CountCorrect4thR*TotalConfidenceR4th == 0 && CountCorrect1stR*TotalConfidenceR3rd == 0 && CountCorrect2ndR*TotalConfidenceR2nd == 0 && CountCorrect1stR*TotalConfidenceR == 0)
	{
		if (CountCorrect1stR > CountCorrect2ndR && CountCorrect1stR > CountCorrect3rdR && CountCorrect1stR > CountCorrect4thR)
		{
			TotalSeq = TotalSeq.concat(maxKeyStringR);
			YesNoR5 = 1;
		}

		if (CountCorrect2ndR > CountCorrect1stR && CountCorrect2ndR > CountCorrect3rdR && CountCorrect2ndR > CountCorrect4thR)
		{
			TotalSeq = TotalSeq.concat(maxKey2ndStringR);
			YesNoR5 = 1;
		}

		if (CountCorrect3rdR > CountCorrect1stR && CountCorrect3rdR > CountCorrect2ndR && CountCorrect3rdR > CountCorrect4thR)
		{
			TotalSeq = TotalSeq.concat(maxKey3rdStringR);
			YesNoR5 = 1;
		}

		if (CountCorrect4thR > CountCorrect1stR && CountCorrect4thR > CountCorrect2ndR && CountCorrect4thR > CountCorrect3rdR)
		{
			TotalSeq = TotalSeq.concat(maxKey4thStringR);
			YesNoR5 = 1;
		}
	}
	}

	System.out.println("YesNoR0 = " + YesNoR0);
	System.out.println("YesNoR1 = " + YesNoR1);
	System.out.println("YesNoR2 = " + YesNoR2);
	System.out.println("YesNoR3 = " + YesNoR3);
	System.out.println("YesNoR4 = " + YesNoR4);
	System.out.println("YesNoR5 = " + YesNoR5);

	if (YesNoR0 == 0 && YesNoR1 == 0 && YesNoR2 == 0 && YesNoR3 == 0 && YesNoR4 == 0 && YesNoR5 == 0)
	{
		StopRight = 1;
	}

	if (TotalSeq.endsWith("KKKK") | TotalSeq.endsWith("LLLL") )
	{
		StopRight = 1;
	}
	
	if (TotalSeq.substring(TotalSeq.length()-1).equals(maxKeyStringR)) 
	{
		ConfidenceOutputFreq = ConfidenceOutputFreq.concat(" " + String.valueOf(maxValueR));
		CountCorrectOutput = CountCorrectOutput.concat(" " + String.valueOf(CountCorrect1stR) );
		AvgConfidenceOutput = AvgConfidenceOutput.concat(" " + String.valueOf(AvgConfidenceR.intValue()) );
		TotalConfidenceOutput = TotalConfidenceOutput.concat(" " + String.valueOf(TotalConfidenceR.intValue()) ); 
	}

	if (TotalSeq.substring(TotalSeq.length()-1).equals(maxKey2ndStringR))
	{
		ConfidenceOutputFreq = ConfidenceOutputFreq.concat(" " + String.valueOf(maxValue2ndR));
		CountCorrectOutput = CountCorrectOutput.concat(" " + String.valueOf(CountCorrect2ndR) );
		AvgConfidenceOutput = AvgConfidenceOutput.concat(" " + String.valueOf(AvgConfidenceR2nd.intValue()) );
		TotalConfidenceOutput = TotalConfidenceOutput.concat(" " + String.valueOf(TotalConfidenceR2nd.intValue()) );
	}

	if (TotalSeq.substring(TotalSeq.length()-1).equals(maxKey3rdStringR))
	{
		ConfidenceOutputFreq = ConfidenceOutputFreq.concat(" " + String.valueOf(maxValue3rdR));
		CountCorrectOutput = CountCorrectOutput.concat(" " + String.valueOf(CountCorrect3rdR) );
		AvgConfidenceOutput = AvgConfidenceOutput.concat(" " + String.valueOf(AvgConfidenceR3rd.intValue()) );
		TotalConfidenceOutput = TotalConfidenceOutput.concat(" " + String.valueOf(TotalConfidenceR3rd.intValue()) );
	}

	if (TotalSeq.substring(TotalSeq.length()-1).equals(maxKey4thStringR))
	{
		ConfidenceOutputFreq = ConfidenceOutputFreq.concat(" " + String.valueOf(maxValue4thR));
		CountCorrectOutput = CountCorrectOutput.concat(" " + String.valueOf(CountCorrect4thR) );
		AvgConfidenceOutput = AvgConfidenceOutput.concat(" " + String.valueOf(AvgConfidenceR4th.intValue()) );
		TotalConfidenceOutput = TotalConfidenceOutput.concat(" " + String.valueOf(TotalConfidenceR4th.intValue()) );
	}

System.out.println(TotalSeq);
}

	System.out.println ("StopLeft = " + StopLeft + "    " + "StopRight = " + StopRight);

	if ((ExtensionCountR <= 5 && ExtensionCountL <= 5) | (StopLeft != 0 && StopRight != 0) | (TotalSeq.length() > 100))
	{
	break;
	}

}

//------------------------------------------------------
//Printing the Output
//------------------------------------------------------

	if (TotalSeq.contains("LLLL") == false && TotalSeq.length() > 20)
	{
		System.out.println("TotalSeq: " + TotalSeq);

		p.println("#" + I + " of " + SeedSeqs.length);
		p.println(">" + SeedSeq);
		p.println("Redundancy: " + SeedRedundancy);
		p.println("Seed Sequence Average Residue Confidence:");
		p.println((int)SeedSeqsConfidenceOutput[0] + " " + (int)SeedSeqsConfidenceOutput[1] + " " + (int)SeedSeqsConfidenceOutput[2] + " " + (int)SeedSeqsConfidenceOutput[3] + " " + (int)SeedSeqsConfidenceOutput[4]);
		p.println("Sequence:");
		p.println(TotalSeq);
		p.println("Sequence (Delimited):");

		try{

		for (int i = 0; i <= TotalConfidenceOutput.split(" ").length-1; i++)
		{
			if (TotalSeq.length() > 5)
			{

			if (TotalConfidenceOutput.split(" ")[i].equals("X"))
			{
				if ( !(TotalSeq.split("")[i] + TotalSeq.split("")[i+1] + TotalSeq.split("")[i+2] + TotalSeq.split("")[i+3] + TotalSeq.split("")[i+4]).equals(SeedSeq))
				{
					p.print("- ");
				}
			break;
			} 
			}
		}

		for(int i = 0; i <= TotalSeq.length()-1; i++)
		{
			p.print(TotalSeq.charAt(i) + " ");
		}

		} catch (Exception e) {System.out.println(e);}

		p.println(" ");
		p.println("Average MS Confidence Score:");
		p.println(AvgConfidenceOutput);
		p.println("Total MS Confidence Score:");
		p.println(TotalConfidenceOutput);
		p.println("Frequency of residue occurrence (Window of 4 residues):");
		p.println(ConfidenceOutputFreq);
		p.println("Frequency of residue occurrence (Window of 5 residues):");
		p.println(CountCorrectOutput);
		p.println(" ");
	}
}
}

	p.close();

	}
}