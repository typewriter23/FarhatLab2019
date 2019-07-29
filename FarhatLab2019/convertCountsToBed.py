import pandas as pd
import os
import ntpath

# This script provides functions that will take in a 
    # counts.tsv file (as outputted by featurecounts)
    # and converts it to a bed file (as inputted by fastQTL)
    
countsFolder = "/n/scratch2/jy250/countsFolder"
outputFile = "/n/scratch2/jy250/testCountsToBedOutput.bed"
# Want to create a large bed file containing all the info
    # of each counts file in the folder
    
# Going to have to play with the filename to generate IDs
listOfCounts = []

for filename in os.listdir(countsFolder):
    if filename.endswith(".counts.tsv"): 
        listOfCounts.append(filename)

def genBase(sampleFile):
    """Given a sampleFile, generates the left side of the   
        bedfile to be taken in by fastQTL
        
        samplefile should be of the form:
        # Program:featureCounts v1.6.4; Command:"featureCounts" "-t" "gene" "-g" "locus_tag" "-a" ...
        Geneid  Chr     Start   End     Strand  Length  /n/scratch2/jy250/test100/SRR1202642.cleaned_1.bam
        Rv0001  NC_000962.3     1       1524    +       1524    338
        Rv0002  NC_000962.3     2052    3260    +       1209    184
        Rv0003  NC_000962.3     3280    4437    +       1158    67
        Rv0004  NC_000962.3     4434    4997    +       564     67
        Rv0005  NC_000962.3     5240    7267    +       2028    891
        ...
        
        Returns a dataframe of the form:
       Geneid          Chr    Start      End Strand  Length
0      Rv0001  NC_000962.3        1     1524      +    1524
1      Rv0002  NC_000962.3     2052     3260      +    1209
2      Rv0003  NC_000962.3     3280     4437      +    1158
3      Rv0004  NC_000962.3     4434     4997      +     564
4      Rv0005  NC_000962.3     5240     7267      +    2028
5      Rv0006  NC_000962.3     7302     9818      +    2517
6      Rv0007  NC_000962.3     9914    10828      +     915
7      Rvnt01  NC_000962.3    10887    10960      +      74
8      Rvnt02  NC_000962.3    11112    11184      +      73
        ...
        """
    df = pd.read_csv(sampleFile, delimiter = "\t", header = 1) # Header => ignore the first line
    # print(df.columns)
    df = df.drop("Strand", axis = 1)
    df = df.drop("Length", axis = 1)
    df = df[['Chr', 'Start', 'End', 'Geneid']]
    return df   # .drop(df.columns[-1], axis = 1) #Returns a dataframe with everything but the last column
                                             # ASSUMES LAST COLUMN IS WHERE THE DATA IS
def averageCol(col, sampleID):
    """Returns a 2 column by ROW rows dataframe
        with all values of counts normalized"""
    sum = col[sampleID].sum()
    def normalize(x):
        return x / sum
        
    # print(col)
    if(sum == 0):
        return None
    col[sampleID] = col[sampleID].apply(normalize)
    return col
    
def getCol(countsFile):
    """Gets the column of counts from the specified counts file"""
    df = pd.read_csv(countsFile, delimiter = "\t", header = 1)
    col = df[[df.columns[-1], "Geneid"]]
    sampleID = genSampleID(countsFile)
    col = col.rename(columns = {df.columns[-1]: sampleID})
    return averageCol(col, sampleID)
    # return pd.concat([master, df], axis=1)
    
def genSampleID(path):
    """Generates the sampleID from the filename given
        Code partially taken from StackOverflow:
        https://stackoverflow.com/questions/8384737/extract-file-name-from-path-no-matter-what-the-os-path-format"""
    head, tail = ntpath.split(path)
    result =  tail or ntpath.basename(head)
    return result.split(".")[0] # Gets just the sample name, cleans out the ".cleaned.[EXT]"
    
df = genBase(os.path.join(countsFolder, listOfCounts[0]))        

for countsFile in listOfCounts:
    # print(countsFile)
    countsLoc = os.path.join(countsFolder, countsFile)
    sampleID = genSampleID(countsLoc) # gets the sampleID from the fileName
    col = getCol(countsLoc)
    # print(col)
    if col is not None:
        df = df.merge(col, left_on = "Geneid", right_on = "Geneid", how = "outer") # combine the columns of the tables
#  = df.drop("Geneid", axis = 1)
    # df.loc[sampleID] = getCol(countsLoc) # appends the column of counts to the end of the dataframe
    
# print(df)

df.to_csv(outputFile, sep = "\t", index = False)
    