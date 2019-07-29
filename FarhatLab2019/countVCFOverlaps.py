import os
from JerryPipeToolbox import *
import re
import sys

outputFolder = "/n/data1/hms/dbmi/farhat/Jerry/VCFOverlaps"
def countOverlaps(vcfX, vcfY):
    """Counts the number of overlaps between vcfX and vcfY
        using vcftools"""
    # First run vcftools
    summaryFile = runVCFTools(vcfX, vcfY)
    
    # Next, extract from the output using Regular expressions
        # https://stackoverflow.com/questions/4666973/how-to-extract-the-substring-between-two-markers


    # text = '... Found 4411535 sites common to both files. ...'
    
    
    #String to parse for X union Y
    allString = 'Found (.+?) sites common to both files.'
    
    #String to parse ~X intersect Y
    xxString = "Found (.+?) sites only in main file."
    
    # String to parse for X intersect ~Y
    yyString = "Found (.+?) sites only in second file."
    
    # String ro parse for X intersect Y
    xyString = "Found (.+?) non-matching overlapping sites."
    
    summaryString = ""
    with open(summaryFile, 'r') as f:
        summaryString = f.read()
    
    xy = int(parseSummary(summaryString, xyString), base = 10)
    xx = int(parseSummary(summaryString, xxString), base = 10)
    yy = int(parseSummary(summaryString, yyString), base = 10)
    all = int(parseSummary(summaryString, allString), base = 10) # Sites common to X, Y and Reference
    return {"both": xy,
            "onlyX": xx,
            "onlyY": yy,
            "all": all}
            
def countOverlapsVCFGZVCF(vcf = "", gzvcf = ""):
    """Counts the number of overlaps between vcf and gvcf
        using vcftools"""
    # First run vcftools
    summaryFile = runGZVCFTools(vcf, gzvcf)
    
    # Next, extract from the output using Regular expressions
        # https://stackoverflow.com/questions/4666973/how-to-extract-the-substring-between-two-markers


    # text = '... Found 4411535 sites common to both files. ...'
    
    
    #String to parse for X union Y
    allString = 'Found (.+?) sites common to both files.'
    
    #String to parse ~X intersect Y
    xxString = "Found (.+?) sites only in main file."
    
    # String to parse for X intersect ~Y
    yyString = "Found (.+?) sites only in second file."
    
    # String ro parse for X intersect Y
    xyString = "Found (.+?) non-matching overlapping sites."
    
    summaryString = ""
    with open(summaryFile, 'r') as f:
        summaryString = f.read()
    
    xy = int(parseSummary(summaryString, xyString), base = 10)
    xx = int(parseSummary(summaryString, xxString), base = 10)
    yy = int(parseSummary(summaryString, yyString), base = 10)
    all = int(parseSummary(summaryString, allString), base = 10) # Sites common to X, Y and Reference
    return {"both": xy,
            "onlyX": xx,
            "onlyY": yy,
            "all": all}
            
def genOutputName(vcfX, vcfY):
    xBaseName = genBaseName(os.path.basename(vcfX))
    yBaseName = genBaseName(os.path.basename(vcfY))
    return os.path.join(outputFolder, "{0}And{1}VCFIntersects".format(xBaseName, yBaseName))
    
    
def runGZVCFTools(vcf, gzvcf, outputFileName = None):
    """ (Adaptation of runVCFTools that takes in one
            .gzvcf)
        Given the filenames (in the form of absolute paths) to 
            vcf and gzvcf
       - Runs VCFtools
       - and returns the output file to be parsed"""
    if outputFileName == None:
        outputFileName = genOutputName(vcf, gzvcf)
    vcfCommand = "vcftools --vcf {vcf} --gzdiff {gzvcf} --diff-site --out {outputFile}".format(
                            vcf = vcf, gzvcf = gzvcf, outputFile = outputFileName)
    print(vcfCommand)
    os.system(vcfCommand)
    return outputFileName + ".log" #TODO: Unclean!!
    
def runVCFTools(vcfX, vcfY, outputFileName = None):
    """Given the filenames (in the form of absolute paths) to 
            vcfX and vcfY
       - Runs VCFtools
       - and returns the output file to be parsed"""
    if outputFileName == None:
        outputFileName = genOutputName(vcfX, vcfY)
    vcfCommand = "vcftools --vcf {vcfX} --diff {vcfY} --diff-site --out {outputFile}".format(
                            vcfX = vcfX, vcfY = vcfY, outputFile = outputFileName)
    print(vcfCommand)
    os.system(vcfCommand)
    return outputFileName+ ".log" #TODO: Unclean!!
    
def parseSummary(summaryStr, targetStr):
    """Parses a long string for a particular amount of text"""
    parsed = re.search(targetStr, summaryStr)
    A = ""
    if parsed:
        A = parsed.group(1)
    return A
    
# Test command: 
# x = countOverlaps("/n/data1/hms/dbmi/farhat/Jerry/3SampleFastqs/SRR6176439_1.cleaned.vcf", 
    # "/n/data1/hms/dbmi/farhat/Jerry/3SampleFastqs/SRR6176436_1.cleaned.vcf")
# print(x)