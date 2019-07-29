

def runETQL(vcfFile, pheFile, outputName = None):
    """Function that runs fastEQTL on the given vcf and phenotype file
    E.g. command: 
        'fastQTL --vcf genotypes.vcf.gz --bed phenotypes.bed.gz --region \
        22:17000000-18000000 --out nominals.default.txt.gz' """
    if outputName is not None:
        os.system("fastQTL --vcf {vcfFile} --bed phenotypes.bed.gz \
            --region 22:17000000-18000000 --out {outName}".format(
            vcfFile = vcfFile, pheFile = pheFile, outName = outputName)


def getVCF(id):
    """Gets the name of the vcf file
        corresponding to the 
        > getVCF("SRR1202488")
        "SRR1202488.cleaned_1.vcf"
        """
    return id + ".cleaned_1.vcf"
def getCounts(id):
    """"""
    return id + ".phe" # TODO: change with output you generated

def getOutputName(id):
    return id + ".etql.results"
    
for id in idList:
    runETQL(getVCF(id), getCounts(id), outputName = getOutputName(id))
    