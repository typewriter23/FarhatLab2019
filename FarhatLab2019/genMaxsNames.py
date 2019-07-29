import os

inputFolder = "/n/scratch2/jy250/fastqFiles"
dict = {}
for filename in os.listdir(inputFolder):
    if filename.endswith(".vcf"): 
        sample = os.path.splitext(filename)[0]
        dict[sample] = os.path.join(inputFolder, filename)
print(dict)