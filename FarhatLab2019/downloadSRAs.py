# Theoretically, this script should be able to be modified to run any command line exe file

import sys, os
# # ONLY FOR COMPUTERS THAT ARE {1.) JERRY'S, 2.) WINDOWS} 
# bashDir = "C:/Users/typew/Documents/GitHub/FarhatLab/sratoolkit.2.9.6-1-win64/bin/{0}" # 0 is placeholder for exe
# exe = "fastq-dump.exe {0}" # 0 is placeholder for SRA run id
# # idListLocation = "C:/Users/typew/Documents/GitHub/FarhatLab/3SamplePaperRunIDs.txt"
# idListLocation = "C:/Users/typew/Documents/GitHub/FarhatLab/test.txt"
# idList = []
# # ---------- END Windows information -------------------------------

# idListLocation = "/n/data1/hms/dbmi/farhat/Jerry/3SamplePaperRunIDs.txt" # Location of list of SRA IDs
idListLocation = "1000SampleRunIDs.txt" # Location of list of SRA IDs
# ONLY FOR COMPUTERS THAT ARE {1.) Linux, 2.) O2} 
cwd = os.getcwd() # Gets the current directory in O2 (under the assumption that my data folder will be relatively
                        # close to my current folder
                        

#outputFolder = "/n/data1/hms/dbmi/farhat/Jerry/inputs/"
outputFolder = "/n/scratch2/jy250/1000SampleTest"
    # "/n/data1/hms/dbmi/farhat/Jerry/3SampleFastqs" # Send to much larger scratch folder
bashDir = "/n/data1/hms/dbmi/farhat/Jerry/FarhatLab2019/sratoolkit.2.9.6-1-centos_linux64/bin/{0}" # Location of sratoolkit's fastq-dump.exe command
                                                                                                    # {0} is placeholder for exe file
slurmCommand = 'sbatch -p short -c 1 -t 0-11:59 --mem=1G --mail-user=Jerry_Yang@hms.harvard.edu --wrap="{0}"'
exe = "fastq-dump.2.9.6 --split-files --outdir '{1}' {0} " 
        # {0} is placeholder for SRA run id
        # {1} is a placeholder for output directory
# idListLocation = "C:/Users/typew/Documents/GitHub/FarhatLab/3SamplePaperRunIDs.txt"

idList = []

# Import list of SRA ids
# IO info taken from here:
# https://docs.python.org/3/tutorial/inputoutput.html
with open(idListLocation) as idFile:
    for line in idFile:
        idList.append(line)

# idList is now a list of IDs

# Now, want to loop thorugh every id and run SRA toolkit/fastq-dump.exe on all of them 

for id in idList:
    bashCommand = bashDir.format(exe).format(id, outputFolder)
    sCommand = slurmCommand.format(bashCommand)
    print(sCommand)
    os.system(sCommand)