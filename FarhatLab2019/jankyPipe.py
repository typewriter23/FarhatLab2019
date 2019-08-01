from JerryPipeToolbox import *


def Launch_JankyPipe(fqf1 , fqf2 , tag , output_dir , scratch_dir , O2_SLURM_logs_dir):

    
    '''
    This script launches a job to call variants for the input fastq files against H37Rv
    using a number of packages. The important output (VCF, lineage info files, quality report)
    is stored in the output directory while the intermediary files (SAMs, trimmed fastqs, BAM, etc)
    are stored in a scratch directory.
    '''
    
    #store all commands in a list
    commands_list = []
    
    #change directory to scratch
    commands_list.append( 'cd ' + scratch_dir )

    ###################################
    ### Load Necessary Modules ########
    ###################################

    
    # load gcc
    commands_list.append( 'module load gcc/6.2.0' )
    
    #load perl
    commands_list.append( 'module load perl/5.24.0' )

    #load java
    commands_list.append( 'module load java/jdk-1.8u112' )

    #load BWA
    commands_list.append( 'module load bwa/0.7.15' )

    #load Samtools
    commands_list.append( 'module load samtools/1.3.1' )

    #load BCFtools
    commands_list.append( 'module load bcftools/1.3.1' )

    #load Picard
    commands_list.append( 'module load picard/2.8.0' )

    
    
    
    #Create Index files for Reference Genome
    commands_list.append( 'mkdir RefGen' )

    #copy reference genome over to RefGen folder
    commands_list.append( 'cp /n/data1/hms/dbmi/farhat/mm774/References/GCF_000195955.2_ASM19595v2_genomic.fna RefGen/TBRefGen.fasta' )

    #change directory to RefGen folder
    commands_list.append( 'cd RefGen' )

    ###################################
    ### Create Index Files for H37Rv ##
    ###################################
    commands_list.append( 'samtools faidx TBRefGen.fasta' )
    commands_list.append( 'bwa index TBRefGen.fasta' )

    RefGen = scratch_dir + '/RefGen/TBRefGen.fasta' #H37Rv reference

    #go back to parent directory
    commands_list.append( 'cd ..' )

    ###################################
    ### UnZip FastQ files #############
    ###################################
    fqf1_base_name = fqf1.split('/')[-1][0:-9]
    fqf2_base_name = fqf2.split('/')[-1][0:-9]

    #work with the unzipped files for the rest of the pipeline (after unzipping them)
    fqf1_unzipped = scratch_dir + '/{}'.format(fqf1_base_name) + '.fastq'
    fqf2_unzipped = scratch_dir + '/{}'.format(fqf2_base_name) + '.fastq'

    commands_list.append( 'zcat {0} > {1}'.format(fqf1, fqf1_unzipped) )
    commands_list.append( 'zcat {0} > {1}'.format(fqf2, fqf2_unzipped) )

    #use the unzipped fastq files now
    fqf1 = fqf1_unzipped
    fqf2 = fqf2_unzipped
    
    ###################################
    ### Clean FastQ read names ########
    ###################################
    
    #delete any weird caracters from the read names in the FastQ files
    commands_list.append( 'python /n/data1/hms/dbmi/farhat/lfreschi/repos/megapipe/bin/megapipe-correct-names-reads.py {}'.format(fqf1) )
    commands_list.append( 'python /n/data1/hms/dbmi/farhat/lfreschi/repos/megapipe/bin/megapipe-correct-names-reads.py {}'.format(fqf2) )

    ####################################
    ### PRINSEQ (trim reads) ##########
    ###################################

    #create directory for prinseq in output directory
    commands_list.append( 'mkdir ' + output_dir + '/prinseq' )

    commands_list.append( 'perl /n/data1/hms/dbmi/farhat/bin/prinseq-lite-0.20.4/prinseq-lite.pl -fastq {0} -fastq2 {1} -out_format 3 -out_good {2}/{3}-trimmed -out_bad null -log {4}/{3}-trimmed.log -min_qual_mean 20 -verbose'.format(fqf1, fqf2, scratch_dir, tag , output_dir+'/prinseq') )

    #use newly trimmed fastq files now
    fqf1 = scratch_dir + '/{}'.format(tag) + '-trimmed_1.fastq'
    fqf2 = scratch_dir + '/{}'.format(tag) + '-trimmed_2.fastq'

    ######################################
    ### BWA (align reads to reference) ###
    ######################################

    #create SAM file
    samfile = scratch_dir + '/{}.sam'.format(tag)

    #run BWA
    commands_list.append( 'bwa mem -M {3} {0} {1} > {2}'.format(fqf1 , fqf2 , samfile , RefGen) )

    #####################################
    ### PICARD (sort & convert to BAM) ##
    #####################################

    #create BAM file
    bamfile = scratch_dir + '/{0}.sorted.bam'.format(tag)

    commands_list.append( 'java -Xmx16G -jar /n/data1/hms/dbmi/farhat/bin/picard/picard/build/libs/picard.jar SortSam INPUT={0} OUTPUT={1} SORT_ORDER=coordinate'.format(samfile, bamfile) )

    ####################################
    ### PICARD (remove duplicates) ####
    ###################################

    #create BAM file with removed duplicates
    drbamfile = bamfile # .replace(".bam", ".duprem.bam")

    #remove duplicates from BAM file
    # commands_list.append( "java -Xmx32G -jar /n/data1/hms/dbmi/farhat/bin/picard/picard/build/libs/picard.jar MarkDuplicates I={0} O={1} REMOVE_DUPLICATES=true M={2} ASSUME_SORT_ORDER=coordinate".format(bamfile, drbamfile, drbamfile[:-4]+'.metrics') )

    ####################################
    ### SAMTOOLS (to index BAM file) ###
    ####################################
    
    commands_list.append( "samtools index {0}".format(drbamfile) )
    
    ######################################
    ### QUALIMAP (quality of BAM file) ###
    ######################################
    
    #store quality report, pilon VCF & lineage call information all in Output directory
    commands_list.append( 'cd ' + output_dir )
    commands_list.append( 'mkdir QualiMap' ) #make a folder for pilon output in output directory
    commands_list.append( 'unset DISPLAY' ) #unset JAVA virtual machine variable [http://qualimap.bioinfo.cipf.es/doc_html/faq.html]
    commands_list.append( "/n/data1/hms/dbmi/farhat/bin/qualimap_v2.2.1/qualimap bamqc -bam {0} --outdir {1} --outfile {2}.pdf --outformat PDF".format(drbamfile, output_dir+'/QualiMap', tag+'_stats') )

    ###################################
    ### PILON (call variants) #########
    ###################################
    
    #store quality report, pilon VCF & lineage call information all in Output directory
    commands_list.append( 'mkdir pilon' ) #make a folder for pilon output in output directory
    out_pilon_dir = output_dir + '/pilon/' #variable for pilon output path

    commands_list.append( 'java -Xmx32G -jar /n/data1/hms/dbmi/farhat/bin/pilon/pilon-1.22.jar --genome {0} --bam {1} --output {2} --outdir {3} --variant'.format(RefGen, drbamfile, tag, out_pilon_dir) )

    #####################################
    ### Luca's LINEAGE CALLING script ###
    #####################################

    #create directory 
    commands_list.append( 'mkdir ' + scratch_dir + '/fast-lineage-caller/' )#make a folder for lineage call in output directory
    commands_list.append( 'mkdir ' + output_dir + '/fast-lineage-caller/' )#make a folder for lineage call in scratch directory

    #create VRT file
    vrtfile = scratch_dir + '/fast-lineage-caller/{}.vrt'.format(tag)

    commands_list.append( 'cd ' + scratch_dir + '/fast-lineage-caller' )#change directory to store output in scratch

    #convert VCF to VRT
    commands_list.append( '/n/data1/hms/dbmi/farhat/lfreschi/repos/vrt-tools/bin/vrtTools-vcf2vrt.py {0} {1} 1'.format(out_pilon_dir+tag+'.vcf', vrtfile) )

    #call lineage with SNP database an VRT file
    commands_list.append( 'cd ' + output_dir + '/fast-lineage-caller' )#change directory to store output in VCF output

    commands_list.append( '/n/data1/hms/dbmi/farhat/lfreschi/repos/fast-lineage-caller/bin/FastLineageCaller-assign2lineage.py --lin_snps /home/rv76/Bio_Pipelines/fast-lineage-caller-master/example/db_snps.tsv ' + vrtfile + ' &> ' + 'lineage_call.txt' )

    countReads(commands_list, drbamfile)
    ###############################################################################################################
    ######################################## SUBMIT as a job to O2 ################################################
    ###############################################################################################################
    
    #append all commands in a single string to be submitted as a job
    JankyPipe_job = ''
    for command_i in commands_list:
        JankyPipe_job = JankyPipe_job + '\n' + command_i
        
        #print(command_i)
        #print(' ')
    
    
    #directory where you want output + error files
    os.chdir(O2_SLURM_logs_dir)

    job_name = tag

    # s = Slurm(job_name , {'partition':'short' , 'n':'1' , 't':'0-10:00:00' , 'mem-per-cpu':'48G' , 'mail-type':'FAIL' , 'mail-user':'Jerry_Yang@hms.harvard.edu'})
    submitSlurmScript(commands_list, O2_SLURM_logs_dir)
    #submits the job
    # job_id = s.run(JankyPipe_job)

    # print(job_name  + ' : ' +  str(job_id))

# Launch_JankyPipe("/n/scratch2/jy250/SRR6176429_1.fastq" , fqf2 , "Jerry" , "./Output" , "/n/scratch2/jy250/" , "./Logs")