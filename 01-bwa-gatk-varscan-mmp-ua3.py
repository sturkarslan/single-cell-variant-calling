import glob, sys, os, string

# Process fastQ Files
def runPipeline():
    fastqFoldersAll = glob.glob('%s/*/' %(dataDir))

    # remove trimmed directory
    fastqFolders=[element for element in fastqFoldersAll if element not in ('DvH_single_cell_amplified_genome_strain_UA3_152_09-30575559/trimmed/', 'DvH_single_cell_amplified_genome_strain_UA3_152_09-30575559/DvH_09_1*',
    'DvH_single_cell_amplified_genome_strain_UA3_152_09-30575559/DvH_09_2*'
    'DvH_single_cell_amplified_genome_strain_UA3_152_09-30575559/DvH_09_3*'
    'DvH_single_cell_amplified_genome_strain_UA3_152_09-30575559/DvH_09_4*'
    'DvH_single_cell_amplified_genome_strain_UA3_152_09-30575559/DvH_09_5*')]
    print
    print 'FASTQ Folders: %s' %len(fastqFolders)
    #fastqFolders.remove('DvH_single_cell_amplified_genome_strain_UA3_152_03-30586562/trimmed/')



    folderCount = 1
    for folder in fastqFolders:
        sampleFolder = folder.split("/")[1] # Folder containing reads
        print '\033[33mProcessing Folder: %s of %s (%s)\033[0m' %(folderCount, len(fastqFolders), sampleFolder )
        fastqFilesFirst = glob.glob('%s/*R1*.fastq.gz' %(folder)) # 1st file the folder

        for fastqfile in fastqFilesFirst:
            i = fastqfile.split("/")[2]
            print '%s' %(i)

        fileCount = 1
        for file in fastqFilesFirst:
            i = file.split("/")[2]
            print
            print '\033[32m Processing File: %s of %s (%s)\033[0m' %(fileCount, len(fastqFilesFirst), i )

            # Collect Sample attributes
            fileName = file.split("_R1")[0]
            sampleResultsDir = resultsDir+ '/'+organism+'/'+fileName.split("/")[1]
            sampleTitle = fileName.split("/")[2]
            lane = fileName.split("/")[2].split("_")[2]
            sampleName = fileName.split("/")[2]
            #sampleId = fileName.split("/")[2].split("_")[0].split("-")[3] # uncomment for actual run
            sampleId = fileName.split("/")[2].split("_")[0] # uncomment for test run
            # create readgroup info
            RGId = fileName.split("/")[2].split("_")[1] + "_" + lane
            RGSm = sampleId
            RGLb = sampleId
            RGPu = lane
            print "RG: ID: %s SM: %s LB: %s PU: %s" %(RGId, RGSm, RGLb, RGPu)

            # create read pairs
            firstPair = fileName + "_R1_001.fastq.gz"
            secondPair = fileName + "_R2_001.fastq.gz"
            print
            print "First Pair : %s" %(firstPair.split("/")[2])
            print "Second Pair: %s" %(secondPair.split("/")[2])
            print

            # Run pipeline commands
            # Run Fastqc
            #runQC(firstPair, secondPair)

            # Run trimmomatic
            firstPairTrimmedPaired, secondPairTrimmedPaired = runTrim(firstPair, secondPair)

            # Run bwa
            runBWA(firstPairTrimmedPaired, secondPairTrimmedPaired,lane,sampleName,sampleId,sampleResultsDir,sampleTitle, RGId, RGSm, RGLb, RGPu)

            # Run samtools fixmate
            runSamtoolsFixmate(sampleResultsDir,sampleTitle)

            # Run samtools sort
            runSamtoolsSort(sampleResultsDir,sampleTitle)

            # Run GATK Indel realignment
            runGATK(sampleResultsDir,sampleTitle)

            fileCount = fileCount +1



        # Run markduplicates
        libraryName = runMarkDuplicates(sampleResultsDir,sampleTitle,sampleFolder)


        # Samtools flagstat
        flagstats(sampleResultsDir, libraryName)

        # Run Samtools Variant calling
        variantCalling(sampleResultsDir,sampleTitle,libraryName)

        # Run variant calling with Varscan
        varscan(sampleResultsDir,sampleTitle,libraryName)


        # Variant calling with GATK HaploTypeCaller
        haplotypeCaller(sampleResultsDir,libraryName)



        folderCount = folderCount + 1
        #sys.exit()

    return firstPair, secondPair, libraryName, sampleResultsDir


# Quality control
def runQC(firstPair, secondPair):
    print
    print "\033[34m Running FastQC Quality Control \033[0m"

    # create results folder
    # if not os.path.exists('%s' %(fastqcDir)):
#         #print '\033[31m %s directory does NOT exists.Creating one! \033[0m' %(fastqcDir)
#         os.makedirs('%s' %(fastqcDir))
#     else:
        #print '\033[31m %s directory exists. Not creating. \033[0m' %(fastqcDir)
    # run Command and write output into both screen and logfile with 2>&1 | tee -a %s
    cmd = '%s -t 4 -o %s %s %s' %(fastqc, fastqcDir, firstPair, secondPair)
    print
    print 'FastQC Command:', cmd
    #os.system(cmd)





# Trim reads
def runTrim(firstPair, secondPair):
    print "\033[34m Trimmed Files... \033[0m"
    # Program Parameters
    illuminaClip = "ILLUMINACLIP:/users/sturkars/Trimmomatic-0.35/adapters/TruSeq3-PE-2.fa:2:30:10" #Remove adapters
    leading = "LEADING:3" #Remove leading low quality or N bases
    trailing = "TRAILING:3" #Remove trailing low quality or N bases
    slidingWindow = "SLIDINGWINDOW:4:20" #Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 20
    minLen = "MINLEN:36"
    paired = "PE" # or SE for single end
    threads = "-phred33 -threads 8"

    # define result files
    filesFolder = firstPair.split('/')[0]
    firstPairTrimmedPaired = filesFolder+"/trimmed/"+firstPair.split('/')[2].split(".fastq.gz")[0] + "_paired_trimmed.fastq.gz"
    secondPairTrimmedPaired = filesFolder+"/trimmed/"+secondPair.split('/')[2].split(".fastq.gz")[0] + "_paired_trimmed.fastq.gz"
    firstPairTrimmedUnpaired = filesFolder+"/trimmed/"+firstPair.split('/')[2].split(".fastq.gz")[0] + "_unpaired_trimmed.fastq.gz"
    secondPairTrimmedUnpaired = filesFolder+"/trimmed/"+secondPair.split('/')[2].split(".fastq.gz")[0] + "_unpaired_trimmed.fastq.gz"
    trimDir = filesFolder+"/trimmed/"

    # # create trim folder
    # if not os.path.exists('%s' %(trimDir)):
    #     #print '\033[31m %s directory does NOT exists.Creating one! \033[0m' %(trimDir)
    #     os.makedirs('%s' %(trimDir))
    # else:
    #     #print '\033[31m %s directory exists. Not creating. \033[0m' %(trimDir)

    # define command
    cmd = '/users/sturkars/java/bin/java -Xmx128m -jar %s %s %s %s %s %s %s %s %s %s %s %s %s %s 2>&1 | tee -a %s' %(trimmomaticPath, paired, threads, firstPair, secondPair, firstPairTrimmedPaired, firstPairTrimmedUnpaired, secondPairTrimmedPaired, secondPairTrimmedUnpaired, illuminaClip, leading, trailing, slidingWindow, minLen, pipelineLog)
    #print "Trimmomatic Command: ", cmd
    #os.system(cmd)
    print firstPairTrimmedPaired
    print secondPairTrimmedPaired
    print
    return firstPairTrimmedPaired, secondPairTrimmedPaired





def runBWA(firstPairTrimmedPaired, secondPairTrimmedPaired,lane,sampleName,sampleId,sampleResultsDir,sampleTitle, RGId, RGSm, RGLb, RGPu):
    print "\033[34m Running BWA alignment... \033[0m"

    # create results folder
    if not os.path.exists('%s' %(sampleResultsDir)):
        print '%s directory does NOT exists.Creating one!...' %(sampleResultsDir)
        os.makedirs('%s' %(sampleResultsDir))
    else:
        print '%s directory exists. Not creating...' %(sampleResultsDir)

    # modify read group information
    readGroup = "@RG\\tID:%s\\tPL:ILLUMINA\\tSM:%s\\tLB:%s\\tPU:%s" %(RGId, RGSm, RGLb, RGPu)
    # bwa run command
    cmd = "%s mem -R '%s' %s %s %s > %s/%s.sam" %(bwaPath, readGroup, genomeFasta, firstPairTrimmedPaired, secondPairTrimmedPaired, sampleResultsDir, sampleTitle)
    print "Run BWA Command", cmd
    print
    #os.system(cmd)




def runSamtoolsFixmate(sampleResultsDir,sampleTitle):
    print
    print "\033[34m Running SAMtools fixmate... \033[0m"
    # fixmate run command
    cmd = '%s fixmate -O bam %s/%s.sam %s/%s_fixmate.bam' %(samtoolsPath, sampleResultsDir, sampleTitle, sampleResultsDir, sampleTitle)
    print "Samtools Fixmate Command: ", cmd
    #os.system(cmd)




def runSamtoolsSort(sampleResultsDir,sampleTitle):
    print
    print "\033[34m Running SAMtools sort.. \033[0m"
    # run command
    cmd = '%s sort -O bam -o %s/%s_sorted.bam -T tmp/%s_temp %s/%s_fixmate.bam' %(samtoolsPath, sampleResultsDir, sampleTitle, sampleTitle, sampleResultsDir, sampleTitle)
    print "Samtools Sort Command: ", cmd
    # index bam file
    cmd2 = '%s index %s/%s_sorted.bam' %(samtoolsPath, sampleResultsDir, sampleTitle)
    os.system(cmd)
    os.system(cmd2)


   #

def runGATK(sampleResultsDir,sampleTitle):
    print
    print "\033[34m Running GATK Realigner.. \033[0m"
    #run Target Interbal Creater command
    cmd1 = '%s -Xmx128m -jar %s -T RealignerTargetCreator -R %s -I %s/%s_sorted.bam -o %s/%s.intervals' %(javaPath, gatkPath, genomeFasta, sampleResultsDir, sampleTitle, sampleResultsDir, sampleTitle)
    #run Indel Realigner command
    cmd2 = '%s -Xmx4G -jar %s -T IndelRealigner -R %s -I %s/%s_sorted.bam -targetIntervals %s/%s.intervals -o %s/%s_realigned.bam' %(javaPath, gatkPath, genomeFasta, sampleResultsDir, sampleTitle, sampleResultsDir, sampleTitle, sampleResultsDir, sampleTitle)
    # index bam file
    cmd3 = '%s index %s/%s_realigned.bam' %(samtoolsPath, sampleResultsDir, sampleTitle)
    # Detect covariates
    cmd4 = '%s -Xmx4G -jar %s -T BaseRecalibrator -R %s -I %s/%s_realigned.bam -knownSites %s-variants-compiled_sorted.vcf -o %s/%s_recal.table' %(javaPath, gatkPath, genomeFasta, sampleResultsDir, sampleTitle, organism, sampleResultsDir, sampleTitle)
    # Adjust quality scores
    cmd5 = '%s -Xmx4G -jar %s -T PrintReads -R %s -I %s/%s_realigned.bam -BQSR %s/%s_recal.table -o %s/%s_recal.bam' %(javaPath, gatkPath, genomeFasta, sampleResultsDir, sampleTitle, sampleResultsDir, sampleTitle, sampleResultsDir, sampleTitle)

    print "Command GATK Interval Creater: ", #cmd1
    print
    os.system(cmd1)
    print
    print "Command GATK Realigner: ", #cmd2
    print
    os.system(cmd2)
    print
    print "Command GATK index BAM: ", #cmd3
    print
    os.system(cmd3)
    print
    print "Command GATK BaseRecalibrator: ", #cmd4
    print
    os.system(cmd4)
    print
    print "Command GATK PrintReads: ", #cmd5
    print
    os.system(cmd5)
    print






def runMarkDuplicates(sampleResultsDir,sampleTitle,sampleFolder):
    print
    print "\033[34m Running Mark Duplicates.. \033[0m"
    libraryName = sampleTitle.split("_L")[0]
    print "Library Name:" + libraryName
    print
    # collect list of recalibrated bam files
    recalBams = glob.glob('results-09/%s/%s/*_recal.bam' %(organism,sampleFolder))
    #print sampleFolder, recalBams
    metricsFile = sampleResultsDir+"/"+sampleTitle+'.metrics'
    # creaate command line parameter for each file
    print 'Input Recalibrated BAM Files...'
    bamList = []
    for i in recalBams:
        inputAdd = 'INPUT=%s' %(i)
        bamList.append(inputAdd)
        print i
    bamListJoined = " ".join(bamList)
    print
    print 'Output BAM File...'
    print '%s/%s_marked.bam' %(sampleResultsDir, libraryName)
    print bamListJoined
    # MArk Duplicates
    cmd = '%s -Xmx4G -jar %s MarkDuplicates %s VALIDATION_STRINGENCY=LENIENT METRICS_FILE=%s OUTPUT=%s/%s_marked.bam' %(javaPath, piccardPath, bamListJoined, metricsFile, sampleResultsDir, libraryName)
    print
    print "Mark Duplicated Command:... ", #cmd
    print
    os.system(cmd)

    # index bam file
    cmd2 = '%s index %s/%s_marked.bam 2>&1 | tee -a %s' %(samtoolsPath, sampleResultsDir, libraryName, pipelineLog)
    print
    print "Index bamfile Command:... ", #cmd2
    os.system(cmd2)
    return libraryName





def variantCalling(sampleResultsDir,sampleTitle,libraryName): # With samtools
    print
    print "\033[34m Running SAMtools Variant Calling.. \033[0m"
    # Produce BCF file with all locations in the genome
    cmd = '%s mpileup -go %s/%s.bcf -f %s %s/%s_marked.bam 2>&1 | tee -a %s' %(samtoolsPath, sampleResultsDir, libraryName, genomeFasta, sampleResultsDir, libraryName, pipelineLog)
    # Reduce list of sites
    cmd2 = '%s call -vmO z -o %s/%s.vcf.gz %s/%s.bcf 2>&1 | tee -a %s' %(bcftoolsPath, sampleResultsDir, libraryName, sampleResultsDir, libraryName, pipelineLog)
    # Prepare vcf file for querying
    cmd3 = '%s -p vcf %s/%s.vcf.gz 2>&1 | tee -a %s' %(tabixPath, sampleResultsDir, libraryName, pipelineLog)
    #Filtering
    percentageString = "%"
    cmd4 = "%s filter -O z -o %s/%s.filtered.vcf.gz -s LOWQUAL -i '%sQUAL>30' %s/%s.vcf.gz 2>&1 | tee -a %s" %(bcftoolsPath, sampleResultsDir, libraryName, percentageString, sampleResultsDir, libraryName, pipelineLog)
    print "Variant Calling mpileup: ", #cmd
    os.system(cmd)
    print
    print "Variant Calling Reduce: ", #cmd2
    os.system(cmd2)
    print
    print "Variant Calling tabix: ", #cmd3
    os.system(cmd3)
    print
    print "Variant Calling filtering: ", #cmd4
    os.system(cmd4)






def varscan(sampleResultsDir,sampleTitle,libraryName): # with varscan
    print
    print "\033[34m Running Varscan.. \033[0m"
    # varscan for snps
    cmd0 = '%s mpileup -B -f %s -o %s/%s.pileup %s/%s_marked.bam' %(samtoolsPath, genomeFasta, sampleResultsDir, libraryName, sampleResultsDir, libraryName)
    # varscan for snps
    cmd = '%s -Xmx128m -jar %s mpileup2snp %s/%s.pileup --output-vcf 1 --min-coverage 10 --min-reads2 3 --p-value 0.01 > %s/%s_varscan_snp.vcf' %(javaPath, varscanPath, sampleResultsDir, libraryName, sampleResultsDir, libraryName)
    # varscan for indels
    cmd2 = '%s -Xmx128m -jar %s mpileup2indel %s/%s.pileup --output-vcf 1 --min-coverage 10 --min-reads2 3 --p-value 0.01 > %s/%s_varscan_indel.vcf' %(javaPath, varscanPath, sampleResultsDir, libraryName, sampleResultsDir, libraryName)
    print "samtools Mpileup: ", cmd0
    print
    #os.system(cmd0)
    print "Varscan for SNPs: ", cmd
    print
    #os.system(cmd)
    print
    print "Varscan for INDELS: ", cmd2
    #os.system(cmd2)



def bamStats(sampleResultsDir, sampleTitle):
    print
    print "\033[34m Running BamStats.. \033[0m"
    libraryName = sampleTitle.split("_L")[0]
    # MArk Duplicates
    cmd = '%s -Xmx4g -jar %s -i %s/%s_marked.bam -f %s/%s_genomic.bed -o %s/%s_stats.html -v html 2>&1 | tee -a %s' %(javaPath,bamstatsPath,sampleResultsDir, libraryName, genomeDir, organism, sampleResultsDir, libraryName, pipelineLog)
    print "BamStats Command: ", cmd
    print
    #os.system(cmd)



def flagstats(sampleResultsDir, libraryName):
    print
    print "\033[34m Running Samtools Flagstat.. \033[0m"
    # Samtools flagstat
    cmd = '%s flagstat %s/%s_marked.bam > %s/%s.flagstats.txt' %(samtoolsPath, sampleResultsDir, libraryName, sampleResultsDir, libraryName)
    print "Flagstat Command: ", cmd
    print
    os.system(cmd)


def haplotypeCaller(sampleResultsDir,libraryName): #with GATK HaploTypeCaller
    print
    print "\033[34m Running GATK Haplotype Variant Caller.. \033[0m"
    # haplotype command
    cmd1 = '%s -Xmx4G -jar %s -T HaplotypeCaller -R %s -I %s/%s_marked.bam --genotyping_mode DISCOVERY -ploidy 1 -stand_emit_conf 10 -stand_call_conf 30 -o %s/%s_gatk-variants-raw.vcf' %(javaPath, gatkPath, genomeFasta, sampleResultsDir, libraryName, sampleResultsDir, libraryName)
    # Select snp variants
    cmd2 = '%s -Xmx128m -jar %s -T SelectVariants -R %s -V %s/%s_gatk-variants-raw.vcf -selectType SNP -o %s/%s_gatk-variants-raw-snps.vcf' %(javaPath, gatkPath, genomeFasta, sampleResultsDir, libraryName, sampleResultsDir, libraryName)
    # Apply filters to SNPs
    cmd3 = "%s -Xmx128m -jar %s -T VariantFiltration -R %s -V %s/%s_gatk-variants-raw-snps.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' --filterName 'my_snp_filter' -o %s/%s_gatk-variants-filtered-snps.vcf" %(javaPath, gatkPath, genomeFasta, sampleResultsDir, libraryName, sampleResultsDir, libraryName)
    # Select indel variants
    cmd4 = '%s -Xmx128m -jar %s -T SelectVariants -R %s -V %s/%s_gatk-variants-raw.vcf -selectType INDEL -o %s/%s_gatk-variants-raw-indels.vcf' %(javaPath, gatkPath, genomeFasta, sampleResultsDir, libraryName, sampleResultsDir, libraryName)
    # Apply filters to indels
    cmd5 = "%s -Xmx128m -jar %s -T VariantFiltration -R %s -V %s/%s_gatk-variants-raw-indels.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0' --filterName 'my_indel_filter' -o %s/%s_gatk-variants-filtered-indels.vcf" %(javaPath, gatkPath, genomeFasta, sampleResultsDir, libraryName, sampleResultsDir, libraryName)
    # Merge vcf files
    cmd6 = "%s -Xmx128m -jar %s -T CombineVariants -R %s --variant %s/%s_gatk-variants-filtered-snps.vcf --variant %s/%s_gatk-variants-filtered-indels.vcf -o %s/%s_gatk-variants-filtered-ALL.vcf  -genotypeMergeOptions UNIQUIFY" %(javaPath, gatkPath, genomeFasta, sampleResultsDir, libraryName, sampleResultsDir, libraryName, sampleResultsDir, libraryName)

    print "GATK HaplotypeCaller Comnand: ", cmd1
    #os.system(cmd1)
    print "Select SNP Variants: ", cmd2
    #os.system(cmd2)
    print
    print "Applying filters for SNPs: ", cmd3
    #os.system(cmd3)
    print "Select IndelVariants: ", cmd4
    #os.system(cmd4)
    print
    print "Applying filters for Indelss: ", cmd5
    #os.system(cmd5)
    print
    #print "Merging vcf files: ", cmd6
    #os.system(cmd6)



def genomeIndexes():
    print
    print "\033[34m Creating Genome indexes... \033[0m"
    # indexing genomes
    cmd1=bwaPath+' index '+genomeFasta
    cmd2=samtoolsPath+' faidx '+genomeFasta
    dictInput = genomeFasta.replace('.fasta', '')
    cmd3='%s -Xmx4G -jar %s CreateSequenceDictionary R=%s O=%s.dict 2>&1 | tee -a %s' %(javaPath,piccardPath,genomeFasta,dictInput, pipelineLog)
    print "BWA Genome index command: ", cmd1
    print
    print "Samtools Genome index Command: ", cmd2
    print
    print "GATK CreateDictionary Command: ", cmd3

    #os.system(cmd1)
    #os.system(cmd2)
    #os.system(cmd3)


################# Main #################
# Input files
organism = "mmp"
dataDir = "DvH_single_cell_amplified_genome_strain_UA3_152_09-30575559"
genomeDir = "reference"
genomeGff = '%s/%s.GCA_000195755.1.30.gtf' %(genomeDir, organism)
resultsDir = "results-09"
fastqcDir = '%s/fastqc' %(resultsDir)
pipelineLog = '%s/pipelineLog.txt' %(resultsDir)
knownSites = '%s-variants-compiled.vcf' %(organism)
# snpEff databases
if organism == "mmp":
    snpEffDatabase = "Methanococcus_maripaludis_S2_uid58035"
    genomeFasta = '%s/Methanococcus_maripaludis_s2.GCA_000011585.1.30.dna.genome.fasta'  %(genomeDir)
if organism == "dvh":
    snpEffDatabase = "dvh-genome"
    genomeFasta = '%s/Desulfovibrio_vulgaris_str_hildenborough.GCA_000195755.1.30.dna.genome.fasta' %(genomeDir)

### Create sequnce dictioonary
#/usr/bin/java -Xmx4G -jar /users/sturkars/picard-tools-1.139/picard.jar CreateSequenceDictionary REFERENCE=Desulfovibrio_vulgaris_str_hildenborough.GCA_000195755.1.30.dna.genome.fasta OUTPUT=Desulfovibrio_vulgaris_str_hildenborough.GCA_000195755.1.30.dna.genome.fasta.dict

# Programs
bwaPath = "/users/sturkars/bwa-0.7.12/bwa" # path to STAR executable
samtoolsPath = "/users/sturkars/samtools-1.2/bin/samtools"   # path to Trimmomatic executable
bcftoolsPath = "/users/sturkars/bcftools-1.1/bcftools"
fastqc = "/users/sturkars/FastQC/fastqc" # path to fastqc executable
javaPath = "/usr/bin/java"
piccardPath = "/users/sturkars/picard-tools-1.139/picard.jar"
trimmomaticPath = "/users/sturkars/Trimmomatic-0.35/trimmomatic-0.35.jar"
gatkPath = "/users/sturkars/gatk/GenomeAnalysisTK.jar"
tabixPath = "/users/sturkars/htslib-1.1/tabix"
varscanPath = "/users/sturkars/VarScan.v2.3.9.jar"
snpEffPath = "/users/sturkars/snpEff/snpEff.jar"

# Run functions
#genomeIndexes()
firstPair, secondPair, libraryName, sampleResultsDir = runPipeline()
