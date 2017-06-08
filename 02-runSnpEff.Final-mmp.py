import glob, sys, os, string

organism = "mmp"
snpEffPath = "/users/sturkars/snpEffv43/snpEff/snpEff.jar"
javaPath = "/usr/bin/java"
bcftoolsPath = "/users/sturkars/bcftools-1.1/bcftools"
# snpEff databases
if organism == "mmp":
    snpEffDatabase = "Methanococcus_maripaludis_s2"
if organism == "dvh":
    snpEffDatabase = "Desulfovibrio_vulgaris_str_hildenborough"
targetFolder = "results-UA3-152-09"
vcfFolders = glob.glob("%s/%s/*/"%(targetFolder,organism))
#print vcfFolders
#sys.exit()

folderCount = 1
for vcfFolder in vcfFolders:
    sampleFolder = vcfFolder.split("/")[2] # Folder containing reads
    print '\033[33mProcessing Folder: %s of %s (%s)\033[0m' %(folderCount, len(vcfFolders), sampleFolder )


    gatkIndel_vcf = glob.glob(vcfFolder+"*_gatk-variants-filtered-indels.vcf")
    gatkSnp_vcf = glob.glob(vcfFolder+"*_gatk-variants-filtered-snps.vcf")
    varscanIndel_vcf = glob.glob(vcfFolder+"*varscan_indel.vcf")
    varscanSnp_vcf = glob.glob(vcfFolder+"*varscan_snp.vcf")
    samtools_zipped = glob.glob(vcfFolder+"*filtered.vcf.gz")

    dvhVcfFiles = []
    dvhVcfFiles.extend(gatkIndel_vcf)
    dvhVcfFiles.extend(gatkSnp_vcf)
    dvhVcfFiles.extend(varscanIndel_vcf)
    dvhVcfFiles.extend(varscanSnp_vcf)
    dvhVcfFiles.extend(samtools_zipped)
    print dvhVcfFiles

    fileCount = 1
    for vcfFile in dvhVcfFiles:
        sampleName = vcfFile.split("/")[3].split("_")[0] + "_" + vcfFile.split("/")[3].split("_")[1]
        fileName = vcfFile.split("/")[3]
        print
        print '\033[32m Processing File: %s of %s (%s)\033[0m' %(fileCount, len(dvhVcfFiles), fileName )
        #print vcfFile
        #print sampleName
        #print fileName


        #print "....." + vcfFolder
        outVcfFile = vcfFile.split('.vcf')[0]+"-snpeff.vcf"
        outStatsFile = vcfFile.split('.vcf')[0]+"-snpeff.stats.txt"
        cmd1 ='%s -Xmx2g -jar %s -classic -csvStats %s -geneId -lof -v -formatEff -o gatk -ud 0 %s %s > %s'  %(javaPath, snpEffPath, outStatsFile, snpEffDatabase, vcfFile, outVcfFile)
        print cmd1
        os.system(cmd1)

        path2script=snpEffPath.split('/snpEff.jar')[0]
        outFilteredVcf = outVcfFile.split('-snpeff.vcf')[0]+"-filteredvariants.vcf"

        if vcfFile in (vcfFolder+sampleName+"_varscan_snp.vcf", vcfFolder+sampleName+"_varscan_indel.vcf"):
            print "I found a varscan file..."
            print vcfFile

            if os.stat(outVcfFile).st_size != 0:
                print ' File is not EMPTY'
                cmd2 = 'cat %s | %s -Xmx128m -jar %s/SnpSift.jar filter "(FILTER = \'PASS\')" > %s' %(outVcfFile, javaPath, path2script, outFilteredVcf)
                #cmd2 = 'cp %s %s' %(outVcfFile, outFilteredVcf)
                print cmd2
                os.system(cmd2)

            else:
                print "File %s is empty" %outVcfFile
                cmd2 = 'mv %s %s.empty' %(outFilteredVcf, outFilteredVcf)
            #print vcfFolder+sampleName+"_varscan_snp.vcf"

                #cmd2 = 'cat %s | %s -Xmx128m -jar %s/SnpSift.jar filter "(FILTER = \'PASS\')" > %s' %(outVcfFile, javaPath, path2script, outFilteredVcf)

                print cmd2
                os.system(cmd2)
        else:
                cmd2 = 'cat %s | %s -Xmx128m -jar %s/SnpSift.jar filter "(FILTER = \'PASS\')" > %s' %(outVcfFile, javaPath, path2script, outFilteredVcf)
                #cmd2 = 'cp %s %s' %(outVcfFile, outFilteredVcf)
                print cmd2
                os.system(cmd2)

    	finalVariants = outFilteredVcf.split('-filteredvariants.vcf')[0]+"-finalvariants.txt"
    	cmd3 ='cat %s | perl %s/scripts/vcfEffOnePerLine.pl | %s -Xmx128m -jar %s/SnpSift.jar extractFields - CHROM POS REF ALT AF AC DP MQ "(FILTER = \'PASS\')" "EFF[*].EFFECT" "EFF[*].IMPACT" "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" "EFF[*].GENE" "EFF[*].CODING" "EFF[*].RANK" "EFF[*].DISTANCE"> %s'%(outFilteredVcf, path2script, javaPath, path2script, finalVariants)
        vcfInput = ' '.join(dvhVcfFiles)
        print "Printing my line"
        print cmd3
        os.system(cmd3)

        fileCount = fileCount + 1

    folderCount = folderCount + 1
    #sys.exit()
