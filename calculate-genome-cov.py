

import glob, sys, os
# get list of bam files
bamFiles = glob.glob('%s/%s/*/*_marked.bam' %(resultsDir,organism))

count = 0
for bamFile in bamFiles:
    print bamFile
    count = count + 1
print "I found total of %s bamfiles" %s

resultsDir = "/Volumes/omics4tb/sturkarslan/singleCell-UA3-mmp/results-UA3-152-09"
organism = "mmp"
