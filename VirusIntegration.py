import os
import subprocess
import sys
import re
from collections import defaultdict

BWA = "bwa"
SAMTOOLS = "samtools"
CDHIT = "cd-hit"
BLASTN = "blastn"
TEMPDIR = "."

def checkDependencies():
    cmd = [BWA]
    try:
        p = subprocess.Popen(cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    except:
        sys.exit("Error: bwa not installed!")

    cmd = [SAMTOOLS, "--version"]
    try:
        p = subprocess.Popen(cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    except:
        sys.exit("Error: samtools not installed!")

    cmd = [CDHIT]
    try:
        p = subprocess.Popen(cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    except:
        sys.exit("Error: CD-HIT not installed!")

    cmd = [BLASTN]
    try:
        p = subprocess.Popen(cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    except:
        sys.exit("Error: Blastn not installed!")

def preprocess():
    os.chdir(TEMPDIR)

def readSoftClips(file):
    os.system("rm *.fa")
    # read softclip file
    f= open(file,'r')
    #leftDict stores the sequences that matched sequence is on right side and softclip is on left side
    leftDict_softclip=defaultdict(list)
    leftDict_virusSeq=defaultdict(list)
    #rightDict stores the sequences that matched sequence is on left side and softclip is on right side
    rightDict_softclip=defaultdict(list)
    rightDict_virusSeq=defaultdict(list)
    for line in open(file):
        line = f.readline()
        subStr= line.split('\t')
        match1 = re.search(r'^(\d+)M(\d+)S$', subStr[3])
        match2 = re.search(r'^(\d+)S(\d+)M$', subStr[3])
        if match1:
            pos1 = int(match1.group(1))
            pos2 = int(match1.group(2))
            if pos1>=15 and pos2>=15:
                softclip = subStr[4][pos1:]
                virusSeq = subStr[4][0:pos1]
                breakpoint = int(subStr[1])+pos1
                leftDict_softclip[breakpoint].append(softclip)
                if breakpoint in leftDict_virusSeq:
                    if len(virusSeq) > leftDict_virusSeq[breakpoint]:
                        leftDict_virusSeq[breakpoint] = virusSeq
                else:
                    leftDict_virusSeq[breakpoint].append(virusSeq)
        if match2:
            pos1 = int(match2.group(1))
            pos2 = int(match2.group(2))
            if pos1>=15 and pos2>=15:
                softclip = subStr[4][0:pos1]
                virusSeq = subStr[4][pos1:]
                breakpoint = int(subStr[1])
                rightDict_softclip[breakpoint].append(softclip)
                if breakpoint in rightDict_virusSeq:
                    if len(virusSeq) > rightDict_virusSeq[breakpoint]:
                        rightDict_virusSeq[breakpoint] = virusSeq
                else:
                    rightDict_virusSeq[breakpoint].append(virusSeq)
    f.close()

    # write breakpoints into fa files
    for keyname, seqList in rightDict_softclip.items():
        f2=open("right_breakpoint"+str(keyname)+".fa", 'w+')
        count2=1
        for seq in seqList:
            f2.writelines("> seq"+str(count2)+'\n')
            f2.writelines(seq+'\n')
            count2=count2+1
        f2.close()

    for keyname, seqList in leftDict_softclip.items():
        f2=open("left_breakpoint"+str(keyname)+".fa", 'w+')
        count2=1
        for seq in seqList:
            f2.writelines("> seq"+str(count2)+'\n')
            f2.writelines(seq+'\n')
            count2=count2+1
        f2.close()

def main():
    checkDependencies()
    preprocess()
    readSoftClips(file='HPV16.softclip.txt')

if __name__ == "__main__":
	main()