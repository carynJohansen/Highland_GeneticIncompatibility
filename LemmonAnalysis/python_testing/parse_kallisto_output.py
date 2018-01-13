# Python script to parse the Kallisto output to assess how many reads were aligned.

###############################
#           Env               #

import sys
import re

###############################
#          Methods            #

# open file
def open_file( filename ):
	with open(filename, 'r') as f:
		log_data = f.read()
	log_split = log_data.split("\n")
	return log_split

# function to open file containing file list, and to make that a list
def make_file_list( masterfile ):
	fileList = open_file(masterfile)
	return fileList

# function to parse the file and extract the numbers
def parse_and_extract( masterFile ):
	logList = open_file(masterFile)
	runInfo=[]
	for file in logList:
		file_content = open_file(file)
		quantLine = file_content[20]
		processedReads = quantLine.split(" ")[2]
		processedReads = processedReads.replace(",","")
		pseudoalignedReads = quantLine.split(" ")[4]
		pseudoalignedReads = pseudoalignedReads.replace(",","")
		srr = file_content[8]
		info=[srr, processedReads, pseudoalignedReads]
		runInfo.append(info)
	return runInfo

# math function to turn the number strings into numbers and to calculate the percentage of reads pseudoaligned

# function to return the information in a file
	# this will need the file name, or at least the SRR number
	# tab delim format, to be analyzed and visualized in R
def write_to_file( runInfo, newFile ):
	with open(newFile, 'w') as nf:
		for i in runInfo:
			for j in i:
				nf.write(j+"\t")
			nf.write("\n")


# main function for the eventual looping through multiple SRR kallisto runs
def main(masterFile):
	runInfo = parse_and_extract(masterFile)
	write_to_file(runInfo, "alignedReads.txt")

###############################
#            Main             #

if __name__ == '__main__':
	file = sys.argv[1]
	print file
	main( file )
