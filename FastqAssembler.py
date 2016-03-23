#Daniel Ence
#November 18, 2015

import os

class FastqAssembler(object):

	def __init__ (self,fastqList,outputDirectory):
		fastqList = fastqList
		outputDirectory = outputDirectory

	def assembleFastqs(self, assembler):
		if(assembler == "velvet"):
			assembleFastqsWithVelvet()
		else:
			print "unknown assembler specified!!!"
	
	def assembleFastqsWithVelvet(self):
		for fastq in self.fastqList:
			velvetDir = self.outputDirectory + "/" fastq + "_velvet"
			
			#Check whether we already made this directory
			if not os.path.exists(velvetDir):
				os.mkdir(velvetDir)
		
			#Check whether we already did this velvet run
			if not os.path.exists(velvetDir + "/contigs.fa"):
				#assume that velvet is in the path
				velvetHCommand = "velveth " + velvetDir + " 29 -fastq " + fastq
				subprocess.call(velvetH_command,shell=True)

				velvetGCommand = "velvetg " + velvetDir + " -cov_cutoff auto -exp_cov auto "
				subprocess.call(velvetGCommand,shell=True) 
		
