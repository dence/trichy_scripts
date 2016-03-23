#Daniel Ence
#November 17, 2015


class TaxonomerRunner(object):
	
	def __init__(self,kmerLength=31,FastqFile,outputDirectory,cpus=5,dbPrefx,protein=0):
		#self.TaxonomerPath = "path to taxonomer"
		self.kmerLength = kmerLength
		self.FastqFile = FastqFile
		self.outputDirectory = outputDirectory
	def makeConfigFiles(self):
		taxonomer_config = 
		self.outputDirectory + "/" + "taxonomer." + self.kmerLength + ".config"
		config_handle = open(taxnomer_config,'w')
		#print custom lines
		config_handle.write("db_prefix:dbXX")
		config_handle.write("sti_prefix:stiXX")
		config_handle.write("tri_prefix:triXX")
		config_handle.write("kanalyze_input:kcXX")
		config_handle.write("input_fasta:inXX")
		config_handle.write("read_file:" + self.FastqFile)
		config_handle.write("output_file:outXX")
	
		config_handle.write("afterburner:0")	
		config_handle.write("kmer_classification_cutoff:1")
		config_handle.write("output_read_seq:0")
			
	def runTaxonomer(self):
		#Assume that taxonomer is in the path
		print "Running Taxonomer!"




