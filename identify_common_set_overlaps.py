#Daniel Ence
#February 1, 2016

import sys
import argparse
import re
#read in a blastx report (tabbed output)
#identify the reports that contain blast hits from the common set
#print out a list of those reports (the loci that the reports are based on?)


blast_reports = sys.argv[1:]

def update_no_gould(no_gould,line):

	if(re.search("\.RNA\.",line)):
		return no_gould
	
	if(re.search("Pentatrich.hominis",line)):
		no_gould["Pentatrich"] = 1
	if(re.search("Trich.tenax",line)):
		no_gould["tenax"] = 1
	if(re.search("Tritrich.foetus",line)):
		no_gould["foetus"] = 1
#	if(re.search("Tvag.30235",line)):
#		no_gould["30235"] = 1
#	if(re.search("Tvag.30238",line)):
#		no_gould["30238"] = 1
#	if(re.search("Tvag.50143",line)):
#		no_gould["50143"] = 1
#	if(re.search("Tvag.B7RC2",line)):
#		no_gould["B7RC2"] = 1
#	if(re.search("Tvag.G3",line)):
#		no_gould["G3"] = 1
#	if(re.search("Tvag.JRSTV41",line)):
#		no_gould["JRSTV41"] = 1
#	if(re.search("Tvag.T1",line)):
#		no_gould["T1"] = 1
#	if(re.search("Tvag.t016",line)):
#		no_gould["t016"] = 1

	if(re.search("Tvag",line)):
		no_gould["Tvag"] = 1

	return no_gould

def update_gould(with_gould, line):
	if(re.search("Pentatrich",line)):
		with_gould["Pentatrich"] = 1
	if(re.search("Tetratrich",line)):
		with_gould["Tetratrich"] = 1
	if(re.search("tenax",line)):
		with_gould["tenax"] = 1
	if(re.search("foetus",line)):
		with_gould["foetus"] = 1
	if(re.search("Trich.gall",line)):
		with_gould["Trich.gall"] = 1
	if(re.search("Trichomitus.batr",line)):
		with_gould["Trichomitus.batr"] = 1
#	if(re.search("Tvag.30235",line)):
#		with_gould["3025"] = 1
#	if(re.search("Tvag.30238",line)):
#		with_gould["30238"] = 1
#	if(re.search("Tvag.50143",line)):
#		with_gould["50143"] = 1
#	if(re.search("Tvag.B7RC2",line)):
#		with_gould["B7RC2"] = 1
#	if(re.search("Tvag.G3",line)):
#		with_gould["G3"] = 1
#	if(re.search("Tvag.JRSTV41",line)):
#		with_gould["JRSTV41"] = 1
#	if(re.search("Tvag.T1",line)):
#		with_gould["T1"] = 1
#	if(re.search("Tvag.t016",line)):
#		with_gould["t016"] = 1
	if(re.search("Tvag",line)):
		with_gould["Tvag"] = 1	
	return with_gould


#set up list of things to look for 
#report a list of common set with the Gould samples and a separate list without the Gould samples
with_gould_handle = open("with_gould_common_set.txt",'w')
without_gould_handle = open("without_gould_common_set.txt",'w')

for report in blast_reports:
	#hash of booleans to keep track of what we found for this blast_report
	#target_common_set_gould = {'Pentatrich':0,'Tetratrich':0,'tenax':0,'foetus':0,'Trich_gall':0,'Trichomitus_batr':0
	#			,'30235':0,'30238':0,'50143':0,'B7RC2':0,'G3':0,'JRSTV41':0,'T1':0,'t016':0};	
	#target_common_set_no_gould = {'Pentatrich':0,'Tetratrich':0,'tenax':0,'foetus':0
	#				,'Tvag.30235':0,'Tvag.30238':0,'Tvag.50143':0,'Tvag.B7RC2':0,'Tvag.G3':0,'Tvag.JRSTV41':0,'Tvag.T1':0,'Tvag.t016':0}


	
	target_common_set_gould = {'Pentatrich':0,'Tetratrich':0,'tenax':0,'foetus':0,'Trich.gall':0,'Trichomitus.batr':0,"Tvag":0}
	target_common_set_no_gould = {'Pentatrich':0,'tenax':0,'foetus':0,"Tvag":0}

	curr_file = open(report,'read')
	#hash of booleans to keep track of what we found for this blast_report
	for line in curr_file:
		target_common_set_no_gould = update_no_gould(target_common_set_no_gould,line)
		target_common_set_gould = update_gould(target_common_set_gould,line)
	curr_file.close()


	report_name_parts = report.split(".velvet_contigs")
	report_name = report_name_parts[0]

	tmp_list = ["Pentatrich","tenax","foetus","Tvag"]

	without_gould_line = "\n\n\n" + report_name + "\twithout_gould:\t"
	with_gould_line = report_name + "\twith_gould:\t"
	for tmp in tmp_list:
		wo_val = target_common_set_no_gould[tmp]
		w_val = target_common_set_gould[tmp]
		
		without_gould_line = without_gould_line + tmp + ":" + str(wo_val) + "\t"
		with_gould_line = with_gould_line + tmp + ":" + str(w_val) + "\t"

	#print without_gould_line + "\n"
	#print with_gould_line + "\n"

	with_gould_all_present = 1
	#with_gould_line = report_name + "\t"
	for key,value in target_common_set_gould.items():
		#print key + "\t" + str(value) + "\n"
		if(value == 0):
			with_gould_all_present = 0
	#		with_gould_line = with_gould_line + key + ":0\t"
	#	else:
	#		with_gould_line = with_gould_line + key + ":1\t"	

	without_gould_all_present = 1
	#without_gould_line = report_name + "\t"
	for key,value in target_common_set_no_gould.items():
	#	print key + "\t" + str(value) + "\n"
		if(value == 0):
	#		print "unsetting without_gould_all_present:\t" + key + "\t" + str(value) + "\n"
			without_gould_all_present = 0
	#		without_gould_line = without_gould_line + key + ":0\t"
	#	else:
	#		without_gould_line = without_gould_line + key + ":1\t"

	
	if(with_gould_all_present):
		with_gould_handle.write(report_name + "\n")
	#with_gould_handle.write(with_gould_line + "\n")
	
	if(without_gould_all_present):
		without_gould_handle.write(report_name + "\n")
	#without_gould_handle.write(without_gould_line + "\n")
with_gould_handle.close()
without_gould_handle.close()
