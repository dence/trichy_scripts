#Daniel Ence
#07/09/2015


import sys
import re
from Bio import SeqIO

#Takes in a fasta file of intron sequences
#look at first and last 20 bps in all the 
#sequences

#makes a positions weight matrix for the first and last 20bps in the introns
def make_nucl_matrix(matrix_length):
	A = []
	G = []
	C = []
	T = []
	
	for i in range(0,matrix_length):
		A.append(0)
		G.append(0)
		C.append(0)
		T.append(0)
	nucl_matrix = {}
	nucl_matrix['A'] = A
	nucl_matrix['G'] = G
	nucl_matrix['C'] = C
	nucl_matrix['T'] = T

	return nucl_matrix	
			
def update_nucl_matrix(nucl_matrix, string):
	
	for i in range(0,len(string)):
		char = string[i]
		nucl_matrix[string[i]][i] = nucl_matrix[string[i]][i] + 1
	return nucl_matrix

def print_nucl_matrix(nucl_matrix):
	A_line = ""
	G_line = ""
	C_line = ""
	T_line = ""
	#A_line = "A:\t"
	#G_line = "G:\t"
	#C_line = "C:\t"
	#T_line = "T:\t"
	for i in range(0,matrix_length):
		A_percent = float(first_bps_matrix['A'][i]) / float(total_seqs)
		G_percent = float(first_bps_matrix['G'][i]) / float(total_seqs)
		T_percent = float(first_bps_matrix['T'][i]) / float(total_seqs)
		C_percent = float(first_bps_matrix['C'][i]) / float(total_seqs)

		A_line = A_line + '{:.2}'.format(A_percent) + "\t"
		T_line = T_line + '{:.2}'.format(T_percent) + "\t"
		G_line = G_line + '{:.2}'.format(G_percent) + "\t"
		C_line = C_line + '{:.2}'.format(C_percent) + "\t"

	print A_line
	print C_line
	print G_line
	print T_line

	return A_line + "\n" + C_line + "\n" + G_line + "\n" + T_line + "\n"

INTRONS = SeqIO.parse(sys.argv[1],"fasta")
matrix_length=20
first_bps_matrix=make_nucl_matrix(matrix_length)
last_bps_matrix=make_nucl_matrix(matrix_length)

total_seqs = 0
for record in INTRONS:
	#print record.seq	
	first_bps = record.seq[0:20]
	first_bps_matrix = update_nucl_matrix(first_bps_matrix,first_bps)
	last_bps = record.seq[(len(record.seq) - matrix_length):len(record.seq)]
	last_bps_matrix = update_nucl_matrix(last_bps_matrix,last_bps)
	total_seqs = total_seqs + 1
	#print first_10bps
	#print second_10bps

firstbps_file = open("first_bps.pwm.txt",'w')
lastbps_file = open("last_bps.pwm.txt",'w')
firstbps_file.write(print_nucl_matrix(first_bps_matrix))
lastbps_file.write(print_nucl_matrix(last_bps_matrix))

firstbps_file.close()
lastbps_file.close()



