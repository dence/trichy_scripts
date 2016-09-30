#Daniel Ence
#09/30/2016

import re
import os
import argparse

def main(percent_cov_report):
			


if __name__ == "__main__":
	parser = argparser(description="Go through a percent coverage report and output a table of percent of genes fitting cutoffs")
	parser.add_argument("--percent_cov_report",type=str,help="The percent_coverage_report to analyze")
	parser.add_argument("--percent_cov_report2",type=str,help="Another percent coverage report to compare with the first one")
	parser.add_argument("--samples",type=str,help="comma-delimited list of samples to include in the table. Should be in booth reports if two reports supplied")
	args=parser.parse_args()
	main(args.percent_cov_report) 


