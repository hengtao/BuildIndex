#!/usr/bin/env python3

import re,time,os,sys
import traceback
import argparse
import json
import subprocess
import configparser
from functools import reduce
import logging
from logging.handlers import RotatingFileHandler
import glob
from datetime import datetime
import gzip

class HelpFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
	pass

def Traceback(traceback_recordfile):
	with open(traceback_recordfile, 'w+') as f:
		traceback.print_exc(file = f)
		f.flush()
		
def LogRecord(logger, loginfo):
	now_time = datetime.now()
	logger.info(loginfo)
	return now_time, logger
	
def Cmd_and_Time(cmd, logger):
	starttime, logger = LogRecord(logger, cmd)
	p = subprocess.Popen(cmd, shell = True)
	#p.wait()
	if(p.returncode != 0):
		logger.exception("COMMAND fail:----\n{cmd}".format(cmd = cmd))
		sys.exit(1)
	else:
		spend_time = datetime.now() - starttime
		days	   = spend_time.days
		totalseconds = spend_time.seconds
		hours	  = int(int(totalseconds)/(60*60))
		minites	= int((int(totalseconds)%3600)/60)
		seconds	= int(totalseconds)%60
		logger.info("Total time for this step is: {days} days, {hours} hours, {minites} minites, {seconds} seconds".format(days = days, hours = hours, minites = minites, seconds = seconds))
		return spend_time
		
def Check_software(logger, software_path):
	if os.path.exists(software_path):
		logger.debug("Choose software:" + software_path + "!")
	else:
		output = os.popen('which ' + software_path).read()
		if output:
			software_temp = output.split("\n")[0]
			if os.path.exists(software_temp):
				software_path = software_temp
				logger.debug("Choose software:" + software_path + "!")
		else:
			location = os.popen("locate " + software_path).read()
			print(location)
			if location:
				software_path = location.split("\n")[0]
				print(software_path)
				if not os.path.exists(software_path):
					logger.error("Can't locate the " + software_path + "!")
					logger.error("Please install this software: " + software_path + "!")
					exit(1)
	return software_path

def gtf2bed(gtf, bed):
	GTF_in_file = open(gtf, "r")
	BED_out_file = open(bed, "w")
	GTFs = csv.reader(GTF_in_file, dialect = "excel-tab")
	BEDs = csv.writer(BED_out_file, dialect = "excel-tab")
	transcripts = {}
	for GTF in GTFs:
		if not GTF[0].startswith("#"):
			if GTF[2] == "exon":
				transcript_id = re.findall('transcript_id "([^;]*)"', GTF[8])[0]
				if transcript_id in transcripts:
					transcripts[transcript_id][6].append([int(GTF[3]), int(GTF[4])])
				else:
					transcripts[transcript_id] = []
					transcripts[transcript_id].append(GTF[0])
					transcripts[transcript_id].append(GTF[3])
					transcripts[transcript_id].append(GTF[4])
					transcripts[transcript_id].append(GTF[5])
					transcripts[transcript_id].append(GTF[6])
					transcripts[transcript_id].append(transcript_id)
					transcripts[transcript_id].append([[int(GTF[3]), int(GTF[4])]])

	for transcript in transcripts:
		transcripts[transcript][6].sort()
		transcript_start = transcripts[transcript][6][0][0]
		transcript_end = transcripts[transcript][6][len(transcripts[transcript][6]) - 1][1]
		exon_sizes = ""
		exon_starts = ""
		for exon in transcripts[transcript][6]:
			exon_size = exon[1] - exon[0] + 1
			exon_start = exon[0] - transcript_start
			exon_sizes = exon_sizes + str(exon_size) + ","
			exon_starts = exon_starts + str(exon_start) + ","
		BED = []
		BED.append(transcripts[transcript][0])
		BED.append(transcript_start - 1)
		BED.append(transcript_end)
		BED.append(transcripts[transcript][5])
		BED.append(transcripts[transcript][3])
		BED.append(transcripts[transcript][4])
		BED.append(transcript_start - 1)
		BED.append(transcript_end)
		BED.append("0, 0, 0")
		BED.append(len(transcripts[transcript][6]))
		BED.append(exon_sizes)
		BED.append(exon_starts)
		BEDs.writerow(BED)
	GTF_in_file.close()
	BED_out_file.close()

class BuildIndex():
	def __init__(self, args, genomefile, gtf, samtools, picard, hisat2, hisplice, hiexons, bwasoft):
		self.genomefile = genomefile
		self.gtf  = gtf
		self.fai  = args.faindex
		self.dict = args.dictindex
		self.bwa  = args.bwaindex
		self.ht2  = args.hisat2index
		##TOOLS
		self.samtools = samtools
		self.picard   = picard
		self.hisat2   = hisat2
		self.hisplice = hisplice
		self.hiexons  = hiexons
		self.bwasoft  = bwasoft
	def buildindex(self):
		if self.fai == "yes":
			cmd = "{samtools} faidx {genome}".format(samtools = self.samtools, genome = self.genomefile)
			Cmd_and_Time(cmd, logger)
		if self.dict == "yes":
			dictfile = self.genomefile.replace("fa", "dict")
			cmd = "{samtools} dict -o {dictfile} {genome}".format(samtools = self.samtools, dictfile = dictfile, genome = self.genomefile)
			Cmd_and_Time(cmd,logger)
			dirname = os.path.dirname(self.genome)
			bedfile = self.gtf.replace("gtf", "bed")
			gtf2bed(self.gtf, bedfile)
			intervallist = self.gtf.replace("gtf", "interval_list")
			cmd = "java -Xmx5g -XX:ParallelGCThreads=2 -Djava.io.tmpdir={dirname} -jar {picard} I={bedfile} SORT=true SD={dict} O={intervallist}".format(dirname = dirname, picard = self.picard, bedfile = bedfile , dict = dictfile, intervallist = intervallist)
			Cmd_and_Time(cmd, logger)
		if self.ht2 == "yes":
			pattern = re.sub(r".fa|.fasta|.fna|.fa.gz|.fasta.gz|.fna.gz")
			splicefile = self.genomefile.replace("fa", "splicesite.txt")
			exonfile   = self.genomefile.replace("fa", "exons.txt")
			refbase = re.sub(pattern, "", os.path.basename(self.genome))
			cmd = "{hisplice} {gtf} > {splicefile} \n".format(hisplice = self.hisplice, gtf = self.gtf, splicefile = splicefile)
			cmd += "{hiexons} {gtf} > {exonfile} \n".format(hiexons = self.hiexons, gtf = self.gtf, exonfile = exonfile)
			cmd += "{hisat2_build} --exon {exons} --ss {splicesite} {ref} {refbase}".format(hisat2_build = self.hisat2, exons = exonfile, splicesite = splicefile, ref = self.genome, refbase = refbase)
			Cmd_and_Time(cmd, logger)
		if self.bwa == "yes":
			cmd = "{bwa} index {ref}".format(bwa = self.bwasoft, ref = self.genome)
			Cmd_and_Time(cmd, logger)
			
			
def main(args):
	genomefile  = args.genomefile
	gtf = args.gtf
	logfile	= args.logfile

	## Import logger, 获取logger实例
	logger	= logging.getLogger(__name__)
	## 指定logger输出格式
	formatter = logging.Formatter("%(asctime)s %(levelname)s: %(message)s")
	#logger.setLevel(level = logging.INFO)
	handler = logging.FileHandler(logfile)
	handler.formatter = formatter
	## 控制台日志
	console_handler = logging.StreamHandler(sys.stdout)
	console_handler.formatter = formatter
	## 为logger添加日志处理器
	logger.addHandler(handler)
	logger.addHandler(console_handler)
	## 指定日志的最低输出级别，默认为WARN
	logger.setLevel(level = logging.INFO)
	
	samtools = Check_software(logger, "samtools")
	picard   = Check_software(logger, "picard.jar")
	hisat2   = Check_software(logger, "hisat2-build")
	hisplice = Check_software(logger, "extract_splice_sites.py")
	hiexons  = Check_software(logger, "extract_exons.py")
	bwasoft  = Check_software(logger, "bwa")
	buildindex = BuildIndex(args, genomefile, gtf, samtools, picard, hisat2, hisplice, hiexons, bwasoft)
	buildindex.buildindex()

if __name__ == "__main__":
	scriptpath = os.path.split(os.path.realpath(__file__))[0]
	parse = argparse.ArgumentParser(formatter_class = HelpFormatter, description = '''
Usage:

python3 {scriptpath}/BuildIndex.py <args> <args> ... <genome.fa>

required arguments:
<genome.fa>     Reference genome file in fa/fasta/fna format.
<genome.gtf>    Annotation file for genome in gtf format.
<logfile>       Logfile to record running log.


'''.format(scriptpath = scriptpath))
	
	parse.add_argument('-log', '--blogfile', required = True, dest = "logfile", help = "logfile to record running logs", type = str, nargs = '?')
	parse.add_argument('-gtf', '--gtf', required = True, dest = "gtf", help = "annotation file in gtf format", type = str, nargs = '?')

	parse.add_argument('-ht2', '--hisat2index', required = False, dest = "hisat2index", help = "build hisat2 index or not", default = "yes", choices = ["yes", "no"], type = str, nargs = '?')
	parse.add_argument('-fai', '--faindex', required = False, dest = "faindex", help = "build samtools fai index or not", default = "yes", choices = ["yes", "no"], type = str, nargs = '?')
	parse.add_argument('-dict', '--dictindex', required = False, dest = "dictindex", help = "build dict index or not", default = "yes", choices = ["yes", "no"], type = str, nargs = '?')
	parse.add_argument('-bwa', '--bwaindex', required = False, dest = "bwaindex", help = "build bwa index or not", default = "yes", choices = ["yes", "no"], type = str, nargs = '?')
	
	parse.add_argument('genomefile', help = "reference genome file in fa/fasta/fna format", type = str, nargs = '?')

	args = parse.parse_args()
	
	main(args)
