
import glob, os, sys, gzip, tabix, re, random
from multiprocessing import Process, Manager, Pool
from functools import partial
import tempfile

result_dic = {}
def vcf_compare(ngs,sampleID,bms):
	result= ''
	chip_open= open(bms)
	out= open(bms+'.out','w')
	ex_chr=["chrX","chrY","chrM"]
	#header	
	result='Position\tREF\tChip\tNGS\tFILTER\tGenotype\tTrue/False\n'

	#chip	
	allele_x_f, allele_y_f = '', ''
	done= ''
	while 1:
		chip = {}
		line= chip_open.readline()
		if not line: break
		c_fields= line[:-1].split(',') #remove newline charactor
		sample_id=c_fields[0] #sample id 
		if sampleID != sample_id: # If sample id with script argument
			print "[ERROR] Sample ID is not the same."  
			print sampleID
			print sample_id
			sys.exit(1)
		if c_fields[1] in ex_chr:
			continue
		else:
			key= c_fields[1]+':'+c_fields[2] # chr#& position
			if '/' in c_fields[6]: # Chip array allele - G:G, A:C, No Call, Invalid ..
				CGT = c_fields[6].split(':')[0]
				GT1 = CGT.split('/')[0] # G:G, A:C
				if GT1 == "." :
					allele_x = "No Call"
				elif int(GT1) == 0: #reference genotype
					allele_x = c_fields[3]
				else: #variant
					allele_x = c_fields[4]
	
				GT2 = CGT.split('/')[1]
				if GT2 == ".":
					allele_y = "No Call"
				elif int(GT2) == 0:#reference genotype
					allele_y = c_fields[3]
				else: #variant
					allele_y = c_fields[4]
			else:
				allele_x = "No Call" # No Call, Invalid ..
				allele_y = "No Call"

			chip[key]='\t'+c_fields[3]+'\t'+allele_x+'/'+allele_y
			done= ''
		
		gvcf= tabix.open(ngs)
		gvcf_tmp=gvcf.query(c_fields[1],int(c_fields[2])-30000,int(c_fields[2])+30000) #chr, start, end

		#GVCF&COMPARE
		for fields in gvcf_tmp:
			GT= fields[9].split(':')[0] # Genotype
			check= 0
			if 'END' not in fields[7]:#variant
				key= fields[0]+':'+fields[1] #vcf position
				if key in chip:#WGS has variant at same postion
					vcf= fields[3]+"/"+fields[4]
					result += key+chip[key]+"\t"
					chip_s= chip[key].split('\t')[-1].split('/')
					vcf_s= vcf.split('/')
					chip_ref= chip[key].split('\t')[1].split('/')[0] # chip allele, exist '.' in gVCF 
					if chip_s[0] == "A" or chip_s[0] == "T" or chip_s[0] == "G" or chip_s[0] == "C":
						if GT == "1/1" or GT == "1": #alt homo
							for i in chip_s:
								if i == vcf_s[1]:
									check += 1
							if check == 2:
								result += fields[4]+"/"+fields[4]+"\t"+fields[6]+"\t"+GT+"\tT\n"
							if check != 2:
								result += fields[4]+"/"+fields[4]+"\t"+fields[6]+"\t"+GT+"\tF\n"
						elif GT == "0/1": #hetero
							if chip_s[0] != chip_s[1]:
								for i in chip_s:
									for j in vcf_s:
										if i == j:
											check += 1
								if check == 2:
									result += vcf+"\t"+fields[6]+"\t"+GT+"\tT\n"
								if check != 2:
									result += vcf+"\t"+fields[6]+"\t"+GT+"\tF\n"
							elif chip_s[0] == chip_s[1]:
								result += vcf+"\t"+fields[6]+"\t"+GT+"\tF\n"
						elif GT == "0/0" or GT == "0": #ref homo, does not exist 'END' string in INFO column
							for i in chip_s:
								if i == chip_ref:
									check += 1
							if check == 2:
								result += chip_ref+"/"+chip_ref+"\t"+fields[6]+"\t"+GT+"\tT\n"
							if check != 2:
								result += chip_ref+"/"+chip_ref+"\t"+fields[6]+"\t"+GT+"\tF\n"
						elif GT == ".": # does not sequencing
							result += ".\t"+fields[6]+"\t"+GT+"\tN\n"
						else:
							result += fields[3]+"/"+fields[4]+"\t"+fields[6]+"\t"+GT+"\tF\n"
					else: # chip results are 'UND,NOAMP,INV ...'
						if GT == "1/1" or GT == "1":
							result += fields[4]+"/"+fields[4]+"\t"+fields[6]+"\t"+GT+"\tN\n"
						elif GT == "0/1":
							result += vcf+"\t"+fields[6]+"\t"+GT+"\tN\n"
						elif GT == "0/0" or GT == "0":
							result += chip_ref+"/"+chip_ref+"\t"+fields[6]+"\t"+GT+"\tN\n"
						elif GT == ".":
							result += ".\t"+fields[6]+"\t"+GT+"\tN\n"
						else:
							result += chip_ref+"/"+chip_ref+"\t"+fields[6]+"\t"+GT+"\tF\n"
					
			elif 'END' in fields[7]: # ref homo or not sequencing
				INFO_items= fields[7].split(';')
				for item in INFO_items:
					if item.startswith('END'):
						value= item.split('=')[1]
						for id in range(int(fields[1]),int(value)+1,1): # from variant position to 'END' position
							key= fields[0]+':'+str(id)
							if key in chip:
								#result_dic[key] = 1
								chip_ref= chip[key].split('\t')[1].split('/')[0] # dbSNP allele
								result += key+chip[key]+"\t"
								chip_s= chip[key].split('\t')[-1].split('/')
								if chip_s[0] == "A" or chip_s[0] == "T" or chip_s[0] == "G" or chip_s[0] == "C":
									if GT == '.': # does not sequencing
										result += ".\t"+fields[6]+"\t"+GT+"\t"+"N\n"
									elif GT == '0/0' or GT == '0':
										for x in chip_s: #ref homo
											if x == chip_ref:
												check += 1
			
										if check == 2:
											result += chip_ref+"/"+chip_ref+"\t"+fields[6]+"\t"+GT+"\t"+"T\n"
										if check != 2:
											result += chip_ref+"/"+chip_ref+"\t"+fields[6]+"\t"+GT+"\t"+"F\n"
								else: # chip results are 'UND,NOAMP,INV ...'
									if GT == ".":
										result += ".\t"+fields[6]+"\t"+GT+"\t"+"N\n"
									elif GT == '0/0' or GT == '0':
										result += chip_ref+"/"+chip_ref+"\t"+fields[6]+"\t"+GT+"\t"+"N\n"
			check = 0
							
	out.write(result)
	out.close()

def calculation(sampleID):
	true_check = []
	true, false, notDefined, match_total = 0, 0, 0, 0
	result_count= open(sampleID+'.out','r')
	while 1:
		line= result_count.readline()
		if not line: break
		if line.startswith("chr"):
			try: result_dic[line[:-1].split('\t')[0]] +=1
			except KeyError:
				result_dic[line[:-1].split('\t')[0]] = 1
	result_count.seek(0)
	while 1:
		line= result_count.readline()
		if not line: break
		match= line[:-1].split('\t')[6]
		if line[:-1].split('\t')[0] not in result_dic:
			continue
		else:
			if result_dic[line[:-1].split('\t')[0]] == 1:
				if match == "T": true+= 1
				elif match == "F": false+= 1
				elif match == "N": notDefined+= 1

			else:
				true_check.append(match)
				number = len(true_check)
				if number == result_dic[line[:-1].split('\t')[0]]:
					if "T" in true_check:
						true+= 1
					elif "F" in true_check:
						false+= 1
					else: notDefined+= 1
					true_check = []

	match_total= int(true)+int(false)+int(notDefined)
	print "total:\t%d" %match_total
	print "true:\t%d, %s%%" %(true, str(round(true/float(match_total)*100, 2)))
	print "false:\t%d, %s%%" %(false, str(round(false/float(match_total)*100, 2)))
	print "not defined:\t%d, %s%%\n" %(notDefined, str(round(notDefined/float(match_total)*100,2)))

def main():
	tempfile.tempdir = os.getcwd()
	ngs = sys.argv[2]
	sampleID = sys.argv[3]
	bms_list = sys.argv[1]
	bms_list_l = []
	with open(bms_list) as f:
		for line in f:
			bms_file = line.rstrip()
			bms_list_l.append(bms_file)

	pool = Pool(20)
	func = partial(vcf_compare, ngs, sampleID)
	pool.map(func, bms_list_l)
	pool.close()
	pool.join()
	cwd = os.getcwd()
	merge_cmd = "cat " + cwd +"/*.??.out > " + cwd +"/" +sampleID +".out && touch " + cwd+ "/" + sampleID + ".ok"
	os.system(merge_cmd)
	while 1:
		check =  cwd + "/" +sampleID + ".ok"
		if os.path.isfile(check):
			break
	calculation(sampleID)

if __name__ == '__main__':
	if len(sys.argv) != 4:
		sys.exit('python python.py [chip_result_list] [gVCF from iSAAC-IVC] [sampleID]\n')
	main()
