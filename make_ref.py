import os, sys
#def main(chip, chip_ID, sample_pos, sampleID, out_dir):
def main(chip_vcf, chip_pos, sample_ID, out_dir):
	ref_dic= {}
	result_dic ={}
	sample = sample_ID
	chip_result = chip_vcf
	out_file = out_dir + sample + "_chip_result.txt"
	out = open(out_file,"w")
	with open(chip_pos,'rb') as f:
		for line in f:
			i = line[:-1].split('\t')
			chr = i[0]
			pos = i[1]
			ref = i[3]
			probe = ":".join(i[2].split(","))
			key =  chr +","+pos+","+ref
			#print key
			try: ref_dic[key] += ":" + probe
			except KeyError:
				ref_dic[key] = key+","+probe

	with open(chip_result,'rb') as f:
		for line in f:
			if line.startswith("#"):
				continue
			else:
				i = line.rstrip('\n')
				#i = line[:-1].split('\t')
				info = i.split('\t')
				ID = ";".join(info[2].split(','))
				alt = ";".join(info[4].split(','))
				chr = "chr" + info[0]
				pos = info[1]
				ref = info[3]
				data = ','.join(info[0:2])+","+ ref +","+ alt +","+','.join(info[8:])
				key = chr +","+pos+","+ref
				#print key
				result_dic[key] = "chr" + data
						
	KeyList = ref_dic.keys()

	for key in KeyList:
		try:  out.write(str(sample.rstrip('\n')) + ","+ result_dic[key] + "\n")
		except KeyError:
			out.write(str(sample.rstrip('\n')) + "," + ','.join(ref_dic[key].split(",")[0:3]) + ",.,.,.\n")
	out.close
	cmd = "split -d -l 37000 " + out_file + " " + out_file +"."
	os.system(cmd)
	cmd2 = "cd " + out_dir + ";" + "ls " + sample + "_chip_result.txt.?? > result.list"
	os.system(cmd2) 

if __name__ == '__main__':
	if len(sys.argv) != 5:
		sys.exit('python python.py [chip_vcf] [chip_position] [sampleID] [output_dir]\n')
	main(sys.argv[1], sys.argv[2],sys.argv[3], sys.argv[4])
