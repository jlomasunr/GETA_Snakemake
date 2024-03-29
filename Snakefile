import os
import yaml
import logging
import Bio
from rich.logging import Console, RichHandler
from termcolor import colored

# Set up logging format
logging.basicConfig(level=logging.INFO,
		format="%(message)s",
		datefmt="[%X]",
		handlers=[RichHandler(console=Console(stderr=True))])

# Load config file
with open('config.yaml', 'r') as config_file:
	try:
		config_dict=yaml.safe_load(config_file)
	except yaml.YAMLError as e:
		logging.error("Could not load config file! Check config.yaml ... see error below")
		sys.exit(e)

# Check config file for valid augustus set-up
if config_dict["Input"]["augustus_species"] == "" and 
config_dict["Input"]["use_augustus"] == "":
	logging.error(f"No value provided for augustus. Check config.yaml...")
		sys.exit()
elif config_dict["Input"]["use_augustus"] != "":
	#TODO: check augustus install for species, error and quit if not present use config_dict["Input"]["AUGUSTUS_CONFIG_PATH"]
	#      Set a variable for use augustus
	#      Set variable to hold agustus species
elif config_dict["Input"]["augustus_species"] != "":
	#TODO: Set augustus species variable

# Check config file for other inputs
for input in config_dict["Input"]:
	if config_dict[input] == "" and 
	input != "augustus_species" and 
	input != "use_augustus":
		logging.error(f"No value provided for {input}. Check config.yaml...")
		sys.exit()
logging.info("config.yaml loaded successfully!")

# Store inputs
outDir = config_dict["out_prefix"]
if not os.path.exists(outDir):
	os.makedirs(outDir)

#TODO: Check software dependencies 
# RepeatModeler, RepeatMasker, HPCgridrunner, java, hisat2, samtools, hmmer
software_list = ["RepeatMasker","RepeatModeler","samtools","java","blastn","hisat2","hmmscan"]
check_software_existence(software_list)

#check_software_existence("RepeatMasker")
#check_software_existence("samtools")
#def check_software_existence(software_name:str)->None:

def check_software_existence(software_name:list)->None:
	for software in software_list:
    	if shutil.which(software_name) is not None:
    		print(f"{software_name} is installed")
		else:
			print(colored(f"{software_name} is not installed","red"))

# Setup output folders
os.makedirs(f"{outDir}/RepeatMasker/repeatMasker/", exist_ok=True)
os.makedirs(f"{outDir}/RepeatMasker/repeatModeler/", exist_ok=True)
os.makedirs(f"{outDir}/TrimGalore/", exist_ok=True)
os.makedirs(f"{outDir}/HiSat2/", exist_ok=True)
os.makedirs(f"{outDir}/transcript/", exist_ok=True)
os.makedirs(f"{outDir}/Augustus/", exist_ok=True)
if config_dict["Input"]["augustus_species"] != "":
	os.makedirs(f"{outDir}/Augustus/training/", exist_ok=True)



# Copy over genome and protein file
# TODO: Make more robust?
if not os.path.exists(f"{outDir}/genome.fasta"):
	os.system(f"cp {config_dict["genome"]} genome.fasta")
if not os.path.exists(f"{outDir}/homolog.fasta"):
	os.system(f"cp {config_dict["protein"]} homolog.fasta")

# Store RNAseq files
PE_reads = {}
SE_reads = {}

if PE1 and PE2:
    PE1 = config_dict["Input"]["rnaseq_first"].split(",")
    PE2 = config_dict["Input"]["rnaseq_second"].split(",")
    if len(PE1) != len(PE2):
        raise ValueError("the input file number of -1 was not equal to -2.")
# RNAseq_dir = # TODO: Get absoulte path to RNAseq data directory from first file name
	RNAseq_dir = abs_path(PE[0])

def abs_path(filename:str)->str:
	path=os.path.dirname(os.path.abspath(filename))
	path += "/"
	return(path)

#    for i in range(len(PE1)):
#        PE1[i] = abs_path(PE1[i])
#        PE2[i] = abs_path(PE2[i])
#        pe_file = PE1[i] + "\t" + PE2[i]
#        pe_reads[pe_file] = 1

#if single_end:
#    SE = single_end.split(",")
#   for i in range(len(SE)):
#        SE[i] = abs_path(SE[i])
#        SE_reads[SE[i]] = 1


def allInput():
	input = [f"{outDir}/RepeatMasker/repeatMasker/genome.fasta.out",		# RepeatMasker_species
			 f"{outDir}/RepeatMasker/repeatModeler/genome.fasta.out",		# RepeatMasker_custom
			 f"{outDir}/RepeatMasker/genome.masked.fasta",					# RepeatMasker_merge
			 f"{outdir}/HiSat2/hisat2.sam",
			 f"{outdir}/HiSat2/hisat2.sorted.sam",
			 f"{outdir}/HiSat2/genome.1.ht2",
			 f"{outdir}/HiSat2/splited_sam_files.list"
			]
	for fq_gz in RNAseq:
		input.append(f"TrimGalore/{re.sub("^.*\/", "", fq_gz)}")			# TrimGalore

	if config_dict["Input"]["RM_lib"] == "None":
		input.append(f"{outDir}/RepeatMasker/repeatModeler/species-families.fa")	# RepeatModeler
		input.append(f"{outDir}/RepeatMasker/repeatModeler/species-families.stk")	# RepeatModeler

	if config_dict["Input"]["rnaseq_single"] != "":
		input.append(f"{outdir}/TrimGalore/readSE.fq.gz")				# combineSingle
	if config_dict["Input"]["rnaseq_first"] != "":
		input.append(f"{outdir}/TrimGalore/read1.fq.gz")				# combinePaired
		input.append(f"{outdir}/TrimGalore/read2.fq.gz")				# combinePaired

	return(input)

rule all:
	input: allInput
		
rule RepeatMasker_species:
	output: f"{outDir}/RepeatMasker/repeatMasker/genome.fasta.out"
	threads: config_dict["Threads"]["RepeatMasker"]
	run:
		shell(f"RepeatMasker {config_dict["Internal"]["RepeatMasker"]} \
			-pa {threads} \
			-species {config_dict["Input"]["RM_species"]} \
			-dir {outDir}/RepeatMasker/repeatMasker/ \
			genome.fasta")

rule RepeatModeler:
	output: 
		f"{outDir}/RepeatMasker/repeatModeler/species-families.fa",
		f"{outDir}/RepeatMasker/repeatModeler/species-families.stk"
	threads: config_dict["Threads"]["RepeatModeler"]
	run:
		if not os.path.exists(f"{outDir}/RepeatMasker/repeatModeler/BuildDatabase.ok"):
			shell(f"BuildDatabase -name species -engine ncbi {outDir}/genome.fasta")
			shell("touch BuildDatabase.ok")
		
		shell(f"RepeatModeler \
			-pa {threads} \
			-database species \
			-LTRStruct")

# Variable input required for RepeatMasker_custom based on value of 'RM_lib'
def repeatMaskerCustom_input():
	if config_dict["Input"]["RM_lib"] == "None":
		input = f"{outDir}/RepeatMasker/repeatModeler/species-families.fa"
	else: 
		input = config_dict["Input"]["RM_lib"]
	return(input)

rule RepeatMasker_custom:
	input: repeatMaskerCustom_input
	output: f"{outDir}/RepeatMasker/repeatModeler/genome.fasta.out"
	threads: config_dict["Threads"]["RepeatMasker_custom"]
	run:
		shell(f"RepeatMasker {config_dict["Internal"]["RepeatMasker"]} \
			-pa {threads} \
			-lib {input} \
			-dir {outDir}/RepeatMasker/repeatModeler \
			{outDir}/genome.fasta")

rule RepeatMasker_merge:
	input:
		repeatmasker = f"{outDir}/RepeatMasker/repeatMasker/genome.fasta.out",
		repeatmodeler = f"{outDir}/RepeatMasker/repeatModeler/genome.fasta.out"
	output: f"{outDir}/RepeatMasker/genome.masked.fasta"
	threads: config_dict["Threads"]["RepeatMasker_merge"]
	run:
		shell(f"scripts/merge_repeatMasker_out.pl \
			{outDir}/genome.fasta \
			{input.repeatmasker} \
			{input.repeatmodeler} \
			> {outDir}/RepeatMasker/genome.repeat.stats")
		shell(f"scripts/maskedByGff.pl \
			genome.repeat.gff3 \
			{outDir}/genome.fasta \
			> {outDir}/RepeatMasker/genome.masked.fasta")

rule TrimGalore_paired:
	input: 
		fwd = RNAseq_dir + "/{run}_1.fq.gz",
		rev = RNAseq_dir + "/{run}_2.fq.gz"
	output: 
		fwd = outdir + "/TrimGalore/{run}_1.fq.gz",
		rev = outdir + "/TrimGalore/{run}_2.fq.gz"
	threads: config_dict["Threads"]["TrimGalore_paired"]
	run:
		shell("trim_galore --paired \
			--three_prime_clip_R1 5 \
			--three_prime_clip_R2 5 \
            --cores {threads} \
			--max_n 40 \
			--gzip \
			-o " + outdir + "/TrimGalore \
			{input.fwd} {input.rev}")

def combine_input(pair, wildcards):
	input = expand(outdir + "/TrimGalore/{run}" + pair + ".fq.gz", run=wildcards.run)
	return(input)

rule combinePaired: 
	input:
		fwd = combinePaired_input("_1", wildcards),
		rev = combinePaired_input("_2", wildcards)
	output:
		fwd = outdir + "/TrimGalore/read1.fq.gz",
		rev = outdir + "/TrimGalore/read2.fq.gz"
	run:
		shell("cat {input.fwd} > {output.fwd}")
		shell("cat {input.rev} > {output.rev}")

rule TrimGalore_single:
	input: RNAseq_dir + "/{run}.fq.gz"
	output: outdir + "/TrimGalore/{run}.fq.gz"
	threads: config_dict["Threads"]["TrimGalore_single"]
	run:
		shell("trim_galore \
			--three_prime_clip 5 \
            --cores {threads} \
			--max_n 40 \
			--gzip \
			-o " + outdir + "/TrimGalore \
			{input}")

rule combinePaired: 
	input: combine_input("", wildcards)
	output: outdir + "/TrimGalore/readSE.fq.gz"
	run:
		shell("cat {input} > {output}")

rule HiSat2_build:
	input: f"{outDir}/RepeatMasker/genome.masked.fasta"
	output: outdir + "/HiSat2/genome.1.ht2"
	run:
		shell(f"hisat2-build {config_dict["Internal"]["hisat2-build"]}" + " {input} " + outdir + "/HiSat2/genome")

def HiSat2_input(wildcards):
	input = [built = "HiSat2/genome"]
	if config_dict["Input"]["rnaseq_single"] != "":
		input.append([se = outdir + "/TrimGalore/readSE.fq.gz"][0])
	if config_dict["Input"]["rnaseq_first"] != "":
		input.append([fwd = outdir + "/TrimGalore/read1.fq.gz"][0])
		input.append([rev = outdir + "/TrimGalore/read2.fq.gz"][0])

rule HiSat2:
	input: HiSat2_input
	output: outdir + "/HiSat2/hisat2.sam"
	threads: config_dict["threads"]["HiSat2"]
	run:
		cmd = ""
		if config_dict["Input"]["rnaseq_single"] != "":
			cmd += " -U {input.se}"
		if config_dict["Input"]["rnaseq_first"] != "":
			cmd += " -1 {input.fwd} -2 {input.rev}"
		if config_dict["Input"]["strand_specific"] == "True":
			cmd += " --rna-strandness RF"
		
		shell("hisat2 -x {input.built} \
			-p {threads} \
			" + cmd + " \
			-S " + outdir + "Hisat2/hisat2.sam \
			" + config["Internal"]["hisat2"])

rule HiSat2_sort:
	input: outdir + "/HiSat2/hisat2.sam"
	output: 
		outdir + "/HiSat2/hisat2.sorted.sam"
		outdir + "/HiSat2/splited_sam_files.list"
	run:
		shell("samtools sort  -o hisat2.sorted.bam -O BAM hisat2.sam")
		shell("samtools view -h hisat2.sorted.bam > hisat2.sorted.sam")
		# Split the SAM for use in transfrag and transdecoder
		shell("scripts/split_sam_from_non_aligned_region \
			" + outdir + "/HiSat2/hisat2.sorted.sam \
			splited_sam_out 10 > " + outdir + "/HiSat2/splited_sam_files.list")

rule Sam2Transfrag:
	input:
	output:
	run:
		shell("$dirname/bin/sam2transfrag $config{'sam2transfrag'} $no_strand_specific --intron_info_out $_.intron $_.sam > $_.gtf\n";

#TODO: Step 4 (homolog) relies on input from step 3

