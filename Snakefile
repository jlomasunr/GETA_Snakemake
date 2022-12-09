import yaml
import logging
import Bio
from rich.logging import Console, RichHandler

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
	#TODO: check augustus install for species, error and quit if not present 
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

#TODO: Check software dependencies 
# RepeatModeler, RepeatMasker, HPCgridrunner, java, hisat2, samtools, hmmer

# Setup output folders
os.makedirs(f"{outdir}/RepeatMasker/repeatMasker/", exist_ok=True)
os.makedirs(f"{outdir}/RepeatMasker/repeatModeler/", exist_ok=True)
os.makedirs(f"{outdir}/trimmomatic/", exist_ok=True)

# Copy over genome and protein file
# TODO: Make more robust?
if not os.path.exists(f"{outDir}/genome.fasta"):
	os.system(f"cp {config_dict["genome"]} genome.fasta")
if not os.path.exists(f"{outDir}/homolog.fasta"):
	os.system(f"cp {config_dict["protein"]} homolog.fasta")

def allInput():
	input = [f"{outDir}/RepeatMasker/repeatMasker/genome.fasta.out",		# RepeatMasker_species
			 f"{outDir}/RepeatMasker/repeatModeler/genome.fasta.out",		# RepeatMasker_custom
			 f"{outDir}/RepeatMasker/genome.masked.fasta",				# RepeatMasker_merge

			]
	if config_dict["Input"]["RM_lib"] == "None":
		input.append(f"{outDir}/RepeatMasker/repeatModeler/species-families.fa")		# RepeatModeler
		input.append(f"{outDir}/RepeatMasker/repeatModeler/species-families.stk")		# RepeatModeler
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
	threads:
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

rule trimmomatic:
	input:
	output:
	threads:
	run:
	use trim galore