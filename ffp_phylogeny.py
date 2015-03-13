#!/usr/bin/python
import optparse
import time
import os
import tempfile
import sys
import shlex, subprocess
from string import maketrans

VERSION_NUMBER = "0.1.00"

class MyParser(optparse.OptionParser):
	"""
	 From http://stackoverflow.com/questions/1857346/python-optparse-how-to-include-additional-info-in-usage-output
	 Provides a better class for displaying formatted help info in epilog() portion of optParse; allows for carriage returns.
	"""
	def format_epilog(self, formatter):
		return self.epilog


def stop_err( msg ):
    sys.stderr.write("%s\n" % msg)
    sys.exit(1)

def getTaxonomyNames(type, multiple, abbreviate, filepaths, filenames):
	"""
	Returns a taxonomic list of names corresponding to each file being analyzed by ffp.
	This may also include names for each fasta sequence found within a file if the
	"-m" multiple option is provided. 	Default is to use the file names rather than fasta id's inside the files.
	NOTE: THIS DOES NOT (MUST NOT) REORDER NAMES IN NAME ARRAY. 
	EACH NAME ENTRY IS TRIMMED AND MADE UNIQUE
	
	@param type string ['text','amino','nucleotide']
	@param multiple boolean Flag indicates to look within files for labels
	@param abbreviate boolean Flag indicates to shorten labels	
	@filenames array original input file names as user selected them
	@filepaths array resulting galaxy dataset file .dat paths
	
	"""
	# Take off prefix/suffix whitespace/comma :
	taxonomy = filenames.strip().strip(',').split(',')
	translations = maketrans(' .-	','____')
	names=[]
	ptr = 0

	for file in filepaths:
		# First, convert space, period to underscore in file names.	  ffprwn IS VERY SENSITIVE ABOUT THIS.
		# Also trim labels to 50 characters.  Turns out ffpjsd is kneecapping a taxonomy label to 10 characters if it is greater than 50 chars.
		taxonomyitem = taxonomy[ptr].strip().translate(translations)[:50]
		# print taxonomyitem
		if not type in 'text' and multiple:
			#Must read each fasta file, looking for all lines beginning ">"
			with open(file) as fastafile:
				lineptr = 0
				for line in fastafile:
					if line[0] == '>':
						name = line[1:].split(None,1)[0].strip()[:50]
						# Odd case where no fasta description found
						if name == '': name = taxonomyitem + '.' + str(lineptr)
						names.append(name)
						lineptr += 1
		else:
			names.append(taxonomyitem)
		
		ptr += 1

	if abbreviate:
		names = trimCommonPrefixes(names)
		names = trimCommonPrefixes(names, True) # reverse = Suffixes.

	return names
	
def trimCommonPrefixes(names, reverse=False):
	"""
	Examines sorted array of names.  Trims off prefix of each subsequent pair.
	
	@param names array of textual labels (file names or fasta taxonomy ids)
	@param reverse boolean whether to reverse array strings before doing prefix trimming.
	"""
	wordybits = '|.0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'

	if reverse:
		names = map(lambda name: name[::-1], names) #reverses characters in names
	
	sortednames = sorted(names)
	ptr = 0
	sortedlen = len(sortednames)
	oldprefixlen=0
	prefixlen=0
	for name in sortednames:
		ptr += 1

		#If we're not at the very last item, reevaluate prefixlen
		if ptr < sortedlen:

			# Skip first item in an any duplicate pair.  Leave duplicate name in full.
			if name == sortednames[ptr]:
				if reverse:
					continue
				else:
					names[names.index(name)] = 'DupLabel-' + name
					continue

			# See http://stackoverflow.com/questions/9114402/regexp-finding-longest-common-prefix-of-two-strings
			prefixlen = len( name[:([x[0]==x[1] for x in zip(name, sortednames[ptr])]+[0]).index(0)] )
				
		if prefixlen <= oldprefixlen:
			newprefix = name[:oldprefixlen]
		else:
			newprefix = name[:prefixlen]
		# Expands label to include any preceeding characters that were probably part of it.
		newprefix = newprefix.rstrip(wordybits)
		newname = name[len(newprefix):]
		# Some tree visualizers don't show numeric labels?!?!
		if not reverse and newname.replace('.','',1).isdigit():
			newname = 'id_' + newname 
		names[names.index(name)] = newname #extract name after prefix part; has nl in it
		oldprefixlen = prefixlen

	if reverse:
		names = map(lambda name: name[::-1], names) #now back to original direction
	
	return names

def getTaxonomyFile(names):
	"""
	FFP's ffpjsd -p [taxon file of labels] option creates a phylip tree with 
	given taxon labels
	
	@param names array of datafile names or fasta sequence ids
	"""

	try:
		temp = tempfile.NamedTemporaryFile(mode='w+t',delete=False)
		taxonomyTempFile = temp.name
		temp.writelines(name + '\n' for name in names)

	except: 
		stop_err("Galaxy configuration error for ffp_phylogeny tool. Unable to write taxonomy file " + taxonomyTempFile)

	finally:
		temp.close()
	
	return taxonomyTempFile


def check_output(command):
	"""
	Execute a command line containing a series of pipes; and handle error cases by exiting at first error case.  This is a substitute for Python 2.7 subprocess.check_output() - allowing piped commands without shell=True call .  Based on Python subprocess docs 17.1.4.2

	ISSUE: warnings on stderr are given with no exit code 0:
		ffpry: Warning: No keys of length 6 found.
		ffpcol: (null): Not a key valued FFP.

	Can't use communicate() because this closes processes' stdout 
	file handle even without errors because of read to end of stdout:
	(stdoutdata, stderrdata) = processes[ptr-1].communicate()

	"""
	commands = command.split("|")
	processes = []
	ptr = 0
	for command_line in commands:
		print 'testing,' , command_line.strip()
		args = shlex.split(command_line.strip())
		if ptr == 0:
			proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			processes.append(proc)
		else:

			#this has to come before error processing?
			newProcess = subprocess.Popen(args, stdin=processes[ptr-1].stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			
			# It seems the act of reading standard error output is enough to trigger
			# error code signal for that process, i.e. so that retcode returns a code.
			retcode = processes[ptr-1].poll()
			stderrdata = processes[ptr-1].stderr.read()
			if retcode or len(stderrdata) > 0:
				stop_err(stderrdata)			

			processes.append(newProcess)			
			processes[ptr-1].stdout.close() # Allow prev. process to receive a SIGPIPE if current process exits.
		
		ptr += 1

	retcode = processes[ptr-1].poll()
	(stdoutdata, stderrdata) = processes[ptr-1].communicate()
	if retcode or len(stderrdata) > 0:
		stop_err(stderrdata)
	
	return stdoutdata
	

class ReportEngine(object):

	def __init__(self): pass

	def __main__(self):


		## *************************** Parse Command Line *****************************
		parser = MyParser(
			description = 'FFP (Feature frequency profile) is an alignment free comparison tool',
			usage = 'python ffp_phylogeny.py [input_files] [output file] [options]',
			epilog="""Details:

			FFP (Feature frequency profile) is an alignment free comparison tool for phylogenetic analysis and text comparison. It can be applied to nucleotide sequences, complete genomes, proteomes and even used for text comparison.
		
		""")

		parser.set_defaults(row_limit=0)
		# Don't use "-h" , it is reserved for --help!

		parser.add_option('-t', '--type', type='choice', dest='type', default='text', 
			choices=['amino','nucleotide','text'],
			help='Choice of Amino acid, nucleotide or plain text sequences to find features in')

		parser.add_option('-l', '--length', type='int', dest='length', default=6, 
			help='Features (any string of valid characters found in data) of this length will be counted.  Synonyms: l-mer, k-mer, n-gram, k-tuple')

		#parser.add_option('-n', '--normalize', dest='normalize', default=True, action='store_true', 
		#	help='Normalize counts into relative frequency')
			
		parser.add_option('-m', '--multiple', dest='multiple', default=False, action='store_true', 
			help='By default all sequences in a fasta file be treated as 1 sequence to profile.  This option enables each sequence found in a fasta file to have its own profile.')

		parser.add_option('-M', '--metric', type='string', dest='metric', 
			help='Various metrics to measure count distances by.')

		parser.add_option('-x', '--taxonomy', type='string', dest='taxonomy', 
			help='Taxanomic label for each profile/sequence.')

		parser.add_option('-d', '--disable', dest='disable', default=False, action='store_true', 
			help='By default amino acid and nucleotide characters are grouped by functional category (protein or purine/pyrimidine group) before being counted.  Disable this to treat individual characters as distinct.')

		parser.add_option('-a', '--abbreviate', dest='abbreviate', default=False, action='store_true', 
			help='Shorten tree taxonomy labels as much as possible.')
		
		parser.add_option('-s', '--similarity', dest='similarity', default=False, action='store_true', 
			help='Enables pearson correlation coefficient matrix and any of the binary distance measures to be turned into similarity matrixes.')
		
		parser.add_option('-f', '--filter', type='choice', dest='filter', default='none',
			choices=['none','f','n','e','freq','norm','evd'],
			help='Choice of [f=raw frequency|n=normal|e=extreme value (Gumbel)] distribution: Features are trimmed from the data based on lower/upper cutoff points according to the given distribution.')

		parser.add_option('-L', '--lower', type='float', dest='lower', 
			help='Filter lower bound is a 0.00 percentages')
		
		parser.add_option('-U', '--upper', type='float', dest='upper',
			help='Filter upper bound is a 0.00 percentages')

		parser.add_option('-o', '--output', type='string', dest='output', 
			help='Path of output file to create')

		parser.add_option('-T', '--tree', dest='tree', default=False, action='store_true', help='Generate Phylogenetic Tree output file')

		parser.add_option('-v', '--version', dest='version', default=False, action='store_true', help='Version number')

		# Could also have -D INT decimal precision included for ffprwn .
			
		options, args = parser.parse_args()

		if options.version:
			print VERSION_NUMBER
			return
		
		import time
		time_start = time.time()

		try:
			in_files = args[:]
		
		except:
			stop_err("Expecting at least 1 input data file.")

		
		#ffptxt / ffpaa / ffpry
		if options.type in 'text':
			command = 'ffptxt'
			
		else:
			if options.type == 'amino':
				command = 'ffpaa'
			else:
				command = 'ffpry'
				
			if options.disable:
				command += ' -d'
				
			if options.multiple:
				command += ' -m'
		
		command += ' -l ' + str(options.length)

		if len(in_files): #Note: app isn't really suited to stdio
			command += ' "' + '" "'.join(in_files) + '"'
				
		#ffpcol / ffpfilt
		if options.filter != 'none':		
			command += ' | ffpfilt'
			if options.filter != 'count':
				command += ' -' + options.filter
			if options.lower > 0:
				command += ' --lower ' + str(options.lower)
			if options.upper > 0:
				command += ' --upper ' + str(options.upper)
			
		else:
			command += ' | ffpcol'

		if options.type in 'text':
			command += ' -t'
			
		else:

			if options.type == 'amino':
				command += ' -a'
			
			if options.disable:
				command += ' -d'
			
		#if options.normalize:
		command += ' | ffprwn'

		#Now create a taxonomy label file, ensuring a name exists for each profile.
		taxonomyNames = getTaxonomyNames(options.type, options.multiple, options.abbreviate, in_files, options.taxonomy)
		taxonomyTempFile = getTaxonomyFile(taxonomyNames)
		
		# -p = Include phylip format 'infile' of the taxon names to use.  Very simple, just a list of fasta identifier names.
		command += ' | ffpjsd -p ' + taxonomyTempFile

		if options.metric and len(options.metric) >0 :
			command += ' --' + options.metric
			if options.similarity:
				command += ' -s'

		# Generate Newick (.nhx) formatted tree if we have at least 3 taxonomy items:
		if options.tree:
			if len(taxonomyNames) > 2:
				command += ' | ffptree -q' 
			else:
				stop_err("For a phylogenetic tree display, one must have at least 3 ffp profiles.")

		print command
		
		result = check_output(command)
		with open(options.output,'w') as fw:
			fw.writelines(result)
		os.remove(taxonomyTempFile)

if __name__ == '__main__':

	time_start = time.time()

	reportEngine = ReportEngine()
	reportEngine.__main__()
	
	print('Execution time (seconds): ' + str(int(time.time()-time_start)))
	
