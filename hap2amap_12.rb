#!/usr/bin/ruby 
$VERBOSE = true

require 'getoptlong'

# == Synopsis
#
# converts a set of hapmap files to ancestry map format.  For each file, it:
#																					 Creates three files: my_SNPs.c1.ancestrymapgeno 
#  																															my_SNPs.c1.snp
#  											   																			my_SNPs.c1.ind
#																					(assuming the input file was my_SNPs.c1.hmp.txt)
# 
#
# The genotype file (.ancestrymapgeno) contains 1 line per valid genotype.  There are 3 columns:
#  1st column is SNP name
#  2nd column is sample ID
#  3rd column is number of reference alleles (0 or 1 or 2)
# Missing genotypes are encoded by the absence of an entry in the genotype file.
#
# The snp file contains 1 line per SNP.  There are 6 columns (last 2 optional):
#   1st column is SNP name
#   2nd column is chromosome.  X chromosome is encoded as 23.
#     Also, Y is encoded as 24, mtDNA is encoded as 90, and XY is encoded as 91.
#     Note: SNPs with illegal chromosome values, such as 0, will be removed
#   3rd column is genetic position (in Morgans).  If unknown, ok to set to 0.0.
#   4th column is physical position (in bases)
#   Optional 5th and 6th columns are reference and variant alleles.
#     For monomorphic SNPs, the variant allele can be encoded as X (unknown).
#
# The indiv file contains 1 line per individual.  There are 3 columns:
#   1st column is sample ID.  Length is limited to 39 characters, including
#     the family name if that will be concatenated.
#   2nd column is gender (M or F).  If unknown, ok to set to U for Unknown.
#   3rd column is a label which might refer to Case or Control status, or
#     might be a population group label.  If this entry is set to "Ignore", 
#     then that individual and all genotype data from that individual will be
#     removed from the data set in all convertf output.
#
# == Usage
# 
#   ./hap2amap.rb -i inlist.txt 
# 
# == Author
# Ned Young, Caicedo Lab, U Mass Biology
# 4/04/12
#

####################################################################################################

######################
# METHOD DEFINITIONS #
######################

# =ref_alleles
#
# Gives the number of alleles that match the reference genome sequence (0, 1, or 2).  Assumes that 
# the first allele listed in field 1 (second field) is the reference genome's allele.
# Gets passed the base and the "A/G" type string that identifies the SNPs alleles.   
# Returns 0, 1, or 2 (single integer) 

def ref_alleles( geno, allele_string ) 
	alleles = allele_string.split(/\//)
	ref = alleles[0]
	num_ref_alleles = ""
	# MRWSYK
	if ref == "A"
		if geno == "A" 
			num_ref_alleles = 2
		elsif (geno == "M") || (geno == "R") ||(geno == "W")
			num_ref_alleles = 1
		elsif geno == "N"
			num_ref_alleles = ""
		else 
			num_ref_alleles = 0
		end
	elsif ref == "C"
		if geno == "C" 
			num_ref_alleles = 2
		elsif (geno == "M") || (geno == "S") ||(geno == "Y")
			num_ref_alleles = 1
		elsif geno == "N"
			num_ref_alleles = ""
		else 
			num_ref_alleles = 0
		end
	elsif ref == "G"
		if geno == "G" 
			num_ref_alleles = 2
		elsif (geno == "R") || (geno == "S") ||(geno == "K")
			num_ref_alleles = 1
		elsif geno == "N"
			num_ref_alleles = ""
		else 
			num_ref_alleles = 0
		end
	elsif ref == "T"
		if geno == "T" 
			num_ref_alleles = 2
		elsif (geno == "W") || (geno == "Y") ||(geno == "K")
			num_ref_alleles = 1
		elsif geno == "N"
			num_ref_alleles = ""
		else 
			num_ref_alleles = 0
		end
	end
	return num_ref_alleles
end

####################################################################################################

# Main program

# First check that arg(s) are given on the command line
fail "\nUsage: #$0 -i inlist.txt\n" if ARGV.size < 1

opts = GetoptLong.new(
	[ "-i",				GetoptLong::REQUIRED_ARGUMENT ]
)
inlist_fn = nil
opts.each do |opt, arg|
	case opt
		when '-i'
			inlist_fn = arg
	end
end											 

# Keep track of Ns for each indiv
# [0] is just those with zero Ns.  They are also in [1], which is for frequencies from 0-0.1; [2] is for freqs. > 0.1, <= 0.2, etc.
# [11] is also special, recording the number of indivs with all Ns.
# ncount = Hash.new  
# ncount[0] = 0
# ncount[1] = 0
# num_indivs = 0

# Start data structure (hash) to store every genotype in
name_list = Array.new
big_struct = Hash.new
num_snps = 0
amap_handle = nil
snp_handle = nil
ind_handle = nil
ids = Array.new
genos = Array.new
first_file = true
inlist_fh = File.new( inlist_fn )
inlist_fh.each do |line| 
	next if line =~ /^\s*$/  # checks for blank line
	# Get hapmap file name from line
	hmp_fn = line.chomp.strip  
	puts "processing file: #{hmp_fn}\n"
	
	if first_file == true
		#create 3 filehandles based on the input filename
		in_fn_chunks = hmp_fn.split(/\./)
		in_fn_chunks.pop
		in_fn_chunks.push("ancestrymapgeno")
		amap_fn = in_fn_chunks.join(".")
		in_fn_chunks.pop
		in_fn_chunks.push("snp")
		snp_fn = in_fn_chunks.join(".")
		in_fn_chunks.pop
		in_fn_chunks.push("ind")
		ind_fn = in_fn_chunks.join(".")
		amap_handle = File.new( amap_fn, "w" )
		snp_handle = File.new( snp_fn, "w" )
		ind_handle = File.new( ind_fn, "w" )
	end
	
  # Make an array to hold the names and genos for each snp.
  snp_genos = Array.new
  # input the hapmap file, which needs to be in the same folder as this program
  hmp_fh = File.new( hmp_fn )
  hmp_fh.each do |line| 
  	next if line =~ /^\s*$/  # checks for blank line
  	raw = line.chomp.split(/\s+/)  # space defines field limits
  	snp_name = raw[ 0 ]
  	# Remove first eleven columns; just keep the indiv data (names on the first line, SNPs thereafter)
  	front_matter = raw.slice!(0..10)
  	if line =~ /^rs#/  # header line
  		if first_file == true
				# Put the ids into an array with zero-based numbering
				raw.each do |id|
					chunks = id.split(/:/)
  				short_id = chunks[0]
					ids.push( short_id )
					ind_handle.puts "#{short_id}\tU\tX"
				end
			end
  	else
			allele_string = front_matter[1]
			alleles = allele_string.split(/\//)
			ref = alleles[0]
			other = alleles[1]
			chrom_num = front_matter[2]
			pos = front_matter[3]
			# Next, need to loop thru the genotypes, convert to num_ref_alleles and save each one to file 
			# to print later
			genos = []
			raw.each do |geno|
  			genos.push( geno )
  		end
  	

  	end  # Finished a line of data
  	septuplet = [ snp_name, allele_string, ref, other, chrom_num, pos, genos ]
  	snp_genos.push( septuplet ) unless snp_name == "rs#"
  end  # Finished processing a single hapmap file
  # Print out the SNPs and genos from one chromosome
  snp_genos.each do |septuplet|
  	(snp_name, allele_string, ref, other, chrom_num, pos, genos) = septuplet
  	#genos_string = genos.join("\t")
  	num_genos = genos.size
  	last_i = num_genos - 1
  	for i in 0..last_i
  		geno = genos[ i ]
  		id = ids[ i ]
  		num_ref_alleles = ref_alleles( geno, allele_string )
  	  unless !num_ref_alleles || num_ref_alleles == ""
  	  	amap_handle.puts "#{snp_name}\t#{id}\t#{num_ref_alleles}" 
  	  end
  	end
  	snp_handle.puts "#{snp_name}\t#{chrom_num}\t0.0\t#{pos}\t#{ref}\t#{other}"
  end
  first_file = false
end # Finished processing all of the hapmap files

exit
