#!/usr/bin/ruby 
$VERBOSE = true

require 'getoptlong'

# == Synopsis
#
# converts a set of hapmap files to PHYLIP sequential format.  It creates one file, 
# regardless of how many are input.
# (num taxa and chars on first line)
# 
# == Usage
# 
#   ./hap2phy.rb -i lyc_inlist.txt -o lyciumSNPs.phy
#
# 
# == Author
# Ned Young, Miller Lab, Amherst College
# 06/12/14
# based on hap2struct.rb, by Ned Young, Caicedo Lab, U Mass Biology, 06/05/12
# based on two_char.rb, 05/25/12
# based on hap2amap.rb, 4/04/12
#

####################################################################################################

# Main program

# First check that arg(s) are given on the command line
fail "\nUsage: #$0 -i lyc_inlist.txt -o lyciumSNPs.phy\n" if ARGV.size < 2

opts = GetoptLong.new(
	[ "-i",				GetoptLong::REQUIRED_ARGUMENT ],
	[ "-o",				GetoptLong::REQUIRED_ARGUMENT ]
)
inlist_fn = nil
output_fn = nil
opts.each do |opt, arg|
	case opt
		when '-i'
			inlist_fn = arg
		when '-o'
			output_fn = arg
	end
end											 

# Start data structure (hash) to store every genotype in
name_list = Array.new
big_struct = Hash.new
num_snps = 0
num_taxa = 0
num_indivs = 0
phy_fh = File.new( output_fn, "w" )
ind_ids = Array.new
snps = Array.new
old_chrom_num = 0
first_file = true
inlist_fh = File.new( inlist_fn )
inlist_fh.each do |line| 
	next if line =~ /^\s*$/  # checks for blank line
	# Get hapmap file name from line
	hmp_fn = line.chomp.strip  
	puts "processing file: #{hmp_fn}\n"
  # Make an array to hold the names and genos for each snp.
  snp_genos = Array.new
  # input the hapmap file, which needs to be in the same folder as this program
  hmp_fh = File.new( hmp_fn )
  hmp_fh.each do |line| 
  	next if line =~ /^\s*$/  # checks for blank line
  	raw = line.chomp.split(/\s+/)  # space defines field limits
  	snp_name = raw[ 0 ]
  	snp_pos = raw[ 3 ].to_i
  	chrom_num =  raw[ 2 ].to_i

  	if chrom_num != old_chrom_num    # Looks to see when new chromosome
  		old_chrom_num = chrom_num
  	end
  	# Remove first eleven columns; just keep the indiv data (names on the first line, SNPs thereafter)
  	raw.slice!(0..10)
  	if line =~ /^rs#/  # header line
  		if first_file == true
				# Put the indiv ids into an array with zero-based numbering
				raw.each do |ind_id|
					chunks = ind_id.split(/:/)
  				short_ind_id = chunks[0]
					ind_ids.push( short_ind_id )
					# initialize an array to store genotype info in
					big_struct[ short_ind_id ] = Array.new
				end
			end
  	else
			snps.push( snp_name )
			num_snps += 1
			# Next, need to loop thru the genotypes, convert to array of codes and save each one in 	 
			# big_struct hash to print later.
			# Use integers to find the individual's name (from ids) and the data (from raw).
			num_indivs = raw.size
			last_i = num_indivs - 1
			for i in 0..last_i
				ind_id = ind_ids[ i ]
				geno = raw[ i ]
  			doublet = [ geno, snp_name ]
  			big_struct[ ind_id ].push( doublet )			
			end
			old_pos = snp_pos
  	end  # Finished a line of data

  end  # Finished processing a single hapmap file
  first_file = false
end # Finished processing all of the hapmap files

# Print header line
puts "Based on the snp names array, it looks like there are #{snps.size} SNPs."  									
num_taxa = ind_ids.size
snp_string = snps.join("\t")
phy_fh.puts "\t#{num_taxa}\t#{num_snps}" # 
# Print out the indivs and genos from all chromosomes to a single file
ind_ids.each do |ind_id|
	genos = big_struct[ ind_id ]
	puts "Based on the genos array, it looks like there are #{genos.size} SNPs." if ind_id == "103830" 		#temp
	# Create back part for both lines
	back_part = Array.new
	genos.each do | doublet |
		geno, snp_name = doublet
		code = geno 
		back_part.push( code )			  				
	end
	#grp = "1"
	back_part_string = back_part.join("")
	phy_fh.puts "#{ind_id}\t#{back_part_string}"
end																																							#temp

exit
