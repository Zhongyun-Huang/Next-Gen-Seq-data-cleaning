#!/usr/bin/ruby 
$VERBOSE = true

require 'getoptlong'

# == Synopsis
#
# converts a set of hapmap files to structure format.  It creates one file, regardless of how many 
# are input.
# 
# Note: the groups must be numbered 1, 2, 3, etc.  Can't skip a number.
#
#
# == Usage
# 
#   ./hap2struct.rb -i inlist.txt -o for_struct_Olsen2_n10000.stru -m test_membership.txt -p populations.txt -n 0
#
#  -i is followed by the input file.
#  -o is followed by the output file.
#  -m is followed by the indiv. pop. membership file.
#  -p is followed by the population codes file.
#  -n is the minimum interval (distance) between SNPs, used to reduce the number of SNPs and 
#  distribute them more evenly.
# 
# == Author
# Ned Young, Caicedo Lab, U Mass Biology
# 06/05/12
# based on two_char.rb, 05/25/12
# based on hap2amap.rb, 4/04/12
#

####################################################################################################

######################
# METHOD DEFINITIONS #
######################

# =get_members
#
# Pull in the taxon numbers and population names from a file.
# incoming argument is a filename.
# returns two hashes; one for locales and one for the broader regions

def get_members( pm_fn )
	members = Hash.new
	reg_members = Hash.new
	pm_fh = File.new( pm_fn )
	pm_fh.each do |line| 
		next if (!line || line =~ /^(\s)*$/)   # skips blank lines								
		next if (!line || line =~ /^Taxon/)   # skips header line								
		raw = line.chomp.split(/\t+/)  # tab defines field limits
		num = raw[0].strip
		gen_loc = raw[1]
		spec_loc = raw[2]
		members[ num ] = spec_loc.strip
		reg_members[ num ] = gen_loc.strip
	end
	return members, reg_members
end

# =get_pop_codes
#
# Pull in the population names and codes from a file.
# incoming argument is a filename.
# returns a hash; keys are pop names; values are numbers

def get_pop_codes( pop_fn )
	pop_codes = Hash.new
	pop_fh = File.new( pop_fn )
	pop_fh.each do |line| 
		next if (!line || line =~ /^(\s)*$/)   # skips blank lines								
		raw = line.chomp.split(/\t+/)  # tab defines field limits
		gen_loc = raw[0]
		code_num = raw[1].strip.to_i
		pop_codes[ gen_loc ] = code_num
	end
	return pop_codes
end

####################################################################################################

# Main program

# First check that arg(s) are given on the command line
fail "\nUsage: #$0 -i inlist.txt -o for_struct_Olsen2_n10000.stru -m test_membership.txt -p populations.txt -n 0\n" if ARGV.size < 5

opts = GetoptLong.new(
	[ "-i",				GetoptLong::REQUIRED_ARGUMENT ],
	[ "-o",				GetoptLong::REQUIRED_ARGUMENT ],
	[ "-m",				GetoptLong::REQUIRED_ARGUMENT ],
	[ "-p",				GetoptLong::REQUIRED_ARGUMENT ],
	[ "-n",				GetoptLong::REQUIRED_ARGUMENT ]
)
inlist_fn = output_fn = pop_fn = pm_fn = nil
interval = 5000
opts.each do |opt, arg|
	case opt
		when '-i'
			inlist_fn = arg
		when '-o'
			output_fn = arg
		when '-m'
			pm_fn = arg
		when '-p'
			pop_fn = arg
		when '-n'
			interval = arg.to_i
	end
end											 

members, reg_members = get_members( pm_fn )
krm = reg_members.keys
pop_codes = get_pop_codes( pop_fn )

# Start data structure (hash) to store every genotype in
ambig_codes = [ "M", "R", "W", "S", "K", "Y" ]
het_counter = 0
name_list = Array.new
big_struct = Hash.new
num_snps = 0
stru_handle = File.new( output_fn, "w" )
he_rpt = File.new( "hets_etc.rpt", "w" )																												#temp
hes_struct = Hash.new																																									
hei_struct = Hash.new
ind_ids = Array.new
snps = Array.new
kept_count = 0
too_close = 0
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
  old_pos =  0 - interval
  # input the hapmap file, which needs to be in the same folder as this program
  hmp_fh = File.new( hmp_fn )
  hmp_fh.each do |line| 
  	next if line =~ /^\s*$/  # checks for blank line
  	raw = line.chomp.split(/\s+/)  # space defines field limits
  	snp_name = raw[ 0 ]
  	snp_pos = raw[ 3 ].to_i
  	chrom_num =  raw[ 2 ].to_i

  	if chrom_num != old_chrom_num
  		old_pos =  0 - interval   # Looks to see when new chromosome and resets old_pos
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
  	elsif snp_pos - old_pos < interval
  		too_close += 1
  		next
  	else
			snps.push( snp_name )
			# Next, need to loop thru the genotypes, convert to array of codes and save each one in 	 
			# big_struct hash to print later.
			# Use integers to find the individual's name (from ids) and the data (from raw).
			num_indivs = raw.size
			last_i = num_indivs - 1
			for i in 0..last_i
				ind_id = ind_ids[ i ]
				if ind_ids.include?( ind_id )							#temp
				else
					#puts "snp_name: #{snp_name}; -#{ind_id}- is not included in ids"													#temp
					#puts "i: #{i}; ids[i + 1]: -#{ids[i + 1]}-"																								#temp
					next
				end																																													#temp
				geno = raw[ i ]
  			doublet = [ geno, snp_name ]
  			big_struct[ ind_id ].push( doublet )			
			end
			#stuff = "keep\t#{snp_name}\t#{snp_pos}\t#{old_pos}\n" 																				#temp
  		#he_list.push( stuff )																																					#temp
			old_pos = snp_pos
  	end  # Finished a line of data

  end  # Finished processing a single hapmap file
  first_file = false
end # Finished processing all of the hapmap files

# Print header line
puts "Based on the snp names array, it looks like there are #{snps.size} SNPs."  									#temp
snp_string = snps.join("\t")
# stru_handle.puts "#individual\tgroup\t#{snp_string}" # Proper STRUCTURE format has no 
# header line.  You can put this back in if it i useful for troubleshooting, etc.
# Print out the indivs and genos from all chromosomes to a single file
g_codes = { "A" => 	[ "1", "1" ],
						"C" => 	[ "2", "2" ],
						"G" => 	[ "3", "3" ],
						"T" => 	[ "4", "4" ],
						"R" => 	[ "1", "3" ],
						"Y" => 	[ "2", "4" ],
						"S" => 	[ "3", "2" ],
						"W" => 	[ "1", "4" ],
						"K" => 	[ "3", "4" ],
						"M" => 	[ "1", "2" ],
						"B" => 	[ "-9", "-9" ],
						"D" => 	[ "-9", "-9" ],
						"H" => 	[ "-9", "-9" ],
						"V" => 	[ "-9", "-9" ],
						"N" => 	[ "-9", "-9" ] }
ind_ids.each do |ind_id|
	genos = big_struct[ ind_id ]
	puts "Based on the genos array, it looks like there are #{genos.size} SNPs." if ind_id == "103830" 		#temp
	# Create back part for both lines
	back_part_of_top_line = Array.new
	back_part_of_bot_line = Array.new
	genos.each do | doublet |
		geno, snp_name = doublet
		small_array = g_codes[ geno ]
		top_code, bot_code = small_array
		puts "top_code: #{top_code}; bot_code: #{bot_code} when geno == W" if geno == "W" #temp
		back_part_of_top_line.push( top_code )			  				
		back_part_of_bot_line.push( bot_code )			  				
		kept_count += 1	if ind_id == "103830" 		#temp			
	end
	region = reg_members[ ind_id ]
	pop_codes.each_pair do |k,v|																														#temp
		puts k, v																																								#temp
	end																																											#temp
	grp = pop_codes[ region ]
	back_part_of_top_string = back_part_of_top_line.join("\t")
	back_part_of_bot_string = back_part_of_bot_line.join("\t")
	stru_handle.puts "#{ind_id}\t#{grp}\t#{back_part_of_top_string}"
	stru_handle.puts "#{ind_id}\t#{grp}\t#{back_part_of_bot_string}"
	puts "Based on the back_part array +2, it looks like there are #{back_part_of_top_line.size + 2} fields in the line." if ind_id == "103830" 		#temp
end
puts "SNPs kept: #{kept_count}"

list_fh = File.new( "list.snp", "w" )																									#temp
snps_kept = Array.new
snp_list_string	= snps.join("\n")																															#temp
list_fh.puts snp_list_string																																				#temp
he_rpt.puts "Total number of ambiguity codes changes to N:        #{het_counter}"
het_ind_ct = 0
het_snp_ct = 0
hei_struct.each_pair do |ind, ct|
	het_ind_ct += 1 if ct > 0
end
he_rpt.puts "Number of individuals with at least one het changed: #{het_ind_ct}"																																									#temp
he_rpt.puts "Total number of individuals:                         #{ind_ids.size}"																																									#temp
het_ind_ct = 0
hes_struct.each_pair do |snp, ct|
	het_snp_ct += 1 if ct > 0
end
he_rpt.puts "Number of SNPs with at least one het changed: #{het_snp_ct}"																																									#temp
he_rpt.puts "Total number of SNPs:                         #{snps.size}"																																									#temp

exit
