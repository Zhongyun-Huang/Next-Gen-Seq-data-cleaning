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
#   ./hap2struct.rb -i inlist.txt -o for_struct_Olsen2_n10000.stru -n 10000
#
# -n is the minimum interval (distance) between SNPs, used to reduce the number of SNPs and 
# distribute them more evenly.
# 
# == Author
# Ned Young, Caicedo Lab, U Mass Biology
# 06/05/12
# based on two_char.rb, 05/25/12
# based on hap2amap.rb, 4/04/12
#

####################################################################################################

# Main program

# First check that arg(s) are given on the command line
fail "\nUsage: #$0 -i inlist.txt -o for_struct_Olsen2_n10000.stru -n 10000\n" if ARGV.size < 3

opts = GetoptLong.new(
	[ "-i",				GetoptLong::REQUIRED_ARGUMENT ],
	[ "-o",				GetoptLong::REQUIRED_ARGUMENT ],
	[ "-n",				GetoptLong::REQUIRED_ARGUMENT ]
)
inlist_fn = nil
output_fn = nil
interval = 5000
opts.each do |opt, arg|
	case opt
		when '-i'
			inlist_fn = arg
		when '-o'
			output_fn = arg
		when '-n'
			interval = arg.to_i
	end
end											 

# Start data structure (hash) to store every genotype in
ambig_codes = [ "M", "R", "W", "S", "K", "Y" ]
name_list = Array.new
big_struct = Hash.new
num_snps = 0
stru_handle = File.new( output_fn, "w" )
kd_rpt = File.new( "kept-dropped.rpt", "w" )																												#temp
kd_list = Array.new																																									#temp
stuff = "d-k\tsnp_name\tsnp_pos\told_pos\n"																													#temp
kd_list.push( stuff )																																								#temp

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
  	#puts "chrom_num:  #{chrom_num}, old_chrom_num: #{old_chrom_num}" 															#temp

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
  		stuff = "drop\t#{snp_name}\t#{snp_pos}\t#{old_pos}\n" 																				#temp
  		kd_list.push( stuff )																																					#temp
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
  			if ambig_codes.include?( geno )   # Changes ambiguity codes to Ns
  				geno = "N"
  			end
  			doublet = [ geno, snp_name ]
  			big_struct[ ind_id ].push( doublet )			
			end
			stuff = "keep\t#{snp_name}\t#{snp_pos}\t#{old_pos}\n" 																				#temp
  		kd_list.push( stuff )																																					#temp
			old_pos = snp_pos
  	end  # Finished a line of data

  end  # Finished processing a single hapmap file
  first_file = false
end # Finished processing all of the hapmap files

# Print header line
puts "Based on the snp names array, it looks like there are #{snps.size} SNPs."  									#temp
snp_string = snps.join("\t")
stru_handle.puts "#individual\tgroup\t#{snp_string}" 
# Proper STRUCTURE format has no 
# header line.  You can put this back in if it i useful for troubleshooting, etc.
# Print out the indivs and genos from all chromosomes to a single file
g_codes = { "A" => 	"1",
						"C" => 	"2",
						"G" => 	"3",
						"T" => 	"4",
						"N" => 	"-9" }
ind_ids.each do |ind_id|
	genos = big_struct[ ind_id ]
	puts "Based on the genos array, it looks like there are #{genos.size} SNPs." if ind_id == "103830" 		#temp
	# Create back part for both lines
	back_part = Array.new
	genos.each do | doublet |
		geno, snp_name = doublet
		code = g_codes[ geno ]
		back_part.push( code )			  				
		kept_count += 1	if ind_id == "103830" 		#temp			
	end
	grp = "1"
	back_part_string = back_part.join("\t")
	stru_handle.puts "#{ind_id}\t#{grp}\t#{back_part_string}"
	puts "Based on the back_part array +2, it looks like there are #{back_part.size + 2} fields in the line." if ind_id == "103830" 		#temp
end
puts "Dropped 'cause too close: #{too_close}"
puts "Dropped 'cause contained het(s): 0"
puts "SNPs kept: #{kept_count}"

list_fh = File.new( "list.snp", "w" )																									#temp
snps_kept = Array.new
snp_list_string	= snps.join("\n")																															#temp
list_fh.puts snp_list_string																																				#temp
kd_rpt.puts kd_list																																									#temp

exit
