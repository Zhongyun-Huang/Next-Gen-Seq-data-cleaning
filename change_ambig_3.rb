#!/usr/bin/ruby 
$VERBOSE = true

require 'getoptlong'

# == Synopsis
#
# changes ambiguity codes to Ns. Reports the number of such SNPs changed.  Makes a new 
# hapmap file, but includes SNPs from all chromosomes.

# 
# == Usage
# 
#   ./change_ambig.rb -i inlist.txt -o ZJW_ambigs_chgd
#
# 
# == Author
# Ned Young, Caicedo Lab, U Mass Biology
# 6/8/12
# based on filt_snp.rb, 5/16/12
# based on assess_hmp_brief.rb, 4/04/12
#

####################################################################################################

# Main program

# First check that arg(s) are given on the command line
fail "\nUsage: #$0 -i inlist.txt -o ZJW_ambigs_chgd\n" if ARGV.size < 2

opts = GetoptLong.new(
	[ "-i",				GetoptLong::REQUIRED_ARGUMENT ],
	[ "-o", 			GetoptLong::REQUIRED_ARGUMENT ]
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

out_handle = File.new( output_fn, "w" )
report_fn = output_fn + "drops"
report = File.new( report_fn, "w" )
kept = {  1 => 	Array.new,
					2 => 	Array.new,
					3 => 	Array.new,
					4 => 	Array.new,
					5 => 	Array.new,
					6 => 	Array.new,
					7 => 	Array.new,
					8 => 	Array.new,
					9 => 	Array.new,
					10 => Array.new,
					11 => Array.new,
					12 => Array.new}

inlist_fh = File.new( inlist_fn )
lines_examined = 0
ambig_snps_changed = 0
a_all = 0
first_chrom = true
inlist_fh.each do |line| 
	next if line =~ /^\s*$/  # checks for blank line
	# Get hapmap file name from line
	hmp_fn = line.chomp.strip  
	puts "processing file: #{hmp_fn}"                                                                    
  # input the hapmap file, which needs to be in the same folder as this program
  hmp_fh = File.new( hmp_fn )
  hmp_fh.each do |line| 
  	next if line =~ /^\s*$/ || line == nil  # checks for blank line
  	# Get position from line
		raw = line.chomp.split(/\s+/)  # space defines field limits
		snp_name = raw[ 0 ]																															#temp
		c_num = raw[ 2 ].to_i
		# Remove first eleven columns; just keep the SNP data
		front_matter = raw.slice!(0..10)
  	if line =~ /^rs#/  # header line
  		out_handle.puts line if first_chrom == true
  	else # not a header line
			n = 0
			a = 0
			raw_pos_num = 0
			raw.each do |pos|
				n += 1 if pos == "N"
				unless ['A', 'C', 'G', 'T', 'N'].include?(pos)
					a += 1 
					raw[ raw_pos_num ] = "N"
				end
				raw_pos_num += 1
			end
			if a > 0
				ambig_snps_changed += 1
				a_all += a
			end
			# Evaluate frequency 
			freq = n.to_f/raw.size.to_f
			column_set = front_matter + raw
			reconst_line = column_set.join( "\t" )
			kept[ c_num ].push( reconst_line ) 
		end # Finished analyzing info from a row
		lines_examined += 1
  end  # Finished processing a single hapmap file
  first_chrom = false
end # Finished processing all of the hapmap files
[ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 ].each do |c_num|
	kept_set = kept[ c_num ]
	kept_set.each do | line |
		out_handle.puts line
	end
end
report.puts "lines_examined: #{lines_examined}"																		
report.puts "SNPs changed because they contained ambiguity code(s): #{ambig_snps_changed}"																						
report.puts "Ambiguity codes: #{a_all}"																						

exit
