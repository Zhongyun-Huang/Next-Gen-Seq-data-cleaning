#!/usr/bin/ruby 
$VERBOSE = true

require 'getoptlong'

# == Synopsis
#
# counts number of distinct (ACGorT) bases observed in a given SNP (in a hapmap file),
# keeping those SNPs that have > 1.
# Makes a new hapmap file.
# 
# == Usage
# 
#   ./drop_mono.rb -i infile.hmp -o outfile.hmp 
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
fail "\nUsage: #$0 -i infile.hmp -o outfile.hmp\n" if ARGV.size < 2

opts = GetoptLong.new(
	[ "-i",				GetoptLong::REQUIRED_ARGUMENT ],
	[ "-o", 			GetoptLong::REQUIRED_ARGUMENT ]
)
infile_fn = nil
output_fn = nil
opts.each do |opt, arg|
	case opt
		when '-i'
			infile_fn = arg
		when '-o'
			output_fn = arg
	end
end												 

out_handle = File.new( output_fn, "w" )
o_parts = output_fn.split(/\./)
pref = o_parts[0] + "." + o_parts[1]
rpt_fn = pref + ".dm_report"
report = File.new( rpt_fn, "w" )

infile_fh = File.new( infile_fn )

mono_drops = 0
lines_examined = 0

infile_fh.each do |line| 
	next if line =~ /^\s*$/  # checks for blank line
	# Get position from line
	raw = line.chomp.split(/\s+/)  # space defines field limits
	# Remove first eleven columns; just keep the SNP data
	raw.slice!(0..10)
	if line =~ /^rs#/  # header line
		out_handle.puts line
	else # not a header line
		counts = { "A" => 	0,
					 "C" => 	0,
					 "G" => 	0,
					 "T" => 	0,
					 "R" => 	0,
					 "Y" => 	0,
					 "S" => 	0,
					 "W" => 	0,
					 "K" => 	0,
					 "M" => 	0,
					 "B" => 	0,
					 "D" => 	0,
					 "H" => 	0,
					 "V" => 	0,
					 "N" => 	0}
		# Drop SNP unless 2 dif bases (ACGorT) are found
		raw.each do |pos|
			counts[ pos ] += 1
		end
		# Evaluate number of unique codes found 
		uniq_codes = 0
		counts.each_pair do |code, count|
			if ['A', 'C', 'G', 'T'].include?(code)
				if count > 0
					uniq_codes += 1
				end
			end
		end
		# Print out if > 1
		if uniq_codes > 1
			out_handle.puts line 
		else
			mono_drops += 1
		end
	end # Finished analyzing info from a row
	lines_examined += 1
end  # Finished processing hapmap file

report.puts "File examined:\n"
report.puts "#{infile_fn}"
report.puts "#{Time.now}\n"
report.puts "Program version:#{$0}\n"

report.puts "Dropped because monomorphic: #{mono_drops}"
report.puts "lines_examined: #{lines_examined - 1}"																						

exit
