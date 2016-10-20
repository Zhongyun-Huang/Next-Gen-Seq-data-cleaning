#!/usr/bin/ruby 
$VERBOSE = true

require 'getoptlong'

# == Synopsis
#
# takes a hapmap file and converts the indiv. ids from long form to short
# 
# == Usage
# 
#   ./convert_id.rb -i input.txt -o outfile.txt
# 
#  -i is the input file in hapmap format. 
#	 -o is a name for a file with output data in hapmap format.
#		
#
# == Author
# Ned Young, Caicedo Lab, U Mass Biology
# 08/06/12
# based on split_indiv.rb, 7/16/12
# based on limit_bad_indiv.rb, 4/04/12
#

####################################################################################################

# Main program

# First check that arg(s) are given on the command line
fail "\nUsage: #$0 -i input.txt -o outfile.txt\n" if ARGV.size < 2

opts = GetoptLong.new(
	[ "-i",				GetoptLong::REQUIRED_ARGUMENT ],
	[ "-o",				GetoptLong::REQUIRED_ARGUMENT ]
)
input_fn = nil
output_fn = nil
opts.each do |opt, arg|
	case opt
		when '-i'
			input_fn = arg
		when '-o'
			output_fn = arg
	end
end											 
# make new fh 
out_handle = File.new( output_fn, "w" )

name_list = Array.new
puts "processing file: #{input_fn}"  
input_fh = File.new( input_fn )
input_fh.each do |line| 
	next if line =~ /^\s*$/  # checks for blank line
	raw = line.chomp.split(/\s+/)  # space defines field limits
	# Separate out first eleven columns; just keep the indiv data (names on the first line, SNPs thereafter)
	front_matter = raw.slice!(0..10)
	front_part = front_matter.join("\t")
	if line =~ /^rs#/  # header line
		raw.each do |id|
			portions = id.split(/:/)
      short_id = portions[0]
			name = short_id
			name_list.push(name)  # First name is name_list[0], second is name_list[1], etc.
		end
		back_part = name_list.join("\t")
		out_handle.puts "#{front_part}\t#{back_part}"
	else
		# non-header line
		back_part = raw.join("\t")
		out_handle.puts "#{front_part}\t#{back_part}"
	end  # Finished a line of data
end  # Finished processing the hapmap file

exit
