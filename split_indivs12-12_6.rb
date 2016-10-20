#!/usr/bin/ruby 
$VERBOSE = true

require 'getoptlong'

# == Synopsis
#
# reduces the SNP data set to the indivs needed, then recreates the hapmap files
# 
# == Usage
# 
#   ./split_indivs.rb -i inlist.txt -s subset.txt -c W
# 
#  -i is the list of input files in hapmap format. 
#	 -s is a single file with the subset of indiv. names you want to keep
#	 -c is a code that goes in each of the 12 output file names
#
# == Author
# Ned Young, Caicedo Lab, U Mass Biology
# 7/16/12
# based on limit_bad_indiv.rb, 4/04/12
#

####################################################################################################

######################
# METHOD DEFINITIONS #
######################

# =get_out_fn
#
# Get a filename for the reconstituted output file.
# Gets passed the fn of the original file.
# Returns fn for output 

def get_out_fn( ori_fn, code ) 
	ori_fn_chunks = ori_fn.split(/\./)
	last = ori_fn_chunks.pop
	next_to_last = ori_fn_chunks.pop
	ori_fn_chunks.push(code)
	ori_fn_chunks.push( next_to_last )
	ori_fn_chunks.push( last )
	new_fn = ori_fn_chunks.join(".")
	return new_fn
end

####################################################################################################

# Main program

# First check that arg(s) are given on the command line
fail "\nUsage: #$0 -i inlist.txt -s subset.txt -c W\n" if ARGV.size < 3

opts = GetoptLong.new(
	[ "-i",				GetoptLong::REQUIRED_ARGUMENT ],
	[ "-s",				GetoptLong::REQUIRED_ARGUMENT ],
	[ "-c",				GetoptLong::REQUIRED_ARGUMENT ]
)
inlist_fn = nil
subset_fn = nil
code = nil
opts.each do |opt, arg|
	case opt
		when '-i'
			inlist_fn = arg
		when '-s'
			subset_fn = arg
		when '-c'
			code = arg
	end
end											 

# Get indiv. name list from file
name_list = Array.new
keepers = Array.new
subset_fh = File.new( subset_fn )
subset_fh.each do |line| 
	next if line =~ /^\s*$/  # checks for blank line
	id = line.chomp  
	# add id to array
 	keepers.push( id )
end
file_counter = 0
inlist_fh = File.new( inlist_fn )
inlist_fh.each do |line| 
	next if line =~ /^\s*$/  # checks for blank line
	# Get hapmap file name from line
	hmp_fn = line.chomp.strip  
	puts "processing file: #{hmp_fn}"  
	file_counter +=1
	
  # input the hapmap file, which needs to be in the same folder as this program
  hmp_fh = File.new( hmp_fn )
  hmp_fh.each do |line| 
  	next if file_counter > 1
  	next if line =~ /^\s*$/  # checks for blank line
  	raw = line.chomp.split(/\s+/)  # space defines field limits
  	# Separate out first eleven columns; just keep the indiv data (names on the first line, SNPs thereafter)
		raw.slice!(0..10)
  	if line =~ /^rs#/  # header line
  	  raw.each do |id|
        portions = id.split(/:/)
        name = portions[0]
        name_list.push(name)  # First name is name_list[0], second is name_list[1], etc.
    	end
    end
  end  # Finished processing a single hapmap file
end # Finished processing all of the hapmap files
inlist_fh.close
outlist_fh = File.new( "out_list.txt", "w" )
# Print out new version of hapmap files
inlist_fh = File.new( inlist_fn )
inlist_fh.each do |line| 
	next if line =~ /^\s*$/  # checks for blank line
	# Get hapmap file name from line
	hmp_fn = line.chomp.strip  

	# Get new fn
	output_fn = get_out_fn( hmp_fn, code )
	outlist_fh.puts output_fn 
	# make new fh 
	out_handle = File.new( output_fn, "w" )
	puts "writing file: #{output_fn}"  
	
  # input the hapmap file, which needs to be in the same folder as this program
  hmp_fh = File.new( hmp_fn )
  hmp_fh.each do |line| 
  	next if line =~ /^\s*$/  # checks for blank line
  	raw = line.chomp.split(/\s+/)  # space defines field limits
  	# Separate out first eleven columns; just keep the indiv data (names on the first line, SNPs thereafter)
  	front_matter = raw.slice!(0..10)
  	front_part = front_matter.join("\t")
  	if line =~ /^rs#/  # header line
  	  back_matter = Array.new
  	  name_counter = 0
  	  raw.each do |id|
				name = name_list[ name_counter ]
				if keepers.include?( name )
        	back_matter.push(id)  # adds the indiv. to the array for later print-out
        end
        name_counter += 1
    	end
    	back_part = back_matter.join("\t")
    	out_handle.puts "#{front_part}\t#{back_part}"
  	else
    	# non-header line
			back_matter = Array.new
  	  name_counter = 0
			raw.each do |id|
				name = name_list[ name_counter ]
				if keepers.include?( name )
					back_matter.push(id)  # adds the indiv. to the array for later print-out
				end
        name_counter += 1
			end
    	back_part = back_matter.join("\t")
    	out_handle.puts "#{front_part}\t#{back_part}"
    end  # Finished a line of data
  end  # Finished processing a single hapmap file
end # Finished processing all of the hapmap files
inlist_fh.close

exit
