#!/usr/bin/ruby 
$VERBOSE = true

require 'getoptlong'

# == Synopsis
#
# limits indivs according to how common Ns are, then recreates the hapmap files
# 
# == Usage
# 
#   ./limit_bad_indiv.rb -i inlist.txt -c 0.95
# 
#  -c is the cutoff. 0.95 means that any indiv. with 95% or more Ns gets tossed. Default is 0.95.
#
# == Author
# Ned Young, Caicedo Lab, U Mass Biology
# 4/04/12
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

def get_out_fn( ori_fn ) 
	ori_fn_chunks = ori_fn.split(/\./)
	last = ori_fn_chunks.pop
	next_to_last = ori_fn_chunks.pop
	ori_fn_chunks.push("lbo")
	ori_fn_chunks.push( next_to_last )
	ori_fn_chunks.push( last )
	new_fn = ori_fn_chunks.join(".")
	return new_fn
end

####################################################################################################

# Main program

# First check that arg(s) are given on the command line
fail "\nUsage: #$0 -i inlist.txt\n" if ARGV.size < 1

opts = GetoptLong.new(
	[ "-i",				GetoptLong::REQUIRED_ARGUMENT ],
	[ "-c",				GetoptLong::REQUIRED_ARGUMENT ]
)
inlist_fn = nil
cutoff = 0.95
opts.each do |opt, arg|
	case opt
		when '-i'
			inlist_fn = arg
		when '-c'
			cutoff = arg.to_f
	end
end											 

# Start data structure (hash) to store every genotype in
name_list = Array.new
big_struct = Hash.new
out_fn_list = Hash.new # Key: original name, Value: name of altered file.
num_snps = 0
# get a fileh
inlist_fh = File.new( inlist_fn )
inlist_fh.each do |line| 
	next if line =~ /^\s*$/  # checks for blank line
	# Get hapmap file name from line
	hmp_fn = line.chomp.strip  
	puts "processing file: #{hmp_fn}"  
	# Get new fn
	output_fn = get_out_fn( hmp_fn )
	out_fn_list[ hmp_fn ] = output_fn 

  # input the hapmap file, which needs to be in the same folder as this program
  hmp_fh = File.new( hmp_fn )
  hmp_fh.each do |line| 
  	next if line =~ /^\s*$/  # checks for blank line
  	raw = line.chomp.split(/\s+/)  # space defines field limits
  	# Separate out first eleven columns; just keep the indiv data (names on the first line, SNPs thereafter)
		raw.slice!(0..10)
  	if line =~ /^rs#/  # skips header line
  	  raw.each do |id|
        name_list.push(id)  # First name is name_list[0], second is name_list[1], etc.
    	end
  	else
    	# Get genotype from line
    	counter = 0
    	raw.each do |pos|
    	  id = name_list[ counter ]
    	  big_struct[ id ] = 0 unless big_struct[ id ]
    	  big_struct[ id ] += 1 if pos == "N"
    	  counter += 1
    	end # Finished the row
    	num_snps += 1
    end
  end  # Finished processing a single hapmap file
end # Finished processing all of the hapmap files
inlist_fh.close

tossers = Array.new
big_struct.each_pair do |id, ens|
  # Evaluate frequency and plan to drop, if too many Ns
  freq = ens.to_f/num_snps.to_f
  if freq >= cutoff
    # add id to array
    tossers.push( id )
  end  
end

# Print out new version of hapmap files
inlist_fh = File.new( inlist_fn )
inlist_fh.each do |line| 
	next if line =~ /^\s*$/  # checks for blank line
	# Get hapmap file name from line
	hmp_fn = line.chomp.strip  
	output_fn = out_fn_list[ hmp_fn ]
	# make new fh 
	out_handle = File.new( output_fn, "w" )
	toss_fn = output_fn + "_t"
	toss_fh = File.new( toss_fn, "w" )
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
  	  toss_matter = Array.new
  	  name_counter = 0
  	  raw.each do |id|
				name = name_list[ name_counter ]
				unless tossers.include?( name )
        	back_matter.push(id)  # adds the indiv. to the array for later print-out
        else
          toss_matter.push(id)  # adds the indiv. to the toss array for later print-out
        end
        name_counter += 1
    	end
    	
    	back_part = back_matter.join("\t")
    	out_handle.puts "#{front_part}\t#{back_part}"
    	toss_pt_2 = toss_matter.join("\t")
    	toss_fh.puts "#{front_part}\t#{toss_pt_2}"
  	else
    	# non-header line
			back_matter = Array.new
  	  toss_matter = Array.new
  	  name_counter = 0
			raw.each do |id|
				name = name_list[ name_counter ]
				unless tossers.include?( name )
					back_matter.push(id)  # adds the indiv. to the array for later print-out
				else
          toss_matter.push(id)  # adds the indiv. to the toss array for later print-out					
				end
        name_counter += 1				
			end

    	back_part = back_matter.join("\t")
    	out_handle.puts "#{front_part}\t#{back_part}"
    	toss_pt_2 = toss_matter.join("\t")
    	toss_fh.puts "#{front_part}\t#{toss_pt_2}"
    end  # Finished a line of data
  end  # Finished processing a single hapmap file
end # Finished processing all of the hapmap files
inlist_fh.close

exit
