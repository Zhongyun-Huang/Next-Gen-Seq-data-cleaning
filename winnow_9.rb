#!/usr/bin/ruby 
$VERBOSE = true

require 'getoptlong'

# == Synopsis
#
# counts Ns in hapmap file to calculate frequency of Ns in a given SNP,
# keeping those SNPs that pass a data quality threshhold.
# Also looks up the x bp before and after the SNP. If either is mononuc, the SNP is removed from the 
# list. Reports the number of such SNPs removed.  Makes a new hapmap file, but includes SNPs from 
# all chromosomes.

# 
# == Usage
# 
#   ./winnow.rb -i infile -o winnowed_snp_060512.rpt -q 0.1 -m 5
#
# -q is the quality cutoff.  SNps with a freq. of Ns >q are dropped. Default is 0.1.
# 
# -m is the maximum allowable mononuc. run length; runs of a single nuc. (or a single nuc. and Ns)
# greater than m in length cannot contribute their SNPs to the SNP set. Default is 5.
# 
# == Author
# Ned Young, Caicedo Lab, U Mass Biology
# 6/8/12
# based on filt_snp.rb, 5/16/12
# based on assess_hmp_brief.rb, 4/04/12
#

####################################################################################################

######################
# METHOD DEFINITIONS #
######################

# =get_chrom_as_string
#
# Extract sequence from fasta file.  Unite the lines into string.
# Argument is a fn
# Returns string

def get_chrom_as_string( seq_fn )
	contigs = Array.new
	pieces = Array.new
	IO.foreach(seq_fn) do |line|
		unless line =~ /^>/
			pieces.push(line.chomp)
		end
	end
	whole_seq = pieces.join
	return whole_seq
end

# =has_long_run
#
# Examines sequence to see if it is inside a long mononuc. run or a mix of Ns and a single nuc.
# Arguments are a string consisting of the SNPped base +/- n bp; snp_name; n.
# Returns true or false and a msg for trouble-shooting

def has_long_run( frag, snp_name, n )
	has_run = false
	first_time = true
	snp_pos = n + 1 # Zero-based numbering
	snp_base = frag[ snp_pos, 1 ]
	first_time = true
	run_len = 1
	alt_early_run_ct = 0
	alt_late_run_ct = 0
	alt_early_run_base = ""
	alt_late_run_base = ""
	run_start_found = false
	run_end_found = false	
	early_run_start_found = false
	late_run_start_found = false
	msg = Array.new
	for i in 1..( n + 1 )
		early_pos = snp_pos - i
		late_pos = snp_pos + i 
		early = frag[ early_pos, 1 ]
		late = frag[ late_pos, 1 ]
		if early == snp_base
			if run_start_found == false
				run_len += 1
			end
		else 
			run_start_found = true
			if first_time == true
				alt_early_run_base = early
				alt_early_run_ct += 1
			else
				if ( early == alt_early_run_base ) && early_run_start_found == false
					alt_early_run_ct += 1
				else
					early_run_start_found = true
				end
			end
		end
		if late == snp_base 
			if run_end_found == false
				run_len += 1
			end
		else
			run_end_found = true
			if first_time == true
				alt_late_run_base = late
				alt_late_run_ct += 1
			else
				if ( late == alt_late_run_base ) && late_run_start_found == false
					alt_late_run_ct += 1
				else
					late_run_start_found = true
				end
			end
		end
		first_time = false
	end
	if run_len > n
		has_run = true 
	elsif alt_early_run_ct > n
		has_run = true 
	elsif alt_late_run_ct > n
		has_run = true 
	end
	return has_run, msg
end

####################################################################################################

# Main program

# First check that arg(s) are given on the command line
fail "\nUsage: #$0 -i new_inlist.txt -o wintest_snp_0605data_0628.rpt -q 0.1 -m 5\n" if ARGV.size < 2

opts = GetoptLong.new(
	[ "-i",				GetoptLong::REQUIRED_ARGUMENT ],
	[ "-o", 			GetoptLong::REQUIRED_ARGUMENT ],
	[ "-q", 			GetoptLong::OPTIONAL_ARGUMENT ],
	[ "-m", 			GetoptLong::OPTIONAL_ARGUMENT ]
)
inlist_fn = nil
output_fn = nil
q_cutoff = 0.1
m_cutoff = 5
opts.each do |opt, arg|
	case opt
		when '-i'
			inlist_fn = arg
		when '-o'
			output_fn = arg
		when '-q'
			q_cutoff = arg.to_f
		when '-m'
			m_cutoff = arg.to_i
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
q_drops = 0
mono_drops = 0
lines_examined = 0
drops = Array.new																																										#temp
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
		raw.slice!(0..10)
  	if line =~ /^rs#/  # header line
  		out_handle.puts line if first_chrom == true
  	else # not a header line
			n = 0
			raw.each do |pos|
				n += 1 if pos == "N"
			end
			# Evaluate frequency 
			freq = n.to_f/raw.size.to_f
			# Print out if < "cutoff" Ns- THIS IS WHERE THE QUALITY FILTERING HAPPENS
			if freq < q_cutoff
				kept[ c_num ].push( line ) 
			else
				q_drops += 1
			end
		end # Finished analyzing info from a row
		lines_examined += 1
  end  # Finished processing a single hapmap file
  first_chrom = false
end # Finished processing all of the hapmap files
report.puts "Dropped because of quality cutoff: #{q_drops}"
[ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 ].each do |c_num|
	zero = ""
	zero = "0" if c_num < 10
	seq_fn = "../MSU6.0/chr" + zero + c_num.to_s + ".con"
	# Open chrom. seq. file; bring the seq. in.
	c_string = get_chrom_as_string( seq_fn )
	
	kept_set = kept[ c_num ]
	kept_set.each do | line |
		info = line.chomp.split(/\s+/)  # space defines field limits
		snp_pos = info[ 3 ].to_i
		snp_name = info[ 0 ]
		if snp_pos < m_cutoff
			report.puts "#{snp_name} dropped, less than m_cutoff."
		else
			first = snp_pos - ( m_cutoff + 2 ) 
			final = snp_pos + ( m_cutoff )
			frag = c_string.slice(first..final)
			before = frag.slice(0..(m_cutoff))
			herself = frag.slice((m_cutoff + 1)..(m_cutoff + 1))   
			after = frag.slice((m_cutoff + 2)..-1)
			items = has_long_run( frag, snp_name, m_cutoff )
			report.puts items[ 1 ]
			if items[ 0 ] == false 
				out_handle.puts line
			else
				mono_drops += 1
				drops.push( "#{snp_name} dropped- before: #{before}, herself: #{herself}, after: #{after}.\n" )																						 												#temp
			end
		end
	end
end
report.puts "lines_examined: #{lines_examined}"																						
report.puts "dropped because of mononuc. repeats: #{mono_drops}"
report.puts "#{drops}"																									 												#temp

exit
