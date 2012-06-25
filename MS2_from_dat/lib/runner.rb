# USAGE: ruby runner.rb ../data/peptides/peps_by_rank_product_005_cutoff.csv ../results/transitions_005.csv

# Conditions of selection to pick up the 3 most significant transitions from MS2 spectrum
# - for the first 50 more significant peptides 
# - pick up the replicate with the highest peptide score (keep the query, modification, protein, parent mass, calculated mass, delta)
# - parse the query section from the .dat file (get the spectrum id, time, charge, ions1 table with mass:intensity peaks/pairs)
# - split masses and intensities in seperate tables
# - feed masses table to assign method 
# - pick up the peaks with the 3 highest intensities

# TODO : Run the 2 other cutoffs with ACE_PEPTIDES routine. Then run MS2_from_dat after getting the missing dat file.


require 'rubygems'
require 'fastercsv'
require 'pep'
require 'mascot/dat'


def find_modification_positions_in_peptide(pep, modification)
    return (0 .. pep.length - 1).find_all { |i| pep[(i-1), 1] == modification }
end

csvfile = ARGV[0]
outfile = ARGV[1]
fieldnames = []
replicate_names = []
line_counter = 0
mod_positions = []
spectrum_hash = {}

FasterCSV.open(outfile,'w') do |csv|
	csv << %w{ QUERY_NUM REPLICATE PROT_ACCESSION PROT_DESCRIPTION PEPTIDE MOD_PEPTIDE SPECTRUM_ID CHARGE RET_TIME_MINS PARENT_MZ PEP_CALC_MASS PEP_DELTA MS2_MZ MS2_INT y_INT_RANK y_INDEX y_ION_SEQ y-1 y y+2 }
	FasterCSV.foreach(csvfile) do |row|
		line_counter += 1
		if line_counter < 50
			if fieldnames.empty?
				fieldnames = row
				replicate_names = row[5..16]
			elsif !fieldnames.empty?
				peptide = row[4].to_s
				highest_scored_protein_accno = row[1].to_s
				highest_scored_protein_desc = row[2].to_s
				modification = row[20].split(":")[0].to_s
				query_no = row[0].to_i
				parent_mass = row[28].to_s
				calc_mass = row[29].to_s
				delta = row[30].to_s
				scores = row[5..16]

				highest_score = scores.max
				highest_score_index = scores.index(highest_score)
				repl_with_highest_score = replicate_names[highest_score_index].split("_")[1].to_s

				mod_positions = find_modification_positions_in_peptide(peptide, modification)

				# take the ions table from dat file
				filename = '../data/dats/' + repl_with_highest_score +  '.dat'
				if !filename.include? 'F001516'
					puts 'Working for ' + filename
					dat = Mascot::DAT.open(filename, true)
					spectrum_hash = dat.query(query_no)
					#FIELDS: peaks, int_max, num_vals, charge, rtinseconds, title, num_used1, mass_min, mass_max, name, index, int_min
					title = spectrum_hash['title'.to_sym]
					# int_max = spectrum_hash['int_max'.to_sym]
					# mass_max = spectrum_hash['mass_max'.to_sym]
					charge = spectrum_hash['charge'.to_sym]
					rtinseconds = spectrum_hash['rtinseconds'.to_sym]
					ions1 = spectrum_hash['peaks'.to_sym]
					mzs = ions1[0]
					intensities = ions1[1]
					pep = Pep.new(peptide, mod_positions)

					y = pep.assign(mzs, pep.yions)
					ranked_idx = []
					# grab legitimate y ions and rank them in intensity desc
					idxset = []
					yset = []
					mzset = []
					intset=[]
					y.each_index do |i|
						if y[i] && y[i][0] && !y[i][1].nil? && y[i][1] > 0
							idxset.push(i)
							yset.push(y[i][0])
							mzset.push(mzs[i])
							intset.push(intensities[i])
						end
					end
					intset.each do
						maxi = intset.index(intset.max())
						ranked_idx.push( idxset[maxi] )
						intset[maxi] = 0
					end

					ranked_idx.each_with_index do |i,ii|
						if(  y[i] &&
							!y[i][0].nil? && 1
							y[i][0] > 0 && ii < 3 )
							yidx = y[i][0] - 1
							csv <<  [query_no, 
							repl_with_highest_score,	
							highest_scored_protein_accno,
							highest_scored_protein_desc,
							peptide,
							pep.to_s,
							title,
							charge,
							rtinseconds.to_f / 60 ,
							parent_mass,
							calc_mass,
							delta,
							mzs[i],
							intensities[i],
							ii, # the mz intensity rank
							y[i][0],
							pep.yions[yidx][0], # y ion sequence
							y[i][1] - H ,
							y[i][1],
							y[i][1] + (H * 2)
							]
						end
					end
				elsif filename.include? 'F001516'
					puts "!!! Needing #{filename} ..."
				end
			end
		end
	end
end
