# USAGE: ruby runner.rb ../data/peptides/peps_by_rank_product_005_cutoff.csv ../results/transitions_005.csv

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
# REMEMBER TO uncomment the correct dataset (foldername)
foldername = '../data/3H_Ace/'
# foldername = '../data/Endogenous_Ace/'
fieldnames = []
replicate_names = []
line_counter = 0
modification = []
mod_positions = []
mod_positions_str = nil
spectrum_hash = {}
FasterCSV.open(outfile,'w') do |csv|
	csv << %w{ QUERY_NUM REPLICATE PROT_ACCESSION PROT_DESCRIPTION PEPTIDE MODIFICATION MOD_PEPTIDE SPECTRUM_ID CHARGE RET_TIME_MINS HIGHEST_SCORE CUTOFF PARENT_MZ PEP_CALC_MASS PEP_DELTA MS2_MZ MS2_INT y_INT_RANK y_INDEX y_ION_SEQ y-1 y y+2 }
	FasterCSV.foreach(csvfile) do |row|
		line_counter += 1
		if line_counter <= 51
			if fieldnames.empty?
				fieldnames = row
				replicate_names = row[5..16]
			elsif !fieldnames.empty?
				peptide = row[4].to_s
				highest_scored_protein_accno = row[1].to_s
				highest_scored_protein_desc = row[2].to_s
				modification.clear
				if foldername.include? '3H_Ace'
					row[32].scan(/Acetyl:.+\(\d\)\s\((\w)\)/).each do |i| #Acetyl:2H(3) (K); 2 Acetyl:2H(3) (S)
						modification << $1
					end
					#get modification positions from string
					mod_positions_str = row[33].split('')
					mod_positions.clear
					mod_positions_str.each_index do |i|	
						if !mod_positions_str[i].include?('0') && modification.include?(peptide.split('')[i])
							mod_positions << i + 1 #1-based
						end
					end
				elsif foldername.include? 'Endogenous_Ace'
					row[32].scan(/Acetyl\s\((\w)\)/).each do |i| # Acetyl (K)
						modification << $1
						# get modification positions from string
						mod_positions_str = row[33].split('')
						mod_positions.clear
						mod_positions_str.each_index do |i|	
							if !mod_positions_str[i].include?('0') && modification.include?(peptide.split('')[i])
								mod_positions << i + 1 #1-based
							end
						end
					end
				end
				
				query_no = row[0].to_i
				parent_mass = row[28].to_s
				calc_mass = row[29].to_s
				delta = row[30].to_s
				cutoff = row[31].to_s
				scores = row[5..16]
				for i in 0..scores.length-1
				   scores[i] = scores[i].to_f
				end
				highest_score = scores.max
				highest_score_index = scores.index(highest_score)
				repl_with_highest_score = replicate_names[highest_score_index].split("_")[1].to_s
				puts "~~~~~~ #{repl_with_highest_score} > #{query_no} > #{peptide} > #{highest_score} ~~~~~~~~~~"

				# take the ions table from dat file
				filename = foldername + 'dats/' + repl_with_highest_score +  '.dat'
				# puts 'Working for ' + filename
				dat = Mascot::DAT.open(filename, true)
				spectrum_hash = dat.query(query_no)
				title = spectrum_hash['title'.to_sym]
				charge = spectrum_hash['charge'.to_sym]
				rtinseconds = spectrum_hash['rtinseconds'.to_sym]
				ions1 = spectrum_hash['peaks'.to_sym]
				mzs = ions1[0]
				intensities = ions1[1]
# 				puts "Intensities: #{intensities.join(',')} \n\n" #.sort {|x,y| y <=> x }
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
						# puts "yions: #{y[i][0]}, #{y[i][1]},\tm/z: #{mzs[i]}"
						idxset.push(i)
						yset.push(y[i][0])
						mzset.push(mzs[i])
						intset.push(intensities[i])
					end
				end
				if query_no == 151595
					puts '#################################'
					puts y.inspect
					puts '#################################'
					# puts idxset.inspect
					# puts yset.inspect
					# puts mzset.inspect
					# puts intset.inspect
				end
				# puts 'MZs: ' + mzset.join(', ')
				# puts 'INTs: ' + intset.join(', ')
				intset.each do
					maxi = intset.index(intset.max())
					# puts 'MAXI: ' + maxi.to_s + ' > ' + idxset[maxi].to_s
					ranked_idx.push( idxset[maxi] )
					intset[maxi] = 0
				end
				# puts 'Ranked index: ' + ranked_idx.join(', ')

				ranked_idx.each_with_index do |i,ii|
					if( y[i] &&
						!y[i][0].nil? && 1
						y[i][0] > 0 && ii < 3 )
						yidx = y[i][0] - 1
						csv <<  [query_no, 
						repl_with_highest_score,	
						highest_scored_protein_accno,
						highest_scored_protein_desc,
						peptide,
						modification,
						pep.to_s,
						title,
						charge,
						rtinseconds.to_f / 60 ,
						highest_score,
						cutoff,
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
			end
		end
	end
end
