# USAGE: ruby runner.rb peps_by_rank_product_005_cutoff.csv transitions_005.csv ionstables_005.csv
# ruby runner.rb peps_by_rank_product_005_cutoff.csv transitions_005.csv ionstables_005.csv

require 'rubygems'
require 'fastercsv'
require 'pep'
require 'mascot/dat'
require 'gnuplot'

csvfile = ARGV[0]
outfile = ARGV[1]
ionstables_file = ARGV[2]
# REMEMBER TO uncomment the correct dataset (foldername)
# foldername = '../data/3H_Ace/'
foldername = '../data/Endogenous_Ace/'
# resultsfolder = '../results/3H_Ace/'
resultsfolder = '../results/Endogenous_Ace/'
fieldnames = []
replicate_names = []
line_counter = 0
modification = []
mod_positions = []
mod_positions_str = nil
spectrum = {}
itf = File.open(resultsfolder + 'ionstables/' + ionstables_file,'w')

# plot the spectra
def gplot(peptide, title, assigned_ions, mzs, intensities, figure_filename)
	Gnuplot.open do |gp|
		Gnuplot::Plot.new( gp ) do |plot|
			plot.output figure_filename
			plot.terminal 'svg'
			plot.title  "Query title:#{title} of #{peptide}"
			plot.ylabel 'intensity'
			plot.xlabel 'm/z'
			x_vals = assigned_ions.collect{|mass| mass.first}
			y_vals = assigned_ions.collect{|arr| arr[1]}
			l_vals = assigned_ions.collect{|idx| idx[2]}
			
			plot.data << Gnuplot::DataSet.new( [mzs, intensities] ) do |ds|
				ds.with = 'impulses linecolor rgb "blue"'
				ds.linewidth = 1
				ds.notitle
			end
			
			plot.data << Gnuplot::DataSet.new( [x_vals, y_vals] ) do |ds|
				ds.with = 'impulses linecolor rgb "red"'
				ds.linewidth = 1.5
				ds.notitle
			end

			max_y_val = y_vals.max
			label_y_vals = y_vals.map{|value| value + max_y_val*0.05}
			plot.data << Gnuplot::DataSet.new( [x_vals, label_y_vals, l_vals] ) do |ds|
				ds.with = 'labels textcolor lt 1 rotate left'
				ds.notitle
			end
		end
	end
end

FasterCSV.open(resultsfolder + outfile,'w') do |csv|
	csv << %w{ QUERY_NUM REPLICATE PROT_ACCESSION PROT_DESCRIPTION PEPTIDE MODIFICATION MOD_PEPTIDE SPECTRUM_ID CHARGE RET_TIME_MINS HIGHEST_SCORE CUTOFF PARENT_MZ PEP_CALC_MASS PEP_DELTA MS2_MZ MS2_INT y_INT_RANK y_INDEX y_ION_SEQ y-1 y y+2 }
	FasterCSV.foreach(foldername + csvfile) do |row|
		line_counter += 1
		assigned_ions = Array.new()
		if true || line_counter <= 51 # for our purposes we need the first 50 peptides in the list
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
					end
					# get modification positions from string
					mod_positions_str = row[33].split('')
					mod_positions.clear
					mod_positions_str.each_index do |i|
						if !mod_positions_str[i].include?('0') && modification.include?(peptide.split('')[i])
							mod_positions << i + 1 #1-based
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
				puts "#{repl_with_highest_score} > #{query_no} > #{peptide} > #{highest_score}"

				# take the ions table from dat file
				filename = foldername + 'dats/' + repl_with_highest_score +  '.dat'
				dat = Mascot::DAT.open(filename, true)
				spectrum = dat.query(query_no)
				title = spectrum.title
				charge = spectrum.charge
				rtinseconds = spectrum.rtinseconds
				ions1 = spectrum.peaks
				mzs = spectrum.mz
				intensities = spectrum.intensity

				pep = Pep.new(peptide, mzs, intensities, mod_positions)
				assigned_yions = pep.assigned_yions
				
				# Pick up the 3 most significant transitions from MS2 spectrum for each peptide & print transitions file
				ranked_idx = pep.ranked_yions_intensities_idx
				ranked_idx.each_with_index do |i,ii|
					if( assigned_yions[i] && !assigned_yions[i][0].nil? && 1
						assigned_yions[i][0] > 0 && ii < 3)
						yidx = assigned_yions[i][0] - 1
						csv << [
							query_no,
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
							assigned_yions[i][0],
							pep.yions[yidx][0], # y ion sequence
							assigned_yions[i][1] - H ,
							assigned_yions[i][1],
							assigned_yions[i][1] + (H * 2)
						]
					end
				end

# 				# print the ionstable file
# 				itf.puts "#{repl_with_highest_score}, #{query_no}, #{peptide}"
# 				itf.puts(%w{ mass intensity y-index }.join(","))
# 				pep.print_all_yionstable_to_filehandle(itf)
# 
# 				# plot the spectra
# 				assigned_yionstable = pep.assigned_yionstable
# 				figure_filename = "#{resultsfolder}/figures/figure_#{repl_with_highest_score}_#{query_no}_#{peptide}.svg"
# 				gplot(peptide, title, assigned_yionstable, mzs, intensities, figure_filename)
			end
		end #ends if for counting the first 50 entries
	end
end

