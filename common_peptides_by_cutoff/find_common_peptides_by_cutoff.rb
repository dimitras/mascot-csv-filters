# USAGE: ruby find_common_peptides_by_cutoff.rb 1 results/common_peptides_no_cutoff.csv results/common_peptides_arrays_no_cutoff.csv
# ruby find_common_peptides_by_cutoff.rb 0.05 results/common_peptides_005_cutoff.csv results/common_peptides_arrays_005_cutoff.csv
# ruby find_common_peptides_by_cutoff.rb 0.5 results/common_peptides_05_cutoff.csv results/common_peptides_arrays_05_cutoff.csv

require 'csv_parser_4protein_hits'
require 'protein_hit'

class RunCSV

	cutoff = ARGV[0]
	infile_list = {}
	# input files
	Dir["data/*.csv"].each do |infile|
		infile_list[infile] = CSVParser4ProteinHits.open(infile, cutoff)
	end

	# output files
	common_peptides_out_file = ARGV[1]
	common_peptides_arrays_out_file = ARGV[2]

	
	# initialize arguments
	common_peptides_out = File.open(common_peptides_out_file, "w")
	common_peptides_arrays_out = File.open(common_peptides_arrays_out_file, "w")
	
	common_peptides = {}
	infile_list.each_value do |csvp|
		csvp.each_peptide do |peptide|
			if !common_peptides.has_key?(peptide)
				common_peptides[peptide] = 1
			else
				common_peptides[peptide] += 1
			end
		end
	end

	common_peptides_in_files = {}
	peptides_score_counts = {}
	infile_list.each_value do |csvp|
		index = csvp.filename.split("/")[1].split("-")[0].to_i - 1
		csvp.each_peptide do |peptide|
			if !common_peptides_in_files.has_key?(peptide)
				common_peptides_in_files[peptide] = infile_list.keys.map{|x| ''}
				peptides_score_counts[peptide] = 0
			end
			common_peptides_in_files[peptide][index] = 1

			if index >= 0 && index <= 6 # interested to count only the 7 first files
				peptides_score_counts[peptide] += 1
			end
		end
	end


	common_peptides_out.puts "prot_hit_num, prot_acc, prot_desc, prot_score, prot_mass, prot_matches, prot_matches_sig, prot_sequences, prot_sequences_sig, pep_query, pep_rank, pep_isbold, pep_isunique, pep_exp_mz, pep_exp_mr, pep_exp_z, pep_calc_mr, pep_delta, pep_miss, pep_score, pep_expect, pep_res_before, pep_seq, pep_res_after, pep_var_mod, pep_var_mod_pos, pep_scan_title, filename, files_found_pep"

	common_peptides.sort { |a,b| b[1] <=> a[1] }.each do |peptide,count|
		infile_list.each_value do |csvp|
			if csvp.has_peptide(peptide)
				csvp.protein_hits(peptide).each do |hit|
					common_peptides_out.puts hit.to_csv + ",\"" + common_peptides[peptide].to_s + "\""
				end
			end
		end
	end
	common_peptides_out.close


	common_peptides_arrays_out.puts "PROT_ACC , PROT_DESC , PEPTIDE , R1 , R2 , R3 , R4 , R5 , R6 , R7 , R8 , R9 , R10 , SCORE"
	peptides_score_counts.sort { |a,b| b[1] <=> a[1] }.each do |peptide,score|
		infile_list.each_value do |csvp|
			if csvp.has_peptide(peptide)
				csvp.protein_hits(peptide).each do |hit|
					common_peptides_arrays_out.puts '"' + hit.prot_acc.to_s + '","' + hit.prot_desc.to_s + '","' + peptide.to_s + '","' + common_peptides_in_files[peptide].join('","') + '","' + peptides_score_counts[peptide].to_s + '"'
				end
			end
		end
	end
	common_peptides_arrays_out.close
end







			