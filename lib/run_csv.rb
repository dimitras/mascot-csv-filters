# USAGE: ruby run_csv.rb ../data/xuanwen_filters/F001504_regular_ASA_trypsin_regular_KST____1.csv ../data/xuanwen_filters/F001493_regular_ASA_Trypsin_3HKST_____4.csv ../data/xuanwen_filters/F001517_1st_3HASA_T_regular_KST_____2.csv ../data/xuanwen_filters/F001516_1st_3HASA_T_3HKST_____5.csv ../data/xuanwen_filters/common_peptides.csv

require 'csv_parser'
require 'protein_hit'

class RunCSV
	
	count = 0
	expA_fileA = ARGV[0] #1
	expA_fileB = ARGV[1] #4
	expB_fileA = ARGV[2] #2
	expB_fileB = ARGV[3] #5
	common_peptides_out_file = ARGV[4]
	
	expA_fileA_csvp = CSVParser.open(expA_fileA)
	expA_fileB_csvp = CSVParser.open(expA_fileB)
	expB_fileA_csvp = CSVParser.open(expB_fileA)
	expB_fileB_csvp = CSVParser.open(expB_fileB)
	common_peptides_out_fp = File.open(common_peptides_out_file, "w")
	
	common_peptides = {}
	expA_fileA_csvp.each_peptide do |peptide|
		if expB_fileA_csvp.has_peptide(peptide) || expB_fileB_csvp.has_peptide(peptide)
			common_peptides[peptide] = nil
		end
	end
	
	expA_fileB_csvp.each_peptide do |peptide|
		if expB_fileA_csvp.has_peptide(peptide) || expB_fileB_csvp.has_peptide(peptide)
			common_peptides[peptide] = nil
		end
	end
	
	common_peptides_out_fp.puts "prot_hit_num, prot_acc, prot_desc, prot_score, prot_mass, prot_matches, prot_matches_sig, prot_sequences, prot_sequences_sig, pep_query, pep_rank, pep_isbold, pep_isunique, pep_exp_mz, pep_exp_mr, pep_exp_z, pep_calc_mr, pep_delta, pep_miss, pep_score, pep_expect, pep_res_before, pep_seq, pep_res_after, pep_var_mod, pep_var_mod_pos, pep_scan_title, filename"
	
	common_peptides.each_key do |peptide|
		if expA_fileA_csvp.has_peptide(peptide)
			expA_fileA_csvp.protein_hits(peptide).each do |hit|
				common_peptides_out_fp.puts hit.to_csv
			end
		end
		if expA_fileB_csvp.has_peptide(peptide)
			expA_fileB_csvp.protein_hits(peptide).each do |hit|
				common_peptides_out_fp.puts hit.to_csv
			end
		end
		if expB_fileA_csvp.has_peptide(peptide)
			expB_fileA_csvp.protein_hits(peptide).each do |hit|
				common_peptides_out_fp.puts hit.to_csv
			end
		end
		if expB_fileB_csvp.has_peptide(peptide)
			expB_fileB_csvp.protein_hits(peptide).each do |hit|
				common_peptides_out_fp.puts hit.to_csv
			end
		end
		count = count + 1
	end
	puts count
end
			