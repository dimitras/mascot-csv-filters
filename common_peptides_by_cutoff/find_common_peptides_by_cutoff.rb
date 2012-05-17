# USAGE
# ruby find_common_peptides_by_cutoff.rb 10 raw/uniprot_human_db_04182012.fasta results/joined_peps_no_cutoff.csv results/peps_per_existence_no_cutoff.csv results/peps_by_rank_product_no_cutoff.csv
# ruby find_common_peptides_by_cutoff.rb 0.05 raw/uniprot_human_db_04182012.fasta results/joined_peps_005_cutoff.csv results/peps_per_existence_005_cutoff.csv results/peps_by_rank_product_005_cutoff.csv
# ruby find_common_peptides_by_cutoff.rb 0.5 raw/uniprot_human_db_04182012.fasta results/joined_peps_05_cutoff.csv results/peps_per_existence_05_cutoff.csv results/peps_by_rank_product_05_cutoff.csv

require 'csv_parser_for_protein_hits'
require 'protein_hit'
require 'fasta_parser'
require 'fasta_entry'

cutoff = ARGV[0]
proteome_db_fasta_file = ARGV[1]
infile_list = {}
# input files
Dir["data/*.csv"].each do |infile|
	file_num = infile.split("/")[1].split("-")[0].to_i - 1
	if file_num >= 0 && file_num <= 6
		infile_list[infile] = CSVParserForProteinHits.open(infile, cutoff, false, file_num)
	else
		infile_list[infile] = CSVParserForProteinHits.open(infile, cutoff, true, file_num)
	end
end

# output files
joined_replicates_per_pep_out_file = ARGV[2]
pep_existance_in_replicates_out_file = ARGV[3]
peps_sorted_by_rank_product_out_file = ARGV[4]


# initialize arguments
@proteome_db_fap = FastaParser.open(proteome_db_fasta_file)
joined_replicates_per_pep_out = File.open(joined_replicates_per_pep_out_file, "w")
pep_existance_in_replicates_out = File.open(pep_existance_in_replicates_out_file, "w")
peps_sorted_by_rank_product_out = File.open(peps_sorted_by_rank_product_out_file, "w")


# replace the uncharacterized proteins with the newest description from uniprot database
def uncharacterized_to_annotated(protein_hit)
	if protein_hit.prot_desc.include? "Uncharacterized"
		fasta_entry = @proteome_db_fap.entry_by_id(protein_hit.prot_acc)
		if fasta_entry != nil
			protein_hit.prot_desc = fasta_entry.desc
		end
	end
end

# 	infile_list.each_value do |csvp|
# 		csvp.each_hit do |protein_hit|
# 			if protein_hit.prot_desc.include? "Uncharacterized"
# 				fasta_entry = proteome_db_fap.entry_by_id(protein_hit.prot_acc)
# 				if fasta_entry != nil
# 					protein_hit.prot_desc = fasta_entry.desc
# 				end
# 			end
# 		end
# 	end

# count total peptide hits in all replicates
total_peptide_hits = {}
infile_list.each_value do |csvp|
	csvp.each_peptide do |peptide|
		if !total_peptide_hits.has_key?(peptide)
			total_peptide_hits[peptide] = 1
		else
			total_peptide_hits[peptide] += 1
		end
	end
end


# locate peptide existance in all replicates & score adding 1 foreach existance for 7 first reps
pep_existance_in_replicates = {}
pep_existance_counts = {}
infile_list.each_value do |csvp|
	csvp.each_peptide do |peptide|
		if !pep_existance_in_replicates.has_key?(peptide)
			pep_existance_in_replicates[peptide] = infile_list.keys.map{|x| ''}
			pep_existance_counts[peptide] = 0
		end
		pep_existance_in_replicates[peptide][csvp.id] = 1

		if csvp.is_descending
			pep_existance_counts[peptide] += 1
		end
	end
end


# get the rank product for the highest scored hit of each peptide, in all replicates
pep_scores = {}
pep_rank_product = Hash.new { |h,k| h[k] = 1 }
total_peptide_hits.each_key do |peptide|
	infile_list.each_value do |csvp|
		if !pep_scores.has_key?(peptide)
			pep_scores[peptide] = infile_list.keys.map{|x| ''}
		end
		if !csvp.has_peptide(peptide)
			pep_scores[peptide][csvp.id] = 0
			pep_rank_product[peptide] *= csvp.worst_rank
		else
			highest_scored_hit = csvp.highest_scored_hit_for_pep(peptide)
			pep_scores[peptide][csvp.id] = highest_scored_hit.pep_score
			pep_rank_product[peptide] *= highest_scored_hit.rank
		end
	end
end


# create the joined file with all peptide hits from all replicates
joined_replicates_per_pep_out.puts "prot_hit_num, prot_acc, prot_desc, prot_score, prot_mass, prot_matches, prot_matches_sig, prot_sequences, prot_sequences_sig, pep_query, pep_rank, pep_isbold, pep_isunique, pep_exp_mz, pep_exp_mr, pep_exp_z, pep_calc_mr, pep_delta, pep_miss, pep_score, pep_expect, pep_res_before, pep_seq, pep_res_after, pep_var_mod, pep_var_mod_pos, pep_scan_title, filename, files_found_pep"

total_peptide_hits.sort { |a,b| b[1] <=> a[1] }.each do |peptide,count|
	infile_list.each_value do |csvp|
		if csvp.has_peptide(peptide)
			csvp.protein_hits(peptide).each do |hit|
				joined_replicates_per_pep_out.puts hit.to_csv + ",\"" + total_peptide_hits[peptide].to_s + "\""
			end
		end
	end
end
joined_replicates_per_pep_out.close


# create array with the existance of each peptide in all replicates & add a score adding 1 foreach existance for 7 first reps
pep_existance_in_replicates_out.puts "PROT_ACC , PROT_DESC , PEPTIDE , R1 , R2 , R3 , R4 , R5 , R6 , R7 , R8 , R9 , R10 , SCORE"
pep_existance_counts.sort { |a,b| b[1] <=> a[1] }.each do |peptide,score|
	infile_list.each_value do |csvp|
		if csvp.has_peptide(peptide)
			csvp.protein_hits(peptide).each do |hit|
				pep_existance_in_replicates_out.puts '"' + hit.prot_acc.to_s + '","' + hit.prot_desc.to_s + '","' + peptide.to_s + '","' + pep_existance_in_replicates[peptide].join('","') + '","' + pep_existance_counts[peptide].to_s + '"'
			end
		end
	end
end
pep_existance_in_replicates_out.close


# create array with the ranks of the highest scored hit of each peptide in all replicates & add a score by calculating the product rank for each peptide, firstly sorting by pep_score (descending for the first 7 reps & ascnding for the 3 last ones)
peps_sorted_by_rank_product_out.puts "PROT_ACC , PROT_DESC , PEPTIDE , R1 , R2 , R3 , R4 , R5 , R6 , R7 , R8 , R9 , R10 , RP"
pep_rank_product.sort_by{|peptide,rp| rp}.each do |peptide,rp|
	protein_acc = nil
	description = nil
	infile_list.each_value do |csvp|
		if csvp.has_peptide(peptide)
			highest_scored_hit = csvp.highest_scored_hit_for_pep(peptide)
			uncharacterized_to_annotated(highest_scored_hit)
			protein_acc = highest_scored_hit.prot_acc.to_s
			description = highest_scored_hit.prot_desc.to_s
			peps_sorted_by_rank_product_out.puts '"' + protein_acc + '","' + description + '","' + peptide.to_s + '","' + pep_scores[peptide].join('","') + '","' + (pep_rank_product[peptide]**(1/infile_list.length.to_f)).to_s + '"'
		end
	end
end
peps_sorted_by_rank_product_out.close




