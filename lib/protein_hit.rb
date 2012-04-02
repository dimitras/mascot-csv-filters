# # MODULE : filter by peptide && modification (keep only Acetyl*) and keep the entry (peptide) with highest score 
# def common_peptide() end # Finding Common Keys in Hashes : common = hash1.keys & hash2.keys
# 	# for keys of hash 1 check if hash2 key exists & hash3 key exists > keep all entries from hashes

# # MODULE : create new file with common peptides keeping all row info plus a column for filename
# def to_filtered_csv() end

class ProteinHit
	attr_accessor :prot_hit_num, :prot_acc, :prot_desc, :prot_score, :prot_mass, :prot_matches, :prot_matches_sig, :prot_sequences, :prot_sequences_sig, :pep_query, :pep_rank, :pep_isbold, :pep_isunique, :pep_exp_mz, :pep_exp_mr, :pep_exp_z, :pep_calc_mr, :pep_delta, :pep_miss, :pep_score, :pep_expect, :pep_res_before, :pep_seq, :pep_res_after, :pep_var_mod, :pep_var_mod_pos, :pep_scan_title, :filename

	def initialize(prot_hit_num, prot_acc, prot_desc, prot_score, prot_mass, prot_matches, prot_matches_sig, prot_sequences, prot_sequences_sig, pep_query, pep_rank, pep_isbold, pep_isunique, pep_exp_mz, pep_exp_mr, pep_exp_z, pep_calc_mr, pep_delta, pep_miss, pep_score, pep_expect, pep_res_before, pep_seq, pep_res_after, pep_var_mod, pep_var_mod_pos, pep_scan_title, filename)
		@prot_hit_num = prot_hit_num
		@prot_acc = prot_acc
		@prot_desc = prot_desc
		@prot_score = prot_score
		@prot_mass = prot_mass
		@prot_matches = prot_matches
		@prot_matches_sig = prot_matches_sig
		@prot_sequences = prot_sequences
		@prot_sequences_sig = prot_sequences_sig
		@pep_query = pep_query
		@pep_rank = pep_rank
		@pep_isbold = pep_isbold
		@pep_isunique = pep_isunique
		@pep_exp_mz = pep_exp_mz
		@pep_exp_mr = pep_exp_mr
		@pep_exp_z = pep_exp_z
		@pep_calc_mr = pep_calc_mr
		@pep_delta = pep_delta
		@pep_miss = pep_miss
		@pep_score = pep_score
		@pep_expect = pep_expect
		@pep_res_before = pep_res_before
		@pep_seq = pep_seq
		@pep_res_after = pep_res_after
		@pep_var_mod = pep_var_mod
		@pep_var_mod_pos = pep_var_mod_pos
		@pep_scan_title = pep_scan_title
		@filename = filename
	end
	
	def to_csv()
		hit = '"' + [@prot_hit_num, @prot_acc, @prot_desc, @prot_score, @prot_mass, @prot_matches, @prot_matches_sig, @prot_sequences, @prot_sequences_sig, @pep_query, @pep_rank, @pep_isbold, @pep_isunique, @pep_exp_mz, @pep_exp_mr, @pep_exp_z, @pep_calc_mr, @pep_delta, @pep_miss, @pep_score, @pep_expect, @pep_res_before, @pep_seq, @pep_res_after, @pep_var_mod, @pep_var_mod_pos, @pep_scan_title, @filename].join('","') + '"'
		return hit
	end
end
	