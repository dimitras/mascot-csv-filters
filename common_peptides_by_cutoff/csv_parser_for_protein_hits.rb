require 'protein_hit'

class CSVParserForProteinHits
	attr_accessor :filename, :cutoff, :is_ascending, :id

	def initialize(filename, cutoff, is_ascending, id)
		@filename = filename
		@cutoff = cutoff
		@is_ascending = is_ascending
		@id = id.to_i
		@filehandle = File.new(filename)
		@pep_index = {}
		@rank_index = {}
		create_pep_index
		create_rank_index
	end

	def self.open(filename, cutoff, is_ascending, id)
		csvp = CSVParserForProteinHits.new(filename,cutoff, is_ascending, id)
		if block_given?
			csvp.each do |hit|
				yield hit
			end
		else
			return csvp
		end
	end

	def create_pep_index()
		@filehandle.readline
		@filehandle.each do |line|
			line_pos = @filehandle.pos - line.length
			hit = line_parse(line)
			if hit.pep_expect < @cutoff
				if !@pep_index.has_key?(hit.pep_seq)
					@pep_index[hit.pep_seq] = []
				end
				@pep_index[hit.pep_seq] << line_pos
			end
		end
	end

	def create_rank_index()
		score_index = []
		@pep_index.each_key do |key|
			@pep_index[key].each do |line_pos|
				@filehandle.pos = line_pos
				hit = line_parse(@filehandle.readline)
				score_index << [line_pos, hit.pep_score]
			end
		end

		if @is_ascending
			score_index.sort!  {|x,y| x[1] <=> y[1]}#.each do |e| 
			# 	puts "#{e[0]} => #{e[1]}"
			# end
		else
			score_index.sort! {|x,y| y[1] <=> x[1]}#.each do |e| 
			# 	puts "#{e[0]} => #{e[1]}"
			# end
		end

		score_index.each_with_index do |(line_pos, pep_score), index|
			@rank_index[line_pos.to_s] = index + 1
		end
	end

	def worst_rank()
		if @is_ascending
			return 1
		else
			return @rank_index.length
		end
	end

	def highest_scored_hit_for_pep(peptide)
		highest_scored_hit = nil
		hits = []
		hits = protein_hits(peptide)
		score = 0
		hits.each do |hit|
			if hit.pep_score.to_f >= score
				score = hit.pep_score
				highest_scored_hit = hit	
			end
		end
		return highest_scored_hit
	end

	def each()
		@pep_index.each_key do |key|
			hits = []
			@pep_index[key].each do |hit_pos|
				hit = hit_from_pos(hit_pos)
				hits << hit
			end
			yield hits
		end
	end
	
	def each_hit()
		@pep_index.each_key do |key|
			@pep_index[key].each do |hit_pos|
				yield hit_from_pos(hit_pos)
			end
		end
	end

	def each_peptide()
		@pep_index.each_key do |key|
			yield key
		end
	end

	def has_peptide(peptide)
		if @pep_index.has_key?(peptide)
			return true
		else
			return false
		end
	end

	def hit_from_pos(hit_pos)
		@filehandle.pos = hit_pos
		hit = line_parse(@filehandle.readline)
		hit.rank = @rank_index[hit_pos.to_s]
		return hit
	end

	def protein_hits(peptide)
		hits = []
		if !@pep_index.has_key?(peptide)
			return hits
		end
		@pep_index[peptide].each do |hit_pos|
			hit = hit_from_pos(hit_pos)
			hits << hit
		end
		return hits
	end

	def is_descending()
		return !@is_ascending
	end
	
	def line_parse(line)
		(prot_hit_num, prot_acc, prot_desc, prot_score, prot_mass, prot_matches, prot_matches_sig, prot_sequences, prot_sequences_sig, pep_query, pep_rank, pep_isbold, pep_isunique, pep_exp_mz, pep_exp_mr, pep_exp_z, pep_calc_mr, pep_delta, pep_miss, pep_score, pep_expect, pep_res_before, pep_seq, pep_res_after, pep_var_mod, pep_var_mod_pos, pep_scan_title) = line.chomp.split("|")
		return ProteinHit.new(prot_hit_num, prot_acc, prot_desc, prot_score, prot_mass, prot_matches, prot_matches_sig, prot_sequences, prot_sequences_sig, pep_query, pep_rank, pep_isbold, pep_isunique, pep_exp_mz, pep_exp_mr, pep_exp_z, pep_calc_mr, pep_delta, pep_miss, pep_score, pep_expect, pep_res_before, pep_seq, pep_res_after, pep_var_mod, pep_var_mod_pos, pep_scan_title, @filename)
	end

end

