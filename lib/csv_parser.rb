# FasterCSV.foreach("data.csv", :quote_char => '"', :col_sep =>'\t', :row_sep =>:auto) do |row|	end

# require 'rubygems'
# require 'faster_csv'
require 'protein_hit'

class CSVParser
	attr_accessor :filename

	def initialize(filename)
		@filename = filename
		@filehandle = File.new(filename)
		@index = {}
		create_index
	end

	def create_index()
		in_protein_hits_table = false
		@filehandle.each do |line|
			line_pos = @filehandle.pos - line.length
			(prot_hit_num, prot_acc, prot_desc, prot_score, prot_mass, prot_matches, prot_matches_sig, prot_sequences, prot_sequences_sig, pep_query, pep_rank, pep_isbold, pep_isunique, pep_exp_mz, pep_exp_mr, pep_exp_z, pep_calc_mr, pep_delta, pep_miss, pep_score, pep_expect, pep_res_before, pep_seq, pep_res_after, pep_var_mod, pep_var_mod_pos, pep_scan_title) = line.chomp!.split("|")
			if prot_hit_num.to_s == "prot_hit_num"
				in_protein_hits_table = true
				next
			end
			if in_protein_hits_table == true
				if pep_var_mod.include? 'Acetyl'
					if !@index.has_key?(pep_seq) then
						@index[pep_seq] = []
					end
					@index[pep_seq] << line_pos
				end
			end
		end
	end
		
	def self.open(filename)
		csvp = CSVParser.new(filename)
		if block_given?
			csvp.each do |hit|
				yield hit
			end
		else
			return csvp
		end
	end

	def each()
		@index.each_key do |key|
			hits = []
			@index[key].each do |hit_pos|
				@filehandle.pos = hit_pos
				hits << line_parse(@filehandle.readline)
			end
			yield hits
		end
	end
	
	def line_parse(line)
		(prot_hit_num, prot_acc, prot_desc, prot_score, prot_mass, prot_matches, prot_matches_sig, prot_sequences, prot_sequences_sig, pep_query, pep_rank, pep_isbold, pep_isunique, pep_exp_mz, pep_exp_mr, pep_exp_z, pep_calc_mr, pep_delta, pep_miss, pep_score, pep_expect, pep_res_before, pep_seq, pep_res_after, pep_var_mod, pep_var_mod_pos, pep_scan_title) = line.chomp.split("|")
		return ProteinHit.new(prot_hit_num, prot_acc, prot_desc, prot_score, prot_mass, prot_matches, prot_matches_sig, prot_sequences, prot_sequences_sig, pep_query, pep_rank, pep_isbold, pep_isunique, pep_exp_mz, pep_exp_mr, pep_exp_z, pep_calc_mr, pep_delta, pep_miss, pep_score, pep_expect, pep_res_before, pep_seq, pep_res_after, pep_var_mod, pep_var_mod_pos, pep_scan_title, @filename)
	end
	
	def each_peptide()
		@index.each_key do |key|
			yield key
		end
	end
	
	def has_peptide(peptide)
		if @index.has_key?(peptide)
			return true
		else 
			return false
		end
	end
	
	def protein_hits(peptide)
		hits = []
		@index[peptide].each do |hit_pos|
			@filehandle.pos = hit_pos
			hits << line_parse(@filehandle.readline)
		end
		return hits
	end
	
end

# csvp = CSVParser.open(ARGV[0])
# # csvp.each do |hit|
# # # 	puts ">>> peptide:\t" + hit.pep_seq.to_s + "\tmodification:\t" + hit.pep_var_mod.to_s + "\tfilename:\t" + hit.filename.to_s
# # end
# 
# csvp.protein_hits("SQLAAAAR").each do |hit|
# 	puts ">>> peptide:\t" + hit.pep_seq.to_s + "\tmodification:\t" + hit.pep_var_mod.to_s + "\tfilename:\t" + hit.filename.to_s
# end	
	
