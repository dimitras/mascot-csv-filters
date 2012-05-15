require 'csv_parser_for_protein_hits'
require 'test/unit'

class TestCSVParserForProteinHits < Test::Unit::TestCase

	def setup
		cutoff = 0.5
		@peptide = 'SNVVEEFQSPFASQSR'
		@infile_list = {}
		Dir["test/data/*.csv"].each do |infile|
			file_num = infile.split("/")[2].split("-")[0].to_i - 1
			if file_num >= 0 && file_num <= 6
				@infile_list[infile] = CSVParserForProteinHits.open(infile, cutoff, false, file_num)
			else
				@infile_list[infile] = CSVParserForProteinHits.open(infile, cutoff, true, file_num)
			end
		end
	end
	
	def test_csvfile
		@infile_list.each_value do |csvp|
			assert_kind_of(CSVParserForProteinHits, csvp)
		end
	end

	def test_csvp
		@infile_list.each_value do |csvp|
			if csvp.id == 8
				assert_equal(true, csvp.is_ascending)
			else
				assert_equal(false, csvp.is_ascending)
			end
		end
	end

	def test_worst_rank
		@infile_list.each_value do |csvp|
			if csvp.id == 8
				assert_equal(1, csvp.worst_rank)
			end
		end
	end

	def test_is_descending
		@infile_list.each_value do |csvp|
			if csvp.id == 0
				assert_equal(true, csvp.is_descending)
			end
		end
	end

	def test_highest_scored_hit
		@infile_list.each_value do |csvp|
			if csvp.id == 8
				assert_equal(33.25, csvp.highest_scored_hit_for_pep(@peptide).pep_score)
			end
		end
	end

	def test_has_peptide
		@infile_list.each_value do |csvp|
			assert_equal(true, csvp.has_peptide(@peptide))
		end
	end

	def test_line_parse
		@hit = '2471|E9PDJ7|Uncharacterized protein OS=Homo sapiens GN=NMU PE=4 SV=1|40|16983|5|0|1|0|335474|1|1|1|634.6417|1900.9033|3|1900.9074|-0.004|0|33.25|0.081|K|SNVVEEFQSPFASQSR|G|2 Acetyl:2H(3) (S)|0.2000000000000020.0|platelet_asa_100611_A3.27427.27427.3'
		@infile_list.each_value do |csvp|
			if csvp.id == 8
				assert_equal(@peptide, csvp.line_parse(@hit).pep_seq)
			end
		end
	end

	def test_protein_hits
		@infile_list.each_value do |csvp|
			if csvp.id == 8
				assert_equal(@peptide, csvp.protein_hits(@peptide).first.pep_seq)
			end
		end
	end

	def test_sort_with_each
		@infile_list.each_value do |csvp|
			if csvp.id == 0
				csvp.each_peptide do |pep|
					# the second peptide should be the SNVVEEFQSPFASQSR
					next
					assert_equal(@peptide, pep)
				end
			end
		end
	end

end
