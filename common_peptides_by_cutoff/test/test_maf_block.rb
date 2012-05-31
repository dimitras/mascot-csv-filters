require 'test/unit'
require 'maf_block'
require 'maf_like_parser'
require 'maf_entry'

class TestMAFEntry < Test::Unit::TestCase

	def setup
		@mafp = MAFlikeParser.open("test/data/aligned_sample.maf")
		@block = @mafp.maf_block_by_id('NM_001004310')
		@subseq = 'WTAVLL'
		@letter = 'L'
		@species_ids = @mafp.species_ids
	end
	
	def test_maf_block
		assert_kind_of(MAFBlock, @block)
	end

	def test_maf_file
		assert_kind_of(MAFlikeParser, @mafp)
	end

	def test_attributes
		assert_equal("hg19", @block.ref_species)
	end

	def test_subseq_contained_in_ref_species?()
		assert_equal(true, @block.subseq_contained_in_ref_species?(@subseq))
	end

	def test_find_seq_positions_in_ref_species
		assert_equal([3], @block.find_seq_positions_in_ref_species(@subseq))
	end

	def test_find_positions_for_seq_letter_in_ref_species
		assert_equal([7,8], @block.find_positions_for_seq_letter_in_ref_species(@subseq, @letter))
	end

	def test_corresponding_letters_for_secondary_species
		assert_equal(["L(7)L(8)", "L(7)L(8)", "L(7)L(8)", "NANA"], @block.corresponding_letters_for_secondary_species(@subseq, @letter, @species_ids))
	end

end
