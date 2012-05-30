require 'test/unit'
require 'maf_like_parser'
 
class TestMAFlikeParser < Test::Unit::TestCase
 	
	def setup
		@mafp = MAFlikeParser.open("test/data/aligned_sample.maf")
		counter = 0
	   	@mafp.each do |block|
			counter += 1
			if counter == 2
		   		@block = block
		   	end
	   	end
	end

	def test_maf_file
		assert_kind_of(MAFlikeParser, @mafp)
		assert_kind_of(MAFBlock, @block)
	end

	def test_block_number
		assert_equal(2, @mafp.count)
	end

	def test_block
		assert_equal("NM_001004310", @mafp.first.accno)
		assert_equal("NM_001093725", @mafp.last.accno)
	end

	def test_each
		assert_equal("NM_001093725", @block.accno) # expects the second one that finds
	end

	def test_maf_block
		assert_equal("NM_001004310", @mafp.maf_block(1).accno)
	end

	def test_maf_block_by_id
		assert_equal("NM_001004310", @mafp.maf_block_by_id('NM_001004310').accno)
	end

	def test_species_ids
		assert_equal(4, @mafp.species_ids.length)
	end

end