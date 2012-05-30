require 'test/unit'
require 'accno_to_refseq_translator'

class TestAccnoToRefseqTranslator < Test::Unit::TestCase
 	
	def setup
		@fp = AccnoToRefseqTranslator.open("test/data/accno2refseq_sample.txt")
		counter = 0
	   	@fp.each do |entry|
	   		if counter == 1
		   		@entry = entry
		   	end
	   		counter +=1
	   	end
	end

	def test_file
		assert_kind_of(AccnoToRefseqTranslator, @fp)
	end

	def test_entries_number
		assert_equal(5, @fp.count)
	end

	def test_entry
		assert_equal("NM_032291", @fp.first.accno)
		assert_equal("NM_013943", @fp.last.accno)
	end	

	def test_each
		assert_equal("NM_001080397", @entry.accno) # expects the second one that finds
	end

	def test_refseq_from_genename
		assert_equal("NM_052998", @fp.refseq_from_genename('ADC').accno) # 1 gn mporei na xei n accs
	end

	def test_entry
		assert_kind_of(Entry, @entry)
	end

	def test_attributes
		assert_equal("NM_001080397", @entry.accno)
		assert_equal("SLC45A1", @entry.genename)
	end

end