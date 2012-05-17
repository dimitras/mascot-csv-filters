require 'test/unit'
require 'fasta_entry'

class TestFastaEntry < Test::Unit::TestCase

	def setup
		@entry = FastaEntry.new("sp","A0A183","LCE6A_HUMAN Late cornified envelope protein 6A OS=Homo sapiens GN=LCE6A PE=2 SV=1","MSQQKQQSWKPPNVPKCSPPQRSNPCLAPYSTPCGAPHSEGCHSSSQRPEVQKPRRARQKLRCLSRGTTYHCKEEECEGD")
	end
	
	def test_fasta_entry
		assert_kind_of(FastaEntry, @entry)
	end

	def test_attributes
		assert_equal("A0A183", @entry.accno)
		assert_equal("MSQQKQQSWKPPNVPKCSPPQRSNPCLAPYSTPCGAPHSEGCHSSSQRPEVQKPRRARQKLRCLSRGTTYHCKEEECEGD", @entry.seq)
	end

end
