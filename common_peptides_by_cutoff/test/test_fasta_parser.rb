require 'fasta_parser'
require "test/unit"
 
class TestFastaParser < Test::Unit::TestCase
 	
	def setup
		@fap = FastaParser.open("test/data/proteome_uniprot_db_sample.fa")
		counter = 0
	   	@fap.each do |entry|
	   		if counter == 1
		   		@entry = entry
		   	end
	   		counter +=1
	   	end
	end

	def test_fapfile
		assert_kind_of(FastaParser, @fap)
	end

	def test_entries_number
		assert_equal(3, @fap.count)
	end

	def test_entry
		assert_equal("A0A183", @fap.first.accno)
		assert_equal("A0AUK0", @fap.last.accno)
		assert_equal("A0A5B9", @fap.entry(2).accno)
		assert_equal("LCE6A_HUMAN Late cornified envelope protein 6A OS=Homo sapiens GN=LCE6A PE=2 SV=1", @fap.entry_by_id("A0A183").desc)
		assert_equal(80, @fap.first.seq.length)
	end	

	def test_each
		assert_equal("A0A5B9", @entry.accno) # expects the second one that finds
	end

end