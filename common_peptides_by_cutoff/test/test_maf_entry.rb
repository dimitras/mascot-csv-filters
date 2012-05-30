require 'test/unit'
require 'maf_entry'

class TestMAFEntry < Test::Unit::TestCase

	def setup
		@entry = MAFEntry.new("NM_001004310","hg19", 435,"chr1", 159772215, 159785451, "+", "MLLWTAVLLFVPCVGKTVWLYLQAWPNPVFEGDALTLRCQGWKNTPLSQVKFYRDGKFLHFSKENQTLSMGAATVQSRGQYSCSGQVMYIPQTFTQTSETAMVQVQELFPPPVLSAIPSPEPREGSLVTLRCQTKLHPLRSALRLLFSFHKDGHTLQDRGPHPELCIPGAKEGDSGLYWCEVAPEGGQVQKQSPQLEVRVQAPVSRPVLTLHHGPADPAVGDMVQLLCEAQRGSPPILYSFYLDEKIVGNHSAPCGGTTSLLFPVKSEQDAGNYSCEAENSVSRERSEPKKLSLKGSQVLFTPASNWLVPWLPASLLGLMVIAAALLVYVRSWRKAGPLPSQIPPTAPGGEQCPLYANVHHQKGKDEGVVYSVVHRTSKRSEARSAEFTVGRKDSSIICAEVRCLQPSEVSSTEVNMRSRTLQEPLSDCEEVLCZ")
	end
	
	def test_maf_entry
		assert_kind_of(MAFEntry, @entry)
	end

	def test_attributes
		assert_equal("NM_001004310", @entry.accno)
		assert_equal("hg19", @entry.species)
	end

	def test_positions_of_subseq()
		subseq = 'WTAVL'
		assert_equal([3], @entry.positions_of_subseq(subseq))
	end

	def test_letter_at_position()
		pos = 5
		assert_equal('A', @entry.letter_at_position(pos))
	end

end
