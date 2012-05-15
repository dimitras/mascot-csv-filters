require 'test/unit'
require 'protein_hit'

class TestProteinHit < Test::Unit::TestCase

	def setup
		@prot_hit = ProteinHit.new("2471","E9PDJ7","Uncharacterized protein OS=Homo sapiens GN=NMU PE=4 SV=1","40","16983","5","0","1","0","335474","1","1","1","634.6417","1900.9033","3","1900.9074","-0.004","0","33.25","0.081","K","SNVVEEFQSPFASQSR","G","2 Acetyl:2H(3) (S)","0.2000000000000020.0","platelet_asa_100611_A3.27427.27427.3","test/data/9-rep_sample.csv")
	end
	
	def test_protein_hit
		assert_kind_of(ProteinHit, @prot_hit)
	end

	def test_prot_hit
		assert_equal("33.25", @prot_hit.pep_score)
	end

	def test_to_csv
		hit = '"2471","E9PDJ7","Uncharacterized protein OS=Homo sapiens GN=NMU PE=4 SV=1","40","16983","5","0","1","0","335474","1","1","1","634.6417","1900.9033","3","1900.9074","-0.004","0","33.25","0.081","K","SNVVEEFQSPFASQSR","G","2 Acetyl:2H(3) (S)","0.2000000000000020.0","platelet_asa_100611_A3.27427.27427.3","test/data/9-rep_sample.csv"'
		assert_equal(hit, @prot_hit.to_csv)
	end

end
