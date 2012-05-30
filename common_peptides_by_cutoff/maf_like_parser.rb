# parse a multiple alignment fasta file with blocks of entries like the following (one block contains the alignments for all species for a protein):
# >NM_001004310_hg19 435 chr1:159772215-159785451+
# MLLWTAVLLFVPCVGKTVWLYLQAWPNPVFEGDALTLRCQGWKNTPLSQVKFYRDGKFLHFSKENQTLSMGAATVQSRGQYSCSGQVMYIPQTFTQTSETAMVQVQELFPPPVLSAIPSPEPREGSLVTLRCQTKLHPLRSALRLLFSFHKDGHTLQDRGPHPELCIPGAKEGDSGLYWCEVAPEGGQVQKQSPQLEVRVQAPVSRPVLTLHHGPADPAVGDMVQLLCEAQRGSPPILYSFYLDEKIVGNHSAPCGGTTSLLFPVKSEQDAGNYSCEAENSVSRERSEPKKLSLKGSQVLFTPASNWLVPWLPASLLGLMVIAAALLVYVRSWRKAGPLPSQIPPTAPGGEQCPLYANVHHQKGKDEGVVYSVVHRTSKRSEARSAEFTVGRKDSSIICAEVRCLQPSEVSSTEVNMRSRTLQEPLSDCEEVLCZ
# >NM_001004310_gorGor1 435 Supercontig_0001193:17388-22022+;Supercontig_0323122:103-168-;Supercontig_0001193:29116-29349+
# MLLWTXVLLFVPCVGKT---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------SQVLSTPXSNWLVPWLPASLLGLMVIAAALLVYLRSWRKXGPLPSQIPPTAPGGEQCPLYAX-------------------------RSAEFTVGRKDSSIICAEVRCLQPXEVSSKEVNXRSRTLQEPLGDCEEVLCZ
# >NM_001004310_mm9 435 chr1:174527708-174532591-
# MLLWMVLLLC-------------------------------------------------------------------------------------------MAEAQELFPNPELTEFTNSETMD----ILKCTIKVDPKNPTLQLFYTFYKNNHVIQDRSPHSVFS--EAKEENSGLYQCMVDTEDGLIQKKSGYLDIQFWTPVSHPVLTLQHEATNLAVGDKVEFLCEAHQGSLPIFYSFYINGEILGKPLAPSGRAASLLASVKAEWSTKNYSCEAKNNISREISELKKFPLVESQVL--WMSNMLPIWLPASLLGGMVIAAVVLMYFKPCKKH-----------------------------------------------RPETPT--ZKDS--LYVSVEIKDINEIPTNDLDSKTRTCQDPLGLZDHA-CZ
# >NM_001004310_rn4 435 chr13:88605484-88611105-
# MLLWMVLLLCV------------------------------------------------------------------------------------------MTEAQELFQDPVLSRLNSSETSD----LLKCTTKVDPNKPASELFYSFYKDNHIIQNRSHNPLFFISEANEENSGLYQCVVDAKDGTIQKKSDYLDIDLCTSVSQPVLTLQHEATNLAEGDKVKFLCETQLGSLPILYSFYMDGEILGEPLAPSGRAASLLISVKAEWSGKNYSCQAENKVSRDISEPKKFPLVESQVL--SMSTTMVIWLPV-LLGGMAMAAAVLIZFKPCKKH-----------------------------------------------RPETPT--QKDS--LYASLTKGEINEVPTNHLVPNTRTYQDPLGLYNDA---
# >NM_001004310_oryCun2 435 chr13:33310470-33323414-
# MLPWISLLLFVPCDGANGWLNIQALPNPVLEGESLTLRCQGWNNMELSQVNFYKDGARLHSSNGKNILSIQTATLKSSGQYSCSGSVTYLPHI-RETSGTTTVQVQELFSPPVLNAIPSPKLCEGSPLTLRCQTRLHPKRPALRLLFSFYKDGRTVQGRGQHPELHILAMEEGDSGLYWCEAAPEGGRVQKQSPQLEIMVQAPVSRPLLTLPHRIISPAVGDLVELFCEVQRGSPPVLYSFYFDGKLLGNLSAPHGGAASLLVPVKSEQDTGNYSCEAENSVSKERSEPRKISLKGSQVLSTPTCNWPVVYLPAGLLGVMVIAAALLGYFKPWKKA-----------------------------------------------RPVQLTSEPKDISIIYAEVRKPGHREISPEELSQNH-ARPGPLGDYDEV-YZ
# >NM_001004310_bosTau4 435 chr3:10843145-10863767-
# MLLWTVALLFVPCVGKIVWLSP-AQPQLVFEGDILILRCRGRKNAALSHVKFYRDGKFLHFSMENQPLLLGTATANSSGWYNCTGRVKYPRNMDSWDSGTAMVQVQELFLPPVLTALPSHELCEGSPVTLRCQTKLHSQKSASRLLFSFHKEGRTLKNRSSRPELRIPAAKEGDAGLYWCKASSEGGQIQKRSPQLELRVWAPVSRPLLTLR--PTSLAVGDEVELLCEVQRGSPPILYSFHLNGDILRNHVAPHGGPASCLFRVTSEQDAGNYSCEAGNRVSRETSEPETLSVDDPQVLSAPTSNWLVPWLPASLLAMMVIAAALLGYCRPWRKNGPLPPRNLPSAPSEEQHPLYVNVYHQNENNEGVIYSEIHTIPREHEARPAQPAQQEKDISVIYAEVRHPQLSKNPDKGLNRRSTIHZVPTSDYKEALCZ

require 'maf_block'
require 'maf_entry'

class MAFlikeParser
	attr_accessor :filename, :count, :entrypos, :species_ids

	def initialize(filename)
		@filename = filename
		@filehandle = File.new(filename)
		@count = 0
		@entrypos = []
		@entrypos_by_id = {}
		create_index
		@species_ids = init_species_ids
	end

	def self.open(filename)
		mp = MAFlikeParser.new(filename)
		if block_given?
			mp.each do |block|
				yield block
			end
		else
			return mp
		end
	end

	def create_index()
		@filehandle.each do |line|
			if line[0..0] == ">"
				line =~ />(NM_.+)_(.+)\s(\d+)\s(.+):(\d+)-(\d+)(\D)/
				accno = $1
				if !@entrypos_by_id.has_key?(accno)
					@entrypos_by_id[accno] = @filehandle.pos - line.length
					@entrypos << @filehandle.pos - line.length
					self.increment
				end
			end
		end
		@filehandle.rewind
	end

	def each()
		@entrypos.each do |pos|
			yield read_maf_block_at_pos(pos)
		end
	end
	
	def maf_block(entryno)
		if entryno >= 1 && entryno <= @entrypos.length
			return read_maf_block_at_pos(@entrypos[entryno-1])
		else
			return nil
		end
	end
	
	def maf_block_by_id(accno)
		if @entrypos_by_id.has_key?(accno)
			return read_maf_block_at_pos(@entrypos_by_id[accno])
		else
			return nil
		end
	end

	def read_maf_block_at_pos(pos)
		maf_entries = []
		block_id = nil
		ref_species = nil
		@filehandle.pos = pos
		@filehandle.each do |line|
			line = line.chomp
			line =~ />(NM_.+)_(.+)\s(\d+)\s(.+):(\d+)-(\d+)(\D)/
			if line[0..0] == ">"
				accno, species, size, chrom, genomic_start, genomic_stop, strand, seq = $1, $2, $3, $4, $5, $6, $7, @filehandle.readline.chomp
				maf_entries << MAFEntry.new(accno, species, size, chrom, genomic_start, genomic_stop, strand, seq)
				if block_id == nil
					block_id = accno
					ref_species = species
				end
			elsif line == ""
				@filehandle.readline
				return MAFBlock.new(maf_entries, ref_species)
			end
		end
	end

	def init_species_ids()
		species_ids = {}
		@filehandle.each do |line|
			if line[0..0] == ">"
				line =~ />(NM_.+)_(.+)\s(\d+)\s(.+):(\d+)-(\d+)(\D)/
				species = $2
				species_ids[species] = nil
			end
		end
		@filehandle.rewind
		return species_ids.keys
	end

	def increment()
		@count += 1
	end
	
	def first()
		return self.maf_block(1)
	end

	def last()
		return self.maf_block(@entrypos.length)
	end

end

