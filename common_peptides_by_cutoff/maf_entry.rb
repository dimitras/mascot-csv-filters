# MAFEntry is an entry of the maf_block, that contains one fasta-like entry with header info and a sequence.
# >NM_001004310_hg19 435 chr1:159772215-159785451+
# MLLWTAVLLFVPCVGKTVWLYLQAWPNPVFEGDALTLRCQGWKNTPLSQVKFYRDGKFLHFSKENQTLSMGAATVQSRGQYSCSGQVMYIPQTFTQTSETAMVQVQELFPPPVLSAIPSPEPREGSLVTLRCQTKLHPLRSALRLLFSFHKDGHTLQDRGPHPELCIPGAKEGDSGLYWCEVAPEGGQVQKQSPQLEVRVQAPVSRPVLTLHHGPADPAVGDMVQLLCEAQRGSPPILYSFYLDEKIVGNHSAPCGGTTSLLFPVKSEQDAGNYSCEAENSVSRERSEPKKLSLKGSQVLFTPASNWLVPWLPASLLGLMVIAAALLVYVRSWRKAGPLPSQIPPTAPGGEQCPLYANVHHQKGKDEGVVYSVVHRTSKRSEARSAEFTVGRKDSSIICAEVRCLQPSEVSSTEVNMRSRTLQEPLSDCEEVLCZ

class MAFEntry
  attr_accessor :accno, :species, :size, :chrom, :genomic_start, :genomic_stop, :strand, :seq

  def initialize(accno, species, size, chrom, genomic_start, genomic_stop, strand, seq)
    @accno = accno
    @species = species
    @size = size.to_i
    @chrom = chrom
    @genomic_start = genomic_start.to_i
    @genomic_stop = genomic_stop.to_i
    @strand = strand
    @seq = seq
  end

  def positions_of_subseq(subseq)
    positions = (0 .. @seq.length - 1).find_all { |i| @seq[i, subseq.length] == subseq }
    return positions
  end

  def letter_at_position(pos)
    return @seq[pos..pos]
  end

end

