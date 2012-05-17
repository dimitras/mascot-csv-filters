# parse a fasta file (uniprot proteome database) with entries like the following:
# >sp|A0A5B9|TRBC2_HUMAN T-cell receptor beta-2 chain C region OS=Homo sapiens GN=TRBC2 PE=1 SV=1
# DLKNVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSWWVNGKEVHSGVSTDPQ
# PLKEQPALNDSRYCLSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIV
# SAEAWGRADCGFTSESYQQGVLSATILYEILLGKATLYAVLVSALVLMAMVKRKDSRG

require 'fasta_entry'

class FastaParser
	attr_accessor :filename, :count, :entrypos

	def initialize(filename)
		@filename = filename
		@filehandle = File.new(filename)
		@count = 0
		@entrypos = []
		@entrypos_by_id = {}
		create_index
	end

	def self.open(filename)
		fp = FastaParser.new(filename)
		if block_given?
			fp.each do |entry|
				yield entry
			end
		else
			return fp
		end
	end
	
	def each()
		tag, accno, desc, seq = '','','',''
		@filehandle.each do |line|
			line = line.chomp
			if ((line[0..0] == ">") && (tag == '') && (accno == '') && (desc == ''))
				line =~ />(.+)\|(.+)\|(.+)/
				tag, accno, desc, seq = $1, $2, $3, ''
			elsif (line[0..0] == ">")
				yield FastaEntry.new(tag, accno, desc, seq)
				line =~ />(.+)\|(.+)\|(.+)/
				tag, accno, desc, seq = $1, $2, $3, ''
			else
				seq += line
			end
		end
		yield FastaEntry.new(tag, accno, desc, seq)
	end

	def create_index()
		@filehandle.each do |line|
			if line[0..0] == ">"
				line =~ />.+\|(.+)\|.+/
				accno = $1
				@entrypos_by_id[accno] = @filehandle.pos - line.length
				@entrypos << @filehandle.pos - line.length
				self.increment
			end
		end
		@filehandle.rewind
	end

	def entry(entryno)
		@filehandle.pos = @entrypos[entryno-1] # using 1-based numbering scheme
		line = @filehandle.readline.chomp
		line =~ />(.+)\|(.+)\|(.+)/
		tag, accno, desc, seq = $1, $2, $3, ''
		@filehandle.each do |line|
			line = line.chomp
			if (line[0..0] == ">")
				break
			end
			seq = seq + line
		end
		return FastaEntry.new(tag, accno, desc, seq)
	end
	
	def entry_by_id(accno)
		if @entrypos_by_id.has_key?(accno)
			@filehandle.pos = @entrypos_by_id[accno]
			line = @filehandle.readline.chomp
			line =~ />(.+)\|.+\|(.+)/
			tag, desc, seq = $1, $2, ''
			@filehandle.readline
			@filehandle.each do |line|
				if line[0..0] == ">"
					break
				end
				seq = seq + line.chomp
			end
			return FastaEntry.new(tag, accno, desc, seq)
		else
			return nil
		end
	end

	def increment()
		@count += 1
	end
	
	def first()
		return self.entry(1)
	end

	def last()
		return self.entry(@entrypos.length)
	end

end

