# Trace back to Entryname from refseq accno (because data is from uniprot and multiple alignments from ucsc).

class AccnoToRefseqTranslator

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
		fp = AccnoToRefseqTranslator.new(filename)
		if block_given?
			fp.each do |entry|
				yield entry
			end
		else
			return fp
		end
	end
	
	def create_index()
		@filehandle.each do |line|
			if line[0..0] != "#"
				line = line.chomp
				line =~ /(.+)\s+(.+)/
				accno, genename = $1, $2
				if genename.include? ","
					genename = genename.split(',')[0]
				end
				@entrypos_by_id[genename] = @filehandle.pos - line.length - 1
				@entrypos << @filehandle.pos - line.length - 1
				self.increment
			end
		end
		@filehandle.rewind
	end

	def each()
		@filehandle.each do |line|
			if line[0..0] != "#"
				yield line_parse(line)
			end
		end
	end
	
	def refseq_from_genename(genename)
		if @entrypos_by_id.has_key?(genename)
			@filehandle.pos = @entrypos_by_id[genename]
			return line_parse(@filehandle.readline)
		else
			return nil
		end
	end

	def line_parse(line)
		accno, genename = '',''
		line = line.chomp
		line =~ /(.+)\s+(.+)/
		accno, genename = $1, $2
		if genename.include? ","
			genename = genename.split(',')[0]
		end
		return Entry.new(accno, genename)
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

class Entry
	attr_accessor :accno, :genename

	def initialize(accno, genename)
		@accno = accno
		@genename = genename
	end
end
