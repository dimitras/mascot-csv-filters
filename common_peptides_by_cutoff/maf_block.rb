# MAFBlock is a block of maf_entries that have the same accno but different species.

class MAFBlock
	attr_accessor :accno, :maf_entries, :ref_species

	def initialize(maf_entries, ref_species)
		@accno = nil
		@maf_entries = {}
		maf_entries.each do |maf_entry|
			@accno = maf_entry.accno
			@maf_entries[maf_entry.species] = maf_entry
		end
		@ref_species = ref_species
	end

	def find_seq_positions_in_ref_species(seq)
		return @maf_entries[@ref_species].positions_of_subseq(seq)
	end

	def find_positions_for_seq_letter_in_ref_species(seq, letter)
		letter_positions = []
		positions = find_seq_positions_in_ref_species(seq)
		letter_positions_in_seq = (0..seq.length - 1).find_all { |i| seq[i, 1] == letter }
		positions.each do |position|
			letter_positions_in_seq.each do |letter_position|
				letter_positions  << position + letter_position
			end
		end
		return letter_positions
	end

	def corresponding_letters_for_secondary_species(seq, letter, species_ids)
		letters = []
		positions = find_positions_for_seq_letter_in_ref_species(seq, letter)
		species_ids.each do |species_id|
			positions.each do |position|
				if !@maf_entries.has_key?(species_id)
					letters << "NA"
				else
					letters << @maf_entries[species_id].letter_at_position(position).to_s + "(" + position.to_s + ")"
				end
			end
		end
		return letters
	end
end