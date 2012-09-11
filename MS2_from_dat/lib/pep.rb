require 'rubygems'

# REMEMBER TO uncomment the correct modifications (MODS)
# for 3H acetylation
# MODS = {
#   'K' => 45.029395,
#   'S' => 45.029395,
#   'T' => 45.029395
# }
# for endogenous
MODS = {
  'K' => 42.010565,
  'S' => 42.010565,
  'T' => 42.010565
}

# masses taken from the dat file
MW = {
  'A' => 71.037114 ,
  'R' => 156.101111,
  'N' => 114.042927,
  'D' => 115.026943,
  'C' => 160.030649,
  'E' => 129.042593,
  'Q' => 128.058578,
  'G' => 57.021464 ,
  'H' => 137.058912,
  'I' => 113.084064,
  'L' => 113.084064,
  'K' => 128.094963,
  'M' => 131.040485,
  'F' => 147.068414,
  'P' => 97.052764 ,
  'S' => 87.032028 ,
  'T' => 101.047679,
  'U' => 150.953630,
  'W' => 186.079313,
  'Y' => 163.063329,
  'V' => 99.068414,
  "K*" => 128.094963 + MODS['K'],
  "S*" => 87.032028 + MODS['S'],
  "T*" => 101.047679 + MODS['T']
}

H = 1.007825
OH = 17.00274
H2O = 18.010565
NH3 = 17.026549

class Pep
	attr_accessor :iontable, :bions, :yions, :seq, :mzs, :intensities

	def initialize(seq, mzs, intensities, mods=[])
		@seq = seq.strip.upcase.split ""
		@mzs = mzs
		@intensities = intensities
		add_mods_to_pep_seq(mods)

		@iontable = []
		@bions = []
		@yions = []
		build_ions_tables()
		@assigned_yionstable = nil
	end

	def build_ions_tables()
		@seq.each_index do |i|
			idx = i + 1

			#Bions
			b = y = [nil]*5
			unless i == @seq.length - 1
				tmp = @seq[0..i]
				b  = calc_ion_mz_table(tmp)
				@bions.push(b)
			end

			# Yions
			unless i == 0
				tmp = @seq[i..-1]
				y  = calc_ion_mz_table(tmp,false)
				@yions.unshift(y)
			end

			@iontable<<[@seq[i], idx, @seq.length - i , b, y]
		end
	end

	def add_mods_to_pep_seq(mods)
		mods.each do |m|
			@seq[m - 1] += "*"
		end
	end

	def calc_mw(seq=[])
		mw = 0.0
		seq.each do |aa|
			mw += MW[aa]
		end
		return mw
	end

	def mw(charge=0)
		return calc_mw(@seq) + H2O + (H * charge)
	end

	def calc_ion_mz_table(seq,nterm=true)
		tmw = 0.0
		if !nterm
			tmw = H + OH
		end
		# ions = [seq, mh, mhh, mh-nh3 , mh-h20]
		ions = [seq.join(""), (calc_mw(seq) + tmw + H )]
		ions << ( ions[1] +  H ) / 2
		# need R,N or Q for immonium ions
		if seq.join("").match(/[RNQ]/)
			ions << ions[1] - NH3
		else
			ions << nil
		end
		# Need S or T for water losses
		if seq.join("").match(/[ST]/)
			ions << ions[1] - H2O
		else
			ions << nil
		end
		return ions
	end

	# Get the peaks of all assigned yions table (with mass-intensity pairs)
	def assigned_yionstable()
		ranked_idx = ranked_yions_intensities_idx()

		assigned_ions = Array.new
		all_labels = []
		ranked_idx.each_with_index do |i,ii|
			if( assigned_yions[i] && !assigned_yions[i][0].nil? && 1
				assigned_yions[i][0] > 0 )
				if all_labels.include? "y(#{assigned_yions[i][0]})"
					next
				end
				assigned_ions << [assigned_yions[i][1],intensities[i],"y(#{assigned_yions[i][0]})"]
				all_labels << "y(#{assigned_yions[i][0]})"
			end
		end
# 		assigned_ions.reject! { |c| c.empty? }
		return assigned_ions
	end

	# Get the peaks of all assigned yions table (with mass-intensity pairs)
	def assigned_yionshash()
		assigned_ions_hash = {}
		assigned_ions = self.assigned_yionstable

		assigned_ions.each do |row|
			assigned_ions_hash[row[1]] = row[2]
		end
		return assigned_ions_hash
	end

	# get all yions (mass-intensity) and the label(index) for each assigned one
	def print_all_yionstable_to_filehandle(filehandle)
		assigned_ions_hash = assigned_yionshash()

		for i in 0..mzs.length-1
			if assigned_ions_hash.has_key?(intensities[i].to_f)
				filehandle.puts [mzs[i], intensities[i], "y(#{assigned_yions[i][0]})"].join(',')
			else
				filehandle.puts [mzs[i], intensities[i], ''].join(',')
			end
		end
	end

	# Grab legitimate y ions and rank them in intensity desc
	def ranked_yions_intensities_idx()
		ranked_idx = []
		idxset = []
		intset=[]
		assigned_yions.each_index do |i|
			if assigned_yions[i] && assigned_yions[i][0] && !assigned_yions[i][1].nil? && assigned_yions[i][1] > 0
				idxset.push(i)
				intset.push(intensities[i])
			end
		end
		intset.each do
			maxi = intset.index(intset.max())
			ranked_idx.push( idxset[maxi] )
			intset[maxi] = 0
		end
		return ranked_idx
	end

	def assigned_yions()
		if @assigned_yions == nil
			assign_yions()
		end
		return @assigned_yions
	end

	def assign_yions()
		@assigned_yions = assign(@yions, 1.0)
		return @assigned_yions
	end

	# Assigns the ion table to a given mass spectrum, walking through the m/z and ion series arrays to assign the ions to a mass peaks, given a tolerance
	# Returns a map of peaks that have ions assigned to them
	# http://174.129.8.134/mascot/help/results_help.html#PEP
	# example: http://174.129.8.134/mascot/cgi/peptide_view.pl?from=100&to=800&_label_all=0&file=..%2Fdata%2F20120323%2FF001507.dat&query=20052&hit=1&section=5&ave_thresh=38&_ignoreionsscorebelow=0&report=0&_sigthreshold=0.05&_msresflags=3138&_msresflags2=10&percolate=0&percolate_rt=0&tick1=100&tick_int=50&range=700&index=B4DHQ2&px=1
	def assign(ions,tol=1.0)
	# ion is [idx,b,b++,b*,b0] or y
	pkmap = Array.new(mzs.length)
	i1 = i2 =  0

	# traverse through ios series looking for matches
	while i1 < mzs.length && i2 < ions.length
	  # don't have a mass to look for. next ion
	  unless ions[i2][1]
	    i2 += 1
	    next
	  end
	  # x = [idx,+,++,*,0]
	  x = [nil]*5
	  dff = mzs[i1] - ions[i2][1]
	  if dff > tol
	    # mass is too large, get next ion
	    pkmap[i1] = x
	    i2 += 1
	    next
	  elsif dff <  0 - tol
	    # mass is too small for the + ion
	    # checking other ions
	    # ++ # must check for these outside of this method
	    # immonium ions
	    if ions[i2][3] && ((mzs[i1] - ions[i2][3]).abs <= tol) && !((mzs[i1] - ions[i2][3]).abs > tol)
	      x[0]=  i2 + 1
	      x[3] = ions[i2][3]
	    end
	    # water loss
	    if ions[i2][4] && !((mzs[i1] - ions[i2][4]).abs <= tol) && !((mzs[i1] - ions[i2][4]).abs > tol)
	      x[0]=  i2 + 1
	      x[4] = ions[i2][4]
	    end
	    # mass is too small, advance
	    pkmap[i1] = x
	    i1 += 1
	    next
	  end
	  # set the index
	  x[0] = i2  + 1
	  x[1] = ions[i2][1]
	# puts "ION:\ti1 = #{i1} , i2 = #{i2} => with mass diff: #{mzs[i1]} - #{ions[i2][1]}"
	#       if pkmap[i1 - 1][0] == x[0]
	#         puts 'here'
	# #         compare intensities and pick higher one
	# #         if other one is higher then
	#         if ints[i1 - 1] > ints[i1]
	#           pkmap[i1] = [nil]*5
	#         # else if this one is higher
	#         else
	#           pkmap[i1 - 1] = [nil]*5
	#           pkmap[i1] = x
	#           puts "#{pkmap[i1].inspect} \twith intensity: #{ints[i1]} "
	#         end
	#       end
	  pkmap[i1] = x
	  i1 += 1
	end
	# recheck  spectra for ++ daughter ions
	i1 = i2 =  0
	while i1 < mzs.length && i2 < ions.length
	  unless ions[i2][2]
	    i2 += 1
	    next
	  end
	  dff = mzs[i1] - ions[i2][2]
	  if dff > tol
	    # mass is too large
	    i2 += 1
	    next
	  elsif (dff.abs() <= tol )
	    pkmap[i1][0] = i2+1
	    pkmap[i1][2] = ions[i2][2]
	  end
	  i1 += 1
	end
	return pkmap
	end

  def to_s
    inspect
  end
  def inspect
    @seq.join("")
  end
end
