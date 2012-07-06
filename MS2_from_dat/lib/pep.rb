require 'rubygems'

# masses taken from the dat file
MODS = {
  'K' => 45.029395,
  'S' => 45.029395,
  'T' => 45.029395
}
# for non modified
# MODS = {
#   'K' => 42.010565,
#   'S' => 42.010565,
#   'T' => 42.010565
# }

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
  def initialize(seq,mods=[])
    @seq = seq.strip.upcase.split ""
    # mods [:pos]
    mods.each do |m|
      @seq[m - 1] += "*"
    end
    @iontable = []
    @bions = []
    @yions = []
    # build up the ions tables
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
  attr_accessor :iontable, :bions, :yions, :seq

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

  # Assigns the ion table to a given mass spectrum, walking through the m/z and ion series arrays to assign the ions to a mass peaks, given a tolerance
  # Returns a map of peaks that have ions assigned to them
  # http://174.129.8.134/mascot/help/results_help.html#PEP
  # example: http://174.129.8.134/mascot/cgi/peptide_view.pl?from=100&to=800&_label_all=0&file=..%2Fdata%2F20120323%2FF001507.dat&query=20052&hit=1&section=5&ave_thresh=38&_ignoreionsscorebelow=0&report=0&_sigthreshold=0.05&_msresflags=3138&_msresflags2=10&percolate=0&percolate_rt=0&tick1=100&tick_int=50&range=700&index=B4DHQ2&px=1
  def assign(mz,ions,tol=1.0)
    # ion is [idx,b,b++,b*,b0] or y 
    pkmap = Array.new(mz.length)
    i1 = i2 =  0
    puts "KNOWN/MEASURED MZs: #{mz.join(',')} (LENGTH:  #{mz.length}), \n\nCALCULATED masses/yIONs: #{ions.join(', ')} (LENGTH: #{ions.length})\n\n"
    # traverse through ios series looking for matches

    while i1 < mz.length && i2 < ions.length
      # don't have a mass to look for. next ion
      unless ions[i2][1]
        i2 += 1
        next
      end
      # x = [idx,+,++,*,0]
      x = [nil]*5
      dff = mz[i1] - ions[i2][1]
      if dff > tol
        # mass is too large, get next ion
        puts "ION:\ti1 = #{i1} , i2 = #{i2} => with mass diff: #{mz[i1]} - #{ions[i2][1]} (mass too large)"
        pkmap[i1] = x
        i2 += 1
        puts "GOTO\ti2 = #{i2}"
        next
      elsif dff <  0 - tol
        # mass is too small for the + ion
        # checking other ions
        # ++ # must check for these outside of this method
        # immonium ions
        if ions[i2][3] && ((mz[i1] - ions[i2][3]).abs <= tol) && !((mz[i1] - ions[i2][3]).abs > tol)
          x[0]=  i2 + 1
          x[3] = ions[i2][3]
        end
        # water loss
        if ions[i2][4] && !((mz[i1] - ions[i2][4]).abs <= tol) && !((mz[i1] - ions[i2][4]).abs > tol)
          x[0]=  i2 + 1
          x[4] = ions[i2][4]
        end
        # mass is too small, advance
        pkmap[i1] = x
        puts "ION:\ti1 = #{i1} , i2 = #{i2} => with mass diff: #{mz[i1]} - #{ions[i2][1]} (mass too small)"
        i1 += 1
        i2 = 0
        puts "GOTO\ti1 = #{i1}"
        next
      end
      # set the index
      x[0] = i2  + 1
      x[1] = ions[i2][1]
      puts "ION:\ti1 = #{i1} , i2 = #{i2} => with mass diff: #{mz[i1]} - #{ions[i2][1]}"
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
      i2 = 0
      puts "PKMAP@#{i1}: #{pkmap[i1]}"
      puts "GOTO\ti1 = #{i1},\ti2 = #{i2}"
    end
    # puts "\n-- Finished assigning yions --\n\n"
    # recheck  spectra for ++ daughter ions
    i1 = i2 =  0
    while i1 < mz.length && i2 < ions.length
      unless ions[i2][2]
        i2 += 1
        next
      end
      dff = mz[i1] - ions[i2][2]
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
    puts "\n\nPKMAP:\t#{pkmap.join(', ')} \t(LENGTH: #{pkmap.length})\n\n"
    return pkmap
  end

  def to_s 
    inspect
  end
  def inspect
    @seq.join("")
  end
end
