require 'rubygems'
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

AVGMASS = {
  'A'  => 71.08 ,
  'R'  => 156.19,
  'N'  => 114.1 ,
  'D'  => 115.09,
  'C'  => 103.14,
  'E'  => 129.12,
  'Q'  => 128.13,
  'G'  => 57.05 ,
  'H'  => 137.14,
  'I'  => 113.16,
  'L'  => 113.16,
  'K'  => 128.17,
  'M'  => 131.19,
  'F'  => 147.18,
  'P'  => 97.12 ,
  'S'  => 87.08 ,
  'T'  => 101.1 ,
  'U'  => 150.03,
  'W'  => 186.21,
  'Y'  => 163.18,
  'V'  => 99.13 ,
  "K*" => 128.17 + MODS['K'],
  "S*" => 87.08  + MODS['S'],
  "T*" => 101.1  + MODS['T']
}
MW = {
  'A' => 71.03712 ,
  'R' => 156.10112,
  'N' => 114.04293,
  'D' => 115.02695,
  'C' => 103.00919,
  'E' => 129.0426 ,
  'Q' => 128.05858,
  'G' => 57.02147 ,
  'H' => 137.05891,
  'I' => 113.08407,
  'L' => 113.08407,
  'K' => 128.09497,
  'M' => 131.04049,
  'F' => 147.06842,
  'P' => 97.05277 ,
  'S' => 87.03203 ,
  'T' => 101.04768,
  'U' => 150.95364,
  'W' => 186.07932,
  'Y' => 163.06333,
  'V' => 99.06842,
  "K*" => 128.09497 + MODS['K'],
  "S*" => 87.03203 + MODS['S'],
  "T*" => 101.04768 + MODS['T']
}

H = 1.00794
OH = 17.00706
H2O = 18.015
NH3 = 17.0305


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

  def assign(mz,ions,tol=1.0)
    # ion is [idx,b,b++,b*,b0] or y 
    pkmap = Array.new(mz.length)
    i1 = i2 =  0
    # traverse through ios series looking for matches

    while i1 < mz.length && i2 < ions.length
      # don't have a mass to look for. next ion
      ## BIG BUG, will never pick up ++ daughter ions! Please fix 
      unless ions[i2][1]
        i2 += 1
        next
      end
      # x = [idx,+,++,*,0]
      x = [nil]*5
      dff = mz[i1] - ions[i2][1]
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
        i1 += 1
        next
      end
      # set the index
      x[0] = i2  + 1
      x[1] = ions[i2][1]
      pkmap[i1] = x
      i1 += 1
    end
    # recheck  spectra for ++ daughter ions
    i1 = i2 =  0
    while i1 < mz.length && i2 < ions.length
      unless ions[i2][2]
        i2 += 1
        next
      end
      dff = mz[i1] - ions[i2][2]
#       puts([dff,mz[i1],ions[i2][2]].join(", "))
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

class MrmPicker
  require 'fastercsv'
  def self.mrm(pepindex,output_file="mrms.csv")
    pepindex = FasterCSV.read(pepindex,:headers=>true)
    xmls = {}
    # headers
    count = 0
    FasterCSV.open(output_file,'w') do |csv|
      csv << %w{ accession prot_descr peptide mod_peptide spectrum_id charge ret_time_minutes precursor_mz precursor_int ms2_mz ms2_int y_intensity_rank y_index y_ion_seq y-1 y y+2 l_precursor_mz l_y-1 l_y l_y+2  }
      pepindex.each do |r|
        count += 1 
        STDERR.puts count if count % 100 == 0 
        mods = r['mods'].split(":").map {|m| m.to_i }
        pep = Peptide.new(r["Gpeptide"],mods)
        # create the light peptide, only taking out K,L mods
        pep_light = Peptide.new(r["Gpeptide"],[])
        mods.each do |mpos|
          if pep_light.seq[mpos-1] == "C"
            pep_light.seq[mpos-1] = "C*"
          end
        end
        # puts pep,pep_light
        mrmions = []
        r['spectrum'] =~ /^(.+)\.(\d+)\.\d+\.(\d)$/
        frac = $1
        snum = $2.to_i
        chrg = $3.to_i
        unless(xmls[frac])
          xmls[frac] = Rampy::RampFile.open("#{frac}.mzXML")
        end
        x = xmls[frac]
        s = Rampy::Scan.new(x.scan(snum))
        # first check if we have suitable Y ions for the MRM 
        y = pep.assign(s.mz, pep.yions)
        # while I am at it grab the light Y ions table
        y_light = pep_light.yions
        # grab the MRM transitions for the heavy ions
        ranked_idx = [] 
        # grab legitimate y ions and rank them in intensity desc
        idxset = []
        yset = []
        mzset = []
        intset=[]
        y.each_index do |i|
          if y[i] && y[i][0] && !y[i][1].nil? && y[i][1] > 0
            idxset.push(i)
            yset.push(y[i][0])
            mzset.push(s.mz[i])
            intset.push(s.i[i])
          end
        end
        intset.each do 
          maxi = intset.index(intset.max())
          ranked_idx.push( idxset[maxi] )
          intset[maxi] = 0
        end
        # puts r['Gpeptide']
        # puts "HEAVY"
        # pep.yions.map{|ylion| puts ylion.join(", ")}
        # puts "LIGHT"
        # y_light.map{|ylion| puts ylion.join(", ")}

        ranked_idx.each_with_index do |i,ii|
          if(  y[i] &&
               !y[i][0].nil? && 1
               y[i][0] > 0  )
            yidx = y[i][0] - 1
            # puts pep.yions[yidx].join(", ")
            # puts y_light[yidx].join(", ")
            csv <<  [ r['Gprotein'],
              r["Gprotein_description"],
              r['Gpeptide'],
              pep.to_s,
              r['spectrum'], 
              r['charge'],
              r['ret_time_sec'].to_f / 60 ,
              r['ms1_mz'],
              r['ms1_int'],
              s.mz[i],
              s.i[i],
              ii, # the mz intensity rank
              y[i][0],
              pep.yions[yidx][0], # y ion sequence
              y[i][1] - H ,
              y[i][1],
              y[i][1] + (H * 2),
              pep_light.mw(2) / 2,
              y_light[yidx][1] - H,
              y_light[yidx][1],
              y_light[yidx][1] + (H * 2) 
              ]
          end # end if
        end # end each_with_index
      end # end pepindex.each
    end
  end
end
class Plotter
  # require 'gnuplot'
  def self.plot(pep,mz,inten,fn)
    o = File.open(fn,'w')
    b = pep.assign(mz,pep.bions)
    y = pep.assign(mz,pep.yions)
    o.puts(%w{ mz inten b_idx b bb bn b0 y_idx y yy yn y0}.join(","))
    mz.each_index do |i|
      o.puts([mz[i],inten[i],b[i],y[i]].join(","))
    end
  end
end

if __FILE__ == $0 
  require 'yaml'
  require 'rampy'
  # puts p = Peptide.new("VFQQVAQASK",[10])
  puts p = Peptide.new("AITGFDDPFSGK",[12])
  puts p.yions.map {|e| e.join ", "}
  r = Rampy::RampFile.open("END2334A6MS2-40-15007.mzXML")
  s = r.get_peaks(r.scan(5554))
  Plotter.plot(p,s[0],s[1],'plotdata_desmo.csv')
  # r = Rampy::RampFile.open("SCX-A3.mzXML")
  # s =  r.get_peaks(r.scan(2952))
  # Plotter.plot(p,s[0],s[1],'plotdata.csv')
  # test mpicker 
  # MPicker.mrm("combined.csv","mrms.csv")
end