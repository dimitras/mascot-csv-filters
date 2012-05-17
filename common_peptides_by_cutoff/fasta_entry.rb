class FastaEntry
       attr_accessor :tag, :accno, :desc, :seq

       def initialize(tag, accno, desc, seq)
               @tag = tag
               @accno = accno
               @desc = desc
               @seq = seq
       end
end

