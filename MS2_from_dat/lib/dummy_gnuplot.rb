# $LOAD_PATH.unshift(File.expand_path('../../lib', __FILE__))
require 'rubygems'
require 'gnuplot'

x = [10,20,30,40,50]

Gnuplot.open do |gp|
  Gnuplot::Plot.new( gp ) do |plot|

    plot.xrange "[0:3000]"
    plot.title  "Sin Wave Example"
    plot.ylabel "x"
    plot.xlabel "sin(x)"

    plot.data << Gnuplot::DataSet.new('x') do |ds|
      ds.with = "lines"
      ds.linewidth = 4
    end

  end

end