# find files that failed to export

Dir["/mnt/mascot/data/20120*"].each do |dir|
	if dir.split('/')[4].to_i >= 20120319
		dirname = dir.split('/')[4].to_s
		Dir["/mnt/mascot/data/#{dirname}/*.csv"].each do |file|
			puts file
			line = IO.readlines(file)[1]
			puts line # if line is header then is correct
			filename = file.split('/')[5].split('.')[0].to_s
		end
	end
end
