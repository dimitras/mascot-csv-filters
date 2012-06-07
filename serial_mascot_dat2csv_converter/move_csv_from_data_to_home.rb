# move csv files to home

Dir["/mnt/mascot/data/20120*"].each do |dir|
	if dir.split('/')[4].to_i >= 20120319
		dirname = dir.split('/')[4].to_s
		Dir["/mnt/mascot/data/#{dirname}/*.csv"].each do |file|
			filename = file.split('/')[5].split('.')[0].to_s
			system("cp /mnt/mascot/data/#{dirname}/#{filename}.csv ~/mascot_csvs/#{dirname}/#{filename}.csv")
		end
	end
end
