Dir["/mnt/mascot/data/20120*"].each do |dir|
	if dir.split('/')[4].to_i >= 20120319
		dirname = dir.split('/')[4].to_s
		system("mkdir ~/mascot_csvs/#{dirname}")
	end
end