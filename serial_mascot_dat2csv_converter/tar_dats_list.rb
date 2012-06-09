# tar a specific list of dat files

filenames = ['F001526','F001504','F001517','F001583','F001584','F001585','F001586','F001587','F001594','F001515','F001590','F001514','F001525','F001493','F001516','F001506','F001505','F001507','F001509','F001508','F001593','F001521','F001589','F001513']
list = ''

Dir["/mnt/mascot/data/20120*"].each do |dir|
	if dir.split('/')[4].to_i >= 20120319
		dirname = dir.split('/')[4].to_s
		Dir["/mnt/mascot/data/#{dirname}/*.dat"].each do |file|
			filename = file.split('/')[5].split('.')[0].to_s
			if filenames.include?(filename)
				# system("cp /mnt/mascot/data/#{dirname}/#{filename}.dat /mnt/mascot/data/mascot_dats/#{filename}.dat") # too slow
				list << "/mnt/mascot/data/#{dirname}/#{filename}.dat "
			end
		end
	end
end
# system("tar -czf /mnt/mascot/data/mascot_dats/dats.tgz #{list}")
system("tar -cvzf /mnt/mascot/data/mascot_dats/dats.tgz #{list} &> /mnt/mascot/data/mascot_dats/dats.toc")

