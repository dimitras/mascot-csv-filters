# serial_dat_to_csv_converter

Dir["/mnt/mascot/data/20120*"].each do |dir|
	if dir.split('/')[4].to_i >= 20120319
		dirname = dir.split('/')[4].to_s
		system("mkdir ~/#{dirname}")
		Dir["/mnt/mascot/data/" + dirname + "/*.dat"].each do |file|
			filename = file.split('/')[5].split('.')[0].to_s
			cmd = "./export_dat_2.pl file=../data/#{dirname}/#{filename}.dat do_export=1 prot_hit_num=1 prot_acc=1 pep_query=1 pep_rank=1 pep_isbold=1 pep_isunique=1 pep_exp_mz=1 export_format=CSV _sigthreshold=0.05 _ignoreionsscorebelow=25 report=AUTO _server_mudpit_switch=0.000000001 _requireboldred=1 search_master=1 show_header=1 show_decoy=1 show_mods=1 show_params=1 show_format=1 protein_master=1 prot_score=1 prot_desc=1 prot_mass=1 prot_matches=1 peptide_master=1 pep_exp_mr=1 pep_exp_z=1 pep_calc_mr=1 pep_delta=1 pep_miss=1 pep_score=1 pep_expect=1 pep_seq=1 pep_var_mod=1 pep_scan_title=1 > ~/#{dirname}/#{filename}.csv"
			# pid = fork {
			# 	puts "JOB: #{dirname} > "
			#     exec cmd
			# }
			# Process.wait pid
			puts "JOB: #{dirname} > #{filename}"
			system(cmd)
		end
	end
end
