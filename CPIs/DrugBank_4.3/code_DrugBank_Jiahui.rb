File.open("result2.txt","w") do |fout|
	newdrug=0
	atc_code_list=[]
	target_list=[]
	drugbank_id_temp=""
	drug_name=""
	next_line=""
	drug_group=""
	a=0
File.open("drugbank.txt").each_with_index do |fin,index|
	#if index<=300
		#line=fin.split("\t",-1)
		#line.each do |lin|
			#lin.chomp!
		#end
		#fout.puts line.join("\t")
		fin.chomp!

		drugbank_id=/  <drugbank-id primary="true">DB[\d{1}]/.match(fin)
		if  drugbank_id!=nil#如果drugbank_id不nil,说明进入下一个药物，需输出之前药物atc等信息，否则还在同一个药物中，寻找atc
			#clear the old drug info
			if index!=0
				atc_code_list.each do |atc_each|
					fout.puts drugbank_id_temp+"\t"+drug_name.to_s+"\t"+atc_each+"\t"+drug_group+"\t"+target_list.join(",")
				end
				if atc_code_list==[]
					fout.puts drugbank_id_temp+"\t"+drug_name.to_s+"\t\t"+drug_group+"\t"+target_list.join(",")
				end
				#输出包含无atc编码的药物
			#fout.puts drugbank_id_temp+"\t"+atc_code_list.join("\t")
			#一个DB接多个atc code的输出形式
		    end
			drugbank_id_temp="DB0"+drugbank_id.post_match[0..3]
			atc_code_list=[]
			target_list=[]
			drug_name=""
		end
		#atc编码（一个药物对应多个atc code则分多行写）
		atc_code=/atc-code code+/.match(fin)
		if atc_code!=nil
			atc_code_list.push(fin[20..26])
		end

		#targets
		if /  <targets>/.match(fin)!=nil
			a=1
			next
		end
		if a==1
		    targets=/    <gene-name>+/.match(fin)	
		    if targets!=nil
			    target_gene_name=/<\/gene-name>/.match(targets.post_match)
			    if target_gene_name!=nil
			        target_list.push(target_gene_name.pre_match)
		        end
		    end
		    if /  <\/targets>/.match(fin)!=nil
			    a=0
		    end
	    end

	    #drug_name
		if /  <drugbank-id+/.match(fin)!=nil and drug_name==""
			next_line="drug_name"
			next
		end
		if next_line=="drug_name" 
		   next_line=""		
			drug_name_temp=/  <name>+/.match(fin)
		    if drug_name_temp!=nil
			    drug_name_post=/<\/name>/.match(drug_name_temp.post_match)
			    if drug_name_post!=nil
				    drug_name=drug_name_post.pre_match
			    end
		    end
		    
		end

		#drug group(批准/实验/...)
		drug_group_temp=/    <group>+/.match(fin)
		    if drug_group_temp!=nil
			    drug_group_post=/<\/group>/.match(drug_group_temp.post_match)
			    if drug_group_post!=nil
				    drug_group=drug_group_post.pre_match
			    end
		    end
		

	#end
end
end



