import pandas as pd
import numpy as np
import gzip
from cmapPy.pandasGEXpress.parse import parse
from sys import argv
import os
import cmapPy.pandasGEXpress.write_gctx as wg
# load Cell information, and Inst_information
Cell_info=pd.read_csv('/home/jojo9103/data1/2020Year/drug_repositioning/LINCS/GSE70138_Broad_LINCS_cell_info_2017-04-28.txt.gz',sep='\t',compression='gzip')
inst_info=pd.read_csv('/home/jojo9103/data1/2020Year/drug_repositioning/LINCS/GSE70138_Broad_LINCS_inst_info_2017-03-06.txt.gz',sep='\t',compression='gzip')

# Gene info (match rid)
gene_info=pd.read_csv('/home/jojo9103/data1/2020Year/drug_repositioning/LINCS/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz',sep='\t',compression='gzip')
gene_if=np.array(gene_info['pr_gene_symbol'])
gene_if=pd.DataFrame(gene_if,index=gene_info['pr_gene_id'].astype('str'),columns=['Gene_Symbol'])


Cell_df=np.array([Cell_info['cell_id'],Cell_info['primary_site'],Cell_info['sample_type']])
Cell_df=Cell_df.T
Cell_df=pd.DataFrame(data=Cell_df,columns=['cell_id','primary_site','sample_type'])
# Cell id annotate primary_site
Cell_dir={Cell_df['cell_id'][i]:[Cell_df['cell_id'][i],Cell_df['primary_site'][i],Cell_df['sample_type'][i]] for i in range(0,Cell_df.shape[0])}

# inst_matrix 
inst_df=np.array([inst_info['inst_id'],inst_info['cell_id'],inst_info['pert_time'],inst_info['pert_id'],inst_info['pert_iname'],inst_info['pert_type']])
inst_df=inst_df.T
inst_df=pd.DataFrame(data=inst_df,columns=['inst_id','cell_id','pert_time','pert_id','pert_iname','pert_type'])

# inst_df annotate using Cell_dir
Cell_if=[Cell_dir[i] for i in inst_df['cell_id']]
Cell_if=pd.DataFrame(Cell_if,columns=['cell_id_1','primary_site','sample_type'])

# find matched Cell_id and primary_site
Summary_df=pd.concat([inst_df,Cell_if],axis=1)

#### Filter input argv[1]='Cell_id', argv[2]='primary_site', argv[3]='pert_time', argv[4]= 'pert_id'
#if argv[1]=='h' or argv[1]=='help'or argv[1]=='Help'or argv[1]=='HELP':


if argv[1]=='all' or argv[1]=='All' or argv[1]=='ALL' or argv[1]=='a' or argv[1]=='A':
	per_s=pd.crosstab(index=Summary_df['primary_site'],columns='count').index
	for ps in per_s: # primary_site
		Summary_df_1=Summary_df[(Summary_df['primary_site'].values==ps)]
		per_t=pd.crosstab(index=Summary_df_1['pert_time'],columns='count').index
		for pt in per_t: # pert_time
			Summary_df_2=Summary_df_1[(Summary_df_1['pert_time'].values==pt)]
			ce_id=pd.crosstab(index=Summary_df_1['cell_id'],columns='count').index
#			print(Summary_df_1.head)
			for ci in ce_id: # ce_id
				Summary_df_3=Summary_df_2[(Summary_df_2['cell_id'].values==ci)]
				per_id=pd.crosstab(index=Summary_df_2['pert_id'],columns='count').index
#				print(Summary_df_1.head)
				for pid in per_id:
#					f_p=str(ps)+'/'+str(pt)+'/'+pid+'/'+str(ci)+'.txt'
					f_p=str(ps)+'/'+str(pt)+'/'+pid+'_'+str(ci)
					if os.path.exists(f_p+'gctx'):
						continue
					Extract_id=Summary_df_3[(Summary_df_3['pert_id'].values==pid)]
					Extract_id.set_index('inst_id',inplace=True)
					print(f_p)
					# Level 3 data
					my_col_metadata_level3=parse('/home/jojo9103/data1/2020Year/drug_repositioning/LINCS/GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx',cid=Extract_id.index)
					my_col_metadata_level3.col_metadata_df=Extract_id
					my_col_metadata_level3.row_metadata_df=gene_if
					if my_col_metadata_level3.data_df.shape[1]>1:
						if not os.path.exists(ps):
							os.mkdir(ps)
						if not os.path.exists(ps+'/'+str(pt)):
							os.mkdir(str(ps)+'/'+str(pt))
						wg.write(my_col_metadata_level3,f_p)











