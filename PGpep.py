#!/usr/bin/env python
# coding: utf-8

# In[239]:





# In[240]:


import gc
import glob
import os
import subprocess
import inspect
import time
import sys
# import ahocorapy
ahocorapy_install = 'pip3.10 install ahocorapy'
r = subprocess.Popen(ahocorapy_install, shell=True).wait()
if r == 1:
    print("install ahocorapy failed!")
from ahocorapy.keywordtree import KeywordTree

whole_input = sys.argv[1:]
input_param = " ".join(whole_input)


# In[ ]:

def main_one(input_param):
    #mode o p_name patient_name t_id_c 0 g_id_c 3 fpkm_c 9 fpkm_p fpkm_path bam_p bam_path f_s_db_p 1st_2nd_db_path c_db_p composite_db_path u_db_p uniprot_db_path 
    transcript_id_column_num = int(input_param.split('t_id_c ')[1].split(' ')[0])
    gene_id_column_num = int(input_param.split('g_id_c ')[1].split(' ')[0])
    fpkm_column_num = int(input_param.split('fpkm_c ')[1].split(' ')[0])
    fpkm_path = str(input_param.split('fpkm_p ')[1].split(' ')[0])
    bam_path = str(input_param.split('bam_p ')[1].split(' ')[0])
    first_2nd_db_path = str(input_param.split('f_s_db_p ')[1].split(' ')[0])
    composite_db_path = str(input_param.split('c_db_p ')[1].split(' ')[0])
    uniprot_db_path = str(input_param.split('u_db_p ')[1].split(' ')[0])
    patient_name = str(input_param.split('p_name ')[1].split(' ')[0])
    
    
    #1 : peaks 결과 후처리 ( 1. db search 결과와 겹치는 scan 제거 2. 기존 db에 있는 서열들과 겹치는 서열 제거 )
    post_processing_dn()
    
    #2 : fpkm 계산 결과로 환자별 gtf 생성 ^
    make_gtf_every_patients(transcript_id_column_num, gene_id_column_num, fpkm_column_num) 
    
    #3 : 환자별 actg 실행 ^
    gogo_actg_one(patient_name) 
    
    #4 : human gene table 생성
    global fasta_dict 
    fasta_dict = make_fasta_dict() 
    
    #5 : actg 결과와 rna-seq data mapping  ^ (samtools 이용하여 actg로 매핑된 peptide의 유전자 영역에 rna read 확인)
    start_rna_evidence_one(patient_name, bam_path)
    
    #6  ^
    after_rna_mapped_es_one(patient_name)
    #7  ^
    after_rna_mapped_ee_ri_one(patient_name)
    
    #8  ^
    see_junction_read_one(patient_name, bam_path)
    
    #9  ^
    start_novel_transcript_evidence_one(patient_name)
    
    #10  ^
    start_novel_transcript_evidence_result_one(patient_name)
    
    
def main_whole(input_param):
    
    #0 : 명령어 받기 mode(단일, 전체), fpkm 경로, actg 경로, bam 경로
    #mode w t_id_c 0 g_id_c 3 fpkm_c 9
    #ex) python PGpep.py mode w t_id_c 0 g_id_c 3 fpkm_c 9
    transcript_id_column_num = int(input_param.split('t_id_c ')[1].split(' ')[0])
    gene_id_column_num = int(input_param.split('g_id_c ')[1].split(' ')[0])
    fpkm_column_num = int(input_param.split('fpkm_c ')[1].split(' ')[0])
    
    #1 : peaks 결과 후처리 ( 1. db search 결과와 겹치는 scan 제거 2. 기존 db에 있는 서열들과 겹치는 서열 제거 )
    post_processing_dn()
    
    #2 : fpkm 계산 결과로 환자별 gtf 생성
    make_gtf_every_patients(transcript_id_column_num, gene_id_column_num, fpkm_column_num) 
    
    #3 : 환자별 actg 실행
    gogo_actg() 
    
    #4 : human gene table 생성
    global fasta_dict 
    fasta_dict = make_fasta_dict() 
    
    #5 : actg 결과와 rna-seq data mapping ( samtools 이용하여 actg로 매핑된 peptide의 유전자 영역에 rna read 확인)
    start_rna_evidence() 
    
    #6
    after_rna_mapped_es()
    #7
    after_rna_mapped_ee_ri()
    
    #8
    see_junction_read()
    
    #9
    start_novel_transcript_evidence()
    
    #10
    start_novel_transcript_evidence_result()


# In[53]:


plus_stop_codon = ['TGA', 'TAA' ,'TAG']
minus_stop_codon = ['TCA', 'TTA', 'CTA']
plus_donor = 'GT'
plus_acceptor = 'AG'
minus_donor = 'AC'
minus_acceptor = 'CT'


# In[54]:


"""------------------PEAKS 전처리 부분 시작------------------"""


# In[55]:


"""
description:

arguments : None

output : None
"""
def post_processing_dn():
    # peaks 결과 전처리
    dn_dict = {}
    file_list = os.listdir('./DENOVO//')
    replace_list = '[+-1234567890\\(\\).]*'
    dn = open('./DENOVO//'+file_list[0])
    dn_ori = dn.readlines()
    for line in dn_ori[1:]:
        peptide = line.split(",")[3].replace("I","L")
        peptide = peptide.sub(replace_list, "", peptide)
        fraction_scan = line.split(",")[4]
        if dn_dict.get(fraction_scan) == None:
            dn_dict[fraction_scan] = peptide
    db_1st_2nd_dict = {}
    # peaks 결과 db 서치 fraction:scan 겹치는 거 삭제
    file_list = os.listdir('./DB/1st_2nd_DB_search')
    for i in file_list:
        set_num = int(i.split("_")[1].split("set")[0])
        search_result = os.listdir('./DB/1st_2nd_DB_search/'+i)
        folder_path = './DB/1st_2nd_DB_search/'+i+'/'
        for file in search_result:
            file_now = open(folder_path+file)
            file_now_ori = file_now.readlines()
            db_1st_2nd_dict = get_scan_from_1st_2nd(file_now_ori,set_num,db_1st_2nd_dict)
    for i in db_1st_2nd_dict.keys():
        if dn_dict.get(i) != None:
            del dn_dict[i]
    # peaks 결과 uniprot db와 composite db 대조 후 삭제
    db_pep_dict = {}
    db_pep_dict = read_make_CompositeDB_dict(db_pep_dict)
    db_pep_dict = read_make_uniprot_DB_dict(db_pep_dict)
    filter_over7_dn_dict = {}
    filter_over7_dn_dict = make_dict_switch_key_value(dn_dict)
    filter_over7_no_db_dict = {}
    filter_over7_no_db_dict = do_some_ahocorapy(filter_over7_dn_dict, db_pep_dict)
    make_actg_input_pep_list(filter_over7_no_db_dict)
    del dn_dict
    del db_1st_2nd_dict
    del db_pep_dict
    throw_garbage()
    
def post_processing_dn_one(first_2nd_db_path, composite_db_path, uniprot_db_path):
    # peaks 결과 전처리
    dn_dict = {}
    file_list = os.listdir('./DENOVO//')
    replace_list = '[+-1234567890\\(\\).]*'
    dn = open('./DENOVO//'+file_list[0])
    dn_ori = dn.readlines()
    for line in dn_ori[1:]:
        peptide = line.split(",")[3].replace("I","L")
        peptide = peptide.sub(replace_list, "", peptide)
        fraction_scan = line.split(",")[4]
        if dn_dict.get(fraction_scan) == None:
            dn_dict[fraction_scan] = peptide
    db_1st_2nd_dict = {}
    # peaks 결과 db 서치 fraction:scan 겹치는 거 삭제
    file_list = os.listdir('./DB/1st_2nd_DB_search')
    file_now = open(first_2nd_db_path)
    file_now_ori = file_now.readlines()
    db_1st_2nd_dict = get_scan_from_1st_2nd(file_now_ori,set_num,db_1st_2nd_dict)
    for i in db_1st_2nd_dict.keys():
        if dn_dict.get(i) != None:
            del dn_dict[i]
    # peaks 결과 uniprot db와 composite db 대조 후 삭제
    db_pep_dict = {}
    db_pep_dict = read_make_CompositeDB_dict_one(composite_db_path, db_pep_dict)
    db_pep_dict = read_make_uniprot_DB_dict_one(uniprot_db_path, db_pep_dict)
    filter_over7_dn_dict = {}
    filter_over7_dn_dict = make_dict_switch_key_value(dn_dict)
    filter_over7_no_db_dict = {}
    filter_over7_no_db_dict = do_some_ahocorapy(filter_over7_dn_dict, db_pep_dict)
    make_actg_input_pep_list(filter_over7_no_db_dict)
    del dn_dict
    del db_1st_2nd_dict
    del db_pep_dict
    throw_garbage()

# In[56]:


"""
description:

arguments:
    filter_over7_no_db_dict

output : 
"""
def make_actg_input_pep_list(filter_over7_no_db_dict):
    write_list = []
    set_pep_dict_list = []
    
    file_list = os.listdir('./ACTG/')
    file_len =  len(file_list)
    
    for i in range(1,file_len+1):
        set_pep_dict = {}
        set_pep_dict_list.append(set_pep_dict)
        pep_write = open("./ACTG/set"+str(i)+"/"+str(i)+"set_dn_pep_1014.txt","wt")
        write_list.append(pep_write)
    for i in filter_over7_no_db_dict.keys():
        pep = i
        fraction_info = filter_over7_no_db_dict.get(i).split("//")
        for fraction in fraction_info:
            set_info = int(fraction.split(":")[0].split("F")[1])
            now_set = set_info//24+1
            if set_info%24 == 0:
                now_set = now_set-1
            set_pep_dict_list[now_set-1][pep] = 1
    for i in range(len(set_pep_dict_list)):
        for pep in set_pep_dict_list[i].keys():
            write_list[i].write(pep+"\n")
#             write_list[now_set-1].write(pep+":"+fraction+"\n")
    for i in write_list:
        i.close()
        


# In[57]:


"""
description:

arguments:
    filter_over7_dn_dict
    db_pep_dict
    
output : 
    filter_over7_dn_dict
"""
def do_some_ahocorapy(filter_over7_dn_dict, db_pep_dict):
    dn_tree = KeywordTree(case_insensitive=True)
    print('tree make start')
    for pep in filter_over7_dn_dict.keys():
        dn_tree.add(pep)
    dn_tree.finalize()
    print('tree make finish')
    erase = {}
    line_count = 0
    print('search start')
    for i in db_pep_dict.keys():
        results = dn_tree.search_all(i)
        line_count+=1
        if line_count%100000==0:
            print(line_count)
        if results != None:
            for result in results:
                result = str(result).split("'")[1]
                erase[result] = 1
    print('erased : '+str(len(erase)))
    print("before erase dn pep num : "+str(len(filter_over7_dn_dict)))
    for i in erase.keys():
        del filter_over7_dn_dict[i]
    print("after erase dn pep num : "+str(len(filter_over7_dn_dict)))
    del dn_tree
    return filter_over7_dn_dict


# In[58]:


"""
description:

arguments:
    dn_dict

output : 
    reverse_dn_dict
"""
def make_dict_switch_key_value(dn_dict):
    reverse_dn_dict = {}
    print(len(dn_dict))
    under8 = 0
    for key in dn_dict.keys():
        pep = dn_dict.get(key)
        if len(pep)>7:
            if reverse_dn_dict.get(pep) is None:
                reverse_dn_dict[pep] = key
            else:
                reverse_dn_dict[pep] = reverse_dn_dict.get(pep)+"//"+key
        else:
            under8+=1
    print(len(reverse_dn_dict))
    print('under8 : '+str(under8))
    return reverse_dn_dict


# In[59]:


"""
description:

arguments:
    db_pep_dict

output : 
    db_pep_dict
"""
def read_make_uniprot_DB_dict(db_pep_dict):
    file_list = os.listdir('./DB/uniprot')
    for i in file_list:
        if i.split(".")[-1] == 'fasta':
            file_now = open('./DB/uniprot/'+i)
            file_now_ori = file_now.readlines()
            db_pep_dict = make_DB_with_fasta(file_now_ori,db_pep_dict)
    return db_pep_dict
    
def read_make_uniprot_DB_dict_one(uniprot_db_path, db_pep_dict):    
    file_now = open(uniprot_db_path)
    file_now_ori = file_now.readlines()
    db_pep_dict = make_DB_with_fasta(file_now_ori,db_pep_dict)
    return db_pep_dict


# In[60]:


"""
description:

arguments:
    db_pep_dict

output : 
    db_pep_dict
"""
def read_make_CompositeDB_dict(db_pep_dict):
    file_list = os.listdir('./DB/CompositeDB')
    for i in file_list:
        search_result = os.listdir('./DB/CompositeDB/'+i)
        folder_path = './DB/CompositeDB/'+i+'/'
        for file in search_result:
            if file.split(".")[-1] == 'fasta':
                file_now = open(folder_path+file)
                file_now_ori = file_now.readlines()
                db_pep_dict = make_DB_with_fasta(file_now_ori,db_pep_dict)
    return db_pep_dict
    
def read_make_CompositeDB_dict_one(composite_db_path, db_pep_dict):    
    file_now = open(composite_db_path)
    file_now_ori = file_now.readlines()
    db_pep_dict = make_DB_with_fasta(file_now_ori,db_pep_dict)
    return db_pep_dict


# In[61]:


"""
description:

arguments:
    ori_file
    db_pep_dict

output : 
    db_pep_dict
"""
def make_DB_with_fasta(ori_file,db_pep_dict):
#     uni_dict = {}
    head = ''
    pep = []
    count = 0
    for line in ori_file[:]:
        if ">" in line:
            if count == 0:
                count+=1
                head = line.replace("\n","")
            else:
                pep_seq = "".join(pep)
                pep_seq = pep_seq.replace("I","L")
                db_pep_dict[pep_seq] = head
                count+=1
                pep = []
                pep_seq = ''
                head = line.replace("\n","")
        else:
            pep.append(line.replace("\n",""))
            
    pep_seq = "".join(pep)
    pep_seq = pep_seq.replace("I","L")
    db_pep_dict[pep_seq] = head
    print(len(db_pep_dict))
    return db_pep_dict


# In[62]:


"""
description:

arguments:
    ori_file
    set_num
    db_1st_2nd_dict
    
output : 
    db_1st_2nd_dict
"""
def get_scan_from_1st_2nd(ori_file,set_num,db_1st_2nd_dict):
    c = '[+-1234567890.]*'
    for i in ori_file[1:]:
        fraction = i.split("\t")[0].split("_")[-1]
        scan = i.split("\t")[-1].replace("\n","")
        if fraction[0] == '0':
            now_fraction = str(int(fraction[1])+(int(set_num)-1)*24)
            fraction = 'F'+now_fraction+":"
        else:
            now_fraction = str(int(fraction[1])+(int(set_num)-1)*24)
            fraction = 'F'+now_fraction+":"
        psm = fraction+scan
        pep = i.split("\t")[8][2:-2].replace("I","L")
        pep = pep.sub(c, "", pep)
        db_1st_2nd_dict[psm] = pep
    return db_1st_2nd_dict


# In[63]:


"""------------------PEAKS 전처리 부분 끝------------------"""


# In[64]:


"""------------------ACTG 시작------------------"""


# In[65]:


"""
description:

arguments:
    gtf_path
    patient
    const_ori
    
output : 
    ser_name
"""
def make_const_param(gtf_path,patient,const_ori):
    const_temp = open("temp_const_params.xml","wt")
    for i in const_ori:
        if '<Input format="GTF" type="transcriptome">' in i:
#             print("what_happen?")
            const_temp.write('			<Input format="GTF" type="transcriptome">'+gtf_path+'</Input>\n')
        elif '			<Output format="ser" type="graphFile">result_for_0930.ser</Output>' in i:
#             print("what_happen??")
            const_temp.write('			<Output format="ser" type="graphFile">'+patient+'.ser</Output>\n')
        else:
            const_temp.write(i)
    const_temp.close()
#     print('yes')
    ser_name = patient+".ser"
    return ser_name
    


# In[66]:


"""
description:

arguments:
    pep_path
    patient
    mapping_ori
    ser_name
    output_path
    
output : 
"""
def make_mapping_param(pep_path,patient,mapping_ori,ser_name,output_path):
    temp_mapping = open("temp_mapping_params.xml","wt")
    for i in mapping_ori:
        if '			<Input format="list" type="peptideList">uniq_db_pep_for_novel_event_0503.txt</Input>' in i:
            temp_mapping.write('			<Input format="list" type="peptideList">'+pep_path+'</Input>\n')
        elif '			<Input format="ser" type="graphFile">result_for_only_novel.ser</Input>' in i:
            temp_mapping.write('			<Input format="ser" type="graphFile">'+ser_name+'</Input>\n')
        elif '			<Output type="outputPath">output</Output>' in i:
            temp_mapping.write('			<Output type="outputPath">'+output_path+'</Output>\n')
        else:
            temp_mapping.write(i)
    temp_mapping.close()


# In[67]:


"""
description:

arguments:

output : 
"""
def gogo_actg():
    file_list = os.listdir('./ACTG/')
    const = open("const_params.xml")
    const_ori = const.readlines()
    mapping = open("mapping_params.xml")
    mapping_ori = mapping.readlines()
    for set_num in file_list:
        patient_list = os.listdir('./ACTG/'+set_num)
        print(set_num)
        pep_path = ''
        for pep_list in patient_list:
            if len(pep_list.split(".")) == 2:
                pep_path = 'ACTG/'+set_num+"/"+pep_list
        for patient in patient_list:
            if len(patient.split(".")) == 1:
                print(patient)
                gtf_path = 'ACTG/'+set_num+"/"+patient+"/GTF"
                make_dir_for_actg(patient,set_num)
                output_path = 'ACTG/'+set_num+"/"+patient+"/output"
                ser_name = make_const_param(gtf_path,patient,const_ori)
                actg_ser_input = 'java -Xmx8G -Xss16M -jar ACTG_construction.jar temp_const_params.xml'
                r = subprocess.Popen(actg_ser_input, shell=True).wait()
                if r == 1: 
                    print("making ser failed")
                make_mapping_param(pep_path,patient,mapping_ori,ser_name,output_path)
                actg_mapping_input = 'java -Xmx8G -Xss16M -jar ACTG_mapping.jar temp_mapping_params.xml'
                r = subprocess.Popen(actg_mapping_input, shell=True).wait()
                if r == 1: 
                    print("ACTG mapping failed")
                remove_ser_input = 'rm '+ser_name
                r = subprocess.Popen(remove_ser_input, shell=True).wait()
                if r == 1: 
                    print("deleting ser failed")
                
                
def gogo_actg_one(patient_name):
    file_list = os.listdir('./ACTG/')
    const = open("const_params.xml")
    const_ori = const.readlines()
    mapping = open("mapping_params.xml")
    mapping_ori = mapping.readlines()
    for set_num in file_list:
        patient_list = os.listdir('./ACTG/'+set_num)
        print(set_num)
        pep_path = ''
        for pep_list in patient_list:
            if len(pep_list.split(".")) == 2:
                pep_path = 'ACTG/'+set_num+"/"+pep_list
        for patient in patient_list:
            if len(patient.split(".")) == 1:
                if patient == patient_name:
                    gtf_path = 'ACTG/'+set_num+"/"+patient+"/GTF"
                    make_dir_for_actg(patient,set_num)
                    output_path = 'ACTG/'+set_num+"/"+patient+"/output"
                    ser_name = make_const_param(gtf_path,patient,const_ori)
                    actg_ser_input = 'java -Xmx8G -Xss16M -jar ACTG_construction.jar temp_const_params.xml'
                    r = subprocess.Popen(actg_ser_input, shell=True).wait()
                    if r == 1: 
                        print("making ser failed")
                    make_mapping_param(pep_path,patient,mapping_ori,ser_name,output_path)
                    actg_mapping_input = 'java -Xmx8G -Xss16M -jar ACTG_mapping.jar temp_mapping_params.xml'
                    r = subprocess.Popen(actg_mapping_input, shell=True).wait()
                    if r == 1: 
                        print("ACTG mapping failed")
                    remove_ser_input = 'rm '+ser_name
                    r = subprocess.Popen(remove_ser_input, shell=True).wait()
                    if r == 1: 
                        print("deleting ser failed")
                


# In[68]:


"""
description:

arguments:
    patient_num
    set_num

output : 
"""
def make_dir_for_actg(patient_num,set_num):
    make_dir_input = 'mkdir -m 777 -p ACTG/'+set_num+'/'+patient_num+'/output'
    r = subprocess.Popen(make_dir_input, shell=True).wait()
    if r == 1: 
        print("making output folder failed")


# In[69]:


"""
description:

arguments:
    patient_num
    set_num

output : 
"""
def make_dir_gtf(patient_num, set_num):
    make_dir_input = 'mkdir -m 777 -p ACTG/'+set_num+'/'+patient_num+'/GTF'
    r = subprocess.Popen(make_dir_input, shell=True).wait()
    if r == 1: 
        print("making patient folder failed")


# In[70]:


"""
description:

arguments:
    gtf_ori

output : 
"""
def read_make_fpkm(gtf_ori,transcript_id_column_num, gene_id_column_num, fpkm_column_num):
    patient_name_list = []
    file_list = os.listdir('./FPKM/')
    fpkm_dict= {}
    for i in file_list:
        set_num = int(i.split('set')[1])
        print(set_num)
        search_result = os.listdir('./FPKM/'+i)
        folder_path = './FPKM/'+i+'/'
        for folder in search_result:
            now_path = folder_path+'/'+folder
            files = os.listdir(now_path)
            for file in files:
                if file.split(".")[1] == 'isoforms':
                    patient_name = file.split("_RSq")[0]
                    patient_name_list.append(patient_name)
                    fpkm_info = open(now_path+"/"+file)
                    fpkm_info_ori = fpkm_info.readlines()
                    make_dir_gtf(patient_name, i)
                    now_new_gtf = open('./ACTG/'+i+'/'+patient_name+'/GTF/'+patient_name+"_grch37.gtf",'wt')
                    p_trans_dict = {}
                    p_gene_dict = {}
                    for line_fpkm in fpkm_info_ori[1:]:
                        fpkm = float(line_fpkm.split("\t")[fpkm_column_num])
                        if fpkm >0 :
                            trans = line_fpkm.split("\t")[transcript_id_column_num]
                            gene = line_fpkm.split("\t")[gene_id_column_num]
                            p_trans_dict[trans] = 1
                            p_gene_dict[gene] = 1
                    for line_gtf in gtf_ori[:5]:
                        now_new_gtf.write(line_gtf)
                    for line_gtf in gtf_ori[5:]:
                        trans_type = line_gtf.split("\t")[2]
                        trans_info = line_gtf.split("\t")[8]
                        if trans_type == 'gene':
                            gene = line_gtf.split("\t")[8].split("\"")[1]
                            if p_gene_dict.get(gene) != None:
                                now_new_gtf.write(line_gtf)
                        if "; transcript_id \"" in trans_info:
                            trans = line_gtf.split("\t")[8].split("\"")[3]
                            if p_trans_dict.get(trans) != None:
                                now_new_gtf.write(line_gtf)
                    print(now_new_gtf)
                    now_new_gtf.close()
                    
def read_make_fpkm_one(fpkm_path, gtf_ori,transcript_id_column_num, gene_id_column_num, fpkm_column_num):
    patient_name_list = []
    file_list = os.listdir('./FPKM/')
    fpkm_dict= {}
    
    patient_name = fpkm_path.split("_RSq")[0].split("/")[-1]
    patient_name_list.append(patient_name)
    fpkm_info = open(fpkm_path)
    fpkm_info_ori = fpkm_info.readlines()
    make_dir_gtf(patient_name, i)
    now_new_gtf = open('./ACTG/'+i+'/'+patient_name+'/GTF/'+patient_name+"_grch37.gtf",'wt')
    p_trans_dict = {}
    p_gene_dict = {}
    for line_fpkm in fpkm_info_ori[1:]:
        fpkm = float(line_fpkm.split("\t")[fpkm_column_num])
        if fpkm >0 :
            trans = line_fpkm.split("\t")[transcript_id_column_num]
            gene = line_fpkm.split("\t")[gene_id_column_num]
            p_trans_dict[trans] = 1
            p_gene_dict[gene] = 1
    for line_gtf in gtf_ori[:5]:
        now_new_gtf.write(line_gtf)
    for line_gtf in gtf_ori[5:]:
        trans_type = line_gtf.split("\t")[2]
        trans_info = line_gtf.split("\t")[8]
        if trans_type == 'gene':
            gene = line_gtf.split("\t")[8].split("\"")[1]
            if p_gene_dict.get(gene) != None:
                now_new_gtf.write(line_gtf)
        if "; transcript_id \"" in trans_info:
            trans = line_gtf.split("\t")[8].split("\"")[3]
            if p_trans_dict.get(trans) != None:
                now_new_gtf.write(line_gtf)
    print(now_new_gtf)
    now_new_gtf.close()


# In[71]:


"""
description:

arguments:

output : 
"""
def make_gtf_every_patients(transcript_id_column_num, gene_id_column_num, fpkm_column_num):
    file_list = os.listdir('./GTF')
    for i in file_list:
        print(i)
        if i.split(".")[-1] == 'gtf':
            gtf = open("./GTF/"+i)
            gtf_ori = gtf.readlines()
            print(len(gtf_ori))
            read_make_fpkm(gtf_ori,transcript_id_column_num, gene_id_column_num, fpkm_column_num)
            
def make_gtf_every_patients_one(fpkm_path, transcript_id_column_num, gene_id_column_num, fpkm_column_num):
    file_list = os.listdir('./GTF')
    for i in file_list:
        print(i)
        if i.split(".")[-1] == 'gtf':
            gtf = open("./GTF/"+i)
            gtf_ori = gtf.readlines()
            print(len(gtf_ori))
            read_make_fpkm(gtf_ori,transcript_id_column_num, gene_id_column_num, fpkm_column_num)


# In[72]:


"""------------------ACTG 끝------------------"""


# In[73]:


"""------------------Samtools RNA 매핑 시작------------------"""


# In[74]:


"""
description:
    make dictionary with gtf file

arguments:
    gtf_ori : gtf file read by .readlines()

output : 
    gtf_dict : dictionary that has gene information
        key : gene name
        value : all lines about specific gene
"""
def make_gtf_dict(gtf_ori):
    gtf_dict = {}
    for line in gtf_ori[5:]:
        name = line.split("\t")[8].split("\"")[1]
        if gtf_dict.get(name) == None:
            gtf_dict[name] = line
        else:
            gtf_dict[name] = gtf_dict.get(name)+line
    return gtf_dict


# In[75]:


"""
description:
    make dictionary with gtf file

arguments:
    gtf_ori : gtf file read by .readlines()

output : 
    gtf_dict : dictionary that has gene information
        key : gene name
        value : all lines about specific gene
"""
def make_gtf_dict_intron(gtf_ori):
    gtf_dict = {}
    for line in gtf_ori[5:]:
        if line.split("\t")[1] == 'protein_coding':
            name = line.split("\t")[8].split("\"")[1]
            if gtf_dict.get(name) == None:
                if line.split("\t")[2] == 'gene':
                    gtf_dict[name] = line
            else:
                gtf_dict[name] = gtf_dict.get(name)+line
    print(len(gtf_dict))
    return gtf_dict


# In[76]:


"""
description:
    concatenate flat file and gff file from ACTG

arguments:
    flat_ori : flat file from ACTG has peptide_seq and evnet case(exon-extension, exon-skipping, etc)
    ex) 1	EVFLEAGR	GAGGTGTTCTTGGAGGCTGGACGT	ENSG00000243620	exon-extension
    gff_ori : gff file form ACTG has mapping position information of peptide
    ex) chr3	ACTG	exon	146898278	146898301	.	-	2	ID=1

output:
    flat_with_gff_dict : dictionary that has flat and gff information
        key : concatenated one line of flat and gff
        value : 1
"""
def make_flat_with_gff_dict_actg(flat_ori, gff_ori):
    flat_dict = {}
    now_dict = {}
    for k in flat_ori[1:]:
        flat_id = k.split("\t")[0]
        flat_dict[flat_id] = k
    for k in gff_ori[:]:
        gff_id = k.split("\t")[-1].split("=")[1].replace("\n","")
        now_dict[k.replace("\n","")+"\t"+flat_dict.get(gff_id)] = 1
    return now_dict


# In[77]:


"""
description:
    concatenate flat file and gff file from ACTG

arguments:
    flat_ori : flat file from ACTG has peptide_seq and evnet case(exon-extension, exon-skipping, etc)
    ex) 1	EVFLEAGR	GAGGTGTTCTTGGAGGCTGGACGT	ENSG00000243620	exon-extension
    gff_ori : gff file form ACTG has mapping position information of peptide
    ex) chr3	ACTG	exon	146898278	146898301	.	-	2	ID=1

output:
    flat_with_gff_dict : dictionary that has flat and gff information
        key : concatenated one line of flat and gff
        value : 1
"""
def make_flat_with_gff_dict(flat_ori, gff_ori):
    flat_dict = {}
    flat_with_gff_dict = {}
    for k in flat_ori[1:]:
        flat_id = k.split("\t")[0]
        flat_dict[flat_id] = k
    for k in gff_ori[:]:
        gff_id = k.split("\t")[-1].split("=")[1].replace("\n","")
        flat_with_gff_dict[k.replace("\n","")+"\t"+flat_dict.get(gff_id)] = 1
    return flat_with_gff_dict


# In[78]:


"""
description:
    gtf 정보와 ACTG 결과물을 통해 펩타이드가 매핑된 유전자가 인트론인지 엑손인지 구분
    
arguments:
    line : one key of flat_with_gff_dict
    gtf_dict : dictionary that has gtf information
        key : gene name
        value : all lines about specific gene
    
output : 
    map_to_return : list that shows given genomic region is intron(=10) or exon(>10)
    ex) [10, 10, 10, 10, 11, 11, 11, 11]
"""
def how_exon(line,gtf_dict):
    extension_start = float(line.split("\t")[3])
    extension_end = float(line.split("\t")[4])
    name = line.split("\t")[12]
#     print('---')
#     print(name)
#     print('---')
    info = ''
    gene_start = 0
    gene_end = 0
    if gtf_dict.get(name) != None:
        info = gtf_dict.get(name)
        info_list = info.split("\n")
        gene_start = float(info.split("\n")[0].split("\t")[3])
        gene_end = float(info.split("\n")[0].split("\t")[4])
        map_size = gene_end - gene_start + 1
        gene_map = [0 for i in range(int(map_size))]
        for i in info_list[1:-1]:
            exon_start = float(0)
            exon_end = float(0)
            if i.split("\t")[2] == 'exon':
                exon_start = int(float(i.split("\t")[3])-gene_start)
                exon_end = int(float(i.split("\t")[4])-gene_start)
                for j in range(exon_start,exon_end+1):
                    gene_map[j] += 1
        for h in range(int(extension_start-gene_start),int(extension_end-gene_start+1)):
            gene_map[h] += 10
        map_to_return = []
        for h in range(int(extension_start-gene_start),int(extension_end-gene_start+1)):
            map_to_return.append(gene_map[h])
        if len(set(map_to_return)) ==1:
            return map_to_return
        else: 
            transcript_list = []
            transcript_first = 0
            max_index = -1
            max_sum = -1
            gene_map = [0 for i in range(int(map_size))]
            for i in info_list[1:-1]:
                exon_start = float(0)
                exon_end = float(0)
                if i.split("\t")[2] == 'transcript' and transcript_first == 0:
                    transcript_first = 1
                elif i==info_list[-2]:
                    if i.split("\t")[2] == 'exon':
                        exon_start = int(float(i.split("\t")[3])-gene_start)
                        exon_end = int(float(i.split("\t")[4])-gene_start)
                        for j in range(exon_start,exon_end+1):
                            gene_map[j] += 1
                        transcript_map = gene_map
                        gene_map = [0 for i in range(int(map_size))]
                        transcript_list.append(transcript_map)
                    else:
                        if sum(gene_map) != 0:
                            transcript_map = gene_map
                            gene_map = [0 for i in range(int(map_size))]
                            transcript_list.append(transcript_map)
                elif i.split("\t")[2] == 'exon':
                    exon_start = int(float(i.split("\t")[3])-gene_start)
                    exon_end = int(float(i.split("\t")[4])-gene_start)
                    for j in range(exon_start,exon_end+1):
                        gene_map[j] += 1
                elif i.split("\t")[2] == 'transcript' and transcript_first != 0:
                    transcript_map = gene_map
                    gene_map = [0 for i in range(int(map_size))]
                    transcript_list.append(transcript_map)
                    transcript_first += 1
            for i in range(len(transcript_list)):
                for h in range(int(extension_start-gene_start),int(extension_end-gene_start+1)):
                    transcript_list[i][h] += 10
                if sum(transcript_list[i])>max_sum:
                    max_sum = sum(transcript_list[i])
                    max_index = i
            if max_index == -1:
                print(line)
            gene_map = transcript_list[max_index]
            map_to_return = []
            for h in range(int(extension_start-gene_start),int(extension_end-gene_start+1)):
                map_to_return.append(gene_map[h])
            return map_to_return


# In[79]:


"""
description:
    gtf 정보와 ACTG 결과물을 통해 펩타이드가 매핑된 유전자가 인트론인지 엑손인지 구분
    
arguments:
    line : one key of flat_with_gff_dict
    gtf_dict : dictionary that has gtf information
        key : gene name
        value : all lines about specific gene
    
output : 
    map_to_return : list that shows given genomic region is intron(=10) or exon(>10)
    ex) [10, 10, 10, 10, 11, 11, 11, 11]
"""
def how_exon_intron(line,gtf_dict):
    extension_start = float(line.split("\t")[3])
    extension_end = float(line.split("\t")[4])
    name = line.split("\t")[12]
#     print('---')
#     print(name)
#     print('---')
    info = ''
    gene_start = 0
    gene_end = 0
    if gtf_dict.get(name) != None:
        info = gtf_dict.get(name)
        info_list = info.split("\n")
        gene_start = float(info.split("\n")[0].split("\t")[3])
        gene_end = float(info.split("\n")[0].split("\t")[4])
        map_size = gene_end - gene_start + 1
        gene_map = [0 for i in range(int(map_size))]
        for i in info_list[1:-1]:
            exon_start = float(0)
            exon_end = float(0)
            if i.split("\t")[2] == 'exon':
                exon_start = int(float(i.split("\t")[3])-gene_start)
                exon_end = int(float(i.split("\t")[4])-gene_start)
#                 print(info.split("\n")[0])
#                 print(i)
#                 print(exon_start)
#                 print(exon_end)
#                 print(len(gene_map))
                for j in range(exon_start,exon_end+1):
                    gene_map[j] += 1
        for h in range(int(extension_start-gene_start),int(extension_end-gene_start+1)):
            gene_map[h] += 10
        map_to_return = []
        for h in range(int(extension_start-gene_start),int(extension_end-gene_start+1)):
            map_to_return.append(gene_map[h])
        if len(set(map_to_return)) ==1:
            return map_to_return
        else: 
            transcript_list = []
            transcript_first = 0
            max_index = -1
            max_sum = -1
            gene_map = [0 for i in range(int(map_size))]
            for i in info_list[1:-1]:
                exon_start = float(0)
                exon_end = float(0)
                if i.split("\t")[2] == 'transcript' and transcript_first == 0:
                    transcript_first = 1
                elif i==info_list[-2]:
                    if i.split("\t")[2] == 'exon':
                        exon_start = int(float(i.split("\t")[3])-gene_start)
                        exon_end = int(float(i.split("\t")[4])-gene_start)
                        for j in range(exon_start,exon_end+1):
                            gene_map[j] += 1
                        transcript_map = gene_map
                        gene_map = [0 for i in range(int(map_size))]
                        transcript_list.append(transcript_map)
                    else:
                        if sum(gene_map) != 0:
                            transcript_map = gene_map
                            gene_map = [0 for i in range(int(map_size))]
                            transcript_list.append(transcript_map)
                elif i.split("\t")[2] == 'exon':
                    exon_start = int(float(i.split("\t")[3])-gene_start)
                    exon_end = int(float(i.split("\t")[4])-gene_start)
                    for j in range(exon_start,exon_end+1):
                        gene_map[j] += 1
                elif i.split("\t")[2] == 'transcript' and transcript_first != 0:
                    transcript_map = gene_map
                    gene_map = [0 for i in range(int(map_size))]
                    transcript_list.append(transcript_map)
                    transcript_first += 1
            for i in range(len(transcript_list)):
                for h in range(int(extension_start-gene_start),int(extension_end-gene_start+1)):
                    transcript_list[i][h] += 10
                if sum(transcript_list[i])>max_sum:
                    max_sum = sum(transcript_list[i])
                    max_index = i
            if max_index == -1:
                print(line)
            gene_map = transcript_list[max_index]
            map_to_return = []
            for h in range(int(extension_start-gene_start),int(extension_end-gene_start+1)):
                map_to_return.append(gene_map[h])
            return map_to_return
    else:
        return None


# In[80]:


"""
description:
    how_exon으로 계산한 인트론 엑손 정보를 기존 flat_with_gff_dict에 추가

arguments:
    tsv_dict : ACTG ouput(flat_with_gff_dict)
        key : concatenated one line of flat and gff
        value : 1
    gtf_dict : dictionary that has gtf information
        key : gene name
        value : all lines about specific gene

output : 
    tsv_mapped_dict : ACTG ouput(flat_with_gff_dict) + intron, exon genomic region info
"""
def make_mapped_dict(tsv_dict, gtf_dict):
    tsv_mapped_dict = {}
    for i in tsv_dict.keys():
        key = i +"\t"+str(how_exon(i,gtf_dict))
        tsv_mapped_dict[key] = 1
    return tsv_mapped_dict


# In[81]:


"""
description:
    how_exon으로 계산한 인트론 엑손 정보를 기존 flat_with_gff_dict에 추가

arguments:
    tsv_dict : ACTG ouput(flat_with_gff_dict)
        key : concatenated one line of flat and gff
        value : 1
    gtf_dict : dictionary that has gtf information
        key : gene name
        value : all lines about specific gene

output : 
    tsv_mapped_dict : ACTG ouput(flat_with_gff_dict) + intron, exon genomic region info
"""
def make_mapped_dict_intron(tsv_dict, gtf_dict):
    tsv_mapped_dict = {}
    for i in tsv_dict.keys():
        answer = how_exon_intron(i,gtf_dict)
        if answer != None:
            key = i +"\t"+str(answer)
            tsv_mapped_dict[key] = 1
    return tsv_mapped_dict


# In[82]:


"""
description:
    ACTG의 결과물 중 매핑된 펩타이드가 intron과 exon이 어떻게 매핑되었는지 확인
    ex) full_exon, full_intron, Front_exon, Back_exon

arguments:
    tsv_mapped_dict : ACTG ouput(flat_with_gff_dict) + intron, exon genomic region info

output : 
    full_intron_exon_dict : tsv_mapped_dict + mapping case (ex. full_exon, full_intron, Front_exon, Back_exon)
"""
def make_full_intron_exon_dict(tsv_mapped_dict):
    full_intron_exon_dict = {}
    count = 0
    information = 'Null'
    for line in tsv_mapped_dict.keys():
        line = line.replace("\n","")
        start = line.split("\t")[3]
        end = line.split("\t")[4]
        f_b = line.split("\t")[6]
        match_list = line.split("\t")[14].replace("[","").replace("]","").split(", ")
        if len(set(match_list)) ==1:
            if match_list[0] == '10':
                full_intron_exon_dict[line+"\t0//full_intron"] = 1
            else:
                full_intron_exon_dict[line+"\t1//full_exon"] = 1
        elif len(set(match_list)) ==2:
            full = len(match_list)
            if match_list[0] == '10' and match_list[-1] == '10':
                print("what happen!!!!!!!!!")
                for i in match_list:
                    if i == '10':
                        count+=1
                print(line)
                information = "Middle_exon"
                print("what happen!!!!!!!!!")
            elif match_list[0] != '10' and match_list[-1] != '10':
                print("what happen!!!!!!!!!")
                print(line)
                information = 'Both_Back_Front_exon'
                print("what happen!!!!!!!!!")
            else:
                for i in match_list:
                    if i == '10':
                        count+=1
                if match_list[0] != '10':
                    if f_b == '+':
                        information = 'Front_exon'
                    else:
                        information = 'Back_exon'
                else:
                    if f_b == '+':
                        information = 'Back_exon'
                    else:
                        information = 'Front_exon'
            full_intron_exon_dict[line+"\t"+str((full-count)/full)+"//"+information] = 1
            count = 0
            information = 'Null'
        else:
            full = len(match_list)
            for i in match_list:
                if i == '10':
                    count+=1
            if match_list[0] != '10':
                if match_list[-1]!= '10':
                    information = 'Both_Back_Front_exon'
                else:
                    if f_b == '+':
                        information = 'Front_exon'
                    else:
                        information = 'Back_exon'
            else:
                if f_b == '+':
                    information = 'Back_exon'
                else:
                    information = 'Front_exon'
            full_intron_exon_dict[line+"\t"+str((full-count)/full)+"//"+information] = 1
            count = 0
            information = 'Null'
    return full_intron_exon_dict


# In[83]:


"""
description:
    exon, intron 각각 overlap_num 수 만큼 핵산이 매핑 threshold 조건

arguments:
    line : one key of tsv_mapped_dict
    overlap_num : 매핑된 핵산 수 threshold
output : 
    0 : 각 영역에 매핑 핵산 수가 overlap_num 보다 낮은 경우
    1 : 각 영역에 매핑 핵산 수가 overlap_num 보다 같거나 높은 경우
"""
def ee_trash(line,overlap_num):
    match_list = list(map(int, line.split("\t")[14].replace("[","").replace("]","").split(", ")))
    peptide = line.split("\t")[10]
    information = line.split("\t")[-1].split("//")[1]
    intron_overlap_num = 0
    exon_overlap_num = 0
    for num in match_list:
        if num == 10:
            intron_overlap_num += 1
        else:
            exon_overlap_num += 1
    if exon_overlap_num < overlap_num or intron_overlap_num < overlap_num:
        return 0
    else:
        return 1


# In[84]:


"""
description:
    full_intron_exon_dict 중에서 full_intron과 매핑된 핵산 수 threshold를 넘기지 못하는 candidate와 다음 단계로 나아갈 candidate 구분

arguments:
    full_intron_exon_dict : ACTG ouput(flat_with_gff_dict) + intron, exon genomic region info + mapping case 
    overlap_num : 매핑된 핵산 수 threshold
    
output : 
    threshold_overlap_over_4_dict : full_intron_exon_dict 중 매핑된 핵산 수 조건을 충족한 candidate
    full_intron_dict : full_intron_exon_dict 중 매핑된 영역이 전부 intron인 경우
"""
def make_threshold_overlap_over_4_dict_and_full_intron_dict(full_intron_exon_dict,overlap_num):
    threshold_overlap_over_4_dict = {}
    full_intron_dict = {}
    for line in full_intron_exon_dict.keys():
        line = line.replace("\n","")
        peptide = line.split("\t")[10]
        match_list = line.split("\t")[14].replace("[","").replace("]","").split(", ")
        overlap = match_list.count('11')
        infor = line.split("\t")[-1].split("//")[1]
        if infor == 'full_exon':
            threshold_overlap_over_4_dict[line] = 1
        elif infor != 'full_intron' and infor != 'full_exon':
            if len(peptide)>7:
                if infor == 'Back_exon':
                    result = ee_trash(line,overlap_num)
                    if result == 0:
                        pass
                    else:
                        threshold_overlap_over_4_dict[line] = 1
                elif infor == 'Front_exon':
                    result = ee_trash(line,overlap_num)
                    if result == 0:
                        pass
                    else:
                        threshold_overlap_over_4_dict[line] = 1
        elif infor == 'full_intron':
            full_intron_dict[line] = 1
    return threshold_overlap_over_4_dict,full_intron_dict


# In[85]:


"""
description:
    input dict의 key 값들을 출력 

arguments:
    now_dict : input dict
    set_num : set number
    patient_num : patient number

output:
    file (ex. patient_num + dictionary_name + after_34.bed)
"""
def make_bed_dict_with_now_dict(now_dict,set_num,patient_num):
    now = str(retrieve_name(now_dict))
    bed = open('./OUTPUT/'+set_num+"/"+patient_num+"_"+now+"_after_34.bed","wt")
    for i in now_dict.keys():
        bed.write(i)
    bed.close()


# In[86]:


"""
description:
    make bed file and bed format dict

arguments:
    now_dict : input dict
    set_num : set number
    patient_num : patient number

output: 
    file (ex. patient_num + _novel_evnet_candi_before_bam.bed)
    bed_dict : bed format of input dict
"""
def make_bed_and_dict_with_now_full_actg(now_dict,bed_dict, patient_num, set_num):
    for line in now_dict.keys():
        info = line.split("\t")
        chromo = info[0]
        start = str(int(info[3])-1)
        end = info[4]
        bed_dict[chromo+"\t"+start+"\t"+end+"\t"+"_".join(info[6:]).replace("\n","")+"_"+set_num+"_"+patient_num+"\n"] = 1


# In[87]:


"""
description:
    make bed file and bed format dict

arguments:
    now_dict : input dict
    set_num : set number
    patient_num : patient number

output: 
    file (ex. patient_num + _novel_evnet_candi_before_bam.bed)
    bed_dict : bed format of input dict
"""
def make_bed_and_dict_with_now(now_dict, patient_num, set_num):
    bed_dict = {}
    bed = open('./BED/'+set_num+"/"+patient_num+"_novel_evnet_candi_before_bam.bed","wt")
    for line in now_dict.keys():
        info = line.split("\t")
        chromo = info[0]
        start = str(int(info[3])-1)
        end = info[4]
        bed.write(chromo+"\t"+start+"\t"+end+"\t"+"_".join(info[6:])+"\n")
        bed_dict[chromo+"\t"+start+"\t"+end+"\t"+"_".join(info[6:])] = 1
    bed.close()
    return bed_dict


# In[88]:


"""
description:
    make bed file and bed format dict

arguments:
    now_dict : input dict
    set_num : set number
    patient_num : patient number

output: 
    file (ex. patient_num + _novel_evnet_candi_before_bam.bed)
    bed_dict : bed format of input dict
"""
def make_bed_and_dict_with_now_intron(now_dict, patient_num, set_num):
    bed_dict = {}
    bed = open('./BED/'+set_num+"/"+patient_num+"_intron_novel_evnet_candi_before_bam.bed","wt")
    for line in now_dict.keys():
        info = line.split("\t")
        chromo = info[0]
        start = str(int(info[3])-1)
        end = info[4]
        bed.write(chromo+"\t"+start+"\t"+end+"\t"+"_".join(info[6:])+"\n")
        bed_dict[chromo+"\t"+start+"\t"+end+"\t"+"_".join(info[6:])] = 1
    bed.close()
    return bed_dict


# In[89]:


"""
description:
    make bed file and bed format dict

arguments:
    now_dict : input dict
    set_num : set number
    patient_num : patient number

output: 
    file (ex. patient_num + _novel_evnet_candi_before_bam.bed)
    bed_dict : bed format of input dict
"""
def make_bed_and_dict_with_now_whole(now_dict, full_intron_dict, patient_num, set_num):
    bed_dict = {}
    bed = open('./BED/'+set_num+"/"+patient_num+"_novel_evnet_candi_before_bam.bed","wt")
    for line in now_dict.keys():
        info = line.split("\t")
        chromo = info[0]
        start = str(int(info[3])-1)
        end = info[4]
        bed.write(chromo+"\t"+start+"\t"+end+"\t"+"_".join(info[6:])+"\n")
        bed_dict[chromo+"\t"+start+"\t"+end+"\t"+"_".join(info[6:])] = 1
    for line in full_intron_dict.keys():
        info = line.split("\t")
        chromo = info[0]
        start = str(int(info[3])-1)
        end = info[4]
        bed.write(chromo+"\t"+start+"\t"+end+"\t"+"_".join(info[6:])+"\n")
        bed_dict[chromo+"\t"+start+"\t"+end+"\t"+"_".join(info[6:])] = 1
    bed.close()
    return bed_dict


# In[90]:


"""
description:
    samtools depth process for candidate

arguments:
    patient_num : patient number
    set_num : set number

output : 
    file : "OUTPUT/"+set_num+"/"+patient_num+"_candidate_rna_mapped.txt"
"""
def samtools_depth_output(patient_num, set_num):
    
    bam_file_path = "./BAM/"+set_num+"/"
    file_list = os.listdir(bam_file_path)
    bam_file_name = ''
    for file in file_list:
        if patient_num in file:
            if file.split(".")[-1] == 'bam':
                bam_file_name = bam_file_path+file
                
    bed_file_name = ''
    bed_file_path = "./BED/"+set_num+"/"
    file_list = os.listdir(bed_file_path)
    for file in file_list:
        if patient_num in file:
            if file.split(".")[-1] == 'bed':
                bed_file_name = bed_file_path+file
                
    output_file_name = "OUTPUT/"+set_num+"/"+patient_num+"_candidate_rna_mapped.txt"
    
    samtools_input = 'samtools depth -q 0 '+bam_file_name+' -b '+bed_file_name+' -o '+output_file_name
    r = subprocess.Popen(samtools_input, shell=True).wait()
    if r == 1: 
        print("samtools depth process failed")
        
def samtools_depth_output_one(bam_path, patient_num, set_num):
    
    bam_file_path = "./BAM/"+set_num+"/"
    file_list = os.listdir(bam_file_path)
    bam_file_name = bam_path
                
    bed_file_name = ''
    bed_file_path = "./BED/"+set_num+"/"
    file_list = os.listdir(bed_file_path)
    for file in file_list:
        if patient_num in file:
            if file.split(".")[-1] == 'bed':
                bed_file_name = bed_file_path+file
                
    output_file_name = "OUTPUT/"+set_num+"/"+patient_num+"_candidate_rna_mapped.txt"
    
    samtools_input = 'samtools depth -q 0 '+bam_file_name+' -b '+bed_file_name+' -o '+output_file_name
    r = subprocess.Popen(samtools_input, shell=True).wait()
    if r == 1: 
        print("samtools depth process failed")


# In[91]:


"""
description:
    samtools depth process for candidate

arguments:
    patient_num : patient number
    set_num : set number

output : 
    file : "OUTPUT/"+set_num+"/"+patient_num+"_candidate_rna_mapped.txt"
"""
def samtools_depth_output_intron(patient_num, set_num):
    
    bam_file_path = "./BAM/"+set_num+"/"
    file_list = os.listdir(bam_file_path)
    bam_file_name = ''
    for file in file_list:
        if patient_num in file:
            if file.split(".")[-1] == 'bam':
                bam_file_name = bam_file_path+file
                
    bed_file_name = ''
    bed_file_path = "./BED/"+set_num+"/"
    file_list = os.listdir(bed_file_path)
    for file in file_list:
        if patient_num in file:
            if 'intron' in file:
                if file.split(".")[-1] == 'bed':
                    bed_file_name = bed_file_path+file
                
    output_file_name = "OUTPUT/"+set_num+"/"+patient_num+"_candidate_rna_mapped_intron.txt"
    
    samtools_input = 'samtools depth -q 0 '+bam_file_name+' -b '+bed_file_name+' -o '+output_file_name
    r = subprocess.Popen(samtools_input, shell=True).wait()
    if r == 1: 
        print("samtools depth process failed")


# In[92]:


"""
description:
    read output of samtools depth process and make dict

arguments:
    patient_num : patient number
    set_num : set number

output : 
    rna_count_dict : read count at genomic region
        key : chromo+":"+spot
        value : read_count
"""
def make_rna_count_dict_intron(patient_num,set_num):
    rna_count_path = "./OUTPUT/"+set_num+"/"
    file_list = os.listdir(rna_count_path)
    rna_count_dict = {}
    for i in file_list:
        if i.split(".")[-1] == 'txt':
            if patient_num in i:
                if 'intron' in i:
                    rna_read = open(rna_count_path+i)
                    rna_read_ori = rna_read.readlines()
                    for line in rna_read_ori:
                        line = line.replace("\n","")
                        info = line.split("\t")
                        chromo = info[0]
                        spot = info[1]
                        read_count = info[2]
                        rna_count_dict[chromo+":"+spot] = read_count
    return rna_count_dict


# In[93]:


"""
description:
    read output of samtools depth process and make dict

arguments:
    patient_num : patient number
    set_num : set number

output : 
    rna_count_dict : read count at genomic region
        key : chromo+":"+spot
        value : read_count
"""
def make_rna_count_dict(patient_num,set_num):
    rna_count_path = "./OUTPUT/"+set_num+"/"
    file_list = os.listdir(rna_count_path)
    rna_count_dict = {}
    for i in file_list:
        if i.split(".")[-1] == 'txt':
            if patient_num in i:
                rna_read = open(rna_count_path+i)
                rna_read_ori = rna_read.readlines()
                for line in rna_read_ori:
                    line = line.replace("\n","")
                    info = line.split("\t")
                    chromo = info[0]
                    spot = info[1]
                    read_count = info[2]
                    rna_count_dict[chromo+":"+spot] = read_count
    return rna_count_dict


# In[94]:


"""
description:
    concatenate rna read count dict and candidate dict

arguments:
    rna_dict : read count at genomic regcion
        key : chromo+":"+spot
        value : read_count
    bed_dict : candidate dict
        key : candidate
        value : 1

output : 
    candi_with_rna_dict : candidate info + read count
    ex) chr2	979	999	+_2_ID=10956_10956_ALRR_GCCCTGGAGA_ENSG00000168763_exon-extension_[10, 10, 11, 11]_0.5//Back_exon_[2, 2, 2, 2]_2
"""
def make_candi_with_rna_dict_bed_dict(rna_dict, bed_dict):
    candi_with_rna_dict = {}
    candi_without_rna_dict = {}
    no_count = 0
    for i in bed_dict.keys():
        info = i.split("\t")
        start = int(info[1])
        end = int(info[2])
        read_count_list = []
        chromo = info[0]
        total_count = 0
        for spot in range(start+1, end+1):
            key = chromo+":"+str(spot)
            if rna_dict.get(key) is None:
                read_count_list.append(0)
            else:
                read_count_list.append(int(rna_dict.get(key)))
                total_count += int(rna_dict.get(key))
        mean_depth = total_count/(end-start)
        if mean_depth != 0:
            candi_with_rna_dict[i+"_"+str(read_count_list)+"_"+str(mean_depth)] = 1
        else:
            candi_without_rna_dict[i+"_"+str(read_count_list)+"_"+str(mean_depth)] = 1
    print(no_count)
    return candi_with_rna_dict,candi_without_rna_dict


# In[95]:


"""
description:
    모든 intron 부위의 rna read count가 0이면 filter out

arguments:
    candi_with_rna_dict : candidate info + read count

output : 
    filter_intron_rna_dict : 모든 intron 부위의 rna read count가 0이면 해당 candidate가 제외된 candidate dict
"""
def make_filter_intron_rna_none_dict(candi_with_rna_dict):
    filter_intron_rna_dict = {}
    for line in candi_with_rna_dict.keys():
        rna_read = line.split("_")[-2]
        map_exon_intron = line.split("_")[-5]
        rna_read_list = rna_read.replace("[","").replace("]","").split(", ")
        map_exon_intron_list = map_exon_intron.replace("[","").replace("]","").split(", ")
        count_10 = map_exon_intron.count('10')
        if count_10 > 3:
            for one in range(len(map_exon_intron_list)):
                if map_exon_intron_list[one] == '10':
                    if rna_read_list[one] != '0':
                        filter_intron_rna_dict[line] = 1
                        break
        else:
            filter_intron_rna_dict[line] = 1
    return filter_intron_rna_dict


# In[96]:


"""
description:
    여러줄을 한줄로 만들어 줌

arguments:
    ori : 한줄로 만들고 싶은 파일

output : 
    line : 여러줄을 한줄로
"""
def make_one_line(ori):
    chromo_list = []
    for i in ori[1:]:
        chromo_list.append(i.replace("\n",""))
    line = ''.join(chromo_list)
    return line


# In[97]:


"""
description:
    GRch37 fasta(chromosome 25개) 파일을 읽어서 dict으로 만들어줌

arguments:

output : 
    fasta_dict
        key : chormosome name (ex. X,Y,X,1,2,3)
        value : nucleic acid sequence (ex. ATCGATCG)
"""
def make_fasta_dict():
    fasta_dict = {}
    file_list = os.listdir('./FASTA')
    for fasta in file_list:
        sequence = []
        chromo = fasta.split(".")[0]
        chromosome = open("./FASTA/"+fasta)
        chromosome_ori = chromosome.readlines()
        fasta_dict[chromo] = make_one_line(chromosome_ori)
    return fasta_dict


# In[98]:


"""
description:
    samtools depth 프로세스 이후 candidate 중 exon-skipping event

arguments:
    filter_intron_rna_none_dict : 모든 intron 부위의 rna read count가 0이면 해당 candidate가 제외된 candidate dict

output : 
    es_dict : samtools depth 프로세스 이후 candidate 중 exon-skipping event
"""

def make_es_dict(filter_intron_rna_none_dict):
    es_dict = {}
    for i in filter_intron_rna_none_dict.keys():
        if 'exon-skipping' in i.split("_")[-6]:
            es_dict[i] = 1
    return es_dict


# In[99]:


"""
description:
    samtools depth 프로세스 이후 candidate 중 exon-extension event

arguments:
    filter_intron_rna_none_dict : 모든 intron 부위의 rna read count가 0이면 해당 candidate가 제외된 candidate dict

output : 
    ee_dict : samtools depth 프로세스 이후 candidate 중 exon-extension event
"""
def make_ee_dict(filter_intron_rna_none_dict):
    ee_dict = {}
    for i in filter_intron_rna_none_dict.keys():
        if i.split("_")[-6] == 'exon-extension':
            ee_dict[i] = 1
    return ee_dict


# In[100]:


def get_info_dict_for_mean_depth(final_ee_dict_with_info,gtf_dict):
    gene_to_transcript_dict = {}
    transcipt_to_info_dict = {}
    for i in final_ee_dict_with_info.keys():
        gene = i.split("_")[-2]
        if gtf_dict.get(gene) != None:
            info = gtf_dict.get(gene)
            lines = info.split("\n")
            for k in lines[:-1]:
                p_type = k.split("\t")[1]
                if p_type == 'protein_coding' and k.split("\t")[2] == "transcript":
                    t_id = k.split("\"")[3]
                    if transcipt_to_info_dict.get(t_id) is None:
                        transcipt_to_info_dict[t_id] = k
                    if gene_to_transcript_dict.get(gene) is None:
                        gene_to_transcript_dict[gene] = t_id
                    else:
                        gene_to_transcript_dict[gene] = gene_to_transcript_dict.get(gene)+"\n"+t_id
                elif p_type == 'protein_coding' and '; transcript_id ' in k :
                    transcipt_to_info_dict[t_id] = transcipt_to_info_dict.get(t_id)+"\n"+k
    return gene_to_transcript_dict, transcipt_to_info_dict


# In[101]:


"""
description:
    make infomation dict with gtf file

arguments:
    ee_dict : samtools depth 프로세스 이후 candidate 중 exon-extension event
    gtf_dict
        key : gene name
        value : every lines in gtf about specific gene

output : 
    gene_to_transcript_dict
        key : gene name
        value : transcript id name
    transcipt_to_info_dict
        key : transcript id name
        value : every lines in gtf about specific transcript
"""
def get_info_dict(ee_dict,gtf_dict):
    gene_to_transcript_dict = {}
    transcipt_to_info_dict = {}
    for i in ee_dict.keys():
        gene = i.split("_")[-7]
        if gtf_dict.get(gene) != None:
            info = gtf_dict.get(gene)
            lines = info.split("\n")
            for k in lines[:-1]:
                p_type = k.split("\t")[1]
                if p_type == 'protein_coding' and k.split("\t")[2] == "transcript":
                    t_id = k.split("\"")[3]
                    if transcipt_to_info_dict.get(t_id) is None:
                        transcipt_to_info_dict[t_id] = k
                    if gene_to_transcript_dict.get(gene) is None:
                        gene_to_transcript_dict[gene] = t_id
                    else:
                        gene_to_transcript_dict[gene] = gene_to_transcript_dict.get(gene)+"\n"+t_id
                elif p_type == 'protein_coding' and '; transcript_id ' in k :
                    transcipt_to_info_dict[t_id] = transcipt_to_info_dict.get(t_id)+"\n"+k
    return gene_to_transcript_dict, transcipt_to_info_dict


# In[102]:


"""
description:
    어떤 CDS 영역에 candidate가 겹치는지 확인

arguments:
    candi_start : candidate start region
    candi_end : candidate end region
    cds_start : cds start region
    cds_end : cds end region

output : 
    True : CDS 영역에 걸친경우
    False : CDS 영역에 걸치지 않은 경우
"""
def check_region(candi_start, candi_end, cds_start, cds_end):
    if cds_start < candi_start and candi_start < cds_end:
        return True
    elif cds_start < candi_end and candi_end < cds_end:
        return True
    else:
        return False


# In[103]:


"""
description:
    exon-extension 중 candidate의 앞부분이 exon이고 뒤부분이 intron인 경우 frame이 맞는지 확인

arguments:
    CDS_region : CDS 영역 ex) 123:234//345:456
    frame : CDS의 frame ex) 0,1,2
    p_m : candidate mapping 방향 ex) +,-
    candi_start : candidate start region
    candi_end : candidate end region
    
output : 
    True : frame이 맞을 경우
    False : frame이 틀릴 경우
"""
def frame_check_front_novel_event(CDS_region,frame, p_m,candi_start,candi_end):
    if p_m == "+":
        end = int(CDS_region.split("//")[-2].split(":")[1])
        if (end - candi_start + 3 - frame + 1) % 3 == 0:
            return True
        else:
            return False
    else:
        end = int(CDS_region.split("//")[-2].split(":")[0])
        if (candi_end - end + 3 - frame + 1)%3 == 0:
            return True
        else:
            return False


# In[104]:


"""
description:
    exon-extension 중 candidate의 앞부분이 intron이고 뒤부분이 exon인 경우 frame이 맞는지 확인
    
arguments:
    CDS_region : CDS 영역 ex) 123:234//345:456
    frame : CDS의 frame ex) 0,1,2
    p_m : candidate mapping 방향 ex) +,-
    candi_start : candidate start region
    candi_end : candidate end region
    
output : 
    True : frame이 맞을 경우
    False : frame이 틀릴 경우
"""
def frame_check_back_novel_event(CDS_region,frame, p_m,candi_start,candi_end):
    if p_m == "+":
        end = int(CDS_region.split("//")[-2].split(":")[1])
        if (end - candi_start + 3 - frame+1) % 3 == 0:
            return True
        else:
            return False
    else:
        end = int(CDS_region.split("//")[-2].split(":")[0])
        if (candi_end - end + 3 - frame+1)%3 == 0:
            return True
        else:
            return False


# In[105]:


"""
description:
    뒷부분이 exon인 경우 candidate와 가장 가까운 frame이 맞는 stop 코돈 탐색

arguments:
    intron_seq : intron sequence ex) AACCTAG
    intron_position : intron 영역 ex) 123:456
    p_m : candidate mapping 방향 ex) +,-
    candi_start : candidate start region
    candi_end : candidate end region

output : 
    stop_codon_exist : True or False
    stop codon position : stop:end ex) 123:125
"""
def get_far_stop_codon_position_for_back(intron_seq, p_m,intron_position,candi_start,candi_end):
    seq = []
    intron_start = int(intron_position.split(":")[0])
    intron_end = int(intron_position.split(":")[1])
    if p_m == '+':
        temp_end = candi_start-1
        intron_seq = intron_seq[0:temp_end-intron_start+1]
        intron_end = temp_end

        for i in range(len(intron_seq)//3):
            if i == 0:
                seq.append(intron_seq[int((-i*3)-3):])
            else:
                seq.append(intron_seq[int((-i*3)-3):int((-i*3))])
    elif p_m == '-':
        temp_start = candi_end+1
        intron_seq = intron_seq[temp_start-intron_start:]
        intron_start = temp_start

        for i in range(len(intron_seq)//3):
            seq.append(intron_seq[int(i*3):int(i*3+3)])
    stop_codon_exist = False
    now_intron_seq = []
    if p_m == '+':
        for i in range(len(seq)):
            intron_end= intron_end-3
            now_intron_seq.append(seq[i])
            if seq[i] in plus_stop_codon:
                stop_codon_exist = True

                return stop_codon_exist, str(intron_end+1)+":"+str(intron_end+3)
    elif p_m == '-':
        for i in range(len(seq)):
            intron_start += 3
            now_intron_seq.append(seq[i])
            if seq[i] in minus_stop_codon:
                stop_codon_exist = True

                return stop_codon_exist, str(intron_start-3)+":"+str(intron_start-1)
    return stop_codon_exist, None


# In[106]:


def get_close_stop_codon_position_for_front_ri_frame_out(intron_seq, p_m,intron_position, candi_start, candi_end):
    seq = []
    
    intron_start = int(intron_position.split(":")[0])
    intron_end = int(intron_position.split(":")[1])
    
    
    
    if p_m == '+':
        pep_frame = (candi_start- intron_start)%3
        intron_seq = intron_seq[pep_frame:]
        for i in range(len(intron_seq)//3):
            seq.append(intron_seq[int(i*3):int(i*3+3)])
    elif p_m == '-':
        pep_frame = (intron_end- candi_end)%3
        if pep_frame != 0:
            intron_seq = intron_seq[:-pep_frame]
        for i in range(len(intron_seq)//3):
            if i == 0:
                seq.append(intron_seq[int((-i*3)-3):])
            else:
                seq.append(intron_seq[int((-i*3)-3):int((-i*3))])
                
    stop_codon_exist = False
    now_intron_seq = []
    if p_m == '+':
        intron_start = candi_end+1
        for i in range(len(seq)):
            intron_start+=3
            now_intron_seq.append(seq[i])
            if seq[i] in plus_stop_codon:
                stop_codon_exist = True

                return stop_codon_exist, str(intron_start-3)+":"+str(intron_start-1)
    elif p_m == '-':
        intron_end = candi_start-1
        for i in range(len(seq)):
            now_intron_seq.append(seq[i])
            intron_end= intron_end-3
            if seq[i] in minus_stop_codon:
                stop_codon_exist = True

                return stop_codon_exist, str(intron_end+1)+":"+str(intron_end+3)
    return stop_codon_exist, None


# In[107]:


def get_close_stop_codon_position_for_front_frame_out(intron_seq, p_m,intron_position, candi_start, candi_end):
    seq = []
    
    intron_start = int(intron_position.split(":")[0])
    intron_end = int(intron_position.split(":")[1])
    
    if p_m == '+':
        
        intron_seq = intron_seq[candi_end-intron_start+1:]
        for i in range(len(intron_seq)//3):
            seq.append(intron_seq[int(i*3):int(i*3+3)])
    elif p_m == '-':
        
        intron_seq = intron_seq[:candi_start-intron_start]
        for i in range(len(intron_seq)//3):
            if i == 0:
                seq.append(intron_seq[int((-i*3)-3):])
            else:
                seq.append(intron_seq[int((-i*3)-3):int((-i*3))])
                
    stop_codon_exist = False
    now_intron_seq = []
    if p_m == '+':
        intron_start = candi_end+1
        for i in range(len(seq)):
            intron_start+=3
            now_intron_seq.append(seq[i])
            if seq[i] in plus_stop_codon:
                stop_codon_exist = True

                return stop_codon_exist, str(intron_start-3)+":"+str(intron_start-1)
    elif p_m == '-':
        intron_end = candi_start-1
        for i in range(len(seq)):
            now_intron_seq.append(seq[i])
            intron_end= intron_end-3
            if seq[i] in minus_stop_codon:
                stop_codon_exist = True

                return stop_codon_exist, str(intron_end+1)+":"+str(intron_end+3)
    return stop_codon_exist, None


# In[108]:


"""
description:
    앞부분이 exon인 경우 candidate와 가장 가까운 frame이 맞는 stop 코돈 탐색

arguments:
    intron_seq : intron sequence ex) AACCTAG
    intron_position : intron 영역 ex) 123:456
    p_m : candidate mapping 방향 ex) +,-

output : 
    stop_codon_exist : True or False
    stop codon position : stop:end ex) 123:125
"""
def get_close_stop_codon_position_for_front(intron_seq, p_m,intron_position):
    seq = []
    
    intron_start = int(intron_position.split(":")[0])
    intron_end = int(intron_position.split(":")[1])
    
    if p_m == '+':
        for i in range(len(intron_seq)//3):
            seq.append(intron_seq[int(i*3):int(i*3+3)])
    elif p_m == '-':
        for i in range(len(intron_seq)//3):
            if i == 0:
                seq.append(intron_seq[int((-i*3)-3):])
            else:
                seq.append(intron_seq[int((-i*3)-3):int((-i*3))])
                
    stop_codon_exist = False
    now_intron_seq = []
    if p_m == '+':
        for i in range(len(seq)):
            intron_start+=3
            now_intron_seq.append(seq[i])
            if seq[i] in plus_stop_codon:
                stop_codon_exist = True

                return stop_codon_exist, str(intron_start-3)+":"+str(intron_start-1)
    elif p_m == '-':
        for i in range(len(seq)):
            now_intron_seq.append(seq[i])
            intron_end= intron_end-3
            if seq[i] in minus_stop_codon:
                stop_codon_exist = True

                return stop_codon_exist, str(intron_end+1)+":"+str(intron_end+3)
    return stop_codon_exist, None


# In[109]:


"""
description:
    candidate의 upstream cds의 frame을 계산

arguments:
    CDS_region : CDS 영역 ex) 123:234//345:456 (from get_cds_one_transcript_region)

output : 
    frame : CDS의 frame ex) 0,1,2
"""
def get_frame_up_stream(CDS_region):
    CDSs = CDS_region.split("//")
    total_len = 0
    for i in CDSs[:-1]:
        start = int(i.split(":")[0])
        end = int(i.split(":")[1])
        total_len += end-start+1
    frame = total_len%3
    return frame


# In[110]:


"""
description:
    intron 서열과 영역을 계산

arguments:
    CDS_region : CDS 영역 ex) 123:234//345:456
    chromo : chromosome name ex) X,Y,M,1,2,3
    p_m : candidate mapping 방향 ex) +,-
    frame : CDS의 frame ex) 0,1,2
    
output : 
    intron_seq : intron sequence ex) AACCTAG
    intron_position : intron 영역 ex) 123:456
"""
def get_intron(CDS_region, chromo, p_m,frame):
    CDSs = CDS_region.split("//")
    if p_m == "+":
        # [start, end] one-based
        start = int(CDSs[-2].split(":")[1])
        end = int(CDSs[-1].split(":")[0])-1
        key = str(start-frame+1)+":"+str(end)
#         print('frame : '+str(frame))
#         print("intron")
#         print("chr"+chromo+":"+str(start-frame)+"-"+str(end))
        # [start, end) zero-based
        return fasta_dict.get(chromo)[start-frame:end], key
    else:
        # [start, end] one-based
        start = int(CDSs[-1].split(":")[1])
        end = int(CDSs[-2].split(":")[0])-1
        key = str(start+1)+":"+str(end+frame)
#         print('frame : '+str(frame))
#         print("intron")
#         print("chr"+chromo+":"+str(start)+"-"+str(end+frame))
        # [start, end) zero-based
        return fasta_dict.get(chromo)[start:end+frame], key
    #     print("chr"+chromo+":"+str(start)+"-"+str(end))


# In[111]:


def exon_point_calculate(line_ee_dict):
    pep_start = int(line_ee_dict.split("\t")[1])
    pep_end = int(line_ee_dict.split("\t")[2])
    p_m = line_ee_dict.split("\t")[3].split("_")[0]
    overlap_exon_num = line_ee_dict.split("_")[8].count("11")
    b_f = line_ee_dict.split("//")[1].split("_")[0]
    exon_point = 0
    if p_m == "+":
        if b_f == "Back":
            exon_point = pep_end - overlap_exon_num+1
        elif b_f == "Front":
            exon_point = pep_start + overlap_exon_num
    elif p_m == "-":
        if b_f == "Back":
            exon_point = pep_start + overlap_exon_num
        elif b_f == "Front":
            exon_point = pep_end - overlap_exon_num+1
    return exon_point


# In[112]:


def get_cds_all_transcript_region_frame_out(line_ee_dict, gene_to_transcript_dict, transcipt_to_info_dict,over_threshold_34_ee_candi_dict,junction_dict,stop_dict,over_threshold_34_ri_candi_dict):
    front_back = line_ee_dict.split("_")[-4].split("//")[1]
    chromo = line_ee_dict.split("\t")[0].split("chr")[1]
    gene = line_ee_dict.split("_")[-7]
    p_m = line_ee_dict.split("\t")[3].split("_")[0]
    candi_start = int(line_ee_dict.split("\t")[1])+1
    candi_end = int(line_ee_dict.split("\t")[2])
    pep = line_ee_dict.split("\t")[3].split("_")[4]
    ri_info = ''
    if gene_to_transcript_dict.get(gene) == None:
        pass
    else:
        transcripts = gene_to_transcript_dict.get(gene).split("\n")
        now_candidate = False
        retained_intron = False
        for transcript in transcripts:
            transcript_info = transcipt_to_info_dict.get(transcript)
            #ex) cds_region = '123:555//666:777'
            cds_region = get_cds_one_transcript_region(transcript_info,candi_start,candi_end,front_back)
            full_cds_region = get_cds_one_transcript_region_frame_out(transcript_info)
            if cds_region != None:
                stop_codon_exist = False
                frame = get_frame_up_stream(cds_region)
                intron_seq,intron_position = get_intron(cds_region,chromo,p_m,frame)
                
                if front_back == 'Front': #앞쪽이 exon 이라면
                    
                    stop_codon_exist, stop_codon_position = get_close_stop_codon_position_for_front_frame_out(intron_seq, p_m,intron_position, candi_start, candi_end)
                    if stop_codon_exist:
                        start = int(stop_codon_position.split(":")[0])
                        end = int(stop_codon_position.split(":")[1])

                        stop_temp = fasta_dict.get(chromo)[start-1:end]

                        check_donor = get_donor_candi(start, end,chromo, p_m,candi_start, candi_end)
                        if check_donor != None:
                            now_candidate = True
                            stop_dict['chr'+chromo+"\t"+str(start-1)+"\t"+str(end)+"\t"+p_m+'stop_'+line_ee_dict.split("\t")[0]+':'+line_ee_dict.split("\t")[1]+"-"+line_ee_dict.split("\t")[2]+"\n"]=1
                            for donor in check_donor:
                                junction_dict['chr'+chromo+"\t"+donor.split(":")[0]+"\t"+donor.split(":")[1]+"\t"+p_m+"_"+str(frame)+'_donor_'+line_ee_dict.split("\t")[0]+':'+line_ee_dict.split("\t")[1]+"-"+line_ee_dict.split("\t")[2]+"\n"]=1

                    else: #stop 코돈 없어!
                        retained_intron = True
                        ri_info = "_"+p_m+"_"+str(frame)+"_"+intron_position+"\n"
                        check_donor = get_donor_candi_no_stop(intron_position ,chromo, p_m,candi_start, candi_end)
                        if check_donor != None:
                            now_candidate = True # exon-extension
                            for donor in check_donor:
                                junction_dict['chr'+chromo+"\t"+donor.split(":")[0]+"\t"+donor.split(":")[1]+"\t"+p_m+"_"+str(frame)+'_donor_'+line_ee_dict.split("\t")[0]+':'+line_ee_dict.split("\t")[1]+"-"+line_ee_dict.split("\t")[2]+"\n"]=1
                        else: #도너 없어
                            pass
                    
                elif front_back == 'Back': #뒤쪽이 exon 이라면
                    stop_codon_exist, stop_codon_position = get_far_stop_codon_position_for_back(intron_seq, p_m,intron_position,candi_start,candi_end)
                    
                    if stop_codon_exist: #stop 코돈이 있다면
                        start = int(stop_codon_position.split(":")[0])
                        end = int(stop_codon_position.split(":")[1])
                        stop_temp = fasta_dict.get(chromo)[start-1:end]
                        check_accep = get_acceptor_candi(start, end,chromo, p_m,candi_start, candi_end,frame)
                        if check_accep != None: #acceptor가 있다면
                            now_candidate = True
                            stop_dict['chr'+chromo+"\t"+str(start-1)+"\t"+str(end)+"\t"+p_m+'stop_'+line_ee_dict.split("\t")[0]+':'+line_ee_dict.split("\t")[1]+"-"+line_ee_dict.split("\t")[2]+"\n"]=1
                            for accep in check_accep:
                                junction_dict['chr'+chromo+"\t"+accep.split(":")[0]+"\t"+accep.split(":")[1]+"\t"+p_m+"_"+str(frame)+'_acceptor_'+line_ee_dict.split("\t")[0]+':'+line_ee_dict.split("\t")[1]+"-"+line_ee_dict.split("\t")[2]+"\n"]=1
                                pass

                    else: #stop 코돈이 없다면
                        check_accep = get_acceptor_candi_no_stop(intron_position, chromo, p_m,candi_start, candi_end,frame)
                        if check_accep != None: #acceptor가 있다면
                            now_candidate = True
                            for accep in check_accep:
                                junction_dict['chr'+chromo+"\t"+accep.split(":")[0]+"\t"+accep.split(":")[1]+"\t"+p_m+"_"+str(frame)+'_acceptor_'+line_ee_dict.split("\t")[0]+':'+line_ee_dict.split("\t")[1]+"-"+line_ee_dict.split("\t")[2]+"\n"]=1
                            retained_intron = check_back_retained_intron(p_m, intron_position,candi_start,candi_end)
                            if retained_intron:
                                ri_info = "_"+p_m+"_"+str(frame)+"_"+intron_position+"\n"
                        else: #acceptor가 없다면
                            retained_intron = check_back_retained_intron(p_m, intron_position,candi_start,candi_end)
                            if retained_intron:
                                ri_info = "_"+p_m+"_"+str(frame)+"_"+intron_position+"\n"                    

        if now_candidate == True:
            over_threshold_34_ee_candi_dict[line_ee_dict+"\n"] = 1
            
        if retained_intron:
            over_threshold_34_ri_candi_dict[line_ee_dict.replace("\n","")+ri_info] =1
            
    return over_threshold_34_ee_candi_dict, junction_dict, stop_dict, over_threshold_34_ri_candi_dict


# In[113]:


"""
description:
    기존 ACTG exon-extension 결과를 analysis 4번 기준에 따라 retained intron과 exon extension으로 분류

arguments:
    line_ee_dict
    gene_to_transcript_dict
    transcipt_to_info_dict
    over_threshold_34_ee_candi_dict
    junction_dict,stop_dict
    over_threshold_34_ri_candi_dict
    
output : 
    over_threshold_34_ee_candi_dict
    junction_dict
    stop_dict
    over_threshold_34_ri_candi_dict
"""
def get_cds_all_transcript_region(line_ee_dict, gene_to_transcript_dict, transcipt_to_info_dict,over_threshold_34_ee_candi_dict,junction_dict,stop_dict,over_threshold_34_ri_candi_dict):
    front_back = line_ee_dict.split("_")[-4].split("//")[1]
    chromo = line_ee_dict.split("\t")[0].split("chr")[1]
    gene = line_ee_dict.split("_")[-7]
    p_m = line_ee_dict.split("\t")[3].split("_")[0]
    candi_start = int(line_ee_dict.split("\t")[1])+1
    candi_end = int(line_ee_dict.split("\t")[2])
    ri_info = ''
    if gene_to_transcript_dict.get(gene) == None:
        pass
    else:
        transcripts = gene_to_transcript_dict.get(gene).split("\n")
        now_candidate = False
        retained_intron = False
        for transcript in transcripts:
            transcript_info = transcipt_to_info_dict.get(transcript)
            #ex) cds_region = '123:555//666:777'
            cds_region = get_cds_one_transcript_region(transcript_info,candi_start,candi_end,front_back)
            if cds_region != None:
                stop_codon_exist = False
                frame = get_frame_up_stream(cds_region)
                intron_seq,intron_position = get_intron(cds_region,chromo,p_m,frame)
                
                if front_back == 'Front': #앞쪽이 exon 이라면
                    check_frame = frame_check_front_novel_event(cds_region,frame, p_m,candi_start,candi_end)
                    
                    if check_frame: #frame 확인
                        stop_codon_exist, stop_codon_position = get_close_stop_codon_position_for_front(intron_seq, p_m,intron_position)
                        if stop_codon_exist:
                            start = int(stop_codon_position.split(":")[0])
                            end = int(stop_codon_position.split(":")[1])
                            
                            stop_temp = fasta_dict.get(chromo)[start-1:end]

                            check_donor = get_donor_candi(start, end,chromo, p_m,candi_start, candi_end)
                            if check_donor != None:
                                now_candidate = True
                                stop_dict['chr'+chromo+"\t"+str(start-1)+"\t"+str(end)+"\t"+p_m+'stop_'+line_ee_dict.split("\t")[0]+':'+line_ee_dict.split("\t")[1]+"-"+line_ee_dict.split("\t")[2]+"\n"]=1
                                for donor in check_donor:
                                    junction_dict['chr'+chromo+"\t"+donor.split(":")[0]+"\t"+donor.split(":")[1]+"\t"+p_m+"_"+str(frame)+'_donor_'+line_ee_dict.split("\t")[0]+':'+line_ee_dict.split("\t")[1]+"-"+line_ee_dict.split("\t")[2]+"\n"]=1
                                    
                        else: #stop 코돈 없어!
                            retained_intron = True
                            ri_info = "_"+p_m+"_"+str(frame)+"_"+intron_position+"\n"
                            check_donor = get_donor_candi_no_stop(intron_position ,chromo, p_m,candi_start, candi_end)
                            if check_donor != None:
                                now_candidate = True # exon-extension
                                for donor in check_donor:
                                    junction_dict['chr'+chromo+"\t"+donor.split(":")[0]+"\t"+donor.split(":")[1]+"\t"+p_m+"_"+str(frame)+'_donor_'+line_ee_dict.split("\t")[0]+':'+line_ee_dict.split("\t")[1]+"-"+line_ee_dict.split("\t")[2]+"\n"]=1
                            else: #도너 없어
                                pass
                    else: #frame 맞지 않아
                        pass
                else: #뒤쪽이 exon 이라면
                    stop_codon_exist, stop_codon_position = get_far_stop_codon_position_for_back(intron_seq, p_m,intron_position,candi_start,candi_end)
                    
                    if stop_codon_exist: #stop 코돈이 있다면
                        start = int(stop_codon_position.split(":")[0])
                        end = int(stop_codon_position.split(":")[1])
                        stop_temp = fasta_dict.get(chromo)[start-1:end]
                        check_accep = get_acceptor_candi(start, end,chromo, p_m,candi_start, candi_end,frame)
                        if check_accep != None: #acceptor가 있다면
                            now_candidate = True
                            stop_dict['chr'+chromo+"\t"+str(start-1)+"\t"+str(end)+"\t"+p_m+'stop_'+line_ee_dict.split("\t")[0]+':'+line_ee_dict.split("\t")[1]+"-"+line_ee_dict.split("\t")[2]+"\n"]=1
                            for accep in check_accep:
                                junction_dict['chr'+chromo+"\t"+accep.split(":")[0]+"\t"+accep.split(":")[1]+"\t"+p_m+"_"+str(frame)+'_acceptor_'+line_ee_dict.split("\t")[0]+':'+line_ee_dict.split("\t")[1]+"-"+line_ee_dict.split("\t")[2]+"\n"]=1
                                pass

                    else: #stop 코돈이 없다면
                        check_accep = get_acceptor_candi_no_stop(intron_position, chromo, p_m,candi_start, candi_end,frame)
                        if check_accep != None: #acceptor가 있다면
                            now_candidate = True
                            for accep in check_accep:
                                junction_dict['chr'+chromo+"\t"+accep.split(":")[0]+"\t"+accep.split(":")[1]+"\t"+p_m+"_"+str(frame)+'_acceptor_'+line_ee_dict.split("\t")[0]+':'+line_ee_dict.split("\t")[1]+"-"+line_ee_dict.split("\t")[2]+"\n"]=1
                            retained_intron = check_back_retained_intron(p_m, intron_position,candi_start,candi_end)
                            if retained_intron:
                                ri_info = "_"+p_m+"_"+str(frame)+"_"+intron_position+"\n"
                        else: #acceptor가 없다면
                            retained_intron = check_back_retained_intron(p_m, intron_position,candi_start,candi_end)
                            if retained_intron:
                                ri_info = "_"+p_m+"_"+str(frame)+"_"+intron_position+"\n"

        if now_candidate == True:
            over_threshold_34_ee_candi_dict[line_ee_dict+"\n"] = 1
            
        if retained_intron:
            over_threshold_34_ri_candi_dict[line_ee_dict.replace("\n","")+ri_info] =1
            
    return over_threshold_34_ee_candi_dict, junction_dict, stop_dict, over_threshold_34_ri_candi_dict


# In[114]:


"""
description:
    뒷 부분이 exon인 retained intron candidate에 대해서 frame이 맞는지 계산

arguments:
    p_m
    intron_position
    candi_start
    candi_end

output : 
"""
def check_back_retained_intron(p_m, intron_position,candi_start,candi_end):
    if p_m == '+':
        now_len = int(intron_position.split(":")[0])-candi_start
        if now_len%3 == 0:
            return True
        else:
            return False
    else:
        now_len = candi_end - int(intron_position.split(":")[1])
        if now_len%3 == 0:
            return True
        else:
            return False


# In[115]:


"""
description:
    stop 코돈이 없을때의 Acceptor candidate 찾기

arguments:
    intron_position
    chromo
    p_m
    candi_start
    candi_end
    frame

output : 
    return_list
"""
def get_acceptor_candi_no_stop(intron_position, chromo, p_m,candi_start, candi_end,frame):
    return_list = []
    if p_m == "+":
        start = int(intron_position.split(":")[0])-1+frame
        end = candi_start-1
        now_seq = fasta_dict.get(chromo)[start:end]
        for i in range(len(now_seq)-1):
            now_accep = now_seq[i:i+2]
            if now_accep == plus_acceptor:
                after_seq = now_seq[i+2:]
                if (len(after_seq)+frame)%3 ==0:
                    key = str(start)+":"+str(start+2)
                    return_list.append(key)

            start+=1
        if len(return_list)>0:
            return return_list
    else:
        start = candi_end
        end = int(intron_position.split(":")[1])+frame
        now_seq = fasta_dict.get(chromo)[start:end]
        
        for i in range(len(now_seq)-1):
            now_accep = now_seq[i:i+2]
            if now_accep == minus_acceptor:
                after_seq = now_seq[:i]
                if (len(after_seq)+frame)%3 ==0:
                    key = str(start)+":"+str(start+2)
                    return_list.append(key)

            start+=1
        if len(return_list)>0:
            return return_list
    return None


# In[116]:


"""
description:
    stop 코돈이 있을때의 Acceptor candidate 찾기

arguments:
    stop_start
    stop_end
    chromo
    p_m
    candi_start
    candi_end
    frame

output : 
    return_list
"""
def get_acceptor_candi(stop_start, stop_end,chromo, p_m,candi_start, candi_end,frame):
    return_list = []
    if p_m == "+":
        start = stop_start
        end = candi_start-1
        now_seq = fasta_dict.get(chromo)[start:end]
        for i in range(len(now_seq)-1):
            now_accep = now_seq[i:i+2]
            if now_accep == plus_acceptor:
                after_seq = now_seq[i+2:]
                if (len(after_seq)+frame)%3 ==0:
                    key = str(start)+":"+str(start+2)
                    return_list.append(key)

            start+=1
        if len(return_list)>0:
            return return_list
    else:
        start = candi_end
        end = stop_end-1
        now_seq = fasta_dict.get(chromo)[start:end]
        for i in range(len(now_seq)-1):
            now_accep = now_seq[i:i+2]
            if now_accep == minus_acceptor:
                after_seq = now_seq[:i]
                if (len(after_seq)+frame)%3 ==0:
                    key = str(start)+":"+str(start+2)
                    return_list.append(key)

            start+=1
        if len(return_list)>0:
            return return_list
    return None


# In[117]:


"""
description:
    stop 코돈이 없을 때의 donor candidate 찾기

arguments:
    intron_position
    chromo
    p_m
    candi_start
    candi_end

output : 
    return_list
"""
def get_donor_candi_no_stop(intron_position ,chromo, p_m,candi_start, candi_end):
    return_list = []
    if p_m == "+":
        start = candi_end
        end = int(intron_position.split(":")[1])
        now_seq = fasta_dict.get(chromo)[start:end]
        for i in range(len(now_seq)-1):
            now_donor = now_seq[i:i+2]
            if now_donor == plus_donor:
                key = str(start)+":"+str(start+2)
                return_list.append(key)

            start+=1
        if len(return_list)>0:
            return return_list
    else:
        start = int(intron_position.split(":")[0])-1
        end = candi_start-1
        now_seq = fasta_dict.get(chromo)[start:end]
        for i in range(len(now_seq)-1):
            now_donor = now_seq[i:i+2]
            if now_donor == minus_donor:
                key = str(start)+":"+str(start+2)
                return_list.append(key)

            start+=1
        if len(return_list)>0:
            return return_list
    return None


# In[118]:


def get_cds_one_transcript_region_frame_out(transcript_info):
    CDSs = transcript_info.split("\tCDS\t")
    cds_map = []
    
    for CDS in CDSs[1:]:
        cds_start = int(CDS.split("\t")[0])
        cds_end = int(CDS.split("\t")[1])
        key = str(cds_start)+":"+str(cds_end)
        cds_map.append(key)
    if len(cds_map) >1:
        return "//".join(cds_map)
    return None


# In[119]:


"""
description:
    하나의 transcript에서 첫 CDS부터 candidate가 걸친 CDS까지의 CDS 영역들 반환(front 엑손의 경우 그 다음 CDS 영역까지 반환)

arguments:
    transcript_info
    candi_start
    candi_end
    front_back

output : 
    Back:
    하나의 transcript에서 첫 CDS부터 candidate가 걸친 CDS까지의 CDS 영역들 반환
    Front:
    front 엑손의 경우 그 다음 CDS 영역까지 반환
"""
def get_cds_one_transcript_region(transcript_info,candi_start,candi_end,front_back):
    CDSs = transcript_info.split("\tCDS\t")
    cds_map = []
    last_for_front = False
    
    for CDS in CDSs[1:]:
        cds_start = int(CDS.split("\t")[0])
        cds_end = int(CDS.split("\t")[1])
        check = check_region(candi_start, candi_end, cds_start, cds_end)
        if check == False:
            key = str(cds_start)+":"+str(cds_end)
            cds_map.append(key)
            if last_for_front:
                if len(cds_map) >1:
                    return "//".join(cds_map)
        elif check == True:
            if front_back == "Back":
                key = str(cds_start)+":"+str(cds_end)
                cds_map.append(key)
                if len(cds_map) >1:
                    return "//".join(cds_map)
            else:
                key = str(cds_start)+":"+str(cds_end)
                cds_map.append(key)
                last_for_front = True
    return None


# In[120]:


"""
description:
    하나의 transcript에서 첫 CDS부터 candidate가 있는 CDS까지의 영역들을 반환

arguments:
    transcript_info
    candi_region
    p_m
    
output : 
    True # exon-skipping maybe
    junction_region # both junction sites
    skipped_region # [start, end] one-based of skipped CDS
"""
def get_cds_one_transcript_es_region(transcript_info,candi_region,p_m):
    candi_es_list = candi_region.split("//")
    CDSs = transcript_info.split("\tCDS\t")
    cds_map = []
    skipped_region = ''
    
    #exon 하나만 skipping 한거 맞아?
    exon_skipping = 0

    #지금 몇 번째 보고 있어?
    observed = 0
    
    cds_len = 0
    
    for CDS in CDSs[1:]:
        
        cds_start = CDS.split("\t")[0]
        cds_end = CDS.split("\t")[1]
        candi_start = candi_es_list[observed].split(":")[0]
        candi_end = candi_es_list[observed].split(":")[1]
        
        cds_region = cds_start+":"+cds_end
        
        cds_len += int(cds_end) - int(cds_start) + 1
        
        check = check_region_es(candi_start, candi_end, cds_start, cds_end)
        
        if observed == 0:
            if check == True:
                now = int(candi_end) - int(candi_start)+1
#                 if (cds_len-now)%3 != 0: # frame이 맞는지
#                     return False, None, None
                observed += 1
        else:
            if check == False:
                if exon_skipping == 1: # 하나의 exon만 skip 된 것이 맞는지 prevent double-skipping
                    return False, None, None
                elif exon_skipping == 0:
                    exon_skipping = 1
                    if p_m == "+":
                        junction_end = str(int(candi_start)-1)
                        junction_start = candi_es_list[observed-1].split(":")[1]
                    else:
                        junction_end = str(int(candi_es_list[observed-1].split(":")[0])-1)
                        junction_start = candi_end
                    junction_region = junction_start+':'+junction_end
                    skipped_region = cds_start+':'+cds_end
            else:
                if exon_skipping == 1:
                    return True , junction_region, skipped_region
                else:
                    observed += 1
    return False, None, None #


# In[121]:


"""
description:
    지금 CDS 영역이 candidate가 속한 영역이 맞는지 반환

arguments:
    candi_start
    candi_end
    cds_start
    cds_end
    
output : 
"""
def check_region_es(candi_start, candi_end, cds_start, cds_end):
    if cds_start == candi_start or candi_end == cds_end:
        return True
    else:
        return False


# In[122]:


"""
description:
    exon skipping candidate 중 짝이 없이 혼자 있는 candidate 삭제
    
arguments:
    es_dict
    gene_to_transcript_dict
    transcipt_to_info_dict
    set_num
    patient_name

output : 
"""
def es_process(es_dict, gene_to_transcript_dict, transcipt_to_info_dict, set_num, patient_name):
    temp_es_dict = {}
    real_over2_es_dict = {}
    es_junction_dict = {}
    final_es_dict = {}
    for line in es_dict.keys():
        id_now = line.split("\t")[3].split("_")[3]
        if temp_es_dict.get(id_now) == None:
            temp_es_dict[id_now] = 1
        else:
            temp_es_dict[id_now] = temp_es_dict.get(id_now)+1
    for line in es_dict.keys():
        id_now = line.split("\t")[3].split("_")[3]
        if temp_es_dict.get(id_now) > 1:
            if real_over2_es_dict.get(id_now) == None:
                real_over2_es_dict[id_now] = line
            else:
                real_over2_es_dict[id_now] = real_over2_es_dict.get(id_now).replace("\n","")+"%%"+line.replace("\n","")
    for id_now in real_over2_es_dict.keys():
        es = real_over2_es_dict.get(id_now)
        final_es_dict,es_junction_dict = es_deep_process(es, gene_to_transcript_dict, transcipt_to_info_dict,final_es_dict,es_junction_dict)
    make_bed_dict_with_now_dict(final_es_dict,set_num,patient_name)
    make_bed_dict_with_now_dict(es_junction_dict,set_num,patient_name)


# In[123]:


"""
description:
    exon skipping candidate 중 조건 3번에 해당하는 candidate만 filter

arguments:
    es
    gene_to_transcript_dict
    transcipt_to_info_dict
    final_es_dict
    junction_dict

output : 
    final_es_dict
    junction_dict
"""
def es_deep_process(es, gene_to_transcript_dict, transcipt_to_info_dict,final_es_dict,junction_dict):
    es_list = es.split("%%")
    region_list = []
    chromo = es.split("\t")[0].split("chr")[1]
    gene = es.split("_")[6]
    p_m = es.split("\t")[3].split("_")[0]
    for es_candi in es_list:
        candi_start = int(es_candi.split("\t")[1])+1
        candi_end = int(es_candi.split("\t")[2])
        if (candi_end-candi_start)+1<4: #4이상 겹쳐야해
            return final_es_dict, junction_dict
        key = str(candi_start)+":"+str(candi_end)
        region_list.append(key)
    if p_m == '-':
        region_list.reverse()
    candi_region = "//".join(region_list)
    
    if gene_to_transcript_dict.get(gene) == None:
        pass
    else:
        transcripts = gene_to_transcript_dict.get(gene).split("\n")
        now_candidate = False
        for transcript in transcripts:
            transcript_info = transcipt_to_info_dict.get(transcript)
            check , junction_region, skipped_region = get_cds_one_transcript_es_region(transcript_info,candi_region,p_m)
            if check == True:
                now_candidate = True
                final_es_dict[es.replace("%%","_junction:"+junction_region+"_skipped:"+skipped_region+"\n")+"_junction:"+junction_region+"_skipped:"+skipped_region+"\n"] = 1
                junction_dict[es.split("\t")[0]+"\t"+junction_region.split(":")[0]+"\t"+junction_region.split(":")[1]+"\t"+es.replace("\t","_")+"\n"]=1
    return final_es_dict, junction_dict


# In[124]:


def ee_process_frame_out(ee_dict, gene_to_transcript_dict,transcipt_to_info_dict,set_num, patient_name):
    over_threshold_34_ee_candi_dict = {}
    over_threshold_34_ri_candi_dict = {}
    junction_frame_dict = {}
    stop_frame_dict = {}
    for line in ee_dict.keys():
        over_threshold_34_ee_candi_dict,junction_frame_dict, stop_frame_dict,over_threshold_34_ri_candi_dict = get_cds_all_transcript_region_frame_out(line, gene_to_transcript_dict, transcipt_to_info_dict,over_threshold_34_ee_candi_dict,junction_frame_dict, stop_frame_dict,over_threshold_34_ri_candi_dict)
    final_ee_candi = over_threshold_34_ee_candi_dict
    final_ee_junction = junction_frame_dict
    final_ri_candi = over_threshold_34_ri_candi_dict
    make_bed_dict_with_now_dict(stop_frame_dict,set_num,patient_name)
    make_bed_dict_with_now_dict(final_ee_candi,set_num,patient_name)
    make_bed_dict_with_now_dict(final_ee_junction,set_num,patient_name)
    make_bed_dict_with_now_dict(final_ri_candi,set_num,patient_name)


# In[125]:


"""
description:
    exon-extension candidates 중 조건 4를 충족하는 retained intron, exon extension 분리 및 donor, stop position 계산

arguments:
    ee_dict
    gene_to_transcript_dict
    transcipt_to_info_dict
    set_num
    patient_name
    
output : 
"""
def ee_process(ee_dict, gene_to_transcript_dict,transcipt_to_info_dict,set_num, patient_name):
    over_threshold_34_ee_candi_dict = {}
    over_threshold_34_ri_candi_dict = {}
    junction_dict = {}
    stop_dict = {}
    for line in ee_dict.keys():
        over_threshold_34_ee_candi_dict,junction_dict, stop_dict,over_threshold_34_ri_candi_dict = get_cds_all_transcript_region(line, gene_to_transcript_dict, transcipt_to_info_dict,over_threshold_34_ee_candi_dict,junction_dict, stop_dict,over_threshold_34_ri_candi_dict)
    final_ee_candi = over_threshold_34_ee_candi_dict
    final_ee_junction = junction_dict
    final_ri_candi = over_threshold_34_ri_candi_dict
    make_bed_dict_with_now_dict(stop_dict,set_num,patient_name)
    make_bed_dict_with_now_dict(final_ee_candi,set_num,patient_name)
    make_bed_dict_with_now_dict(final_ee_junction,set_num,patient_name)
    make_bed_dict_with_now_dict(final_ri_candi,set_num,patient_name)


# In[126]:


"""------------------Samtools RNA 매핑 끝------------------"""


# In[127]:


"""------------------junction support read 시작------------------"""


# In[128]:


"""
description:

arguments:
    da_file_ori

output : 
    left_dict
    right_dict
    chromo_dict
"""
def preprocess_ee_junction(da_file_ori):
    left_dict = {}
    right_dict = {}
    chromo_dict = {}
    
    if len(da_file_ori)>0:
        for line in da_file_ori:
            chromo = line.split("\t")[0]
            chromo_dict[chromo] = 1
            p_m = line.split("\t")[3].split("_")[0]
            d_a = line.split("\t")[3].split("_")[2]
            if p_m == "-" and d_a == 'acceptor':
                start = int(line.split("\t")[1])+1
                left_dict[line] = 1
            elif p_m == "+" and d_a == 'donor':
                start = int(line.split("\t")[1])+1
                left_dict[line] = 1
            if p_m == "+" and d_a == 'acceptor':
                end = int(line.split("\t")[2])-1
                right_dict[line] = 1
            elif p_m == "-" and d_a == 'donor':
                end = int(line.split("\t")[2])-1
                right_dict[line] = 1
                
    return left_dict, right_dict, chromo_dict


# In[129]:


"""
description:

arguments:
    ri_file_ori

output : 
    chromo_ri_dict,intron_ri_dict
"""
def preprocess_ri(ri_file_ori):
    chromo_ri_dict = {}
    intron_ri_dict = {}
    
    if len(ri_file_ori)>0:
        for line in ri_file_ori:
            chromo = line.split("\t")[0]
            p_m = line.split("_")[-3]
            frame = int(line.split("_")[-2])
            intron_start = int(line.split("_")[-1].split(":")[0])
            intron_end = int(line.split("_")[-1].split(":")[1])
            
            if p_m == '+':
                intron_start = intron_start+frame
            else:
                intron_end = intron_end-frame
            
            intron_region = chromo+":"+str(intron_start)+":"+str(intron_end)
            
            intron_ri_dict[intron_region] = line
            
            if chromo_ri_dict.get(chromo) == None:
                chromo_ri_dict[chromo] = intron_region
            else:
                chromo_ri_dict[chromo] = chromo_ri_dict.get(chromo)+"//"+intron_region
                
    return chromo_ri_dict,intron_ri_dict


# In[130]:


"""
description:

arguments:
    es_file_ori
    
output : 
    chromo_es_dict
    junction_es_dict
"""
def preprocess_es_junction(es_file_ori):
    chromo_es_dict = {}
    junction_es_dict = {}
    
    if len(es_file_ori)>0:
        for line in es_file_ori:
            chromo = line.split("\t")[0]
            junction = line.split("_")[-2].split('junction:')[1]
            skipped = line.split("skipped:")[1]
            j_s = junction.split(":")[0]
            j_e = junction.split(":")[1]
            skipped_s = skipped.split(":")[0]
            key = chromo+":"+j_s+"-"+j_e+":"+skipped_s
            if junction_es_dict.get(key) == None:
                junction_es_dict[key] = line
            else:
                junction_es_dict[key] = junction_es_dict.get(key)+line
            if chromo_es_dict.get(chromo) == None:
                chromo_es_dict[chromo] = junction
            else:
                check_list = chromo_es_dict.get(chromo).split("//")
                if junction not in check_list:
                    chromo_es_dict[chromo] = chromo_es_dict.get(chromo)+"//"+junction
                
    return chromo_es_dict,junction_es_dict


# In[131]:


def how_rna_read_for_ee_l_exact(rna_start,cigar_str, point):
    now_seq = ''
    temp_concat = ''
    last_M_num = 0
    now_position = rna_start
    m_count = cigar_str.count('M')
    now_m_count = 0
    for i in range(len(cigar_str)):
        now_seq = cigar_str[i:i+1]
        if now_seq.isdigit():
            temp_concat += now_seq
        elif now_seq.isalpha():
            if now_seq == 'M':
                now_m_count += 1
                now_position += int(temp_concat)
                if now_position-1 == point:
                    if now_m_count == m_count:
                        return True
                temp_concat = ''
            else:
                if now_seq == 'N':
                    now_position += int(temp_concat)
                    temp_concat = ''
                elif now_seq == 'D':
                    return False
                elif now_seq == 'I':
                    return False
                else:
                    temp_concat = ''
    return False


# In[132]:


def how_rna_read_for_ee_r_exact(rna_start,cigar_str, point):
    now_seq = ''
    temp_concat = ''
    last_M_num = 0
    now_position = rna_start
    check = False
    
    m_count = cigar_str.count('M')
    now_m_count = 0
    now_alpha = 0
    
    for i in range(len(cigar_str)):
        now_seq = cigar_str[i:i+1]
        if now_seq.isdigit():
            temp_concat += now_seq
        elif now_seq.isalpha():
            now_alpha +=1
            if now_seq == 'N':
                now_position += int(temp_concat)
                temp_concat = ''
            else:
                if now_seq == 'M':
                    if now_alpha == 1:
                        if rna_start == point:
                            return True
                    now_position += int(temp_concat)
                    temp_concat = ''

                elif now_seq == 'D':
                    return False
                elif now_seq == 'I':
                    return False
                else:
                    temp_concat = ''
                check = False
    return False


# In[133]:


"""
description:

arguments:
    rna_start
    cigar_str
    point
    
output : 

"""
def how_rna_read_for_ee_l(rna_start,cigar_str, point):
    now_seq = ''
    temp_concat = ''
    last_M_num = 0
    now_position = rna_start
    check = False
     
    for i in range(len(cigar_str)):
        now_seq = cigar_str[i:i+1]
        if now_seq.isdigit():
            temp_concat += now_seq
        elif now_seq.isalpha():
            if now_seq == 'M':
                now_position += int(temp_concat)
                if now_position-1 == point:
                    check = True
                temp_concat = ''
            else:
                if now_seq == 'N':
                    now_position += int(temp_concat)
                    if check:
                        return True
                    temp_concat = ''
                elif now_seq == 'D':
                    return False
                elif now_seq == 'I':
                    return False
                else:
                    temp_concat = ''
                check = False
    return False


# In[134]:


"""
description:

arguments:
    rna_start
    cigar_str
    point
    
output : 
"""
def how_rna_read_for_ee_r(rna_start,cigar_str, point):
    now_seq = ''
    temp_concat = ''
    last_M_num = 0
    now_position = rna_start
    check = False
     
    for i in range(len(cigar_str)):
        now_seq = cigar_str[i:i+1]
        if now_seq.isdigit():
            temp_concat += now_seq
        elif now_seq.isalpha():
            if now_seq == 'N':
                now_position += int(temp_concat)
                if now_position == point:
                    check = True
                temp_concat = ''
            else:
                if now_seq == 'M':
                    if check:
                        return True
                    now_position += int(temp_concat)
                    temp_concat = ''

                elif now_seq == 'D':
                    return False
                elif now_seq == 'I':
                    return False
                else:
                    temp_concat = ''
                check = False
    return False


# In[135]:


"""
description:

arguments:
    rna_start
    cigar_str
    start
    end
    
output : 
"""
def how_rna_read_for_es(rna_start,cigar_str, start, end):
    now_seq = ''
    temp_concat = ''
    last_M_num = 0
    now_position = rna_start
    jnc_len = end-start
    need_jnc = False
    need_mapping = False
    for i in range(len(cigar_str)):
        now_seq = cigar_str[i:i+1]
        if now_seq.isdigit():
            temp_concat += now_seq
        elif now_seq.isalpha():
            if now_seq == 'M': ## TODO: end 도 검사
                now_position += int(temp_concat)
                if now_position-1 == start:
                    need_jnc = True
                if need_mapping == True:
                    return True
                temp_concat = ''
            elif now_seq == 'N':
                if need_jnc:
                    if int(temp_concat) == jnc_len:
                        need_mapping = True
                    else:
                        return False
                now_position += int(temp_concat)
                temp_concat = ''
            elif now_seq == 'D':
                return False
            elif now_seq == 'I':
                return False
            else:
                temp_concat = ''
    return False


# In[136]:


"""
description:

arguments:
    rna_start
    cigar_str
    a1
    a2
    
output : 
"""
def how_rna_read_for_ri(rna_start,cigar_str, a1, a2):
    now_seq = ''
    temp_concat = ''
    last_M_num = 0
    now_position = rna_start
    
    for i in range(len(cigar_str)):
        now_seq = cigar_str[i:i+1]
        if now_seq.isdigit():
            temp_concat += now_seq
        elif now_seq.isalpha():
            if now_seq == 'M':
                last_M_num = int(temp_concat)
                if now_position<=a1 and last_M_num+now_position>a2:
                    return True
                else:
                    now_position += last_M_num
                temp_concat = ''
            elif now_seq == 'N':
                now_position += int(temp_concat)
                temp_concat = ''
            elif now_seq == 'D':
                return False
            elif now_seq == 'I':
                return False
            else:
                temp_concat = ''
    return False


# In[137]:


## 수정 예정 1028_0955
"""
description:

arguments:
    cigar_str
    process_num
output : 
"""
def how_rna_read(cigar_str, process_num):
    now_seq = ''
    temp_concat = ''
    first_M_num = 0
    now_M = 0
    last_M_num = 0
    first_N_num = 0
    for i in range(len(cigar_str)):
        now_seq = cigar_str[i:i+1]
        if now_seq.isdigit():
            temp_concat += now_seq
        elif now_seq.isalpha():
            if now_seq == 'M':
                if now_M == 0:
                    first_M_num = int(temp_concat)
                    last_M_num = int(temp_concat)
                    now_M+=1
                else:
                    last_M_num = last_M_num+int(temp_concat)
                temp_concat = ''
            elif now_seq == 'N':
                first_N_num = int(temp_concat)
                temp_concat = ''
            else:
                temp_concat = ''
    if process_num ==1:
        return first_M_num
    elif process_num ==2:
        return first_M_num+first_N_num
    elif process_num ==3:
        return first_N_num
    elif process_num ==4:
        return last_M_num+first_N_num
        


# In[138]:


## 수정 예정 1028_1038## 
"""
description:

arguments:
    line
    chromo_es_dict
    junction_es_dict
    final_es_after_j

output : 
"""
def line_process_es_j(set_num, patient_num,bam_name, chromo_es_dict, junction_es_dict,final_es_after_j):
                    
    #key = chromo+":"+j_s+"-"+j_e+":"+skipped_s
    for key in junction_es_dict.keys():
        chromo  = key.split(":")[0]
        start = int(key.split(":")[1].split("-")[0])
#         print(key)
        end = int(key.split(":")[1].split("-")[1])
        skipped_point = int(key.split(":")[2])-1
        input_samtools = chromo+":"+str(start)+"-"+str(start)+" -o test.sam"
        make_sam_input = 'samtools view '+bam_name+" "+input_samtools
        r = subprocess.Popen(make_sam_input, shell=True).wait()
        if r == 1: 
            print("making sam failed")
#         else:
#             print("making sam sucess")
        sam = open("test.sam")
        sam_ori = sam.readlines()
        count = 0
        normal_count = 0
        
        line_temp = ''
        for line in sam_ori:
            cigar = line.split("\t")[5]
            start_c = int(line.split("\t")[3])
            if how_rna_read_for_es(start_c,cigar,start,end):
                count+=1
#                 print("what?")
                #line_temp += line
#                 sam_check.write(line)
        if count != 0:
            for line in sam_ori:
                cigar = line.split("\t")[5]
                start_c = int(line.split("\t")[3])
                if how_rna_read_for_es(start_c,cigar,start,skipped_point):
                    normal_count+=1
            input_es = junction_es_dict.get(key).replace("\n",'_normal_jnc_read:'+str(normal_count)+'('+str(round(count/(normal_count+count) , 4))+')\n')
            final_es_after_j[input_es] = count
#             sam_check = open("./SAM/"+key+"_"+set_num+"_"+patient_num+"_es_only_junction.sam","wt")
#             sam_check.write(key+"\n")
#             sam_check.write(junction_es_dict.get(key)+"\n\n")
#             sam_check.write(line_temp)        
#             sam_check.close()
#             print('go')
        
        remove_sam_input = 'rm test.sam'
        r = subprocess.Popen(remove_sam_input, shell=True).wait()
        if r == 1: 
            print("deleting sam failed")
            
            
    
#     chromo = line.split("\t")[2]
    
#     if chromo_es_dict.get(chromo) != None:
        
#         cigar = line.split("\t")[5]
        
#         if "N" in cigar:
#             start = line.split("\t")[3]
#             candidates = chromo_es_dict.get(chromo).split("//")
#             for candidate in candidates:
#                 a = int(candidate.split(":")[0])
#                 b = int(candidate.split(":")[1])
#                 n = b-a
                
#                 s = int(start)
#                 if a-100 < s and s < a:
#                     if str(n)+"N" in cigar:
#                         key = junction_es_dict.get(candidate)
#                         if key != None:
#                             if final_es_after_j.get(key) == None:
#                                 final_es_after_j[key] = 1
#                             else:
#                                 final_es_after_j[key] = final_es_after_j.get(key)+1


# In[139]:


def inside_line_process_ee_r(chromo,exon_point,bam_name):
    input_samtools = chromo+":"+str(exon_point)+"-"+str(exon_point)+" -o test.sam"
    make_sam_input = 'samtools view '+bam_name+" "+input_samtools
    r = subprocess.Popen(make_sam_input, shell=True).wait()
    if r == 1: 
        print("making sam failed")
    sam = open("test.sam")
    sam_ori = sam.readlines()
    count = 0
    for line in sam_ori:
        cigar = line.split("\t")[5]
        start = int(line.split("\t")[3])
        if how_rna_read_for_ee_r(start,cigar,int(exon_point)):
            count+=1
    exact_count = 0
    if count != 0:
        for line in sam_ori:
            cigar = line.split("\t")[5]
            start = int(line.split("\t")[3])
            if how_rna_read_for_ee_r_exact(start,cigar,int(exon_point)):
                exact_count+= 1
        ##exact read 수 새는 코드 추가 
    remove_sam_input = 'rm test.sam'
    r = subprocess.Popen(remove_sam_input, shell=True).wait()
    if r == 1: 
        print("deleting sam failed")
    return count,exact_count


# In[140]:


def inside_line_process_ee_l(chromo,exon_point,bam_name):
    input_samtools = chromo+":"+str(exon_point)+"-"+str(exon_point)+" -o test.sam"
    make_sam_input = 'samtools view '+bam_name+" "+input_samtools
    r = subprocess.Popen(make_sam_input, shell=True).wait()
    if r == 1: 
        print("making sam failed")
    sam = open("test.sam")
    sam_ori = sam.readlines()
    count = 0
    for line in sam_ori:
        cigar = line.split("\t")[5]
        start = int(line.split("\t")[3])
        if how_rna_read_for_ee_l(start,cigar,int(exon_point)):
            count+=1
    exact_count = 0
    if count != 0:
        for line in sam_ori:
            cigar = line.split("\t")[5]
            start = int(line.split("\t")[3])
            if how_rna_read_for_ee_l_exact(start,cigar,int(exon_point)):
                exact_count+= 1
        ##exact read 수 새는 코드 추가 
    remove_sam_input = 'rm test.sam'
    r = subprocess.Popen(remove_sam_input, shell=True).wait()
    if r == 1: 
        print("deleting sam failed")
    return count,exact_count


# In[141]:


## 수정 예정 1028_0955## 
"""
description:

arguments:

output : 
"""
def line_process_ee_j(bam_name, left_dict,right_dict,chromo_dict,final_ee_after_j,set_num,patient_num):
    ee_file = open('./OUTPUT/'+set_num+'/'+patient_num+"_over_threshold_34_ee_candi_dict_after_34.bed")
    ee_file_ori = ee_file.readlines()
    
    
    for key in left_dict.keys():
        exon_point = -1
        search_point = key.split("-")[-1].replace("\n","")
        chromo  = key.split("\t")[0]
        
        original_exon_junction_read_count = 0
        exon_exact_read_count = 0
        
        novel_junction_read_count = 0
        exact_junction_read_count = 0
        
        for line in ee_file_ori[:]:
            if search_point in line:
                if chromo in line:
                    exon_point = exon_point_calculate(line)
        
        
        
        
        
        point = key.split("\t")[1]
        novel_junction_read_count, exact_junction_read_count = inside_line_process_ee_l(chromo,point,bam_name)
        
        
        
        if novel_junction_read_count != 0:
            if exon_point != -1:
                original_exon_junction_read_count, exon_exact_read_count = inside_line_process_ee_l(chromo,exon_point,bam_name)
            ratio = (novel_junction_read_count+exact_junction_read_count)/(original_exon_junction_read_count+exon_exact_read_count+novel_junction_read_count+exact_junction_read_count)
            input_str = key.replace("\n","")+'//ori_jnc:'+str(original_exon_junction_read_count)+"//ori_exact:"+str(exon_exact_read_count)+"//novel_jnc:"+str(novel_junction_read_count)+"//novel_exact:"+str(exact_junction_read_count)+'('+str(round(ratio,4))+')'+"\n"
            
            final_ee_after_j[input_str] = novel_junction_read_count    

        
    for key in right_dict.keys():
        
        exon_point = -1
        search_point = key.split("-")[-1].replace("\n","")
        chromo  = key.split("\t")[0]
        
        original_exon_junction_read_count = 0
        exon_exact_read_count = 0
        
        novel_junction_read_count = 0
        exact_junction_read_count = 0
        
        for line in ee_file_ori[:]:
            if search_point in line:
                if chromo in line:
                    exon_point = exon_point_calculate(line)
        
        point = str(int(key.split("\t")[2])+1)
        novel_junction_read_count, exact_junction_read_count = inside_line_process_ee_r(chromo,point,bam_name)
        
        
        
        if novel_junction_read_count != 0:
            if exon_point != -1:
                original_exon_junction_read_count, exon_exact_read_count = inside_line_process_ee_r(chromo,exon_point,bam_name)
            ratio = (novel_junction_read_count+exact_junction_read_count)/(original_exon_junction_read_count+exon_exact_read_count+novel_junction_read_count+exact_junction_read_count)
            input_str = key.replace("\n","")+'//ori_jnc:'+str(original_exon_junction_read_count)+"//ori_exact:"+str(exon_exact_read_count)+"//novel_jnc:"+str(novel_junction_read_count)+"//novel_exact:"+str(exact_junction_read_count)+'('+str(round(ratio,4))+')'+"\n"
            
            final_ee_after_j[input_str] = novel_junction_read_count 
                    


# In[142]:


## 수정 예정 1028_0955## 
"""
description:

arguments:

output : 
"""
def line_process_ri(bam_name,chromo_ri_dict,intron_ri_dict,final_ri_after_s):
    
    for key in intron_ri_dict.keys():
        chromo  = key.split(":")[0]
        start_point = key.split(":")[1]
        end_point = key.split(":")[2]
        input_samtools = chromo+":"+start_point+"-"+start_point+" -o test.sam"
        make_sam_input = 'samtools view '+bam_name+" "+input_samtools
        r = subprocess.Popen(make_sam_input, shell=True).wait()
        if r == 1: 
            print("making sam failed")
#         else:
#             print("making sam sucess")
        sam = open("test.sam")
        sam_ori = sam.readlines()
        count = 0
        for line in sam_ori:
            cigar = line.split("\t")[5]
            start = int(line.split("\t")[3])
            if how_rna_read_for_ri(start,cigar,int(start_point)-1,int(start_point)):
                count+=1
        now = intron_ri_dict.get(key).replace("\n","_intronPoint:"+start_point)+"\n"
        if count != 0:
            final_ri_after_s[now] = count
        remove_sam_input = 'rm test.sam'
        r = subprocess.Popen(remove_sam_input, shell=True).wait()
        if r == 1: 
            print("deleting sam failed")
            
            
        input_samtools = chromo+":"+end_point+"-"+end_point+" -o test.sam"
        make_sam_input = 'samtools view '+bam_name+" "+input_samtools
        r = subprocess.Popen(make_sam_input, shell=True).wait()
        if r == 1: 
            print("making sam failed")
#         else:
#             print("making sam sucess")
        sam = open("test.sam")
        sam_ori = sam.readlines()
        count = 0
        for line in sam_ori:
            cigar = line.split("\t")[5]
            start = int(line.split("\t")[3])
            if how_rna_read_for_ri(start,cigar,int(end_point),int(end_point)+1):
                count+=1
        now = intron_ri_dict.get(key).replace("\n","_intronPoint:"+end_point)+"\n"
        if count != 0:
            final_ri_after_s[now] = count
        remove_sam_input = 'rm test.sam'
        r = subprocess.Popen(remove_sam_input, shell=True).wait()
        if r == 1: 
            print("deleting sam failed")

                   


# In[143]:


"""
description:

arguments:
    now_dict
    set_num
    patient_num

output : 
"""
def print_after_sam_file_compare(now_dict,set_num,patient_num):
    now_file = open('./FINAL/'+set_num+"/"+patient_num+"_"+retrieve_name(now_dict)+"_after_making_sam.bed",'wt')
    for key in now_dict.keys():
        now_file.write(now_dict.get(key).replace("\n", "_"+(key)+'\n'))
    now_file.close()


# In[144]:


"""
description:

arguments:
    now_dict
    set_num
    patient_num

output : 
"""
def print_after_sam_file(now_dict,set_num,patient_num):
    now_file = open('./FINAL/'+set_num+"/"+patient_num+"_"+retrieve_name(now_dict)+"_after_making_sam.bed",'wt')
    for key in now_dict.keys():
        now_file.write(key.replace("\n","_num:"+str(now_dict.get(key))+"\n"))
    now_file.close()


# In[145]:


"""
description:

arguments:
    da_file_ori
    ri_file_ori
    es_file_ori
    set_num
    patient_num

output : 
"""
def sam_process_intron(ri_file_ori,ri_file_ori2,  set_num, patient_num):
    
    
    chromo_ri_dict,intron_ri_dict = preprocess_ri(ri_file_ori)
    chromo_ri_dict2,intron_ri_dict2 = preprocess_ri(ri_file_ori2)
    
    chromo_ri_dict.update(chromo_ri_dict2)
    intron_ri_dict.update(intron_ri_dict2)
    
    final_ri_after_s_intron_frame = {}
    
    sam_name = patient_num+'.sam'
    bam_name = ''
    file_list = os.listdir('./BAM/'+set_num)
    for i in file_list:
        if patient_num in i:
            if i.split(".")[-1] == 'bam':
                bam_name = i
    bam_name = './BAM/'+set_num+'/'+bam_name
    
    line_process_ri(bam_name,chromo_ri_dict,intron_ri_dict,final_ri_after_s_intron_frame)
    
    final_ri_after_s_p_full_intron = {}
    final_ri_after_s_p_full_intron = postprocess_ri(final_ri_after_s_intron_frame)
    
    if len(final_ri_after_s_p_full_intron) >0:
        make_bed_for_ri_intron_final(final_ri_after_s_p_full_intron,set_num,patient_num)
        
    print_after_sam_file(final_ri_after_s_p_full_intron,set_num,patient_num)
    


        
    throw_garbage()


# In[146]:


"""
description:

arguments:
    da_file_ori
    ri_file_ori
    es_file_ori
    set_num
    patient_num

output : 
"""
def sam_process(da_file_ori, ri_file_ori, es_file_ori, set_num, patient_num):
    
    
    
    left_dict, right_dict, chromo_dict = preprocess_ee_junction(da_file_ori)
    final_ee_after_j = {}
    
    chromo_es_dict,junction_es_dict = preprocess_es_junction(es_file_ori)
    final_es_after_j = {}
    
    chromo_ri_dict,intron_ri_dict = preprocess_ri(ri_file_ori)
    final_ri_after_s = {}
    
    sam_name = patient_num+'.sam'
    bam_name = ''
    file_list = os.listdir('./BAM/'+set_num)
    for i in file_list:
        if patient_num in i:
            if i.split(".")[-1] == 'bam':
                bam_name = i
    bam_name = './BAM/'+set_num+'/'+bam_name
    

    line_process_ee_j(bam_name,left_dict,right_dict,chromo_dict,final_ee_after_j,set_num,patient_num)
    line_process_es_j(set_num, patient_num,bam_name, chromo_es_dict, junction_es_dict,final_es_after_j)
    line_process_ri(bam_name,chromo_ri_dict,intron_ri_dict,final_ri_after_s)
    
    final_ri_after_s_p = {}
    final_ri_after_s_p = postprocess_ri(final_ri_after_s)
    
    if len(final_ri_after_s_p) >0:
        make_bed_for_ri_final(final_ri_after_s_p,set_num,patient_num)
        
    print_after_sam_file(final_ee_after_j,set_num,patient_num)
    print_after_sam_file(final_es_after_j,set_num,patient_num)
    print_after_sam_file(final_ri_after_s_p,set_num,patient_num)
        
    throw_garbage()
    
    
def sam_process_one(bam_path, da_file_ori, ri_file_ori, es_file_ori, set_num, patient_num):
    
    
    
    left_dict, right_dict, chromo_dict = preprocess_ee_junction(da_file_ori)
    final_ee_after_j = {}
    
    chromo_es_dict,junction_es_dict = preprocess_es_junction(es_file_ori)
    final_es_after_j = {}
    
    chromo_ri_dict,intron_ri_dict = preprocess_ri(ri_file_ori)
    final_ri_after_s = {}
    
    sam_name = patient_num+'.sam'
    bam_name = bam_path
    

    line_process_ee_j(bam_name,left_dict,right_dict,chromo_dict,final_ee_after_j,set_num,patient_num)
    line_process_es_j(set_num, patient_num,bam_name, chromo_es_dict, junction_es_dict,final_es_after_j)
    line_process_ri(bam_name,chromo_ri_dict,intron_ri_dict,final_ri_after_s)
    
    final_ri_after_s_p = {}
    final_ri_after_s_p = postprocess_ri(final_ri_after_s)
    
    if len(final_ri_after_s_p) >0:
        make_bed_for_ri_final(final_ri_after_s_p,set_num,patient_num)
        
    print_after_sam_file(final_ee_after_j,set_num,patient_num)
    print_after_sam_file(final_es_after_j,set_num,patient_num)
    print_after_sam_file(final_ri_after_s_p,set_num,patient_num)
        
    throw_garbage()


# In[147]:


"""
description:
    intron 시작 부분과 끝 부분을 bed 형식으로 계산
    
arguments:
    retained intron info
        ex) chrX  153694955  153694979  +_ANVGEMPR_ENSG00000130827_exon-extension_2_153694854:153694965_intronPoint:153694856
        
output : 
    start : intron 시작 부분 (zero-based)
    end : intron 끝 부분 [one-based]
"""

def find_start_end_for_ri(i):
    p_m = i.split("_intronPoint:")[0].split("_")[-3]
    frame = int(i.split("_intronPoint:")[0].split("_")[-2])
    if p_m == '+':
        start = str(int(i.split("_intronPoint:")[0].split("_")[-1].split(":")[0])-1+frame)
        end = str(i.split("_intronPoint:")[0].split("_")[-1].split(":")[1])
    else:
        start = str(int(i.split("_intronPoint:")[0].split("_")[-1].split(":")[0])-1)
        end = str(int(i.split("_intronPoint:")[0].split("_")[-1].split(":")[1])-frame)
    return start,end


# In[148]:


"""
description:
    intron mean depth를 samtools를 통해 계산합니다.
arguments:
    final_ri_after_s_p : retained intron candidate 정보가 있는 dict
        key : retained intron info
            ex) chrX  153694955  153694979  +_ANVGEMPR_ENSG00000130827_exon-extension__2_153694854:153694965_intronPoint:153694856
        value : intron 한쪽의 junction site에 걸쳐있는 rna read 수
            ex) 2
    patient_num
    set_num
output : 
    intron mean depth 결과
"""

def make_bed_for_ri_intron_final(final_ri_after_s_p,set_num,patient_num):
    intron_mean_depth = open('OUTPUT/'+set_num+"/final_ri_intron_"+patient_num+'.bed','wt')
    check_dict = {}
    
    
    #input bed 생성
    for i in final_ri_after_s_p.keys():
        info = i.replace("\n","").split("\t")[-1]
        chromo = i.split("\t")[0]
        start,end = find_start_end_for_ri(i)

        check_str = chromo+"\t"+start+"\t"+end
        line = chromo+"\t"+start+"\t"+end+"\t"+info+"\n"
        if check_dict.get(check_str) == None:
            check_dict[check_str] = 1
            intron_mean_depth.write(line)
    intron_mean_depth.close()
    
    #samtools depth 실행
    bam_file_path = "./BAM/"+set_num+"/"
    file_list = os.listdir(bam_file_path)
    bam_file_name = ''
    for file in file_list:
        if patient_num in file:
            if file.split(".")[-1] == 'bam':
                bam_file_name = bam_file_path+file
    
    bed_file_name = 'OUTPUT/'+set_num+"/final_ri_intron_"+patient_num+'.bed'
    
    output_file_name = 'OUTPUT/'+set_num+"/final_ri_intron_"+patient_num+'_count_for_full_intron.txt'
    
    samtools_input = 'samtools depth -q 0 '+bam_file_name+' -b '+bed_file_name+' -o '+output_file_name
    r = subprocess.Popen(samtools_input, shell=True).wait()
    if r == 1: 
        print("samtools depth process failed")
    
    #samtools depth 실행 결과 적용 데이터 출력
    mapping_ri_dict = make_rna_count_dict_for_final_ri_full_intron(patient_num,set_num)
    for i in check_dict.keys():
#         print(i)
        chromo = i.split("\t")[0]
        start = int(i.split("\t")[1])+1
        end = int(i.split("\t")[2])
        depth = 0
        for spot in range(start,end+1):
            key = chromo+":"+str(spot)
            if mapping_ri_dict.get(key) != None:
                depth += int(mapping_ri_dict.get(key))
        mean = depth/(end-start+1)
        check_dict[i] = '_meanDepth:'+str(round(mean,4))
        
        
    for i in final_ri_after_s_p.keys():
        chromo = i.split("\t")[0]
        start,end = find_start_end_for_ri(i)

        info = chromo+"\t"+start+"\t"+end
        final_ri_after_s_p[i] = str(final_ri_after_s_p.get(i)).replace("\n","")+check_dict.get(info)
    
        
#     remove_output = 'rm '+output_file_name
#     r = subprocess.Popen(remove_output, shell=True).wait()
#     if r == 1: 
#         print("deleting output failed")
        
#     remove_output = 'rm '+bed_file_name
#     r = subprocess.Popen(remove_output, shell=True).wait()
#     if r == 1: 
#         print("deleting bed failed")


# In[149]:




def make_bed_for_ri_final_frame(final_ri_after_s_p,set_num,patient_num):
    intron_mean_depth = open('OUTPUT/'+set_num+"/final_ri_"+patient_num+'.bed','wt')
    check_dict = {}
    
    
    #input bed 생성
    for i in final_ri_after_s_p.keys():
        info = i.replace("\n","").split("\t")[-1]
        chromo = i.split("\t")[0]
        start,end = find_start_end_for_ri(i)

        check_str = chromo+"\t"+start+"\t"+end
        line = chromo+"\t"+start+"\t"+end+"\t"+info+"\n"
        if check_dict.get(check_str) == None:
            check_dict[check_str] = 1
            intron_mean_depth.write(line)
    intron_mean_depth.close()
    
    #samtools depth 실행
    bam_file_path = "./BAM/"+set_num+"/"
    file_list = os.listdir(bam_file_path)
    bam_file_name = ''
    for file in file_list:
        if patient_num in file:
            if file.split(".")[-1] == 'bam':
                bam_file_name = bam_file_path+file
    
    bed_file_name = 'OUTPUT/'+set_num+"/final_ri_frame"+patient_num+'.bed'
    
    output_file_name = 'OUTPUT/'+set_num+"/final_ri_frame"+patient_num+'_count_for_intron_frame.txt'
    
    samtools_input = 'samtools depth -q 0 '+bam_file_name+' -b '+bed_file_name+' -o '+output_file_name
    r = subprocess.Popen(samtools_input, shell=True).wait()
    if r == 1: 
        print("samtools depth process failed")
    
    #samtools depth 실행 결과 적용 데이터 출력
    mapping_ri_dict = make_rna_count_dict_for_final_ri_frame(patient_num,set_num)
    for i in check_dict.keys():
#         print(i)
        chromo = i.split("\t")[0]
        start = int(i.split("\t")[1])+1
        end = int(i.split("\t")[2])
        depth = 0
        for spot in range(start,end+1):
            key = chromo+":"+str(spot)
            if mapping_ri_dict.get(key) != None:
                depth += int(mapping_ri_dict.get(key))
        mean = depth/(end-start+1)
        check_dict[i] = '_meanDepth:'+str(round(mean,4))
        
        
    for i in final_ri_after_s_p.keys():
        chromo = i.split("\t")[0]
        start,end = find_start_end_for_ri(i)

        info = chromo+"\t"+start+"\t"+end
        final_ri_after_s_p[i] = str(final_ri_after_s_p.get(i)).replace("\n","")+check_dict.get(info)
    
        
    remove_output = 'rm '+output_file_name
    r = subprocess.Popen(remove_output, shell=True).wait()
    if r == 1: 
        print("deleting output failed")
        
    remove_output = 'rm '+bed_file_name
    r = subprocess.Popen(remove_output, shell=True).wait()
    if r == 1: 
        print("deleting bed failed")


# In[150]:


"""
description:
    intron mean depth를 samtools를 통해 계산합니다.
arguments:
    final_ri_after_s_p : retained intron candidate 정보가 있는 dict
        key : retained intron info
            ex) chrX  153694955  153694979  +_ANVGEMPR_ENSG00000130827_exon-extension__2_153694854:153694965_intronPoint:153694856
        value : intron 한쪽의 junction site에 걸쳐있는 rna read 수
            ex) 2
    patient_num
    set_num
output : 
    intron mean depth 결과
"""

def make_bed_for_ri_final(final_ri_after_s_p,set_num,patient_num):
    intron_mean_depth = open('OUTPUT/'+set_num+"/final_ri_"+patient_num+'.bed','wt')
    check_dict = {}
    
    
    #input bed 생성
    for i in final_ri_after_s_p.keys():
        info = i.replace("\n","").split("\t")[-1]
        chromo = i.split("\t")[0]
        start,end = find_start_end_for_ri(i)

        check_str = chromo+"\t"+start+"\t"+end
        line = chromo+"\t"+start+"\t"+end+"\t"+info+"\n"
        if check_dict.get(check_str) == None:
            check_dict[check_str] = 1
            intron_mean_depth.write(line)
    intron_mean_depth.close()
    
    #samtools depth 실행
    bam_file_path = "./BAM/"+set_num+"/"
    file_list = os.listdir(bam_file_path)
    bam_file_name = ''
    for file in file_list:
        if patient_num in file:
            if file.split(".")[-1] == 'bam':
                bam_file_name = bam_file_path+file
    
    bed_file_name = 'OUTPUT/'+set_num+"/final_ri_"+patient_num+'.bed'
    
    output_file_name = 'OUTPUT/'+set_num+"/final_ri_"+patient_num+'_count_for_intron.txt'
    
    samtools_input = 'samtools depth -q 0 '+bam_file_name+' -b '+bed_file_name+' -o '+output_file_name
    r = subprocess.Popen(samtools_input, shell=True).wait()
    if r == 1: 
        print("samtools depth process failed")
    
    #samtools depth 실행 결과 적용 데이터 출력
    mapping_ri_dict = make_rna_count_dict_for_final_ri(patient_num,set_num)
    for i in check_dict.keys():
#         print(i)
        chromo = i.split("\t")[0]
        start = int(i.split("\t")[1])+1
        end = int(i.split("\t")[2])
        depth = 0
        for spot in range(start,end+1):
            key = chromo+":"+str(spot)
            if mapping_ri_dict.get(key) != None:
                depth += int(mapping_ri_dict.get(key))
        mean = depth/(end-start+1)
        p_value = calculate_p_value_of_intron(intron_mean_depth_list, mean)
        check_dict[i] = '_meanDepth:'+str(round(mean,4))+'_pValue:'+str(p_value)
        
        
    for i in final_ri_after_s_p.keys():
        chromo = i.split("\t")[0]
        start,end = find_start_end_for_ri(i)

        info = chromo+"\t"+start+"\t"+end
        final_ri_after_s_p[i] = str(final_ri_after_s_p.get(i)).replace("\n","")+check_dict.get(info)
    
        
    remove_output = 'rm '+output_file_name
    r = subprocess.Popen(remove_output, shell=True).wait()
    if r == 1: 
        print("deleting output failed")
        
    remove_output = 'rm '+bed_file_name
    r = subprocess.Popen(remove_output, shell=True).wait()
    if r == 1: 
        print("deleting bed failed")


# In[151]:


def make_rna_count_dict_for_final_mean_depth_compare(patient_num,set_num):
    rna_count_path = "./OUTPUT/"+set_num+"/"
    file_list = os.listdir(rna_count_path)
    rna_count_dict = {}
    for i in file_list:
        if i.split(".")[-1] == 'txt':
            if '_count_for_compare' in i:
                if patient_num in i:
                    rna_read = open(rna_count_path+i)
                    rna_read_ori = rna_read.readlines()
                    for line in rna_read_ori:
                        line = line.replace("\n","")
                        info = line.split("\t")
                        chromo = info[0]
                        spot = info[1]
                        read_count = info[2]
                        rna_count_dict[chromo+":"+spot] = read_count
    return rna_count_dict


# In[152]:


"""
description:
    samtools를 통해 intron 영역의 rna read count 결과를 저장합니다.
arguments:
    patient_num
    set_num
output : 
    rna_count_dict
        key : chromo:spot (chr1:11375689)
        value : read_depth (17)
"""

def make_rna_count_dict_for_final_ri_full_intron(patient_num,set_num):
    rna_count_path = "./OUTPUT/"+set_num+"/"
    file_list = os.listdir(rna_count_path)
    rna_count_dict = {}
    for i in file_list:
        if i.split(".")[-1] == 'txt':
            if '_count_for_full_intron' in i:
                if patient_num in i:
                    rna_read = open(rna_count_path+i)
                    rna_read_ori = rna_read.readlines()
                    for line in rna_read_ori:
                        line = line.replace("\n","")
                        info = line.split("\t")
                        chromo = info[0]
                        spot = info[1]
                        read_count = info[2]
                        rna_count_dict[chromo+":"+spot] = read_count
    return rna_count_dict


# In[153]:


"""
description:
    samtools를 통해 intron 영역의 rna read count 결과를 저장합니다.
arguments:
    patient_num
    set_num
output : 
    rna_count_dict
        key : chromo:spot (chr1:11375689)
        value : read_depth (17)
"""

def make_rna_count_dict_for_final_ri_frame(patient_num,set_num):
    rna_count_path = "./OUTPUT/"+set_num+"/"
    file_list = os.listdir(rna_count_path)
    rna_count_dict = {}
    for i in file_list:
        if i.split(".")[-1] == 'txt':
            if '_count_for_intron_frame' in i:
                if patient_num in i:
                    rna_read = open(rna_count_path+i)
                    rna_read_ori = rna_read.readlines()
                    for line in rna_read_ori:
                        line = line.replace("\n","")
                        info = line.split("\t")
                        chromo = info[0]
                        spot = info[1]
                        read_count = info[2]
                        rna_count_dict[chromo+":"+spot] = read_count
    return rna_count_dict


# In[154]:


"""
description:
    samtools를 통해 intron 영역의 rna read count 결과를 저장합니다.
arguments:
    patient_num
    set_num
output : 
    rna_count_dict
        key : chromo:spot (chr1:11375689)
        value : read_depth (17)
"""

def make_rna_count_dict_for_final_ri(patient_num,set_num):
    rna_count_path = "./OUTPUT/"+set_num+"/"
    file_list = os.listdir(rna_count_path)
    rna_count_dict = {}
    for i in file_list:
        if i.split(".")[-1] == 'txt':
            if '_count_for_intron' in i:
                if patient_num in i:
                    rna_read = open(rna_count_path+i)
                    rna_read_ori = rna_read.readlines()
                    for line in rna_read_ori:
                        line = line.replace("\n","")
                        info = line.split("\t")
                        chromo = info[0]
                        spot = info[1]
                        read_count = info[2]
                        rna_count_dict[chromo+":"+spot] = read_count
    return rna_count_dict


# In[155]:


"""
description:
    인트론 양단 중 한쪽만 support read가 있는 경우 candidate에서 탈락시킵니다.
arguments:
    final_ri_after_s : retained intron candidate 정보가 있는 dict
        key : retained intron info
            ex) chrX  153694955  153694979  +_ANVGEMPR_ENSG00000130827_exon-extension__2_153694854:153694965_intronPoint:153694856
        value : intron 한쪽의 junction site에 걸쳐있는 rna read 수
            ex) 2
output : 
    after_post_process_ri_dict : final_ri_after_s 중에서 양쪽 모두 support read가 있는 candidate
"""
def postprocess_ri(final_ri_after_s):
    after_post_process_ri_dict = {}
    temp_dict = {}
    for i in final_ri_after_s.keys():
        info = "\t".join(i.split("\t")[0:3])
        if temp_dict.get(info) == None:
            temp_dict[info] = i.replace("\n","")+"___"+str(final_ri_after_s.get(i))
        else:
            after_post_process_ri_dict[temp_dict.get(info).split("___")[0]+"\n"] = temp_dict.get(info).split("___")[1]
            after_post_process_ri_dict[i] = str(final_ri_after_s.get(i))
    return after_post_process_ri_dict


# In[156]:



def see_junction_read_intron_frame_out():
    file_list = os.listdir('./ACTG/')
    path = './ACTG/'   
    for set_num in file_list:
        if 'set' in set_num:
            print(set_num)
            set_path = path+set_num+"/"
            patient_list = os.listdir(set_path)



            for patient in patient_list:
#                 if 'IRCR_GBM14_508' in patient:
                print(patient)

                if len(patient.split(".")) == 1:
                    start_time = time.time()

#                     es_file = open('./OUTPUT/'+set_num+'/'+patient+"_final_es_dict_after_34.bed")
                    ri_file = open('./OUTPUT/'+set_num+'/'+patient+"_over_threshold_34_full_intron_ri_frame_candi_dict_after_34.bed")
                    ri_file2 = open('./OUTPUT/'+set_num+'/'+patient+"_over_threshold_34_ri_frame_candi_dict_after_34.bed")
#                     da_file = open('./OUTPUT/'+set_num+'/'+patient+"_junction_dict_after_34.bed")

#                     es_file_ori = es_file.readlines()
                    ri_file_ori = ri_file.readlines()
                    ri_file_ori2 = ri_file2.readlines()
#                     da_file_ori = da_file.readlines()

                    sam_process_intron(ri_file_ori,ri_file_ori2,  set_num, patient)

                    print("one sam whole time :", time.time() - start_time)
            throw_garbage()


# In[157]:


"""
description:

arguments:

output : 
"""
def see_junction_read_intron():
    file_list = os.listdir('./ACTG/')
    path = './ACTG/'   
    for set_num in file_list:
        if 'set' in set_num:
            print(set_num)
            set_path = path+set_num+"/"
            patient_list = os.listdir(set_path)



            for patient in patient_list:
#                 if 'IRCR_GBM14_508' in patient:
                print(patient)

                if len(patient.split(".")) == 1:
                    start_time = time.time()

#                     es_file = open('./OUTPUT/'+set_num+'/'+patient+"_final_es_dict_after_34.bed")
                    ri_file = open('./OUTPUT/'+set_num+'/'+patient+"_over_threshold_34_full_intron_ri_candi_dict_after_34.bed")
#                     da_file = open('./OUTPUT/'+set_num+'/'+patient+"_junction_dict_after_34.bed")

#                     es_file_ori = es_file.readlines()
                    ri_file_ori = ri_file.readlines()
#                     da_file_ori = da_file.readlines()

                    sam_process_intron(ri_file_ori,  set_num, patient)

                    print("one sam whole time :", time.time() - start_time)
            throw_garbage()


# In[158]:


"""
description:

arguments:

output : 
"""
def see_junction_read():
    intron_mean_depth_list = get_intron_mean_depth_info()
    
    file_list = os.listdir('./ACTG/')
    path = './ACTG/'   
    for set_num in file_list:
        if 'set' in set_num:
            print(set_num)
            set_path = path+set_num+"/"
            patient_list = os.listdir(set_path)



            for patient in patient_list:
#                 if 'IRCR_GBM14_508' in patient:
                print(patient)

                if len(patient.split(".")) == 1:
                    start_time = time.time()

                    es_file = open('./OUTPUT/'+set_num+'/'+patient+"_final_es_dict_after_34.bed")
                    ri_file = open('./OUTPUT/'+set_num+'/'+patient+"_over_threshold_34_ri_candi_dict_after_34.bed")
                    da_file = open('./OUTPUT/'+set_num+'/'+patient+"_junction_frame_dict_after_34.bed")

                    es_file_ori = es_file.readlines()
                    ri_file_ori = ri_file.readlines()
                    da_file_ori = da_file.readlines()

                    sam_process(da_file_ori, ri_file_ori, es_file_ori, set_num, patient)

                    print("one sam whole time :", time.time() - start_time)
            throw_garbage()
            
def see_junction_read_one(patient_name, bam_path):
    intron_mean_depth_list = get_intron_mean_depth_info()
    
    file_list = os.listdir('./ACTG/')
    path = './ACTG/'   
    for set_num in file_list:
        if 'set' in set_num:
            print(set_num)
            set_path = path+set_num+"/"
            patient_list = os.listdir(set_path)



            for patient in patient_list:
#                 if 'IRCR_GBM14_508' in patient:

                if len(patient.split(".")) == 1:
                    if patient == patient_name:
                        start_time = time.time()

                        es_file = open('./OUTPUT/'+set_num+'/'+patient+"_final_es_dict_after_34.bed")
                        ri_file = open('./OUTPUT/'+set_num+'/'+patient+"_over_threshold_34_ri_candi_dict_after_34.bed")
                        da_file = open('./OUTPUT/'+set_num+'/'+patient+"_junction_frame_dict_after_34.bed")

                        es_file_ori = es_file.readlines()
                        ri_file_ori = ri_file.readlines()
                        da_file_ori = da_file.readlines()

                        sam_process_one(bam_path, da_file_ori, ri_file_ori, es_file_ori, set_num, patient)

                        print("one sam whole time :", time.time() - start_time)
            throw_garbage()


# In[159]:


"""
description:

arguments:
    gtf_ori
    patient_name
    set_num
    
output : 
"""
def temp_main_after_rna_mapped_es(gtf_ori,patient_name, set_num):
    gtf_dict = make_gtf_dict(gtf_ori)
    filter_intron_rna_none_dict = {}
    filter_intron_rna = open("./OUTPUT/"+set_num+"/"+patient_name+"_whole_candi_after_samtools_RNA_mapped.bed")
    filter_intron_rna_ori = filter_intron_rna.readlines()
    for i in filter_intron_rna_ori[:]:
        filter_intron_rna_none_dict[i.replace("\n","")] = 1
    es_dict = make_es_dict(filter_intron_rna_none_dict)
    gene_to_transcript_dict, transcipt_to_info_dict = get_info_dict(es_dict,gtf_dict)
    es_process(es_dict, gene_to_transcript_dict, transcipt_to_info_dict, set_num, patient_name)


# In[160]:


"""
description:

arguments:

output : 
"""
def after_rna_mapped_es():
    file_list = os.listdir('./ACTG/')
    path = './ACTG/'   
    for set_num in file_list:
        if 'set' in set_num:
            print(set_num)
            set_path = path+set_num+"/"
            patient_list = os.listdir(set_path)
            for patient in patient_list:
                print(patient)
                if len(patient.split(".")) == 1:
                    folder_path = set_path+patient+"/"
                    start_time = time.time()
                    gtf = open(folder_path+"GTF/"+patient+"_grch37.gtf")
                    gtf_ori = gtf.readlines()
                    temp_main_after_rna_mapped_es(gtf_ori,patient, set_num)
                    print("one patient es whole time :", time.time() - start_time)
                throw_garbage()
                
def after_rna_mapped_es_one(patient_name):
    file_list = os.listdir('./ACTG/')
    path = './ACTG/'   
    for set_num in file_list:
        if 'set' in set_num:
            print(set_num)
            set_path = path+set_num+"/"
            patient_list = os.listdir(set_path)
            for patient in patient_list:
                if len(patient.split(".")) == 1:
                    if patient_name == patient:
                        folder_path = set_path+patient+"/"
                        start_time = time.time()
                        gtf = open(folder_path+"GTF/"+patient+"_grch37.gtf")
                        gtf_ori = gtf.readlines()
                        temp_main_after_rna_mapped_es(gtf_ori,patient, set_num)
                        print("one patient es whole time :", time.time() - start_time)
                    throw_garbage()


# In[161]:


"""
description:

arguments:

output : 
"""
def after_rna_mapped_ee_ri():
    file_list = os.listdir('./ACTG/')
    path = './ACTG/'   
    for set_num in file_list:
        if 'set' in set_num:
            print(set_num)
            set_path = path+set_num+"/"
            patient_list = os.listdir(set_path)
            for patient in patient_list:
                print(patient)
                if len(patient.split(".")) == 1:
                    start_time = time.time()
                    folder_path = set_path+patient+"/"
                    gtf = open(folder_path+"GTF/"+patient+"_grch37.gtf")
                    gtf_ori = gtf.readlines()
                    temp_main_after_rna_mapped_frame_out(gtf_ori,patient, set_num)
                    print("one patient frame work whole time :", time.time() - start_time)
                    
def after_rna_mapped_ee_ri_one(patient_name):
    file_list = os.listdir('./ACTG/')
    path = './ACTG/'   
    for set_num in file_list:
        if 'set' in set_num:
            print(set_num)
            set_path = path+set_num+"/"
            patient_list = os.listdir(set_path)
            for patient in patient_list:
                if len(patient.split(".")) == 1:
                    if patient_name == patient:
                        start_time = time.time()
                        folder_path = set_path+patient+"/"
                        gtf = open(folder_path+"GTF/"+patient+"_grch37.gtf")
                        gtf_ori = gtf.readlines()
                        temp_main_after_rna_mapped_frame_out(gtf_ori,patient, set_num)
                        print("one patient frame work whole time :", time.time() - start_time)


# In[162]:


"""
description:

arguments:

output : 
"""
def after_rna_mapped():
    file_list = os.listdir('./ACTG/')
    path = './ACTG/'   
    for set_num in file_list:
        if 'set' in set_num:
            print(set_num)
            set_path = path+set_num+"/"
            patient_list = os.listdir(set_path)
            for patient in patient_list:
                print(patient)
                if len(patient.split(".")) == 1:
                    folder_path = set_path+patient+"/"
                    gtf = open(folder_path+"GTF/"+patient+"_grch37.gtf")
                    gtf_ori = gtf.readlines()
                    temp_main_after_rna_mapped(gtf_ori,patient, set_num)


# In[163]:


"""
description:

arguments:
    gtf_ori
    patient_name
    set_num
    
output : 
"""
def temp_main_after_rna_mapped_frame_out(gtf_ori,patient_name, set_num):
    gtf_dict = make_gtf_dict(gtf_ori)
    filter_intron_rna_none_dict = {}
    filter_intron_rna = open("./OUTPUT/"+set_num+"/"+patient_name+"_whole_candi_after_samtools_RNA_mapped.bed")
    filter_intron_rna_ori = filter_intron_rna.readlines()
    for i in filter_intron_rna_ori[:]:
        filter_intron_rna_none_dict[i.replace("\n","")] = 1
    ee_dict = make_ee_dict(filter_intron_rna_none_dict)
    gene_to_transcript_dict, transcipt_to_info_dict = get_info_dict(ee_dict,gtf_dict)
    ee_process_frame_out(ee_dict, gene_to_transcript_dict,transcipt_to_info_dict,set_num, patient_name)


# In[164]:


"""
description:

arguments:
    gtf_ori
    patient_name
    set_num
    
output : 
"""
def temp_main_after_rna_mapped(gtf_ori,patient_name, set_num):
    gtf_dict = make_gtf_dict(gtf_ori)
    filter_intron_rna_none_dict = {}
    filter_intron_rna = open("./OUTPUT/"+set_num+"/"+patient_name+"_whole_candi_after_samtools_1016.bed")
    filter_intron_rna_ori = filter_intron_rna.readlines()
    for i in filter_intron_rna_ori[:]:
        filter_intron_rna_none_dict[i.replace("\n","")] = 1
    ee_dict = make_ee_dict(filter_intron_rna_none_dict)
    gene_to_transcript_dict, transcipt_to_info_dict = get_info_dict(ee_dict,gtf_dict)
    ee_process(ee_dict, gene_to_transcript_dict,transcipt_to_info_dict,set_num, patient_name)


# In[165]:


"""
description:

arguments:

output : 
"""
def final_candidate():
    file_list = os.listdir('./ACTG/')
    path = './ACTG/'   
    ri_dict_p={}
    es_dict_p={}
    ee_dict_p={}
    ri_dict_s={}
    es_dict_s={}
    ee_dict_s={}
    
    f_ri_dict_p={}
    f_es_dict_p={}
    f_ee_dict_p={}
    f_ri_dict_s={}
    f_es_dict_s={}
    f_ee_dict_s={}
    
    for set_num in file_list:
        if 'set' in set_num:
            print(set_num)
            set_path = path+set_num+"/"
            patient_list = os.listdir(set_path)



            for patient in patient_list:
                print(patient)

                if len(patient.split(".")) == 1:
                    start_time = time.time()

                    es_file = open('./OUTPUT/'+set_num+'/'+patient+"_final_es_dict_after_34.bed")
                    ri_file = open('./OUTPUT/'+set_num+'/'+patient+"_over_threshold_34_ri_candi_dict_after_34.bed")
                    da_file = open('./OUTPUT/'+set_num+'/'+patient+"_junction_dict_after_34.bed")

                    es_file_ori = es_file.readlines()
                    ri_file_ori = ri_file.readlines()
                    da_file_ori = da_file.readlines()
                    
                    ri_dict_p,es_dict_p,ee_dict_p,ri_dict_s,es_dict_s,ee_dict_s=final_candidate_process(da_file_ori, ri_file_ori,
                                                                                                        es_file_ori, set_num, patient,
                                                                                                        ri_dict_p,es_dict_p,ee_dict_p,
                                                                                                        ri_dict_s,es_dict_s,ee_dict_s)
                    
                    final_es_file = open('./FINAL/'+set_num+'/'+patient+"_final_es_after_j_after_making_sam.bed")
                    final_ri_file = open('./FINAL/'+set_num+'/'+patient+"_final_ri_after_s_p_after_making_sam.bed")
                    final_da_file = open('./FINAL/'+set_num+'/'+patient+"_final_ee_after_j_after_making_sam.bed")
                    
                    final_es_file_ori = final_es_file.readlines()
                    final_ri_file_ori = final_ri_file.readlines()
                    final_da_file_ori = final_da_file.readlines()
                    
                    f_ri_dict_p,f_es_dict_p,f_ee_dict_p,f_ri_dict_s,f_es_dict_s,f_ee_dict_s=final_candidate_process_2(final_da_file_ori,
                                                                                                          final_ri_file_ori,
                                                                                                        final_es_file_ori, set_num,
                                                                                                          patient,
                                                                                                        f_ri_dict_p,f_es_dict_p,
                                                                                                          f_ee_dict_p,
                                                                                                        f_ri_dict_s,f_es_dict_s,
                                                                                                          f_ee_dict_s,
                                                                                                        0,0)
                    
                    
                    
    print("________________________________________________________")
    print("before sam")
    get_num(ri_dict_p,es_dict_p,ee_dict_p,ri_dict_s,es_dict_s,ee_dict_s)
    print("________________________________________________________")
    print("after sam")
    get_num_2(f_ri_dict_p,f_es_dict_p,f_ee_dict_p,f_ri_dict_s,f_es_dict_s,f_ee_dict_s)
    
                    
                    


# In[166]:


"""
description:

arguments:

output : 
"""
def final_candidate_2():
    file_list = os.listdir('./ACTG/')
    path = './ACTG/'   
    ri_dict_p={}
    es_dict_p={}
    ee_dict_p={}
    ri_dict_s={}
    es_dict_s={}
    ee_dict_s={}
    
    f_ri_dict_p={}
    f_es_dict_p={}
    f_ee_dict_p={}
    f_ri_dict_s={}
    f_es_dict_s={}
    f_ee_dict_s={}
    
    for set_num in file_list:
        if 'set' in set_num:
            print(set_num)
            set_path = path+set_num+"/"
            patient_list = os.listdir(set_path)



            for patient in patient_list:
                print(patient)

                if len(patient.split(".")) == 1:
                    start_time = time.time()

                    es_file = open('./OUTPUT/'+set_num+'/'+patient+"_final_es_dict_after_34.bed")
                    ri_file = open('./OUTPUT/'+set_num+'/'+patient+"_over_threshold_34_ri_candi_dict_after_34.bed")
                    da_file = open('./OUTPUT/'+set_num+'/'+patient+"_junction_dict_after_34.bed")

                    es_file_ori = es_file.readlines()
                    ri_file_ori = ri_file.readlines()
                    da_file_ori = da_file.readlines()
                    
                    ri_dict_p,es_dict_p,ee_dict_p,ri_dict_s,es_dict_s,ee_dict_s=final_candidate_process(da_file_ori, ri_file_ori,
                                                                                                        es_file_ori, set_num, patient,
                                                                                                        ri_dict_p,es_dict_p,ee_dict_p,
                                                                                                        ri_dict_s,es_dict_s,ee_dict_s)
                    
                    final_es_file = open('./FINAL/'+set_num+'/'+patient+"_final_es_after_j_after_making_sam.bed")
                    final_ri_file = open('./FINAL/'+set_num+'/'+patient+"_final_ri_after_s_p_after_making_sam.bed")
                    final_da_file = open('./FINAL/'+set_num+'/'+patient+"_final_ee_after_j_after_making_sam.bed")
                    
                    final_es_file_ori = final_es_file.readlines()
                    final_ri_file_ori = final_ri_file.readlines()
                    final_da_file_ori = final_da_file.readlines()
                    
                    f_ri_dict_p,f_es_dict_p,f_ee_dict_p,f_ri_dict_s,f_es_dict_s,f_ee_dict_s=final_candidate_process_2(final_da_file_ori,
                                                                                                          final_ri_file_ori,
                                                                                                        final_es_file_ori, set_num,
                                                                                                          patient,
                                                                                                        f_ri_dict_p,f_es_dict_p,
                                                                                                          f_ee_dict_p,
                                                                                                        f_ri_dict_s,f_es_dict_s,
                                                                                                          f_ee_dict_s,
                                                                                                        4,2)
                    
                    
                    
    print("________________________________________________________")
    print("before sam")
    get_num(ri_dict_p,es_dict_p,ee_dict_p,ri_dict_s,es_dict_s,ee_dict_s)
    print("________________________________________________________")
    print("after sam")
    get_num_2(f_ri_dict_p,f_es_dict_p,f_ee_dict_p,f_ri_dict_s,f_es_dict_s,f_ee_dict_s)
    
                    
                    


# In[167]:


"""
description:

arguments:
    ri_dict_p
    es_dict_p
    ee_dict_p
    ri_dict_s
    es_dict_s
    ee_dict_s

output : 
"""
def get_num_2(ri_dict_p,es_dict_p,ee_dict_p,ri_dict_s,es_dict_s,ee_dict_s):
    print('ri : '+str(len(ri_dict_p)))
    print('ee : '+str(len(ee_dict_p)))
    print('es : '+str(len(es_dict_p)))
    ee_p_list = [0 for i in range(50)]
    ee_s_list = [0 for i in range(11)]
    
    es_p_list = [0 for i in range(50)]
    es_s_list = [0 for i in range(11)]
    
    ri_p_list = [0 for i in range(50)]
    ri_s_list = [0 for i in range(11)]
    
    print("final_ri")
    for i in ri_dict_p.keys():
        print(i)
        print(ri_dict_p.get(i))
        print(ri_dict_s.get(i))
        ri_p_list[len(ri_dict_p.get(i).split("//"))-1] += 1
        ri_s_list[len(ri_dict_s.get(i).split("//"))-1] += 1
    
    print("final_es")
    for i in es_dict_p.keys():
        print(i)
        print(es_dict_p.get(i))
        print(es_dict_s.get(i))
        es_p_list[len(es_dict_p.get(i).split("//"))-1] += 1
        es_s_list[len(es_dict_s.get(i).split("//"))-1] += 1
    
    print("final_ee")
    for i in ee_dict_p.keys():
        print(i)
        print(es_dict_p.get(i))
        print(es_dict_s.get(i))
        ee_p_list[len(ee_dict_p.get(i).split("//"))-1] += 1
        ee_s_list[len(ee_dict_s.get(i).split("//"))-1] += 1
    
    print('ee_p_list : '+str(ee_p_list))
    print('ee_s_list : '+str(ee_s_list))
    
    print('ri_p_list : '+str(ri_p_list))
    print('ri_s_list : '+str(ri_s_list))
    
    print('es_p_list : '+str(es_p_list))
    print('es_s_list : '+str(es_s_list))


# In[168]:


"""
description:

arguments:
    ri_dict_p
    es_dict_p
    ee_dict_p
    ri_dict_s
    es_dict_s
    ee_dict_s
output : 
"""
def get_num(ri_dict_p,es_dict_p,ee_dict_p,ri_dict_s,es_dict_s,ee_dict_s):
    print('ri : '+str(len(ri_dict_p)))
    print('ee : '+str(len(ee_dict_p)))
    print('es : '+str(len(es_dict_p)))
    ee_p_list = [0 for i in range(50)]
    ee_s_list = [0 for i in range(11)]
    
    es_p_list = [0 for i in range(50)]
    es_s_list = [0 for i in range(11)]
    
    ri_p_list = [0 for i in range(50)]
    ri_s_list = [0 for i in range(11)]
    
    for i in ri_dict_p.keys():
        ri_p_list[len(ri_dict_p.get(i).split("//"))-1] += 1
        ri_s_list[len(ri_dict_s.get(i).split("//"))-1] += 1
        
    for i in es_dict_p.keys():
        es_p_list[len(es_dict_p.get(i).split("//"))-1] += 1
        es_s_list[len(es_dict_s.get(i).split("//"))-1] += 1
        
    for i in ee_dict_p.keys():
        ee_p_list[len(ee_dict_p.get(i).split("//"))-1] += 1
        ee_s_list[len(ee_dict_s.get(i).split("//"))-1] += 1
    
    print('ee_p_list : '+str(ee_p_list))
    print('ee_s_list : '+str(ee_s_list))
    
    print('ri_p_list : '+str(ri_p_list))
    print('ri_s_list : '+str(ri_s_list))
    
    print('es_p_list : '+str(es_p_list))
    print('es_s_list : '+str(es_s_list))


# In[169]:


"""
description:

arguments:
    final_da_file_ori
    final_ri_file_ori
    final_es_file_ori
    set_num
    patient
    f_ri_dict_p
    f_es_dict_p
    f_ee_dict_p
    f_ri_dict_s
    f_es_dict_s
    f_ee_dict_s
    threshold_1
    threshold_2

output : 
    f_ri_dict_p
    f_es_dict_p
    f_ee_dict_p
    f_ri_dict_s
    f_es_dict_s
    f_ee_dict_s
"""
def final_candidate_process_2(final_da_file_ori, final_ri_file_ori, final_es_file_ori, set_num, patient,
                              f_ri_dict_p,f_es_dict_p,f_ee_dict_p,f_ri_dict_s,f_es_dict_s,f_ee_dict_s,threshold_1,threshold_2):
    
    
    if len(final_da_file_ori)>0:
        for line in final_da_file_ori:
            info = line.split("_")[-2]
            chromo = info.split(":")[0]
            start = info.split(":")[1].split("-")[0]
            end = info.split(":")[1].split("-")[1]
            read_num = int(line.split("_")[-1].split("num:")[1].replace("\n",""))
            key = chromo+"\t"+start+"\t"+end
            
            if read_num>=threshold_1:
            
                if f_ee_dict_p.get(key) == None:
                    f_ee_dict_p[key] = patient
                else:
                    check_list = f_ee_dict_p.get(key).split("//")
                    if patient not in check_list:
                        f_ee_dict_p[key] = f_ee_dict_p.get(key)+"//"+patient

                if f_ee_dict_s.get(key) == None:
                    f_ee_dict_s[key] = set_num
                else:
                    check_list = f_ee_dict_s.get(key).split("//")
                    if set_num not in check_list:
                        f_ee_dict_s[key] = f_ee_dict_s.get(key)+"//"+set_num
            
            

    if len(final_ri_file_ori)>0:
        for line in final_ri_file_ori:
            chromo = line.split("\t")[0]
            start = line.split("\t")[1]
            end = line.split("\t")[2]
            key = chromo+"\t"+start+"\t"+end
            read_num = int(line.split("_")[-1].split("num:")[1].replace("\n",""))
            
            if read_num>=threshold_2:
                if f_ri_dict_p.get(key) == None:
                    f_ri_dict_p[key] = patient
                else:
                    check_list = f_ri_dict_p.get(key).split("//")
                    if patient not in check_list:
                        f_ri_dict_p[key] = f_ri_dict_p.get(key)+"//"+patient

                if f_ri_dict_s.get(key) == None:
                    f_ri_dict_s[key] = set_num
                else:
                    check_list = f_ri_dict_s.get(key).split("//")
                    if set_num not in check_list:
                        f_ri_dict_s[key] = f_ri_dict_s.get(key)+"//"+set_num
    
                        
                        
    es_id_dict = {}     
    if len(final_es_file_ori)>0:
        for line in final_es_file_ori:
            chromo = line.split("\t")[0]
            start = line.split("\t")[1]
            end = line.split("\t")[2]
            now_id = line.split("_")[3]
            key = chromo+"\t"+start+"\t"+end
            read_num = int(line.split("_")[-1].split("num:")[1].replace("\n",""))
            
            if read_num>=threshold_1:
            
                if es_id_dict.get(now_id) == None:
                    es_id_dict[now_id] = key
                else:
                    es_id_dict[now_id] = es_id_dict.get(now_id)+"//"+key

                    key = es_id_dict.get(now_id)

                    if f_es_dict_p.get(key) == None:
                        f_es_dict_p[key] = patient
                    else:
                        check_list = f_es_dict_p.get(key).split("//")
                        if patient not in check_list:
                            f_es_dict_p[key] = f_es_dict_p.get(key)+"//"+patient

                    if f_es_dict_s.get(key) == None:
                        f_es_dict_s[key] = set_num
                    else:
                        check_list = f_es_dict_s.get(key).split("//")
                        if set_num not in check_list:
                            f_es_dict_s[key] = f_es_dict_s.get(key)+"//"+set_num
                    
    return f_ri_dict_p,f_es_dict_p,f_ee_dict_p,f_ri_dict_s,f_es_dict_s,f_ee_dict_s


# In[170]:


"""
description:

arguments:
    da_file_ori
    ri_file_ori
    es_file_ori
    set_num
    patient
    ri_dict_p
    es_dict_p
    ee_dict_p
    ri_dict_s
    es_dict_s
    ee_dict_s
output : 
    ri_dict_p
    es_dict_p
    ee_dict_p
    ri_dict_s
    es_dict_s
    ee_dict_s
"""
def final_candidate_process(da_file_ori, ri_file_ori, es_file_ori, set_num, patient,ri_dict_p,es_dict_p,ee_dict_p
                           ,ri_dict_s,es_dict_s,ee_dict_s):
    
    
    if len(da_file_ori)>0:
        for line in da_file_ori:
            info = line.split("_")[-1]
            chromo = info.split(":")[0]
            start = info.split(":")[1].split("-")[0]
            end = info.split(":")[1].split("-")[1].replace("\n","")
            key = chromo+"\t"+start+"\t"+end
            
            if ee_dict_p.get(key) == None:
                ee_dict_p[key] = patient
            else:
                check_list = ee_dict_p.get(key).split("//")
                if patient not in check_list:
                    ee_dict_p[key] = ee_dict_p.get(key)+"//"+patient
                    
            if ee_dict_s.get(key) == None:
                ee_dict_s[key] = set_num
            else:
                check_list = ee_dict_s.get(key).split("//")
                if set_num not in check_list:
                    ee_dict_s[key] = ee_dict_s.get(key)+"//"+set_num
            
            
            
    if len(ri_file_ori)>0:
        for line in ri_file_ori:
            chromo = line.split("\t")[0]
            start = line.split("\t")[1]
            end = line.split("\t")[2]
            key = chromo+"\t"+start+"\t"+end
            
            if ri_dict_p.get(key) == None:
                ri_dict_p[key] = patient
            else:
                check_list = ri_dict_p.get(key).split("//")
                if patient not in check_list:
                    ri_dict_p[key] = ri_dict_p.get(key)+"//"+patient
                    
            if ri_dict_s.get(key) == None:
                ri_dict_s[key] = set_num
            else:
                check_list = ri_dict_s.get(key).split("//")
                if set_num not in check_list:
                    ri_dict_s[key] = ri_dict_s.get(key)+"//"+set_num
    es_id_dict = {}     
    if len(es_file_ori)>0:
        for line in es_file_ori:
            chromo = line.split("\t")[0]
            start = line.split("\t")[1]
            end = line.split("\t")[2]
            now_id = line.split("_")[3]
            key = chromo+"\t"+start+"\t"+end
            
            if es_id_dict.get(now_id) == None:
                es_id_dict[now_id] = key
            else:
                es_id_dict[now_id] = es_id_dict.get(now_id)+"//"+key
                
                key = es_id_dict.get(now_id)
            
                if es_dict_p.get(key) == None:
                    es_dict_p[key] = patient
                else:
                    check_list = es_dict_p.get(key).split("//")
                    if patient not in check_list:
                        es_dict_p[key] = es_dict_p.get(key)+"//"+patient

                if es_dict_s.get(key) == None:
                    es_dict_s[key] = set_num
                else:
                    check_list = es_dict_s.get(key).split("//")
                    if set_num not in check_list:
                        es_dict_s[key] = es_dict_s.get(key)+"//"+set_num
                    
    return ri_dict_p,es_dict_p,ee_dict_p,ri_dict_s,es_dict_s,ee_dict_s


# In[171]:


"""
description:

arguments:
    patient_num
    set_num

output : 
"""
def after_rna_whole_one(patient_num, set_num):
    bed_dict = {}
    bed_dict = make_bed_dict_with_after_bam(patient_num, set_num)
    rna_mapping_dict = make_rna_count_dict(patient_num, set_num)

    #rna 전부 0이면 날려
    candi_with_rna_dict = make_candi_with_rna_dict_bed_dict(rna_mapping_dict, bed_dict)

    #intron 부분 rna 전부 0이면 날려
    filter_intron_rna_none_dict = make_filter_intron_rna_none_dict(candi_with_rna_dict)
    
    ee_file = open("./OUTPUT/"+set_num+"/"+patient_num+"_whole_candi_after_samtools_1016.bed","wt")
    for i in filter_intron_rna_none_dict.keys():
        ee_file.write(i+"\n")
    ee_file.close()


# In[172]:


"""
description:

arguments:

output : 
"""
def after_rna_whole():
    file_list = os.listdir('./ACTG/')
    for set_num in file_list:
        if 'set' in set_num:
            print(set_num)
            set_path = path+set_num+"/"
            patient_list = os.listdir(set_path)
            for patient in patient_list:
                if len(patient.split(".")) == 1:
                    print(patient)
                    after_rna_whole_one(patient, set_num)


# In[173]:


"""
description:

arguments:

output : 
"""
def get_whole_after_actg():
    file_list = os.listdir('./ACTG/')
    path = './ACTG/'
    
    
    
    flat_with_gff_dict = {}
    bed_dict = {}
    
    for set_num in file_list:
        if 'set' in set_num:
            print(set_num)
            set_path = path+set_num+"/"
            patient_list = os.listdir(set_path)
            for patient in patient_list:
                print(patient)
                if len(patient.split(".")) == 1:
                    folder_path = set_path+patient+"/"
                    
                    output_list = os.listdir(folder_path+"output/")
                    flat_ori = ''
                    gff_ori = ''
                    for output_file in output_list:
                        if output_file.split(".")[1] == 'flat':
                            flat = open(folder_path+"output/"+output_file)
                            flat_ori = flat.readlines()
                        elif output_file.split(".")[1] == 'gff':
                            gff = open(folder_path+"output/"+output_file)
                            gff_ori = gff.readlines()
                    now_dict = make_flat_with_gff_dict_actg(flat_ori, gff_ori)
                    make_bed_and_dict_with_now_full_actg(now_dict,bed_dict, patient, set_num)
    whole_after_actg = open("whole_after_actg.bed","wt")
    for i in bed_dict.keys():
        whole_after_actg.write(i)
    whole_after_actg.close()


# In[174]:


"""
description:

arguments:

output : 
"""
def start_rna_evidence_one(patient_name, bam_path):
    file_list = os.listdir('./ACTG/')
    path = './ACTG/'
    
    #fasta 파일 다 읽어들여와!
    #fasta_dict = make_fasta_dict()
    
    for set_num in file_list:
        if 'set' in set_num:
            print(set_num)
            set_path = path+set_num+"/"
            patient_list = os.listdir(set_path)
            for patient in patient_list:
                if patient == patient_name:
                    if len(patient.split(".")) == 1:
                        folder_path = set_path+patient+"/"
                        gtf = open(folder_path+"GTF/"+patient+"_grch37.gtf")
                        gtf_ori = gtf.readlines()
                        output_list = os.listdir(folder_path+"output/")
                        flat_ori = ''
                        gff_ori = ''
                        for output_file in output_list:
                            if output_file.split(".")[1] == 'flat':
                                flat = open(folder_path+"output/"+output_file)
                                flat_ori = flat.readlines()
                            elif output_file.split(".")[1] == 'gff':
                                gff = open(folder_path+"output/"+output_file)
                                gff_ori = gff.readlines()
                        temp_main_one(bam_path, gtf_ori, gff_ori, flat_ori, patient, set_num)
                        


# In[175]:


"""
description:

arguments:

output : 
"""
def start_rna_evidence_intron():
    file_list = os.listdir('./ACTG/')
    path = './ACTG/'
    
    #fasta 파일 다 읽어들여와!
    #fasta_dict = make_fasta_dict()
    
    for set_num in file_list:
        if 'set' in set_num:
            print(set_num)
            set_path = path+set_num+"/"
            patient_list = os.listdir(set_path)
            for patient in patient_list:
                print(patient)
                if len(patient.split(".")) == 1:
                    folder_path = set_path+patient+"/"
                    gtf = open(folder_path+"GTF/"+patient+"_grch37.gtf")
                    gtf_ori = gtf.readlines()
                    output_list = os.listdir(folder_path+"output/")
                    flat_ori = ''
                    gff_ori = ''
                    for output_file in output_list:
                        if output_file.split(".")[1] == 'flat':
                            flat = open(folder_path+"output/"+output_file)
                            flat_ori = flat.readlines()
                        elif output_file.split(".")[1] == 'gff':
                            gff = open(folder_path+"output/"+output_file)
                            gff_ori = gff.readlines()
                    temp_main_intron(gtf_ori, gff_ori, flat_ori, patient, set_num)
                        


# In[176]:


def get_ri_trans(final_ri_ori):
    final_ri_dict = {}
    final_ri_both_dict = {}
    before_intronPoint = ''
    if len(final_ri_ori)>0:
        for i in final_ri_ori[:]:
            key_list = []
            chromo = i.split("\t")[0]
            key_list.append(chromo)
            start = i.split("\t")[1]
            key_list.append(start)
            end = i.split("\t")[2]
            key_list.append(end)
            p_m = i.split("\t")[3].split("_")[0]
            key_list.append(p_m)
            pep = i.split("\t")[3].split("_")[4]
            key_list.append(pep)
            gene = i.split("\t")[3].split("_")[6]
            key_list.append(gene)
            intronPoint = i.split("\t")[3].split("_")[-3].split("intronPoint:")[1]
#             key_list.append(intronPoint)
            both_info = i.split("intronPoint:")[1]
            
            key = chromo+"\t"+start+"\t"+end
            if final_ri_dict.get(key) == None:
                final_ri_dict[key] = i.replace("\n","")
                before_intronPoint = intronPoint
            else:
                now_value = final_ri_dict.get(key)+"_intronPoint:"+both_info
                if before_intronPoint != '':
                    now_key = "_".join(key_list)+"_"+before_intronPoint+"-"+intronPoint
                    print(now_key)
                    before_intronPoint = ''
                    final_ri_both_dict[now_key] = now_value
    return final_ri_both_dict


# In[177]:


"""
final_ri_both_dict
chromo_start_end_PM_pep_gene_intronPoint-intronPoint
line(intronPoint:123_num:1_meanDepth:1.0_intronPoint:456_num:1_meanDepth:1.0\n)
"""
def get_novel_transcript_ri(gtf_ori,patient_name,set_num, ri_dict):
    gtf_dict = make_gtf_dict(gtf_ori)
    final_ri = open("./FINAL/"+set_num+"/"+patient_name+"_final_ri_after_s_p_after_making_sam.bed")
    final_ri_ori = final_ri.readlines()
    
    both_dict_ri = {}
    both_dict_intron = {}

    both_dict_ri = get_ri_trans(final_ri_ori)
#     both_dict_intron = get_ri_trans(final_ri_intron_ori)
    
#     both_dict_ri.update(both_dict_intron)
#     print(len(both_dict_ri))
    
    gene_to_transcript_dict, transcipt_to_info_dict = get_info_dict_for_mean_depth(both_dict_ri,gtf_dict)
    find_original_exon_point_ri(both_dict_ri,ri_dict, gene_to_transcript_dict,transcipt_to_info_dict,set_num, patient_name)


# In[178]:


"""
ri_dict

"""

def start_novel_transcript_evidence():
    file_list = os.listdir('./ACTG/')
    path = './ACTG/'
    
    
    for set_num in file_list:
        if 'set' in set_num:
            #print(set_num)
            set_path = path+set_num+"/"
            patient_list = os.listdir(set_path)
            for patient in patient_list:
                print(patient)
                start_time = time.time()
                ee_for_compare_dict = {}
                es_for_compare_dict = {}
                ri_for_Novel_Transcript_dict = {}
                if len(patient.split(".")) == 1:
                    
                    folder_path = set_path+patient+"/"
                    gtf = open(folder_path+"GTF/"+patient+"_grch37.gtf")
                    gtf_ori = gtf.readlines()
                    
                    get_novel_transcript_ri(gtf_ori,patient,set_num, ri_for_Novel_Transcript_dict)

                    get_mean_depth_ee(gtf_ori,patient,set_num, ee_for_compare_dict)
                    get_mean_depth_es(gtf_ori,patient,set_num, es_for_compare_dict)
                    
                    #print("one patients samtools calculate points :", time.time() - start_time)
                    
                    start_time = time.time()

                    make_bed_for_mean_depth_ee_es(ee_for_compare_dict, es_for_compare_dict, patient, set_num)
                    
                    #print("one patients samtools calculate mean_depth :", time.time() - start_time)


                    
                    print_after_sam_file_compare(ri_for_Novel_Transcript_dict,set_num,patient)
                    print_after_sam_file_compare(es_for_compare_dict,set_num,patient)
                    print_after_sam_file_compare(ee_for_compare_dict,set_num,patient)

                    throw_garbage()
                    
                    
def start_novel_transcript_evidence_one(patient_name):
    file_list = os.listdir('./ACTG/')
    path = './ACTG/'
    
    
    for set_num in file_list:
        if 'set' in set_num:
            #print(set_num)
            set_path = path+set_num+"/"
            patient_list = os.listdir(set_path)
            for patient in patient_list:
                if patient_name == patient:
                    start_time = time.time()
                    ee_for_compare_dict = {}
                    es_for_compare_dict = {}
                    ri_for_Novel_Transcript_dict = {}
                    if len(patient.split(".")) == 1:
                        
                        folder_path = set_path+patient+"/"
                        gtf = open(folder_path+"GTF/"+patient+"_grch37.gtf")
                        gtf_ori = gtf.readlines()
                        
                        get_novel_transcript_ri(gtf_ori,patient,set_num, ri_for_Novel_Transcript_dict)

                        get_mean_depth_ee(gtf_ori,patient,set_num, ee_for_compare_dict)
                        get_mean_depth_es(gtf_ori,patient,set_num, es_for_compare_dict)
                        
                        #print("one patients samtools calculate points :", time.time() - start_time)
                        
                        start_time = time.time()

                        make_bed_for_mean_depth_ee_es(ee_for_compare_dict, es_for_compare_dict, patient, set_num)
                        
                        #print("one patients samtools calculate mean_depth :", time.time() - start_time)


                        
                        print_after_sam_file_compare(ri_for_Novel_Transcript_dict,set_num,patient)
                        print_after_sam_file_compare(es_for_compare_dict,set_num,patient)
                        print_after_sam_file_compare(ee_for_compare_dict,set_num,patient)

                        throw_garbage()


# In[179]:


"""
ee_dict
key
chromo_start_end_NovelExonPoint_exonPoint_PM_FrontBack_gene_pepSeq_transcript_cdsMap_start-end
cds_map = 123:456//457:789 (모든 엑손부터 순서대로 저장되어 있음 단, novel exon 부분만 수정된채로 집어넣음)

value
final line1

es_dict
key
chromo_start1_end1_start2_end2_skipStart_skipEnd_pM_gene_pepSeq_transcript_cds_map_first(123-456)_second(4654-8798)
cds_map = 123:456//457:789 (모든 엑손부터 순서대로 저장되어 있음 단, skipped 부분만 삭제)

value
final line1+"%%"+final line2

"""

def start_compare_mean_depth_evidence():
    file_list = os.listdir('./ACTG/')
    path = './ACTG/'
    
#     #fasta 파일 다 읽어들여와!
#     fasta_dict = make_fasta_dict()
    
    for set_num in file_list:
        if 'set' in set_num:
            print(set_num)
            set_path = path+set_num+"/"
            patient_list = os.listdir(set_path)
            for patient in patient_list:
                print(patient)
                start_time = time.time()
                ee_for_compare_dict = {}
                es_for_compare_dict = {}
                if len(patient.split(".")) == 1:
                    
                    folder_path = set_path+patient+"/"
                    gtf = open(folder_path+"GTF/"+patient+"_grch37.gtf")
                    gtf_ori = gtf.readlines()

                    get_mean_depth_ee(gtf_ori,patient,set_num, ee_for_compare_dict)
                    get_mean_depth_es(gtf_ori,patient,set_num, es_for_compare_dict)
                    
                    print("one patients samtools calculate points :", time.time() - start_time)
                    
                    start_time = time.time()

                    make_bed_for_mean_depth_ee_es(ee_for_compare_dict, es_for_compare_dict, patient, set_num)
                    
                    print("one patients samtools calculate mean_depth :", time.time() - start_time)

#                     print_after_sam_file(ee_for_compare_dict,set_num,patient)
#                     print_after_sam_file(es_for_compare_dict,set_num,patient)
                    
                    print_after_sam_file_compare(ee_for_compare_dict,set_num,patient)
                    print_after_sam_file_compare(es_for_compare_dict,set_num,patient)
                    

                    throw_garbage()


# In[180]:


"""
ee_dict
key
chromo_start_end_NovelExonPoint_exonPoint_PM_FrontBack_gene_pepSeq_transcript_cdsMap_start-end
cds_map = 123:456//457:789 (모든 엑손부터 순서대로 저장되어 있음 단, novel exon 부분만 수정된채로 집어넣음)

value
final line1

es_dict
key
chromo_start1_end1_start2_end2_skipStart_skipEnd_pM_gene_pepSeq_transcript_cds_map_first(123-456)_second(4654-8798)
cds_map = 123:456//457:789 (모든 엑손부터 순서대로 저장되어 있음 단, skipped 부분만 삭제)

value
final line1+"%%"+final line2

"""
def make_bed_for_mean_depth_ee_es(ee_dict, es_dict, patient, set_num):
    mean_depth = open('OUTPUT/'+set_num+"/mean_depth_"+patient+'.bed','wt')
    check_dict = {}
    
    #input bed 생성
    for i in ee_dict.keys():
        info = i
        chromo = i.split("_")[0]
        start = str(int(i.split("_")[-1].split("-")[0])-1)
        end = i.split("_")[-1].split("-")[1]
        check_str = chromo+"\t"+start+"\t"+end
        line = chromo+"\t"+start+"\t"+end+"\t"+info+"\n"
        if check_dict.get(check_str) == None:
            check_dict[check_str] = 1
            mean_depth.write(line)
            
    for i in es_dict.keys():
        info = i
        chromo = i.split("_")[0]
        start1 = str(int(i.split("_")[-2].split("-")[0])-1)
        end1 = i.split("_")[-2].split("-")[1]
        start2 = str(int(i.split("_")[-1].split("-")[0])-1)
        end2 = i.split("_")[-1].split("-")[1]
        start3 = str(int(i.split("_")[5])-1)
        end3 = i.split("_")[6]
        check_str1 = chromo+"\t"+start1+"\t"+end1
        check_str2 = chromo+"\t"+start2+"\t"+end2
        check_str3 = chromo+"\t"+start3+"\t"+end3
        line1 = chromo+"\t"+start1+"\t"+end1+"\t"+info+"\n"
        line2 = chromo+"\t"+start2+"\t"+end2+"\t"+info+"\n"
        line3 = chromo+"\t"+start3+"\t"+end3+"\t"+info+"\n"
        if check_dict.get(check_str1) == None:
            check_dict[check_str1] = 1
            mean_depth.write(line1)
        if check_dict.get(check_str2) == None:
            check_dict[check_str2] = 1
            mean_depth.write(line2)
        if check_dict.get(check_str3) == None:
            check_dict[check_str3] = 1
            mean_depth.write(line3)
    mean_depth.close()
    
#     #samtools depth 실행
#     bam_file_path = "./BAM/"+set_num+"/"
#     file_list = os.listdir(bam_file_path)
#     bam_file_name = ''
#     for file in file_list:
#         if patient in file:
#             if file.split(".")[-1] == 'bam':
#                 bam_file_name = bam_file_path+file
    
#     bed_file_name = 'OUTPUT/'+set_num+"/mean_depth_"+patient+'.bed'
    
#     output_file_name = 'OUTPUT/'+set_num+"/mean_depth_"+patient+'_count_for_compare.txt'
    
#     samtools_input = 'samtools depth -q 0 '+bam_file_name+' -b '+bed_file_name+' -o '+output_file_name
#     r = subprocess.Popen(samtools_input, shell=True).wait()
#     if r == 1: 
#         print("samtools depth process failed")
    
    
    #samtools depth 실행 결과 적용 데이터 출력
#     mapping_compare_dict = make_rna_count_dict_for_final_mean_depth_compare(patient,set_num)
    for i in check_dict.keys():
        chromo = i.split("\t")[0]
        start = int(i.split("\t")[1])+1
        end = int(i.split("\t")[2])
        depth = 0
        for spot in range(start,end+1):
            key = chromo+":"+str(spot)
#             if mapping_compare_dict.get(key) != None:
#                 depth += int(mapping_compare_dict.get(key))
        mean = depth/(end-start+1)
        check_dict[i] = '_meanDepth:'+str(round(mean,4))
        
        
    for i in ee_dict.keys():
#         print('ee')
#         print('ee_dict.key')
#         print(i)
#         print('ee_dict.value')
#         print(ee_dict.get(i))
#         print('ee_dict.changed_value')
#         print(ee_dict.get(i).replace("\n","")+check_dict.get(check_str)+"\n")
        chromo = i.split("_")[0]
        start = str(int(i.split("_")[-1].split("-")[0])-1)
        end = i.split("_")[-1].split("-")[1]
        check_str = chromo+"\t"+start+"\t"+end

        
        ee_dict[i] = ee_dict.get(i).replace("\n","")+check_dict.get(check_str)+"\n"
#         print('ee')
        
    for i in es_dict.keys():
#         print('es')
#         print('es.key')
#         print(i)
#         print('es.value')
#         print(es_dict.get(i))
#         print('es.changed_value')
        chromo = i.split("_")[0]
        start1 = str(int(i.split("_")[-2].split("-")[0])-1)
        end1 = i.split("_")[-2].split("-")[1]
        start2 = str(int(i.split("_")[-1].split("-")[0])-1)
        end2 = i.split("_")[-1].split("-")[1]
        start3 = str(int(i.split("_")[5])-1)
        end3 = i.split("_")[6]
        check_str1 = chromo+"\t"+start1+"\t"+end1
        check_str2 = chromo+"\t"+start2+"\t"+end2
        check_str3 = chromo+"\t"+start3+"\t"+end3
        es_dict[i] = es_dict.get(i).replace("\n","")+check_dict.get(check_str1)+check_dict.get(check_str2)+check_dict.get(check_str3)+"\n"
#         print('trio 123')
#         print(check_dict.get(check_str1))
#         print(check_dict.get(check_str2))
#         print(check_dict.get(check_str3))
#         print('trio 123')
#         print('es.changed_value')
#         print(es_dict.get(i).replace("\n","")+check_dict.get(check_str1)+check_dict.get(check_str2)+check_dict.get(check_str3)+"\n")
#         print('es')
        
    
        
#     remove_output = 'rm '+output_file_name
#     r = subprocess.Popen(remove_output, shell=True).wait()
#     if r == 1: 
#         print("deleting output failed")
        
#    remove_output = 'rm '+bed_file_name
#    r = subprocess.Popen(remove_output, shell=True).wait()
#    if r == 1: 
#        print("deleting bed failed")
        
            


# In[181]:


"""
final_es_dict_with_info
chromo_start1_end1_start2_end2_skipStart_skipEnd_pM_gene_pepSeq
"""
def find_original_exon_point_es(final_es_dict_with_info, es_dict, gene_to_transcript_dict,transcipt_to_info_dict,set_num, patient_name):
    
    for key in final_es_dict_with_info.keys():
        gene = key.split("_")[-2]
#         front_back = key.split("_")[-3]
#         exon_point = key.split("_")[-5]
        p_m = key.split("_")[-3]
        chromo = key.split("_")[0]
        if p_m == '+':
            one = key.split("_")[2]
            two = str(int(key.split("_")[5]))
            three = key.split("_")[6]
            four = str(int(key.split("_")[3])+1)
        elif p_m == '-':
            one = str(int(key.split("_")[3])+1)
            two = str(int(key.split("_")[5]))
            three = key.split("_")[6]
            four = key.split("_")[2]
#         print("-------------")
#         print(p_m)
#         print(chromo)
#         print(one)
#         print(two)
#         print(three)
#         print(four)
#         print("-------------")

        transcripts = gene_to_transcript_dict.get(gene).split("\n")
        for transcript in transcripts:
            transcript_info = transcipt_to_info_dict.get(transcript)
            first, second, cds_map = deep_original_exon_point_es(transcript_info, p_m, one,two, three, four)

            if second != -1:
                input_key = "_"+transcript+"_"+cds_map+"_"+first+"_"+second
#                 print(input_key)
                es_dict[key+input_key] = final_es_dict_with_info.get(key)


# In[182]:


"""
final_es_dict_with_info
chromo_start1_end1_start2_end2_skipStart_skipEnd_pM_gene_pepSeq
"""
def get_mean_depth_es(gtf_ori,patient_name,set_num, es_dict):
    gtf_dict = make_gtf_dict(gtf_ori)
    final_es = open("./FINAL/"+set_num+"/"+patient_name+"_final_es_after_j_after_making_sam.bed")
    final_es_ori = final_es.readlines()
    
    final_es_id_dict = {}
#     final_es_dict = {}
    final_es_dict_with_info = {}
    
    if len(final_es_ori)>0:
        for line in final_es_ori[:]:
            es_id = line.split("_")[3]
            if final_es_id_dict.get(es_id) == None:
                final_es_id_dict[es_id] = line.replace("\n","")
            else:
                final_es_id_dict[es_id] = final_es_id_dict.get(es_id)+"%%"+line.replace("\n","")
    for key in final_es_id_dict.keys():
        es_list = []
        info_es = final_es_id_dict.get(key)
        chromo = info_es.split("\t")[0]
        es_list.append(chromo)
        start1 = info_es.split("\t")[1]
        es_list.append(start1)
        end1 = info_es.split("\t")[2]
        es_list.append(end1)
        start2 = info_es.split("%%")[1].split("\t")[1]
        es_list.append(start2)
        end2 = info_es.split("%%")[1].split("\t")[2]
        es_list.append(end2)
        skipped = info_es.split("_")[-5]
        s_start = skipped.split(":")[1]
        es_list.append(s_start)
        s_end = skipped.split(":")[2]
        es_list.append(s_end)
        p_m = info_es.split("\t")[3].split("_")[0]
        es_list.append(p_m)
        gene = info_es.split("\t")[3].split("_")[6]
        es_list.append(gene)
        pep = info_es.split("\t")[3].split("_")[4]
        es_list.append(pep)
        input_key = '_'.join(es_list)
        final_es_dict_with_info[input_key] = final_es_id_dict.get(key)
#         print(input_key)
        
    
    gene_to_transcript_dict, transcipt_to_info_dict = get_info_dict_for_mean_depth(final_es_dict_with_info,gtf_dict)
    
    find_original_exon_point_es(final_es_dict_with_info, es_dict, gene_to_transcript_dict,transcipt_to_info_dict,set_num, patient_name)


# In[183]:


def deep_original_exon_point_ri(transcript_info, p_m, intron_start, intron_end):
    CDSs = transcript_info.split("\tCDS\t")
    check = False
    check_first = False
    check_second = False
    novel_start = ''
    novel_end = ''
    cds_map = []
    first = ''
    second = ''
    if p_m == '+':
        first = intron_start
        second = intron_end
    elif p_m == '-':
        first = intron_end
        second = intron_start
    
    for CDS in CDSs[1:]:
        cds_start = CDS.split("\t")[0]
        cds_end = CDS.split("\t")[1]
        if check_first == False:
            if first in CDS:
                if p_m == '+':
                    if cds_end == first:
                        novel_start = cds_start
                        check_first = True
                else:
                    if cds_start == first:
                        novel_end = cds_end
                        check_first = True
            else:
                now_key = cds_start+":"+cds_end
                cds_map.append(now_key)
        else:
            if second in CDS:
                if p_m == '+':
                    if cds_start == second:
                        novel_end = cds_end
                        check_second = True
                    if check_second:
                        now_key = novel_start+":"+novel_end
                        cds_map.append(now_key)
                else:
                    if cds_end == second:
                        novel_start = cds_start
                        check_second = True
                    if check_second:
                        now_key = novel_start+":"+novel_end
                        cds_map.append(now_key)
            else:
                now_key = cds_start+":"+cds_end
                cds_map.append(now_key)
    if check_second:
        return "//".join(cds_map)
    else:
        return None


# In[184]:


"""
final_ri_both_dict
chromo_start_end_PM_pep_gene_intronPoint-intronPoint
line(intronPoint:123_num:1_meanDepth:1.0_intronPoint:456_num:1_meanDepth:1.0\n)
"""

def find_original_exon_point_ri(final_ri_dict_with_info,ri_dict, gene_to_transcript_dict,transcipt_to_info_dict,set_num, patient_name):
#     print(len(final_ee_dict_with_info))
    for key in final_ri_dict_with_info.keys():
#         print(key)
        gene = key.split("_")[-2]

        p_m = key.split("_")[-4]
        intron_start = str(int(key.split("_")[-1].split("-")[0])-1)
        intron_end = str(int(key.split("_")[-1].split("-")[1])+1)
        
        
        
        
        transcripts = gene_to_transcript_dict.get(gene).split("\n")
        for transcript in transcripts:
            transcript_info = transcipt_to_info_dict.get(transcript)
            cds_map = deep_original_exon_point_ri(transcript_info, p_m, intron_start, intron_end)

            if cds_map != None:
#                 print(cds_map)
                input_key = "_"+transcript+"_"+cds_map
                ri_dict[key+input_key] = final_ri_dict_with_info.get(key)


# In[185]:


"""
chromo_start_end_NovelExonPoint_exonPoint_PM_FrontBack_gene_pepSeq
chr17	7150538	7150540	-_0_acceptor_chr17:7150490-7150514//ori_jnc:387//ori_exact:12//novel_jnc:2//novel_exact:0(0.005)_num:2
"""
def find_original_exon_point_ee(final_ee_dict_with_info, ee_dict, gene_to_transcript_dict,transcipt_to_info_dict,set_num, patient_name):
#     print(len(final_ee_dict_with_info))
    for key in final_ee_dict_with_info.keys():
#         print(key)
        gene = key.split("_")[-2]
        front_back = key.split("_")[-3]
        exon_point = key.split("_")[-5]
        p_m = key.split("_")[-4]
        chromo = key.split("_")[0]
        candi_start = str(int(key.split("_")[1])+1)
        candi_end = key.split("_")[2]
        novel_exon_point = key.split("_")[3]
        
        transcripts = gene_to_transcript_dict.get(gene).split("\n")
        for transcript in transcripts:
            transcript_info = transcipt_to_info_dict.get(transcript)
            start, end, cds_map = deep_original_exon_point_ee(transcript_info, p_m, front_back, exon_point, candi_start, candi_end,novel_exon_point)

            if start != -1:
                if front_back == 'Front' and p_m == '+':
                    input_key1 = "_"+transcript+"_"+cds_map+"_"+start+'-'+exon_point
                    input_key2 = "_"+transcript+"_"+cds_map+"_"+str(int(exon_point)+1)+'-'+novel_exon_point
                elif front_back == 'Back' and p_m == '-':
                    input_key1 = "_"+transcript+"_"+cds_map+"_"+start+'-'+exon_point
                    input_key2 = "_"+transcript+"_"+cds_map+"_"+str(int(exon_point)+1)+'-'+novel_exon_point
                elif front_back == 'Front' and p_m == '-':
                    input_key1 = "_"+transcript+"_"+cds_map+"_"+novel_exon_point+'-'+str(int(exon_point)-1)
                    input_key2 = "_"+transcript+"_"+cds_map+"_"+exon_point+'-'+end
                elif front_back == 'Back' and p_m == '+':
                    input_key1 = "_"+transcript+"_"+cds_map+"_"+novel_exon_point+'-'+str(int(exon_point)-1)
                    input_key2 = "_"+transcript+"_"+cds_map+"_"+exon_point+'-'+end
                ee_dict[key+input_key1] = final_ee_dict_with_info.get(key)
                ee_dict[key+input_key2] = final_ee_dict_with_info.get(key)


# In[186]:


def deep_original_exon_point_es(transcript_info, p_m, one,two, three, four):
    first = -1
    second = -1
    cds_map = []
    check_one = False
    check_two = False


#     if one in transcript_info and two in transcript_info and three in transcript_info and four in transcript_info:
    if one in transcript_info and four in transcript_info:
        CDSs = transcript_info.split("\tCDS\t")
        if p_m == '-':
            for CDS in CDSs[1:]:
                cds_start = CDS.split("\t")[0]
                cds_end = CDS.split("\t")[1]
                key = cds_start+":"+cds_end
                if check_one == False:
                    cds_map.append(key)
                    if one == cds_start:
                        first = one+"-"+cds_end
                        check_one = True    
                elif check_one:
                    if check_two == False:
                        check_two = True
                    elif check_two:
                        cds_map.append(key)
                        if four == cds_end:
                            second = cds_start+"-"+four
        else:
            for CDS in CDSs[1:]:
                cds_start = CDS.split("\t")[0]
                cds_end = CDS.split("\t")[1]
                key = cds_start+":"+cds_end
                if check_one == False:
                    cds_map.append(key)
                    if one == cds_end:
                        first = cds_start+"-"+one
                        check_one = True
                elif check_one:
                    if check_two == False:
                        check_two = True
                    elif check_two:
                        cds_map.append(key)
                        if four == cds_start:
                            second = four+"-"+cds_end
    return first, second, "//".join(cds_map)


# In[187]:


def deep_original_exon_point_ee(transcript_info, p_m, front_back, exon_point, candi_start, candi_end,novel_exon_point):
    CDSs = transcript_info.split("\tCDS\t")
    check = False
    cds_map = []
    start = -1
    end = -1
    last_key = ''
    for CDS in CDSs[1:]:
        cds_start = CDS.split("\t")[0]
        cds_end = CDS.split("\t")[1]
        if cds_start == exon_point or cds_end == exon_point:
            check = True
            cds_start = CDS.split("\t")[0]
            cds_end = CDS.split("\t")[1]
#             if cds_start == '187522423':
#                 print("i am here check = exon_point")
#                 print(exon_point)
#                 print(CDS)
#                 print(exon_point)
#                 print("i am here check = exon_point")
            
#             cds_map.append(last_key)
            
            if p_m == '-' and front_back == 'Front':
                if cds_start == exon_point:
                    now_key = novel_exon_point+":"+cds_end
                    cds_map.append(now_key)
                    end = cds_end
                    start = novel_exon_point
            elif p_m == '+' and front_back == 'Back':
                if cds_start == exon_point:
                    now_key = novel_exon_point+":"+cds_end
                    cds_map.append(now_key)
                    end = cds_end
                    start = novel_exon_point
            elif p_m == '+' and front_back == 'Front':
                if cds_end == exon_point:
                    now_key = cds_start+":"+novel_exon_point
                    cds_map.append(now_key)
                    end = novel_exon_point
                    start = cds_start
            elif p_m == '-' and front_back == 'Back':
                if cds_end == exon_point:
                    now_key = cds_start+":"+novel_exon_point
                    cds_map.append(now_key)
                    end = novel_exon_point
                    start = cds_start
        elif check == False:
            cds_start = CDS.split("\t")[0]
            cds_end = CDS.split("\t")[1]
            last_key = cds_start+":"+cds_end
#             if cds_start == '187522423':
#                 print("i am here check = False")
#                 print(last_key)
#                 print("i am here check = False")
            cds_map.append(last_key)
        elif check:
            if start != -1:
                cds_start = CDS.split("\t")[0]
                cds_end = CDS.split("\t")[1]
                now_key = cds_start+":"+cds_end
#                 if cds_start == '187522423':
#                     print("i am here_check=True")
#                     print(now_key)
#                     print("i am here_check=True")
                cds_map.append(now_key)
#     input_key = start+"-"+end
    return start,end, "//".join(cds_map)


# In[188]:


"""
chr17_7150490__0_acceptor_chr17:7150490_7150501_-_Back_ENSG00000175826_HSGSAQVK
chromo_start_end_exonPoint_PM_FrontBack_gene_pepSeq
chr17	7150538	7150540	-_0_acceptor_chr17:7150490-7150514//ori_jnc:387//ori_exact:12//novel_jnc:2//novel_exact:0(0.005)_num:2
"""
def get_mean_depth_ee(gtf_ori,patient_name,set_num, ee_dict):
    gtf_dict = make_gtf_dict(gtf_ori)
    final_ee = open("./FINAL/"+set_num+"/"+patient_name+"_final_ee_after_j_after_making_sam.bed")
    final_ee_ori = final_ee.readlines()
    over34_ee = open("./OUTPUT/"+set_num+"/"+patient_name+"_over_threshold_34_ee_candi_dict_after_34.bed")
    over34_ee_ori = over34_ee.readlines()
    final_ee_dict = {}
    final_ee_dict_with_info = {}
    
    if len(final_ee_ori)>0:
        for i in final_ee_ori[:]:
            chromo = i.split("\t")[0]
            start = i.split(":")[1].split("-")[0]
            end = i.split(":")[1].split("//")[0].split("-")[1]
            a_d = i.split("\t")[3].split("_")[2]
            p_m = i.split("\t")[3].split("_")[0]
            novel_exon_point = ''
            if a_d == 'acceptor' and p_m == '-':
                novel_exon_point = i.split("\t")[1]
            elif a_d == 'donor' and p_m == '+':
                novel_exon_point = i.split("\t")[1]
            if a_d == 'acceptor' and p_m == '+':
                novel_exon_point = str(int(i.split("\t")[2])+1)
            elif a_d == 'donor' and p_m == '-':
                novel_exon_point = str(int(i.split("\t")[2])+1)
            input_key = chromo+"_"+start+"_"+end+"_"+novel_exon_point
            final_ee_dict[input_key] = i
    
    for line in over34_ee_ori[:]:
        for key in final_ee_dict.keys():
            chromo = key.split("_")[0]
            start = key.split("_")[1]
            end = key.split("_")[2]
            if chromo in line:
                if start in line:
                    exon_point = exon_point_calculate(line)
                    p_m = line.split("\t")[3].split("_")[0]
                    front_back = line.split("//")[1].split("_")[0]
                    gene = line.split("\t")[3].split("_")[6]
                    pep = line.split("\t")[3].split("_")[4]
                    input_key = "_"+str(exon_point)+"_"+p_m+"_"+front_back+"_"+gene+"_"+pep
                    final_ee_dict_with_info[key+input_key] = final_ee_dict.get(key)
    
    gene_to_transcript_dict, transcipt_to_info_dict = get_info_dict_for_mean_depth(final_ee_dict_with_info,gtf_dict)
    find_original_exon_point_ee(final_ee_dict_with_info, ee_dict, gene_to_transcript_dict,transcipt_to_info_dict,set_num, patient_name)


# In[189]:


"""------------------junction support read 끝------------------"""


# In[190]:


"""------------------about novel transcript 시작------------------"""


# In[191]:


def make_actg_pep_dict(flat_with_gff_dict):
    actg_pep_dict = {}
    for key in flat_with_gff_dict.keys():
        pep = key.split("\t")[10]
        if actg_pep_dict.get(pep) == None:
            actg_pep_dict[pep] = key
        else:
            actg_pep_dict[pep] = actg_pep_dict.get(pep)+key
    return actg_pep_dict


# In[192]:


def get_file_novel_transcript(set_num, patient):
    ri = open("./FINAL/"+set_num+"/"+patient+"_ri_for_Novel_Transcript_dict_after_making_sam.bed")
    ri_ori = ri.readlines()
    ee = open("./FINAL/"+set_num+"/"+patient+"_ee_for_compare_dict_after_making_sam.bed")
    ee_ori = ee.readlines()
    es = open("./FINAL/"+set_num+"/"+patient+"_es_for_compare_dict_after_making_sam.bed")
    es_ori = es.readlines()
    only_set_num = set_num.split("set")[1]
    flat = open("./ACTG/"+set_num+"/"+patient+"/output/VSG_"+only_set_num+"set_dn_pep_1014.flat")
    flat_ori = flat.readlines()
    gff = open("./ACTG/"+set_num+"/"+patient+"/output/VSG_"+only_set_num+"set_dn_pep_1014.gff")
    gff_ori = gff.readlines()
    flat_with_gff_dict = make_flat_with_gff_dict(flat_ori, gff_ori)
    actg_pep_dict = make_actg_pep_dict(flat_with_gff_dict)
    return ri_ori, ee_ori, es_ori, actg_pep_dict


# In[195]:


def translate_cds_map(p_m, cds_map,chromo):
    new_cds_map = []
    na_seq = ''
    pep = ''
    pep0 = ''
    pep1 = ''
    pep2 = ''
    CDSs = cds_map.split("//")
    for CDS in CDSs:
        cds_start = int(CDS.split(":")[0])
        cds_end = int(CDS.split(":")[1])
        if p_m == '+':
            new_cds_map.append(fasta_dict.get(chromo)[cds_start-1:cds_end])
        else:
            new_cds_map.append(fasta_dict.get(chromo)[cds_start-1:cds_end][::-1])
#             print(fasta_dict.get(chromo)[cds_start-1:cds_end][::-1])
    na_seq = "".join(new_cds_map)
    if p_m == '+':
        pep0 = translate(na_seq)
        pep1 = translate(na_seq[1:])
        pep2 = translate(na_seq[2:])
    elif p_m == '-':
        pep0 = reverse_translate(na_seq)
        pep1 = reverse_translate(na_seq[1:])
        pep2 = reverse_translate(na_seq[2:])
        
#     print(pep)
    return pep0, pep1, pep2


# In[196]:


def frame_check_nt(cds_map, start,p_m):
#     print(start)
    CDSs = cds_map.split("//")
    na_len = 0
    for CDS in CDSs:
        cds_start = int(CDS.split(":")[0])
        cds_end = int(CDS.split(":")[1])
        if p_m == '+':
            if cds_start<= start and start <= cds_end:
#                 print(na_len)
                na_len += start - cds_start
#                 print(na_len)
                if na_len % 3 ==0:
#                     print(CDS)
#                     print('frame true')
                    return True
                else:
#                     print(CDS)
#                     print('frame false')
                    return False
                pass
            else:
                na_len += cds_end-cds_start+1
        else:
            if cds_start<= start and start <= cds_end:
#                 print(na_len)
                na_len += cds_end - start
#                 print(na_len)
                if na_len % 3 ==0:
#                     print(CDS)
#                     print('frame true')
                    return True
                else:
#                     print(CDS)
#                     print('frame false')
                    return False
                pass
            else:
                na_len += cds_end-cds_start+1
    return False


# In[197]:


def get_info_dict_for_mean_depth_ee(gene,gtf_dict):
    gene_to_transcript_dict = {}
    transcipt_to_info_dict = {}
    
    if gtf_dict.get(gene) != None:
        info = gtf_dict.get(gene)
        lines = info.split("\n")
        for k in lines[:-1]:
            p_type = k.split("\t")[1]
            if p_type == 'protein_coding' and k.split("\t")[2] == "transcript":
                t_id = k.split("\"")[3]
                if transcipt_to_info_dict.get(t_id) is None:
                    transcipt_to_info_dict[t_id] = k
                if gene_to_transcript_dict.get(gene) is None:
                    gene_to_transcript_dict[gene] = t_id
                else:
                    gene_to_transcript_dict[gene] = gene_to_transcript_dict.get(gene)+"\n"+t_id
            elif p_type == 'protein_coding' and '; transcript_id ' in k :
                transcipt_to_info_dict[t_id] = transcipt_to_info_dict.get(t_id)+"\n"+k
    return transcipt_to_info_dict


# In[198]:


def ee_novel_t_process(ee_ori):
    novel_protein_ee_dict = {}
    found_dn = {}
#     gtf = open("./GTF/Homo_sapiens_only_chro.GRCh37.75.1006.gtf")
#     gtf_ori = gtf.readlines()
#     gtf_dict = make_gtf_dict(gtf_ori)
    
#     cds_map = deep_original_exon_point_ee(transcript_info, p_m, front_back, exon_point, candi_start, candi_end,novel_exon_point)
#     print(cds_map)
    
    t_dict = {}
    for line in ee_ori:
        
        md1 = -1
        md2 = -1
        fmd = -1
        md_line = ''
        
        chromo = line.split("\t")[0].split("chr")[1]
        start = int(line.split(":")[1].split("-")[0])+1
        end = int(line.split(":")[1].split("//")[0].split("-")[1])
        pep = line.split("_ENSG")[1].split("_")[1]
        found_dn[pep] = 1
        front_back = line.split("_ENSG")[0].split("_")[-1]
        exon_point = line.split("_ENSG")[0].split("_")[-3]
        novel_exon_point = line.split("_ENSG")[0].split("_")[-4]
        candi_start = line.split("_ENSG")[0].split("_")[-6]
        candi_end = line.split("_ENSG")[0].split("_")[-5]
#         print('pep')
#         print(pep)
#         print('pep')
        p_m = line.split("\t")[3].split("_")[0]
        gene = line.split("_")[-5]
        
        transcipt_to_info_dict = get_info_dict_for_mean_depth_ee(gene,gtf_dict)
        t_id = line.split("_")[-3]
        transcript_info = transcipt_to_info_dict.get(t_id)
        a,b,cds_map = deep_original_exon_point_ee(transcript_info, p_m, front_back, exon_point, candi_start, candi_end,novel_exon_point)
#         print('123123123')
#         print(cds_map)
#         print('123123123')
        
#         cds_map = line.split("_")[-2]
        md = float(line.split("_meanDepth:")[1].split("_")[0])
        if t_dict.get(t_id) == None:
            t_dict[t_id] = md
        else:
            md2 = t_dict.get(t_id)
            md1 = md
            if md2>md1:
                fmd = str(round(1, 4))
                md_line = str(md1)+"/"+str(md2)
            else:
                fmd = str(round(1, 4))
                md_line = str(md2)+"/"+str(md1)
        
        
        ori_jnc = int(line.split("ori_jnc:")[1].split("//")[0])
        ori_exact = int(line.split("ori_exact:")[1].split("//")[0])
        novel_jnc = int(line.split("novel_jnc:")[1].split("//")[0])
        novel_exact = int(line.split("novel_exact:")[1].split("(")[0])
        ratio = str(round((novel_jnc+novel_exact)/(ori_jnc+ori_exact),4))
        read_line = "ori_jnc:"+str(ori_jnc)+"+ori_exact:"+str(ori_exact)+"/novel_jnc:"+str(novel_jnc)+"+novel_exact:"+str(novel_exact)
        
        
        
        
        protein = ''
        

#         if frame_check_nt(cds_map, end, p_m):
        protein0, protein1, protein2 = translate_cds_map(p_m, cds_map,chromo)
        protein0 = protein0.replace("I","L")
        protein1 = protein1.replace("I","L")
        protein2 = protein2.replace("I","L")
        
        pep_list = check_pep_in_protein(protein0,protein1,protein2,pep)
        
        
        
#         protein = translate_cds_map(p_m, cds_map,chromo)
#         protein = protein.split("_")[0]
#             print('protein')
#             print(protein)
#             print('protein')
        
#         if len(pep_list) == 1:
#             novel_protein = pep_list[0]
#             novel_protein_ri_dict[novel_protein] = pep+"_"+line
#         elif len(pep_list) > 1:
#             print("oh no!!!")
#             print(line)
#             print(pep_list)
#             print("oh no!!!")
#             novel_protein_ri_dict[pep] = "##Many##_"+pep+"_"+line
#         elif len(pep_list) == 0:
#             novel_protein_ri_dict[pep] = "##No##_"+pep+"_"+line
        
        

        if len(pep_list) == 1:
            novel_protein = pep_list[0]
#             print("ee")
#             print(novel_protein)
#             print(pep)
#             print("ee")
#             if pep in protein:
#                 print(novel_protein)
            if fmd != -1:
#                 input_value = "pep_"+pep+"_nr/or_ratio:"+ratio+"("+read_line+")_meanDepthRatio:"+fmd+"("+md_line+")"+"_"+line
                input_value = "pep_"+pep+"_"+line
#                     print(input_value)
                novel_protein_ee_dict[novel_protein] = input_value
        elif len(pep_list) > 1:
            print("oh no!!!")
            print(line)
            print(pep_list)
            print("oh no!!!")
#             input_value = "nr/or_ratio:"+ratio+"("+read_line+")_meanDepthRatio:"+str(fmd)+"("+md_line+")"+"_"+line
            input_value = "pep_"+pep+"_"+line
            novel_protein_ee_dict[pep] = "##Many##_"+pep+"_"+input_value
        elif len(pep_list) == 0:
            pass
#             novel_protein_ee_dict[pep] = "##No##_"+pep+"_"+input_value
                
#                 novel_protein_ee_dict[novel_protein] = line
    return novel_protein_ee_dict,found_dn


# In[199]:


def meanDepth_es(a_start, a_end, a_md,b_start, b_end, b_md,c_start, c_end, c_md):
    a_len = a_end-a_start+1
    b_len = b_end-b_start+1
    c_len = c_end-c_start+1
    ab_md = (a_len*a_md+b_len*b_md)/(a_len+b_len)
    abc_md = (a_len*a_md+b_len*b_md+c_len*c_md)/(a_len+b_len+c_len)
    if c_md == 0:
        return str(round(1,4))+"(only_novel)", str(round(1,4))
    return str(round(1,4)), str(round(1,4))


# In[200]:


def es_novel_t_process(es_ori):
    found_dn = {}
    novel_protein_es_dict = {}
    for line in es_ori:
        chromo = line.split("\t")[0].split("chr")[1]
        p_m = line.split("\t")[3].split("_")[0]
        pep = line.split("\t")[3].split("_")[4]
        found_dn[pep] = 1
        t_id = line.split("_")[-4]
        cds_map = line.split("_")[-3]
        
        a_start = int(line.split("_")[-2].split("-")[0])
        a_end = int(line.split("_")[-2].split("-")[1])
        b_start = int(line.split("_")[-1].split("-")[0])
        b_end = int(line.split("_")[-1].split("-")[1])
        c_start = int(line.split("_skipped:")[1].split("_")[0].split(":")[0])
        c_end = int(line.split("_skipped:")[1].split("_")[0].split(":")[1])
        a_md = float(line.split("_meanDepth:")[1])
        b_md = float(line.split("_meanDepth:")[2])
        c_md = float(line.split("_meanDepth:")[3].split("_")[0])
        ab_c_md, ab_abc_md = meanDepth_es(a_start, a_end, a_md,b_start, b_end, b_md,c_start, c_end, c_md)
        
        normal_jnc = int(line.split('_jnc_read:')[1].split("(")[0])
        novel_jnc = int(line.split('_num:')[1].split("%%")[0])
        jnc_ratio = ''
        if normal_jnc == 0:
            jnc_ratio = '(only_novel)'
        else:
            jnc_ratio = str(round(novel_jnc/normal_jnc,4))
        
        if p_m == "+":
            start = int(line.split("\t")[1])+1
            end = int(line.split("\t")[2])    
        else:
            start = int(line.split("%%")[1].split("\t")[1])+1
            end = int(line.split("%%")[1].split("\t")[2])
        
        protein = ''
        
        protein0, protein1, protein2 = translate_cds_map(p_m, cds_map,chromo)
        protein0 = protein0.replace("I","L")
        protein1 = protein1.replace("I","L")
        protein2 = protein2.replace("I","L")

#         if frame_check_nt(cds_map, start, p_m):
#             protein = translate_cds_map(p_m, cds_map,chromo)
#             protein = protein.split("_")[0]
        pep_list = check_pep_in_protein(protein0,protein1,protein2,pep)    
    

        if len(pep_list) == 1:
            novel_protein = pep_list[0]
#             print("es")
#             print(novel_protein)
#             print(pep)
#             print("es")
                
#             input_value = "pep_"+pep+"_jnc_ratio:"+jnc_ratio+"("+str(novel_jnc)+"/"+str(normal_jnc)+")_ab/c:"+ab_c_md+"_ab/abc:"+ab_abc_md+"_"+line
            line = line.replace("meanDepth:0.0_meanDepth:0.0_meanDepth:0.0_","evidence1_evidence2_skipped_")
            line = "_".join(line.split("_")[:-2])+"\n"
            input_value = "pep_"+pep+"_"+line
            novel_protein_es_dict[novel_protein] = input_value
        elif len(pep_list) > 1:
            #print("oh no!!!")
            #print(line)
            #print(pep_list)
            #print("oh no!!!")
#             input_value = "jnc_ratio:"+jnc_ratio+"("+str(novel_jnc)+"/"+str(normal_jnc)+")_ab/c:"+ab_c_md+"_ab/abc:"+ab_abc_md+"_"+line
            line = line.replace("meanDepth:0.0_meanDepth:0.0_meanDepth:0.0_","evidence1_evidence2_skipped_")
            line = "_".join(line.split("_")[:-2])+"\n"
            input_value = "pep_"+pep+"_"+line
            novel_protein_es_dict[pep] = "##Many##_pep_"+pep+"_"+input_value
        elif len(pep_list) == 0:
            pass
#             input_value = "jnc_ratio:"+jnc_ratio+"("+str(novel_jnc)+"/"+str(normal_jnc)+")_ab/c:"+ab_c_md+"_ab/abc:"+ab_abc_md+"_"+line
#             novel_protein_es_dict[pep] = "##No##_pep_"+pep+"_"+input_value
    return novel_protein_es_dict,found_dn


# In[201]:


def check_pep_in_protein(protein0,protein1,protein2,pep):
    protein0_list = protein0.split("_")
    protein1_list = protein1.split("_")
    protein2_list = protein2.split("_")
    pep_list = []
    for i in protein0_list:
        if pep in i:
            pep_list.append(i)
    for i in protein1_list:
        if pep in i:
            pep_list.append(i)
    for i in protein2_list:
        if pep in i:
            pep_list.append(i)
    return pep_list


# In[202]:


def ri_novel_t_process(ri_ori):
    novel_protein_ri_dict = {}
    found_dn = {}
    for line in ri_ori:
        chromo = line.split("\t")[0].split("chr")[1]
        start = int(line.split("\t")[1])+1
        end = int(line.split("\t")[2])
        pep = line.split("\t")[3].split("_")[4]
        found_dn[pep] = 1
        p_m = line.split("\t")[3].split("_")[0]
        aa_seq = fasta_dict.get(chromo)[start-1:end]
        t_id = line.split("_")[-2]
        cds_map = line.split("_")[-1]
        protein1 = ''
        protein2 = ''
        protein0 = ''
        

#         if p_m == '+':
#             if frame_check_nt(cds_map, start,p_m):
        protein0, protein1, protein2 = translate_cds_map(p_m, cds_map,chromo)
        protein0 = protein0.replace("I","L")
        protein1 = protein1.replace("I","L")
        protein2 = protein2.replace("I","L")
        
        pep_list = check_pep_in_protein(protein0,protein1,protein2,pep)    

        if len(pep_list) == 1:
            novel_protein = pep_list[0]
            novel_protein_ri_dict[novel_protein] = pep+"_"+line
#             print("ri")
#             print(novel_protein)
#             print(pep)
#             print("ri")
        elif len(pep_list) > 1:
            print("oh no!!!")
            print(line)
            print(pep_list)
            print("oh no!!!")
            novel_protein_ri_dict[pep] = "##Many##_"+pep+"_"+line
        elif len(pep_list) == 0:
            pass
#             novel_protein_ri_dict[pep] = "##No##_"+pep+"_"+line
    return novel_protein_ri_dict, found_dn


# In[203]:


def get_no_actg_pep(actg_pep_dict, dn_pep_ori):
    no_dict = {}
    for line in dn_pep_ori:
        pep = line.replace("\n","")
        if actg_pep_dict.get(pep) == None:
            no_dict[pep] = 1
    return no_dict


# In[204]:


def check_other_dn_pep(novel_dict, found,no_actg_pep_dict,actg_pep_dict):
    novel_pep_dict = {}
    for key in novel_dict.keys():
        for key_1 in no_actg_pep_dict.keys():
            if key_1 in key:
                if found.get(key_1) == None:
                    #print('i found it no_actg')
                    novel_pep_dict[key] = '_novel_pep:'+key+"_dn_no_actg:"+key_1
#                     print(key_1)
#                     print(no_actg_pep_dict.get(key_1))
#                     print(key)
#                     print(novel_dict.get(key))
#                     print('i found it no_actg')
        for key_1 in actg_pep_dict.keys():
            if key_1 in key:
                if found.get(key_1) == None:
                    novel_pep_dict[key] = '_novel_pep:'+key+"_dn_yes_actg:"+key_1
#                     print('i found it yes_actg')
#                     print(key_1)
#                     print(actg_pep_dict.get(key_1))
#                     print(key)
#                     print(novel_dict.get(key))
#                     print('i found it yes_actg')
    for key in novel_dict.keys():
        if novel_pep_dict.get(key) != None:
            novel_dict[key] = novel_dict.get(key)+novel_pep_dict.get(key)


# In[205]:


def make_depth_file(write_file, novel_dict,patient, set_num):
    for i in novel_dict.keys():
        write_file.write(patient+"_"+set_num+"_"+novel_dict.get(i))


# In[206]:


"""


"""

def start_novel_transcript_evidence_result():
    gtf = open("./GTF/Homo_sapiens_only_chro.GRCh37.75.1006.gtf")
    gtf_ori = gtf.readlines()
    gtf_dict = make_gtf_dict(gtf_ori)
    
    file_list = os.listdir('./ACTG/')
    path = './ACTG/'
    
    ri_depth = open("./FINAL/ri_total.txt","wt")
    ee_depth = open("./FINAL/ee_total.txt","wt")
    es_depth = open("./FINAL/es_total.txt","wt")
    
    for set_num in file_list:
        if 'set' in set_num:
            print(set_num)
            set_path = path+set_num+"/"
            patient_list = os.listdir(set_path)
            only_set_num = set_num.split("set")[1]

            dn_pep = open("./ACTG/"+set_num+"/"+only_set_num+"set_dn_pep_1014.txt")
            dn_pep_ori = dn_pep.readlines()

            for patient in patient_list:
                if len(patient.split(".")) == 1:
                    print(patient)
                    start_time = time.time()
                    ri_ori, ee_ori, es_ori, actg_pep_dict = get_file_novel_transcript(set_num, patient)
                    no_actg_pep_dict = get_no_actg_pep(actg_pep_dict, dn_pep_ori)

                    novel_ri_dict, found_ri = ri_novel_t_process(ri_ori)
                    novel_es_dict,found_es = es_novel_t_process(es_ori)
                    novel_ee_dict,found_ee = ee_novel_t_process(ee_ori)
                    
                    check_other_dn_pep(novel_ri_dict, found_ri,no_actg_pep_dict,actg_pep_dict)
                    check_other_dn_pep(novel_ee_dict, found_ee,no_actg_pep_dict,actg_pep_dict)
                    check_other_dn_pep(novel_es_dict, found_es,no_actg_pep_dict,actg_pep_dict)
                    
                    make_depth_file(ri_depth, novel_ri_dict,patient, set_num)
                    make_depth_file(ee_depth, novel_ee_dict,patient, set_num)
                    make_depth_file(es_depth, novel_es_dict,patient, set_num)
                    
                    
                    print("one patients find_other_dn mapping time :", time.time() - start_time)
                throw_garbage()
    ri_depth.close()
    ee_depth.close()
    es_depth.close()
    
    

def start_novel_transcript_evidence_result_one(patient_name):
    gtf = open("./GTF/Homo_sapiens_only_chro.GRCh37.75.1006.gtf")
    gtf_ori = gtf.readlines()
    gtf_dict = make_gtf_dict(gtf_ori)
    
    file_list = os.listdir('./ACTG/')
    path = './ACTG/'
    
    ri_depth = open("./FINAL/ri_total.txt","wt")
    ee_depth = open("./FINAL/ee_total.txt","wt")
    es_depth = open("./FINAL/es_total.txt","wt")
    
    for set_num in file_list:
        if 'set' in set_num:
            set_path = path+set_num+"/"
            patient_list = os.listdir(set_path)
            only_set_num = set_num.split("set")[1]

            dn_pep = open("./ACTG/"+set_num+"/"+only_set_num+"set_dn_pep_1014.txt")
            dn_pep_ori = dn_pep.readlines()

            for patient in patient_list:
                if patient_name == patient:
                    if len(patient.split(".")) == 1:
                        start_time = time.time()
                        ri_ori, ee_ori, es_ori, actg_pep_dict = get_file_novel_transcript(set_num, patient)
                        no_actg_pep_dict = get_no_actg_pep(actg_pep_dict, dn_pep_ori)

                        novel_ri_dict, found_ri = ri_novel_t_process(ri_ori)
                        novel_es_dict,found_es = es_novel_t_process(es_ori)
                        novel_ee_dict,found_ee = ee_novel_t_process(ee_ori)
                        
                        check_other_dn_pep(novel_ri_dict, found_ri,no_actg_pep_dict,actg_pep_dict)
                        check_other_dn_pep(novel_ee_dict, found_ee,no_actg_pep_dict,actg_pep_dict)
                        check_other_dn_pep(novel_es_dict, found_es,no_actg_pep_dict,actg_pep_dict)
                        
                        make_depth_file(ri_depth, novel_ri_dict,patient, set_num)
                        make_depth_file(ee_depth, novel_ee_dict,patient, set_num)
                        make_depth_file(es_depth, novel_es_dict,patient, set_num)
                        
                        
                        print("one patients find_other_dn mapping time :", time.time() - start_time)
                    throw_garbage()
    ri_depth.close()
    ee_depth.close()
    es_depth.close()


# In[207]:


def reverse_translate(seq):
#     reverse_seq = seq[::-1]
    return_seq = []
    for i in seq:
        if i == 'A':
            return_seq.append('T')
        elif i == 'T':
            return_seq.append('A')
        elif i == 'G':
            return_seq.append('C')
        elif i == 'C':
            return_seq.append('G')
    reverse_seq = "".join(return_seq)
    protein = translate(reverse_seq)
    return protein


# In[208]:


def translate(seq):
      
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein =""
#     if len(seq)%3 == 0:
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        if len(codon) == 3:
            protein+= table[codon]
    return protein


# In[209]:


"""------------------about novel transcript 끝------------------"""


# In[210]:


"""------------------기타 시작------------------"""


# In[211]:


def calculate_p_value_of_intron(intron_mean_depth_list, value_from_PGpep):
    intron_mean_depth_list.append(value_from_PGpep)
    intron_mean_depth_list.sort()
    p_value = 1-(intron_mean_depth_list.index(value_from_PGpep)+1)/len(intron_mean_depth_list)
    intron_mean_depth_list.remove(value_from_PGpep)
    p_value = round(p_value,4)
    p_value = int(p_value*1000)/1000
    return p_value


# In[212]:


def get_intron_mean_depth_info():
    intron_depth = open("./ETC/intron_mean_depth_between_two_cds.txt")
    intron_depth_ori = intron_depth.readlines()

    intron_mean_depth_list = []


    for i in intron_depth_ori[:]:
        read = float(i.split("\t")[4])
        start = float(i.split("\t")[1])
        end = float(i.split("\t")[2])
        mean_depth = read/(end-start)
        intron_mean_depth_list.append(mean_depth)

    intron_mean_depth_list.sort()
    return intron_mean_depth_list


# In[213]:


"""
description:

arguments:

output : 
"""
def throw_garbage():
    gc.collect()


# In[214]:


"""
description:
    변수 명을 출력하기 위한 함수

arguments:
    var : 변수

output : 
    입력 변수의 이름
"""
def retrieve_name(var):
        """
        Gets the name of var. Does it from the out most frame inner-wards.
        :param var: variable to get name from.
        :return: string
        """
        for fi in reversed(inspect.stack()):
            names = [var_name for var_name, var_val in fi.frame.f_locals.items() if var_val is var]
            if len(names) > 0:
                return names[0]


# In[215]:


"""
description:

arguments:

output : 
"""
def get_ee_whole():
    file_list = os.listdir('./ACTG/')
    path = './ACTG/'
    for set_num in file_list:
        if 'set' in set_num:
            print(set_num)
            set_path = path+set_num+"/"
            patient_list = os.listdir(set_path)
            for patient in patient_list:
                if len(patient.split(".")) == 1:
                    print(patient)
                    only_ee(patient, set_num)


# In[216]:


"""
description:

arguments:

output : 
"""
def get_ee_es_whole():
    file_list = os.listdir('./ACTG/')
    path = './ACTG/'
    for set_num in file_list:
        if 'set' in set_num:
            print(set_num)
            set_path = path+set_num+"/"
            patient_list = os.listdir(set_path)
            for patient in patient_list:
                if len(patient.split(".")) == 1:
                    print(patient)
                    only_ee_es(patient, set_num)


# In[217]:


"""
description:

arguments:
    patient
    set_num
    
output : 
"""
def only_ee(patient, set_num):
    path = './OUTPUT/'+set_num+"/"
    count = 0
    file_list = os.listdir(path)
    whole_candi = open(path+patient+"_whole_candi_after_samtools_1016.bed")
    only_write = open(path+patient+"_only_ee_after_samtools_1016.bed","wt")
    whole_candi_ori = whole_candi.readlines()
    for line in whole_candi_ori:
        if 'exon-extension' in line.split("_")[-6]:
            only_write.write(line)
            count+=1
    only_write.write("exon-extension #num : "+str(count)+"\n")
    print(count)
    only_write.close()


# In[218]:


"""
description:

arguments:
    patient
    set_num

output : 
"""
def only_ee_es(patient, set_num):
    path = './OUTPUT/'+set_num+"/"
    file_list = os.listdir(path)
    whole_candi = open(path+patient+"_whole_candi_after_samtools_1016.bed")
    only_write = open(path+patient+"_only_ee_es_after_samtools_1016.bed","wt")
    whole_candi_ori = whole_candi.readlines()
    for line in whole_candi_ori:
        if 'exon-extension' in line.split("_")[-6]:
            only_write.write(line)
        elif 'exon-skipping' in line.split("_")[-6]:
            only_write.write(line)
    only_write.close()


# In[219]:


"""
description:

arguments:
    gtf_ori
    gff_ori
    flat_ori
    patient_name
    set_num
output : 
"""
def temp_main(gtf_ori, gff_ori, flat_ori, patient_name, set_num):
    gtf_dict = make_gtf_dict(gtf_ori)
    tsv_dict = make_flat_with_gff_dict(flat_ori, gff_ori)
    tsv_mapped_dict = make_mapped_dict(tsv_dict, gtf_dict)
    full_intron_exon_dict = make_full_intron_exon_dict(tsv_mapped_dict)
    threshold_overlap_over_4_dict,full_intron_dict = make_threshold_overlap_over_4_dict_and_full_intron_dict(full_intron_exon_dict,4)

    bed_dict = make_bed_and_dict_with_now_whole(threshold_overlap_over_4_dict, full_intron_dict, patient_name, set_num)
    
    #do samtools
    samtools_depth_output(patient_name, set_num)

    #samtools 결과 적용
    rna_mapping_dict = make_rna_count_dict(patient_name, set_num)

    #rna 전부 0이면 날려
    candi_with_rna_dict, candi_without_rna_dict = make_candi_with_rna_dict_bed_dict(rna_mapping_dict, bed_dict)

    #intron 부분 rna 전부 0이면 날려
    filter_intron_rna_none_dict = make_filter_intron_rna_none_dict(candi_with_rna_dict)
    
    rna_mapped_file = open("./OUTPUT/"+set_num+"/"+patient_name+"_whole_candi_after_samtools_RNA_mapped.bed","wt")
    for i in filter_intron_rna_none_dict.keys():
        rna_mapped_file.write(i+"\n")
    rna_mapped_file.close()
    

    gc.collect()
    
    
def temp_main_one(bam_path, gtf_ori, gff_ori, flat_ori, patient_name, set_num):
    gtf_dict = make_gtf_dict(gtf_ori)
    tsv_dict = make_flat_with_gff_dict(flat_ori, gff_ori)
    tsv_mapped_dict = make_mapped_dict(tsv_dict, gtf_dict)
    full_intron_exon_dict = make_full_intron_exon_dict(tsv_mapped_dict)
    threshold_overlap_over_4_dict,full_intron_dict = make_threshold_overlap_over_4_dict_and_full_intron_dict(full_intron_exon_dict,4)

    bed_dict = make_bed_and_dict_with_now_whole(threshold_overlap_over_4_dict, full_intron_dict, patient_name, set_num)
    
    #do samtools
    samtools_depth_output_one(bam_path, patient_name, set_num)

    #samtools 결과 적용
    rna_mapping_dict = make_rna_count_dict(patient_name, set_num)

    #rna 전부 0이면 날려
    candi_with_rna_dict, candi_without_rna_dict = make_candi_with_rna_dict_bed_dict(rna_mapping_dict, bed_dict)

    #intron 부분 rna 전부 0이면 날려
    filter_intron_rna_none_dict = make_filter_intron_rna_none_dict(candi_with_rna_dict)
    
    rna_mapped_file = open("./OUTPUT/"+set_num+"/"+patient_name+"_whole_candi_after_samtools_RNA_mapped.bed","wt")
    for i in filter_intron_rna_none_dict.keys():
        rna_mapped_file.write(i+"\n")
    rna_mapped_file.close()
    

    gc.collect()


# In[220]:


"""
description:

arguments:
    gtf_ori
    gff_ori
    flat_ori
    patient_name
    set_num
output : 
"""
def temp_main_intron(gtf_ori, gff_ori, flat_ori, patient_name, set_num):
    gtf_dict = make_gtf_dict_intron(gtf_ori)
    tsv_dict = make_flat_with_gff_dict(flat_ori, gff_ori)
    tsv_mapped_dict = make_mapped_dict_intron(tsv_dict, gtf_dict)
    full_intron_exon_dict = make_full_intron_exon_dict(tsv_mapped_dict)
    threshold_overlap_over_4_dict,full_intron_dict = make_threshold_overlap_over_4_dict_and_full_intron_dict(full_intron_exon_dict,4)
    print(patient_name)
    print(set_num)
    print('full_intron_dict')
    print(len(full_intron_dict))
    print('threshold_overlap_over_4_dict')
    print(len(threshold_overlap_over_4_dict))
    bed_dict = make_bed_and_dict_with_now_intron(full_intron_dict, patient_name, set_num)
    
    
    start_time = time.time()
    
    #do samtools
    samtools_depth_output_intron(patient_name, set_num)
    
    print("one patients samtools mapping time :", time.time() - start_time)

    #samtools 결과 적용
    rna_mapping_dict = make_rna_count_dict_intron(patient_name, set_num)
    

    #rna 전부 0이면 날려
    candi_with_rna_dict,candi_without_rna_dict = make_candi_with_rna_dict_bed_dict(rna_mapping_dict, bed_dict)

    #intron 부분 rna 전부 0이면 날려
    
    ee_file = open("./OUTPUT/"+set_num+"/"+patient_name+"_whole_candi_intron_after_samtools.bed","wt")
    for i in candi_with_rna_dict.keys():
        ee_file.write(i+"\n")
    ee_file.close()
    
    ee_file2 = open("./OUTPUT/"+set_num+"/"+patient_name+"_whole_candi_intron_after_samtools_no_rna.bed","wt")
    for i in candi_without_rna_dict.keys():
        ee_file2.write(i+"\n")
    ee_file2.close()
    


# In[221]:


"""------------------기타 끝------------------"""


# In[222]:


#------------------------ 범죄 수용소 ---------------------------------


# In[223]:


def temp_main_after_rna_mapped_intron_frame_out(gtf_ori,patient_name, set_num):
    gtf_dict = make_gtf_dict(gtf_ori)
    filter_intron_rna_none_dict = {}
    filter_intron_rna = open("./OUTPUT/"+set_num+"/"+patient_name+"_whole_candi_intron_after_samtools.bed")
    filter_intron_rna_ori = filter_intron_rna.readlines()
    for i in filter_intron_rna_ori[:]:
        filter_intron_rna_none_dict[i.replace("\n","")] = 1
    ee_dict = make_ee_dict(filter_intron_rna_none_dict)
    gene_to_transcript_dict, transcipt_to_info_dict = get_info_dict(ee_dict,gtf_dict)
    ee_process_intron_frame_out(ee_dict, gene_to_transcript_dict,transcipt_to_info_dict,set_num, patient_name)


# In[224]:


def temp_main_after_rna_mapped_intron(gtf_ori,patient_name, set_num):
    gtf_dict = make_gtf_dict(gtf_ori)
    filter_intron_rna_none_dict = {}
    filter_intron_rna = open("./OUTPUT/"+set_num+"/"+patient_name+"_whole_candi_intron_after_samtools.bed")
    filter_intron_rna_ori = filter_intron_rna.readlines()
    for i in filter_intron_rna_ori[:]:
        filter_intron_rna_none_dict[i.replace("\n","")] = 1
    ee_dict = make_ee_dict(filter_intron_rna_none_dict)
    gene_to_transcript_dict, transcipt_to_info_dict = get_info_dict(ee_dict,gtf_dict)
    ee_process_intron(ee_dict, gene_to_transcript_dict,transcipt_to_info_dict,set_num, patient_name)


# In[225]:


def ee_process_intron_frame_out(ee_dict, gene_to_transcript_dict,transcipt_to_info_dict,set_num, patient_name):
    over_threshold_34_full_intron_ri_frame_candi_dict = {}
    for line in ee_dict.keys():
        get_cds_all_transcript_region_intron_ri_frame_out(line, gene_to_transcript_dict, transcipt_to_info_dict,over_threshold_34_full_intron_ri_frame_candi_dict)
    final_ri_candi = over_threshold_34_full_intron_ri_frame_candi_dict

    make_bed_dict_with_now_dict(final_ri_candi,set_num,patient_name)


# In[226]:


def ee_process_intron(ee_dict, gene_to_transcript_dict,transcipt_to_info_dict,set_num, patient_name):
    over_threshold_34_full_intron_ri_candi_dict = {}
    for line in ee_dict.keys():
        over_threshold_34_ri_candi_dict = get_cds_all_transcript_region_intron_ri(line, gene_to_transcript_dict, transcipt_to_info_dict,over_threshold_34_full_intron_ri_candi_dict)
    final_ri_candi = over_threshold_34_full_intron_ri_candi_dict
    for i in final_ri_candi.keys():
        print(i)
        print(final_ri_candi.get(i))
    print(len(final_ri_candi))
    make_bed_dict_with_now_dict(final_ri_candi,set_num,patient_name)


# In[227]:


def check_region_intron(candi_start, candi_end, cds_first_point, cds_second_point):

    if cds_first_point < candi_start and candi_start < cds_second_point:
        if cds_first_point < candi_end and candi_end < cds_second_point:

            return True
    elif cds_second_point < candi_start and candi_start < cds_first_point:
        if cds_second_point < candi_end and candi_end < cds_first_point:

            return True
    return False


# In[228]:


def get_cds_one_transcript_region_intron(transcript_info,candi_start,candi_end):
    CDSs = transcript_info.split("\tCDS\t")
    cds_map = []
    cds_first_point = 0
    cds_second_point = 0
    
    for CDS in CDSs[1:]:
        cds_start = int(CDS.split("\t")[0])
        cds_end = int(CDS.split("\t")[1])
        if cds_second_point == 0:
            cds_second_point = int(CDS.split("\t")[0])
            key = str(cds_start)+":"+str(cds_end)
            cds_map.append(key)
        else:
            cds_first_point = cds_second_point
            cds_second_point = int(CDS.split("\t")[0])
            check = check_region_intron(candi_start, candi_end, cds_first_point, cds_second_point) ##123123121242134123451234512345
            if check == False:
                key = str(cds_start)+":"+str(cds_end)
                cds_map.append(key)
            elif check == True:
                key = str(cds_start)+":"+str(cds_end)
                cds_map.append(key)
                if len(cds_map) >1:
#                     print(cds_map)
                    return "//".join(cds_map)
    return None


# In[229]:


def get_cds_all_transcript_region_intron_ri_frame_out(line_ee_dict, gene_to_transcript_dict, transcipt_to_info_dict,over_threshold_34_ri_candi_dict):
    chromo = line_ee_dict.split("\t")[0].split("chr")[1]
    gene = line_ee_dict.split("_")[-7]
    p_m = line_ee_dict.split("\t")[3].split("_")[0]
    candi_start = int(line_ee_dict.split("\t")[1])+1
    candi_end = int(line_ee_dict.split("\t")[2])
    ri_info = ''
    if gene_to_transcript_dict.get(gene) == None:
        pass
    else:
        transcripts = gene_to_transcript_dict.get(gene).split("\n")
        retained_intron = False
        for transcript in transcripts:
            transcript_info = transcipt_to_info_dict.get(transcript)
            #ex) cds_region = '123:555//666:777'
            cds_region = get_cds_one_transcript_region_intron(transcript_info,candi_start,candi_end)
            if cds_region != None:
                stop_codon_exist = False
                frame = 2
                
                intron_seq,intron_position = get_intron(cds_region,chromo,p_m,2)

                stop_codon_exist, stop_codon_position = get_close_stop_codon_position_for_front_ri_frame_out(intron_seq, p_m,intron_position, candi_start, candi_end)
                if stop_codon_exist:
                    pass
                else: #stop 코돈 없어!

                    retained_intron = True
                    ri_info = "_"+p_m+"_"+str(frame)+"_"+intron_position+"\n"

            else: #cds가 없어!
                pass
        if retained_intron:
            over_threshold_34_ri_candi_dict[line_ee_dict.replace("\n","")+ri_info] =1
            
    return over_threshold_34_ri_candi_dict


# In[230]:


def get_cds_all_transcript_region_intron_ri(line_ee_dict, gene_to_transcript_dict, transcipt_to_info_dict,over_threshold_34_ri_candi_dict):
    chromo = line_ee_dict.split("\t")[0].split("chr")[1]
    gene = line_ee_dict.split("_")[-7]
    p_m = line_ee_dict.split("\t")[3].split("_")[0]
    candi_start = int(line_ee_dict.split("\t")[1])+1
    candi_end = int(line_ee_dict.split("\t")[2])
    ri_info = ''
    if gene_to_transcript_dict.get(gene) == None:
        pass
    else:
        transcripts = gene_to_transcript_dict.get(gene).split("\n")
        retained_intron = False
        for transcript in transcripts:
            transcript_info = transcipt_to_info_dict.get(transcript)
            #ex) cds_region = '123:555//666:777'
            cds_region = get_cds_one_transcript_region_intron(transcript_info,candi_start,candi_end)
            if cds_region != None:
                stop_codon_exist = False
                frame = get_frame_up_stream(cds_region)
                
                intron_seq,intron_position = get_intron(cds_region,chromo,p_m,frame)

                
                check_frame = frame_check_back_novel_event(cds_region,frame, p_m,candi_start,candi_end)
                    
                if check_frame: #frame 확인

                    stop_codon_exist, stop_codon_position = get_close_stop_codon_position_for_front(intron_seq, p_m,intron_position)
                    if stop_codon_exist:
                        pass
                    else: #stop 코돈 없어!

                        retained_intron = True
                        ri_info = "_"+p_m+"_"+str(frame)+"_"+intron_position+"\n"
                else: #frame 맞지 않아
                    pass 
            else: #cds가 없어!
                pass
        if retained_intron:
            over_threshold_34_ri_candi_dict[line_ee_dict.replace("\n","")+ri_info] =1
            
    return over_threshold_34_ri_candi_dict


# In[231]:


"""
description:

arguments:

output : 
"""
def after_rna_whole_4():
    file_list = os.listdir('./ACTG/')
    path = './ACTG/'
    for set_num in file_list:
        if 'set4' in set_num:
            print(set_num)
            set_path = path+set_num+"/"
            patient_list = os.listdir(set_path)
            for patient in patient_list:
                if len(patient.split(".")) == 1:
                    print(patient)
                    after_rna_whole_one(patient, set_num)


# In[232]:


def refresh_file_list():
    file_list = os.listdir('./')
    return file_list


# In[233]:


"""
description:

arguments:

output : 
"""
def make_bed_dict_with_after_bam(patient_num, set_num):
    bed_dict = {}
    bed = open('./BED/'+set_num+"/"+patient_num+"_novel_evnet_candi_before_bam.bed","r")
    bed_ori = bed.readlines()
    for line in bed_ori[:]:
        line1 = line.replace("\n","")
        info = line.split("\t")
        chromo = info[0]
        start = info[1]
        end = info[2]
        bed_dict[line1] = 1
    return bed_dict


# In[234]:


"""
description:

arguments:

output : 
"""
def make_bai():
    file_list = os.listdir('./BAM/')
    path = './BAM/'
    go = ['8','9','10','11']
    for set_num in file_list:
        for number in go:
            now = 'set'+number
            if now in set_num:
                print(set_num)
                set_path = path+set_num+"/"
                patient_list = os.listdir(set_path)
                for patient in patient_list:
                    print(patient)
                    now_path = set_path+patient
                    if patient.split(".")[-1] == 'bam':
                        bai_input = 'samtools index '+now_path+' '+set_path+patient+".bai"
                        print(bai_input)
                        r = subprocess.Popen(bai_input, shell=True).wait()
                        if r == 1: 
                            print("making bai failed")


# In[237]:


"""
description:

arguments:

output : 
"""
def before_rna_4():
    file_list = os.listdir('./ACTG/')
    path = './ACTG/'
    
    #fasta 파일 다 읽어들여와!
    #fasta_dict = make_fasta_dict()
    
    for set_num in file_list:
        if 'set4' in set_num:
            print(set_num)
            set_path = path+set_num+"/"
            patient_list = os.listdir(set_path)
            for patient in patient_list:
                print(patient)
                if len(patient.split(".")) == 1:
                    folder_path = set_path+patient+"/"
                    gtf = open(folder_path+"GTF/"+patient+"_grch37.gtf")
                    gtf_ori = gtf.readlines()
                    output_list = os.listdir(folder_path+"output/")
                    flat_ori = ''
                    gff_ori = ''
                    for output_file in output_list:
                        if output_file.split(".")[1] == 'flat':
                            flat = open(folder_path+"output/"+output_file)
                            flat_ori = flat.readlines()
                        elif output_file.split(".")[1] == 'gff':
                            gff = open(folder_path+"output/"+output_file)
                            gff_ori = gff.readlines()
                    temp_main_4(gtf_ori, gff_ori, flat_ori, patient, set_num)


# In[ ]:

mode = str(input_param.split('mode ')[1].split(' ')[0])

if mode == 'w':
    main_whole(input_param)
else:
    main_one(input_param)
    

