PGpep_README

1. Requirement

python version 3.6 이상
Samtools version 1.1 이상
Memory 8G 이상

Linux 환경 추천 (root 계정 이용)
Window 환경 이용시 powershell 을 통한 명령어 입력가능

2. Usage

명령어
a. 모드 whole(환자 전체에 대한 계산, 폴더 구조와 모든 input 파일은 3번 문서를 따라야 함)
python PGpep.py mode w t_id_c 0 g_id_c 3 fpkm_c 9

b. 모드 one(환자 1명에 대한 계산, 각 input 파일의 모든 path를 요구함, 단 결과를 계산하는데 있어 파일이 재사용되어 폴더 구조 및 폴더 이름은 3번 문서를 따라야 함, input 파일은 해당 위치에 없어도 가능함)

python PGpep.py mode o p_name patient_name t_id_c 0 g_id_c 3 fpkm_c 9 fpkm_p fpkm_path bam_p bam_path f_s_db_p 1st_2nd_db_path c_db_p composite_db_path u_db_p uniprot_db_path 
ex) python PGpep.py mode o p_name IRCR_GBM15_774_T01 t_id_c 0 g_id_c 3 fpkm_c 9 fpkm_p /home/bis/KSG/test_0113/FPKM/set1/IRCR_GBM15_774_T01_RSq.isoforms.fpkm_tracking/IRCR_GBM15_774_T01_RSq.isoforms.fpkm_tracking bam_p /home/bis/KSG/test_0113/BAM/set1/IRCR_GBM15_774_T01_RSq_splice.sorted.bam f_s_db_p /home/bis/KSG/test_0113/DB/1st_2nd_DB_search/GBM_1set_rmPosition/TMT_GBM_1set_all_774.id_reconstructed_notFound_Position_Remover.txt c_db_p /home/bis/KSG/test_0113/DB/CompositeDB/1set/CompositeDB_SwissProt_774.fasta u_db_p /home/bis/KSG/test_0113/DB/uniprot/uniprot-filtered-organism__Homo+sapiens+Human+[9606]_.fasta


p_name(환자명)
t_id_c(transcript id column num)
g_id_c(gene id column num)
fpkm_c(fpkm column num)
fpkm_p(fpkm 파일의 이름이 포함된 파일 경로)
bam_p(bam 파일의 이름이 포함된 파일 경로, 해당 경로에 bai 파일도 함께 있어야 함)
f_s_db_p(1st_2nd_DB_search 파일의 이름이 포함된 파일 경로)
c_db_p(composite_db_search 파일의 이름이 포함된 파일 경로)
u_db_p(uniprot_db_search 파일의 이름이 포함된 파일 경로)

3. File 구조(sample의 형식을 참고)

set_name에는 'set'이라는 문자가 포함되어야 함

ACTG
ACTG —> set_name —> patient_name —> GTF, output


BAM
BAM —> set_name —> 파일명에 patient_name이 들어간 bam 파일과 bai 파일

BED
BED —> set_name —> file

DB
DB —> 1st_2nd_DB_search —> set_name —> file
	CompositeDB —> set_name —> file
	uniprot —> file

DENOVO
DENOVO —> peaks result 파일

ETC
ETC —> file

FASTA
FASTA —> file
file명은 염색체 번호여야 함.
Ex) 1.fa, 2.fa, 3.fa, … , 22.fa, MT.fa, X.fa, Y.fa

FINAL
FINAL —> set_name —> file
	

FPKM
FPKM —> set_name —> folder_name —> file_name
folder_name 과 file_name이 같아야 함.
file_name의 형식은 아래와 같음
patient_name+”_RSq.isoforms.fpkm_tracking”
Ex) patient1_RSq.isoforms.fpkm_tracking

GTF
GTF —> file

OUTPUT
OUTPUT —> set_name —> file