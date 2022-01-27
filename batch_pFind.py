import os
import sys
import concurrent.futures as cf

raw_folder = sys.argv[1]
fasta = sys.argv[2]

if len(sys.argv) > 3:
    enzyme_aa = sys.argv[3]
    enzyme_term = sys.argv[4]
    enzyme_specific = sys.argv[5] #0: non-specific, 1 or 2: semi-specific, 3: specific
else:
    enzyme_aa = 'KR'
    enzyme_term = 'C'
    enzyme_specific = 3

pfind_path=r'C:\pFindStudio\pFind3\bin'


used_mod = os.path.abspath('./reduced_mod.ini')
#used_mod = os.path.join(pfind_path,'modification.ini')

n_process = 6


pfind_main='Searcher.exe'
pparse_main='pParse.exe'
enzyme = f'Enzyme {enzyme_aa} _ {enzyme_term}'
digest = enzyme_specific
contaminant=os.path.join(pfind_path, 'contaminant.fasta')
add_contaminant=False


output_folder = os.path.join(raw_folder, 'pFind3')
if not os.path.isdir(output_folder):
    os.makedirs(output_folder)
    
if add_contaminant:
    out_fasta = fasta[:fasta.rfind('.')]+'_con.fasta'
    with open(out_fasta,'w') as outfile:
        with open(fasta) as infile:
            outfile.write(infile.read())
        with open(contaminant) as infile:
            outfile.write(infile.read())
    fasta = out_fasta
    
    
os.chdir(pfind_path)

cfg_template = f'''
# This is a standard pFind configure file
# For help: mail to chihao@ict.ac.cn
# Time: 3/10/2021 11:13:38 PM

[Version]
pFind_Version=EVA.3.0.11

[param]
thread=2
activation_type=HCD-FTMS
mstol=%s
mstolppm=1
msmstol=%s
msmstolppm=1
temppepnum=100
pepnum=10
selectpeak=200
maxprolen=60000000
maxspec=2000000
IeqL=1
npep=2
maxdelta=500
selectmod=
fixmod=
maxmod=3
enzyme={enzyme}
digest={digest}
max_clv_sites=3

[filter]
psm_fdr=0.01
psm_fdr_type=1
mass_lower=600
mass_upper=10000
len_lower=6
len_upper=100
pep_per_pro=1
pro_fdr=0.01

[engine]

### Blind mass search ###
#open=3
#open_tag_len=5
#rest_tag_iteration=0
#rest_tag_len=5
#rest_mod_num=-1
#digest_in_open=3

open=1
open_tag_len=5

rest_tag_iteration=1
rest_tag_len=4
rest_mod_num=10

salvo_iteration=1
salvo_mod_num=5

[file]
modpath={used_mod}
fastapath={fasta}
outputpath=%s
outputname=%s

[datalist]
msmsnum=1
msmspath1=%s
msmstype=%s

[quant]
quant=1|None

[system]
log=LOG_INFO
'''

def run_pFind(pfind_cfg, raw_file):
    if not os.path.isfile(os.path.splitext(raw_file)[0]+'_HCDFT.pf2') and not os.path.isfile(os.path.splitext(raw_file)[0]+'_HCDFT.mgf'):
        if raw_file.lower().endswith('.raw'):
            os.system(f'{pparse_main} -D {raw_file}\n')
        elif raw_file.lower().endswith(".d"):
            from ap_to_pfind import ap_bruker_to_mgf
            from ap_import_bruker import import_bruker_raw
            import_bruker_raw(raw_file)
            ap_bruker_to_mgf(raw_file[:-2]+'.ms_data.hdf')
            
    os.system(f"{pfind_main} {pfind_cfg}\n")
    return raw_file
    
def parallel_run(params_list):
    futures = []
    with cf.ProcessPoolExecutor(n_process) as pp:
        for pfind_cfg, raw_file in params_list:
            futures.append(pp.submit(run_pFind, pfind_cfg, raw_file))
        for future in cf.as_completed(futures):
            print(f"\n[DONE] {os.path.split(future.result())[1]}\n")

folder_list = []
params_list = []

def generate_cfg_one_raw(raw_name, raw_file, ppm=20, fmt='pf2'):
    one_raw_out = os.path.join(output_folder, raw_name)
    folder_list.append(one_raw_out)
    if not os.path.isdir(one_raw_out):
        os.makedirs(one_raw_out)
    with open(os.path.join(one_raw_out, 'pFind.cfg'), 'w') as outfile:
        outfile.write(cfg_template%(ppm, ppm, one_raw_out, raw_name, os.path.join(raw_folder, raw_name+f'_HCDFT.{fmt}'), fmt.upper()))
    pfind_cfg = os.path.join(one_raw_out, 'pFind.cfg')
    raw_file = os.path.join(raw_folder, raw_file)
    params_list.append((pfind_cfg, raw_file))
    
def generate_cfg():
    for raw_file in os.listdir(raw_folder):
        if raw_file.lower().endswith('.raw'):
            print(raw_file)
            raw_name = raw_file[:-4]
            generate_cfg_one_raw(raw_name, raw_file)
        elif raw_file.lower().endswith('.d'):
            print(raw_file)
            raw_name = raw_file[:-2]
            generate_cfg_one_raw(raw_name, raw_file, ppm=20, fmt='mgf')

    print(f'{len(folder_list)} raw files')


if __name__ == '__main__':
    import datetime
    start_time = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    generate_cfg()
    parallel_run(params_list)
    print(f'start time = {start_time}')
    print(f'end time = {datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")}')
