import os
import sys
import concurrent.futures as cf
from alphapept.pyrawfilereader import RawFileReader


raw_folder = sys.argv[1]
if len(sys.argv) > 2:
    enzyme_aa = sys.argv[2] #'_' for non-specific
    enzyme_term = sys.argv[3]
else:
    enzyme_aa = 'KR'
    enzyme_term = 'C'

n_process = 6

pnovo_path=r'C:\pFindStudio\pNovo\bin'

pNovo_main='pNovoplus.exe'
pparse_main='pParse.exe'
enzyme=f'Enzyme {enzyme_aa} {enzyme_term}'
raw_reload = False

modifications={
    #'whatever readable name': 'letter=mass'
    'M+16':'m=147.0354',
    'C+57':'C=160.030654',
}

output_folder = os.path.join(raw_folder, 'pNovo')
if not os.path.isdir(output_folder):
    os.makedirs(output_folder)
    
modlist = '\n'.join(modifications.values())


os.chdir(pnovo_path)

def raw_to_ms1_ms2(raw):
    with open(raw[:-4]+"_HCDFT.mgf", "w", 128*1024*1024) as ms2:
        rawFile = RawFileReader(raw)
        raw_name = os.path.split(raw)[1][:-4]
        def get_one_scan(scan):
            ms_order = rawFile.GetMSOrderForScanNum(scan)
            if ms_order != 2: return
            mass_array, inten_array = rawFile.GetCentroidMassListFromScanNum(scan)
            RTInSeconds = rawFile.RTInSecondsFromScanNum(scan)
            meta_data = rawFile.GetTrailerExtraForScanNum(scan)
            ion_injection_time = meta_data['Ion Injection Time (ms):']
            mass_analyzer_type = rawFile.GetMassAnalyzerTypeForScanNum(scan)
            scan_event = rawFile.GetScanEventStringForScanNum(scan)
            act_sub = scan_event[scan_event.find("@"):].upper()
            if "@HCD" in act_sub and "@ETD" in act_sub:
                activation_type = "ETHCD"
            elif "@CID" in act_sub and "@ETD" in act_sub:
                activation_type = "ETCID"
            elif "@HCD" in act_sub:
                activation_type = "HCD"
            elif "@ETD" in act_sub:
                activation_type = "ETD"
            elif "@CID" in act_sub:
                activation_type = "CID"
            else:
                return
            prec_mz = rawFile.GetPrecursorMassForScanNum(scan)
            mono_mz = float(meta_data['Monoisotopic M/Z:'])
            charge = int(meta_data['Charge State:'])
            if charge == 0 or mono_mz < 1: 
                if charge == 0:
                    charge = 2
                if mono_mz < 1:
                    mono_list = [prec_mz, prec_mz-1.0033/charge]
            else:
                mono_list = [mono_mz]
            for i,mono_mz in enumerate(mono_list):
                ms2.write("BEGIN IONS\n")
                ms2.write(f"TITLE={raw_name}.{scan}.{scan}.{charge}.{i}.dta\n")
                ms2.write(f"SCAN={scan}\n")
                ms2.write(f"ACTIVATION={activation_type}\n")
                ms2.write(f"RTINSECONDS={RTInSeconds}\n")
                ms2.write(f"CHARGE={charge}+\n")
                ms2.write(f"PEPMASS={mono_mz}\n")
                for mass, inten in zip(mass_array, inten_array):
                    ms2.write("%.5f %.1f\n"%(mass, inten))
                ms2.write("END IONS\n")
        for scan in range(rawFile.FirstSpectrumNumber, rawFile.LastSpectrumNumber+1):
            get_one_scan(scan)
            if scan % 1000 == 0:
                print("[RAW] <Loading> {:.2f}%".format(100.0*scan/rawFile.LastSpectrumNumber), flush=True, end='\r')
        print("[RAW] <Loading> 100%", flush=True)

cfg_template = '''
[meta]

#If you want to add a variable modification, 
#please use a letter from (a-z) instead.
#For example, if M+Oxidation is to be added,
#you can add the line below(without '#'), 
#in which 147.0354 = mass(M) + mass(Oxidation)

#n=115.026946
#m=147.0354
#c=160.030654
{modlist}

#N- or C- terminal variable modifications can be added as follows (using 0-9)

#0=42.010565
#1=43.005814


#The lines below show the basic ion types of HCD and ETD data.
HCDIONTYPE=4
HCDIONTYPE1=b	1 1 1 0.0
HCDIONTYPE2=y	1 0 1 18.0105647
HCDIONTYPE3=b	2 1 1 0.0
HCDIONTYPE4=y	2 0 1 18.0105647
ETDIONTYPE=6
ETDIONTYPE1=c 1 1 1 17.026549105
ETDIONTYPE2=z 1 0 1 1.99129206512
ETDIONTYPE3=c-1 1 1 0 16.01872407
ETDIONTYPE4=z+1 1 0 1 2.999665665
ETDIONTYPE5=c 2 1 0 17.026549105
ETDIONTYPE6=z 2 0 0 1.99129206512

#[IMPORTANT]
#An enzyme can be set as: 
#[EnzymeName] [CleavageSites] [N/C] (Cleave at N- or C- terminal)
enzyme={enzyme}

#if you want to use multi-threads, please set the number of threads below (1 ~ 8):
thread=2

#Mass ranges of precursors
#Only the spectra whose precursors are in the specified mass range will be sequenced.
mass_lower_bound=500
mass_upper_bound=3500

#[HCD, CID, ETD]
activation_type=HCD

#[IMPORTANT]
#Tolerance of precursors. 
#If you want to use Daltons, please set 'pep_tol_type_ppm' as 0
pep_tol=7
pep_tol_type_ppm=1

#[IMPORTANT]
#Tolerance of fragment ions. 
#If you want to use Daltons, please set 'frag_tol_type_ppm' as 0
frag_tol=20
frag_tol_type_ppm=1

[file]

#DTA/MS2/MGF are valid options.if DTA is specified, 
#please set the following path(s) as the folder containing the DTA file(s)
spec_type=MGF

#1:means only one activation type, CID/HCD/ETD, is used
#		spec_path1 should be set as the path of the data
#2:(HCD + ETD) is used. In this case, activation_type is ignored.
#		spec_path1 should be set as the path of the HCD data,
#		and spec_path2 should be set as the path of the ETD data.
spec_path=1
spec_path1={mgf}
spec_path2=

#If only one activation type of spectra is used (spec_path=1),
#you can specify a folder containing several MS2 or MGF files.
#Set spec_path1 as the foler,
#and pNovo+ will sequence the MS/MS data files one by one. 
#if folder=no, then the value of 'spec_path1' above must be a MS/MS file path. 
folder=no

#The folder where the result files are to be output
out_path={output_folder}

#The number of peptides reported per spectrum
report_pep=10
'''

def run_pNovo(pnovo_cfg, raw_file):
    # if raw_reload or not os.path.isfile(raw_file+'.done'):
        # raw_to_ms1_ms2(raw_file)
        # with open(raw_file+'.done','w') as f: f.write(raw_file)
    os.system(f"{pNovo_main} {pnovo_cfg}")
    return raw_file
    
def parallel_run(params_list):
    futures = []
    with cf.ProcessPoolExecutor(n_process) as pp:
        futures = []
        for pnovo_cfg, raw_file in params_list:
            futures.append(pp.submit(run_pNovo, pnovo_cfg, raw_file))
        for future in cf.as_completed(futures):
            print(f"\n[DONE] {os.path.split(future.result())[1]}\n")

params_list = []
def generate_cfg():
    for raw_file in os.listdir(raw_folder):
        if raw_file.lower().endswith('.raw'):
            print(raw_file)
            raw_name = raw_file[:-4]
            pNovo_cfg = os.path.join(output_folder, raw_name+'.pNovo.cfg')
            with open(pNovo_cfg, 'w') as outfile:
                outfile.write(cfg_template.format(mgf=os.path.join(raw_folder, raw_name+'_HCDFT.mgf'), enzyme=enzyme, output_folder=output_folder, modlist=modlist))
            raw_file = os.path.join(raw_folder, raw_file)
            params_list.append((pNovo_cfg, raw_file))

    print(f'{len(params_list)} raw files')
    
def merge_results():
    pass

if __name__ == '__main__':
    import datetime
    start_time = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    generate_cfg()
    parallel_run(params_list)
    print(f'start time = {start_time}')
    print(f'end time = {datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")}')
    #merge_results()
