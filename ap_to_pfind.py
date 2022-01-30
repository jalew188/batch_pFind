import sys
import h5py
import numba
import numpy as np
import struct
import os
import tqdm

valid_spectrum_has_at_least_n_peaks = 10
_use_ff = False
def use_ff():
    return _use_ff

def ap_bruker_to_mgf(hdf):
    raw_file = hdf[:-len(".ms_data.hdf")]
    mgf = raw_file+'_HCDFT.mgf'

    def get_list(dataset):
        return h5py.File(hdf).get(dataset)[...]
        
    fmt = 'mgf'
    f = open(mgf, 'w')
    raw = os.path.basename(raw_file)

    if use_ff():
        import alphapept.io
        features = alphapept.io.MS_Data_File(hdf).read(dataset_name='features')
        groups = features.sort_values('query_idx').groupby('query_idx')
        ff_idx_set = set(groups.groups.keys())
    else:
        ff_idx_set = set()
    mass_list = get_list('Raw/MS2_scans/mass_list_ms2')
    inten_list = get_list('Raw/MS2_scans/int_list_ms2')
    scan_list = get_list('Raw/MS2_scans/scan_list_ms2')
    rt_list = get_list('Raw/MS2_scans/rt_list_ms2')
    mobility_list = get_list('Raw/MS2_scans/mobility2')
    idx_list = get_list('Raw/MS2_scans/indices_ms2')
    idx_start = idx_list[:-1]
    idx_end = idx_list[1:]
    precursor_list = get_list('Raw/MS2_scans/mono_mzs2')
    charge_list = get_list('Raw/MS2_scans/charge2')
    
    print(f'{raw}.d contains {len(scan_list)} MS2 spectra')

    if fmt == 'pf2':
        ntitle = len(raw)
        f.write(struct.pack('2i', len(scan_list), ntitle))
        f.write(struct.pack(f'{len(raw)}s', bytes(raw, 'utf-8')))
    for ms2_idx in range(len(scan_list)):
        start, end = idx_start[ms2_idx], idx_end[ms2_idx]
        if end - start < valid_spectrum_has_at_least_n_peaks:
            continue
        masses = mass_list[start:end]
        intens = inten_list[start:end]
        scan = scan_list[ms2_idx]
        if ms2_idx in ff_idx_set:
            spec_prec_mz_list = []
            spec_charge_list = []
            for mono_mz, charge in groups.get_group(ms2_idx)[['mz_matched','charge_matched']].values:
                if np.isnan(mono_mz): continue
                spec_prec_mz_list.append(mono_mz)
                spec_charge_list.append(int(charge))
        else:
            if np.isnan(precursor_list[ms2_idx]): continue
            spec_prec_mz_list = [precursor_list[ms2_idx]]
            spec_charge_list = [int(charge_list[ms2_idx]) if not np.isnan(charge_list[ms2_idx]) else 2]

        if fmt == "pf2":
            peaks = np.zeros(len(masses)*2)
            peaks[np.arange(0,len(peaks),2)] = masses
            peaks[np.arange(1,len(peaks),2)] = intens
            f.write(struct.pack('i', scan))
            f.write(struct.pack('i', len(masses)))
            f.write(struct.pack(f'{len(peaks)}d', *peaks))
            f.write(struct.pack('i', len(spec_prec_mz_list)))
            for mz, charge in zip(spec_prec_mz_list, spec_charge_list):
                f.write(struct.pack('d', mz))
                f.write(struct.pack('i', charge))
        else:
            for i in range(len(spec_prec_mz_list)):
                mz = spec_prec_mz_list[i]
                charge = spec_charge_list[i]
                f.write("BEGIN IONS\n")
                f.write(f"TITLE={raw}.{ms2_idx}.{ms2_idx}.{charge}.{i}.dta\n")
                f.write(f"SCAN={ms2_idx}\n")
                f.write(f"APSCAN={scan}\n")
                f.write(f"CHARGE={charge}+\n")
                f.write(f"RTINSECONDS={rt_list[ms2_idx]*60:.3f}\n")
                f.write(f"MOBILITY={mobility_list[ms2_idx]:.6f}\n")
                f.write(f"PEPMASS={mz:.6f}\n")
                for mass, inten in zip(masses, intens):
                    f.write(f"{mass:.6f}\t{inten:.0f}\n")
                f.write("END IONS\n")
    f.close()
    with open(raw_file+'.done','w') as f:
        f.write(f'Converted from {raw_file}.d by ap_to_pfind')

if __name__ == '__main__':
    hdf = sys.argv[1]
        
    ap_bruker_to_mgf(hdf)
