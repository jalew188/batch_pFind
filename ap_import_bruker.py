from alphapept.io import MS_Data_File
from alphapept.feature_finding import find_features as ap_ff
from ap_to_pfind import use_ff

def import_bruker_raw(raw_path):
    hdf_path = raw_path[:-2]+".ms_data.hdf"
    ms_data_file = MS_Data_File(
        hdf_path,
        is_new_file=True
    )
    ms_data_file.import_raw_DDA_data(
        raw_path,
        n_most_abundant = -1
    )
    if use_ff():
        find_features_bruker(raw_path)
    
def find_features_bruker(raw_path):
    settings = {}
    settings["experiment"] = {}
    settings["experiment"]["file_paths"] = [raw_path]
    settings['workflow'] = {}
    settings['workflow']["find_features"] = True
    settings['features'] = {}
    ap_ff((0,settings))
    
if __name__ == "__main__":
    import sys
    raw_path = sys.argv[1]
    load_bruker(raw_path)