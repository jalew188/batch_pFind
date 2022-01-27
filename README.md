0. Install Python 3.8

1. (Only for Bruker data) Install alphapept to support Bruker data:
Create a new conda env if possible (python=3.8)
`pip install alphapept`  

2. Open 'run_pFind.ps1' in any text editing APP (e.g. notepad.exe):
```
$pfind_bin_folder="C:/pFindStudio/pFind3/bin/"
$fasta="D:/fasta/uniprot_human_reviewed_20210309.fasta"


Copy-Item "reduced_mod.ini" -Destination $pfind_bin_folder"reduced_mod.ini" -Force
Copy-Item "batch_pFind.py" -Destination $pfind_bin_folder"batch_pFind.py" -Force
Copy-Item "ap_import_bruker.py" -Destination $pfind_bin_folder"ap_import_bruker.py" -Force
Copy-Item "ap_to_pfind.py" -Destination $pfind_bin_folder"ap_to_pfind.py" -Force


& "python" "batch_pFind.py" --raw_folder $Args[0] --fasta $fasta --enzyme_aa KR --enzyme_term C --enzyme_specific 3 --pfind_bin_folder $pfind_bin_folder

```
Change $pfind_bin_folder and $fasta to your own folder and fasta paths, save and close.

3. Run to start the search:
```
powershell ./run_pFind.ps1 ...\folder\contains\.d\or\.raw
```

4. The results can be found in pFind3 folder in the raw data folder