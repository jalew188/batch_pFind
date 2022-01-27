$pfind_folder="C:/pFindStudio/pFind3/bin/"
Copy-Item "reduced_mod.ini" -Destination $pfind_folder"reduced_mod.ini" -Force
Copy-Item "batch_pFind.py" -Destination $pfind_folder"batch_pFind.py" -Force
Copy-Item "ap_import_bruker.py" -Destination $pfind_folder"ap_import_bruker.py" -Force
Copy-Item "ap_to_pfind.py" -Destination $pfind_folder"ap_to_pfind.py" -Force


& "python" "batch_pFind.py" $Args[0] "D:/fasta/uniprot_human_reviewed_20210309.fasta" KR C 3