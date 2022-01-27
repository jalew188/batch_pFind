$pfind_folder="C:/pFindStudio/pFind3/bin/"
Copy-Item "reduced_mod.ini" -Destination $pfind_folder"reduced_mod.ini" -Force
Copy-Item "batch_pFind.py" -Destination $pfind_folder"batch_pFind.py" -Force


& "python" "batch_pFind.py" $Args[0] "D:/fasta/uniprot_human_reviewed_20210309.fasta" U C 0