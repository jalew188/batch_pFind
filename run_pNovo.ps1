$pnovo_folder="C:/pFindStudio/pNovo/bin/"
Copy-Item "batch_pNovo.py" -Destination $pnovo_folder"batch_pNovo.py" -Force

& "python" "batch_pNovo.py" $Args[0]