cd build

#SARS-CoV dataset containing 200 sequences, kmerSize (default) =  15, sketchSize (default) = 1000
nsys nvprof --print-gpu-trace ./test-mash-gpu-serial -f ../dataset/sars_200.fa > ../run_logs/serial_sars200.txt

#SARS-CoV dataset containing 2000 sequences, kmerSize (default) =  15, sketchSize (default) = 1000
nsys nvprof --print-gpu-trace ./test-mash-gpu-serial -f ../dataset/sars_2000.fa > ../run_logs/serial_sars2000.txt

#SARS-CoV dataset containing 10000 sequences, kmerSize (default) =  15, sketchSize (default) = 1000
nsys nvprof --print-gpu-trace ./test-mash-gpu-serial -f ../dataset/sars_10000.fa > ../run_logs/serial_sars10000.txt

#SARS-CoV dataset containing 2000 sequences, kmerSize (default) =  15, sketchSize = 500, numBlocks (default) = 128
nsys nvprof --print-gpu-trace ./test-mash-gpu-serial -s 500 -f ../dataset/sars_2000.fa > ../run_logs/serial_sars2000_s500.txt

#SARS-CoV dataset containing 2000 sequences, kmerSize (default) =  15, sketchSize = 2000, numBlocks (default) = 128
nsys nvprof --print-gpu-trace ./test-mash-gpu-serial -s 2000 -f ../dataset/sars_2000.fa > ../run_log/serial_sars2000_s2000.txt