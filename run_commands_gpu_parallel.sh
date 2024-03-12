cd build

#SARS-CoV dataset containing 200 sequences, kmerSize (default) =  15, sketchSize (default) = 1000, numBlocks (default) = 128
nsys nvprof --print-gpu-trace ./test-mash-gpu -f ../dataset/sars_200.fa > ../run_logs/parallel_sars200.txt

#SARS-CoV dataset containing 2000 sequences, kmerSize (default) =  15, sketchSize (default) = 1000, numBlocks (default) = 128
nsys nvprof --print-gpu-trace ./test-mash-gpu -f ../dataset/sars_2000.fa > ../run_logs/parallel_sars2000.txt

#SARS-CoV dataset containing 10000 sequences, kmerSize (default) =  15, sketchSize (default) = 1000, numBlocks (default) = 128
nsys nvprof --print-gpu-trace ./test-mash-gpu -f ../dataset/sars_10000.fa > ../run_logs/parallel_sars10000.txt

#SARS-CoV dataset containing 10000 sequences, kmerSize (default) =  15, sketchSize (default) = 1000, numBlocks = 32
nsys nvprof --print-gpu-trace ./test-mash-gpu -b 32 -f ../dataset/sars_2000.fa > ../run_logs/parallel_sars2000_b32.txt

#SARS-CoV dataset containing 10000 sequences, kmerSize (default) =  15, sketchSize (default) = 1000, numBlocks = 64
nsys nvprof --print-gpu-trace ./test-mash-gpu -b 64 -f ../dataset/sars_2000.fa > ../run_logs/parallel_sars2000_64.txt

#SARS-CoV dataset containing 10000 sequences, kmerSize (default) =  15, sketchSize (default) = 1000, numBlocks = 256
nsys nvprof --print-gpu-trace ./test-mash-gpu -b 256 -f ../dataset/sars_2000.fa > ../run_logs/parallel_sars2000_b256.txt

#SARS-CoV dataset containing 10000 sequences, kmerSize (default) =  15, sketchSize (default) = 1000, numBlocks = 512
nsys nvprof --print-gpu-trace ./test-mash-gpu -b 512 -f ../dataset/sars_2000.fa > ../run_logs/parallel_sars2000_b512.txt

#SARS-CoV dataset containing 2000 sequences, kmerSize (default) =  15, sketchSize = 500, numBlocks (default) = 128
nsys nvprof --print-gpu-trace ./test-mash-gpu -s 500 -f ../dataset/sars_2000.fa > ../run_logs/parallel_sars2000_s500.txt

#SARS-CoV dataset containing 2000 sequences, kmerSize (default) =  15, sketchSize = 2000, numBlocks (default) = 128
nsys nvprof --print-gpu-trace ./test-mash-gpu -s 2000 -f ../dataset/sars_2000.fa > ../run_logs/parallel_sars2000_s2000.txt