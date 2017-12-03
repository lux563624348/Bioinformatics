TO use this function, please copy the following three files to the folder you want to save your results: 
1. gene_promoterGenebody_U_4k_iv_unique.bed
2. read_count.py
3. For_Thomas.sh

" INPUT_A.bed " is the Input peaks file you want to do the overlap and counting. 

!!!EACH TIME CHANGE THE NAME OF YOUR INPUT FILE AS "INPUT_A.bed" and keep it together with the above three files.

And to run it in terminal environment, enter the folder you keep the above three files + INPUT_A.bed, using command:

bash For_Thomas.sh 


Such as: " Xiang@LAPTOP-Q0TSHFKK:/mnt/d/bioproject/Nanping/version.beta/TEST/Example$ bash For_Thomas.sh "

Then the result is shown in the folder Output, with name: "RPKM_INPUT_A.bed.txt" 