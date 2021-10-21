###get the first 2 columns from ELI_scores.xlsx, then we will get test.txt
sort -k2,2 -nr test.txt |head -n 100 |cut -f 1 >up100.txt
sort -k2,2 -nr test.txt |tail -n 100|cut -f 1 >dn100.txt
