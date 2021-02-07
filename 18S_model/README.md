###### Identification of Variable Region Coordinates

Eukaryotic variable region coordinates are based off of the AF258606 18S sequence.
The coordinates of each region are taken from https://academic.oup.com/femsle/article/238/2/455/492077

Steps to obtain the variable region coordinates relative to the RF01960 model:
1. Sequences for each region are extracted from AF258606.fa using the coordinates obtained from the paper referenced above.
The regions are saved to VR.AF258606.fa

2. The variable region sequences are mapped to the covariance model:

    `cmalign -o AF258606_to_RF01960_cmalign.txt RF01960.cm AF258606.fa > VR.AF258606.tab
`      
3. "cm from" and "cm to" fields from the resulting tab file contain the coordinates relative to the model.
Note: some of the regions have a low bit score. To double check the coordinates, the output_hmm_output_hmm_AF258606.out file was used.