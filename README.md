# PANDO
PANDO - MTL With Tree Based Gradient Boosting

This repository contains the **research code** behind the PANDO paper.
As of now, the code will not install the required packages on its own, so you will need to manually install required packages if they are not already installed on your machine.
all the necessary files are in this repository.

The stats from the paper were generated on Ubuntu Linux 16.04 64-bit, R version 3.3.2.

generateData.r --> this contains the data generation procedures for the simulated data section. It also includes the actual experiments for both linear and non linear data, as well as the assocaited plots. The plots will be saved locally to where the code will reside

spam2.r --> this contains the spam dataset experiment. You can run it to get the results shown in the paper

spam/spam.py --> parses the original spam dataset data (the one from ECML/PKDD), and dumps a csv to be used by spam.r. The original spam dataset is also found on this repository

spam.r --> this generates the data for the spam experiment and saves it as CSV

PANDO.r,PANDO2.r ---> the files which implement PANDO and PANDO2 respectively.

Since the code is not efficient, we used multi-core processing. this is part of the TunePando method found in helper.r. you can replace %dopar% with %do% and the code will run in single core if you experience troubles


The results in the paper were obtained on the following specs:

Linux 4.4.0-112-generic #135-Ubuntu SMP x86_64 x86_64 x86_64 GNU/Linux
R version 3.3.2 (2016-10-31) -- "Sincere Pumpkin Patch"
