# Source-Sink

This algorithm aims to calculate the following features from interictal EEG:
- source index
- sink index
- source influence
- sink connectivity

A LTI Model is fit on 500 ms windows of Interictal Data: x(t + 1) = Ax(t) 

The avg absolute values of the A matrix are used to generate the source-sink space and calculate the relevant features!!

It replicates methods outlined in the following paper : Gunnarsdottir KM, Li A, Smith RJ, Kang JY, Korzeniewska A, Crone NE, Rouse AG, Cheng JJ, Kinsman MJ, Landazuri P, Uysal U, Ulloa CM, Cameron N, Cajigas I, Jagid J, Kanner A, Elarjani T, Bicchi MM, Inati S, Zaghloul KA, Boerwinkle VL, Wyckoff S, Barot N, Gonzalez-Martinez J, Sarma SV. Source-sink connectivity: a novel interictal EEG marker for seizure localization. Brain. 2022 Nov 21;145(11):3901-3915. doi: 10.1093/brain/awac300. PMID: 36412516; PMCID: PMC10200292.


  
