# Lumos
R code to test fragmentation parameters on the Lumos
Krista Longnecker
1 November 2017
Update on 11 October 2023: this will be moved to the organization now called KujawinskiLaboratory. All commits prior to this point were by Krista Longnecker, using the account name KujawinskiLaboratory.
 

Testing different parameters on the Lumos to see which set of fragmentation options will serve our needs best. Using several parameters to consider what will be the best choice. How many MS2 scans we get? Is orbitrap or ion trap detection best when using the in silico fragmentation tools? How clean is the peak in the MS1 scan? 

Run the files in this order:
--> run Lumos_MDexp1b_neg_mspurity to get the aligned XCMS data
--> run Lumos_MDexp1b_neg_miner to get the MS2 fragments for each sample (as RDdata files)
--> combineTools_neg

General methods:

Peak picking using the centWave algorithm (Tautenhahn et al., 2008) within XCMS (Smith et al., 2006), followed by the MSpurity program from Lawson et al. (2017) and the compMS2Miner program from Edmands et al. (2017).


 Edmands, W. M. B., L. Petrick, D. K. Barupal, A. Scalbert, M. J. Wilson, J. K. Wickliffe and S. M. Rappaport (2017). "compMS2Miner: An Automatable Metabolite Identification, Visualization, and Data-Sharing R Package for High-Resolution LC–MS Data Sets." Analytical Chemistry 89(7): 3919-3928.
 
 Lawson, T. N., R. J. M. Weber, M. R. Jones, A. J. Chetwynd, G. Rodrı́guez-Blanco, R. Di Guida, M. R. Viant and W. B. Dunn (2017). "msPurity: Automated Evaluation of Precursor Ion Purity for Mass Spectrometry-Based Fragmentation in Metabolomics." Analytical Chemistry 89(4): 2432-2439
 
  Smith, C. A., E. J. Want, G. O'Maille, R. Abagyan and G. Siuzdak (2006). "XCMS: processing mass spectrometry data for metabolite profiling using nonlinear peak alignment, matching, and identification." Analytical Chemistry 78: 779 - 787. 
  
  Tautenhahn, R., C. Bottcher and S. Neumann (2008). "Highly sensitive feature detection for high resolution LC/MS." BMC Bioinformatics 9(1): 504.
