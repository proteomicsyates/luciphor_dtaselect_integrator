# ﻿luciphor\_dtaselect\_integrator  
  
  
Program that integrates results from Luciphor (PTM localization scoring algorithm) into DTASelect output files.

Input:  

- luciphor output file (**-luc** path)
- dtaselect output file (**-dta** path)
- localLFR threshold \[optional\] (**-lflr** double) 
- globalLFR threshold \[optional\] (**-gflr** double)
- remove or not PSMs that don't pass the thresholds (**-rem**)

Output:

- It overrides the sequences of the PSMs on the DTASelect file that pass the thresholds (if available) on the luciphor output. 
- It also adds 6 new columns ('original_sequence', 'luciphor_pep1Score', 'luciphor_pep2Score', 'luciphor_deltaScore', 'globalFLR', 'localFLR') for ONLY the PSMs that pass the luciphor score thresholds.
- If the PTM localization proposed by Luciphor is the same as the original one, the column original_sequence will remain empty. 
- The original DTASelect file is backed up to a new file adding “_original” to the end of its name.
 





