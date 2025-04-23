#!/bin/bash
Rscript --vanilla ./SA_script.R mara 0 2000 &
Rscript --vanilla ./SA_script.R mara 0.25 2000 &
Rscript --vanilla ./SA_script.R mara 0.5 2000 &
Rscript --vanilla ./SA_script.R mara 0.75 2000 &
Rscript --vanilla ./SA_script.R mara 1 2000 

Rscript --vanilla ./SA_script.R sobol 0 2000 &
Rscript --vanilla ./SA_script.R sobol 0.25 2000& 
Rscript --vanilla ./SA_script.R sobol 0.5 2000 
Rscript --vanilla ./SA_script.R sobol 0.75 2000& 
Rscript --vanilla ./SA_script.R sobol 1 2000 

