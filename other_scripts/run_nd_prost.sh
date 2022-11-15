#!/bin/bash
#
# Example shell script to run DPClust on a multi-sample case
#

# cd /pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Results/mutsig_clonality_development/Run_DPClust_on_nchr1_ccf

# DO51953
# number 1 represent the order of DO51953 donor in prost_met_dpinput_summ.txt

/pub6/Temp/Liaojl/miniconda3/bin/R --vanilla --slave -q -f dpclust_pipeline.R --args -r 1 -d Input/Data/ -o Output -i Input/prost_met_dpinput_summ.txt

# DO51954

/pub6/Temp/Liaojl/miniconda3/bin/R --vanilla --slave -q -f dpclust_pipeline.R --args -r 2 -d Input/Data/ -o Output -i Input/prost_met_dpinput_summ.txt

# DO51956

/pub6/Temp/Liaojl/miniconda3/bin/R --vanilla --slave -q -f dpclust_pipeline.R --args -r 3 -d Input/Data/ -o Output -i Input/prost_met_dpinput_summ.txt

# DO51958

/pub6/Temp/Liaojl/miniconda3/bin/R --vanilla --slave -q -f dpclust_pipeline.R --args -r 4 -d Input/Data/ -o Output -i Input/prost_met_dpinput_summ.txt

# DO51959

/pub6/Temp/Liaojl/miniconda3/bin/R --vanilla --slave -q -f dpclust_pipeline.R --args -r 5 -d Input/Data/ -o Output -i Input/prost_met_dpinput_summ.txt

# DO51960

/pub6/Temp/Liaojl/miniconda3/bin/R --vanilla --slave -q -f dpclust_pipeline.R --args -r 6 -d Input/Data/ -o Output -i Input/prost_met_dpinput_summ.txt

# DO51962

/pub6/Temp/Liaojl/miniconda3/bin/R --vanilla --slave -q -f dpclust_pipeline.R --args -r 7 -d Input/Data/ -o Output -i Input/prost_met_dpinput_summ.txt

# DO51964

/pub6/Temp/Liaojl/miniconda3/bin/R --vanilla --slave -q -f dpclust_pipeline.R --args -r 8 -d Input/Data/ -o Output -i Input/prost_met_dpinput_summ.txt

# DO51965

/pub6/Temp/Liaojl/miniconda3/bin/R --vanilla --slave -q -f dpclust_pipeline.R --args -r 9 -d Input/Data/ -o Output -i Input/prost_met_dpinput_summ.txt

