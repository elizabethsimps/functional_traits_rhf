# functional_traits_rhf

This repository includes raw functional trait data, metadata, and analysis for vegetation plots near Right Hand Fork, Logan, UT, USA to reproduce analyses for the paper titled "Variation in near-surface soil temperature drives plant assemblage insurance potential" by Elizabeth G. Simpson, Ian Fraser, Hillary Woolf, and William D. Pearse

This work is licensed under the Creative Commons Attribution-ShareAlike 3.0 Unported License. To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/3.0/ or send a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.

The script "traits_clean.R" (a) uses the stalkless pipeline (https://github.com/willpearse/stalkless) to analyze the leaf area of scanned leaves and (b) processes the other raw trait data collected at the vegetation plots and (c) combines all trait data into a clean data file located in "/clean_data/clean_traits_rhf.csv".

The script "env.R" processes and combines all soil texture and temperature data into a clean microenvironment data file located in "/cleandata/clean_traits_rhf.csv".

To replicate the analysis in the paper, run "fdiv-env-phy_core18.R" and "fdiv-env_life-hist_core18.R" sequentially.

Please reach out to Elizabeth Simpson (elizabeth.simpson@usu.edu) with any questions!
