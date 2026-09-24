# Growth-reproduction trade-offs in woody plants: meta-analysis

### Data and R code 

**This analysis is related to the study:**

"Growth–reproduction trade-offs are common but changing in woody plants: a meta-analysis"
Barczyk et al. 2026

## Summary

Growth and reproduction draw on a common resource pool, yet empirical studies of woody plants report widely differing relationships between seed production and growth. Here we synthesize 685 estimates from 78 studies covering 79 woody species, with study locations concentrated at northern temperate latitudes, to test how growth–reproduction correlations vary across time, species, and environments. Growth and reproduction measured within the same year were negatively correlated, suggesting an immediate cost of reproduction. The strength and direction of growth–reproduction correlations showed at most weak phylogenetic structure among the sampled species and were not systematically related to functional traits and climate. Instead, trade-offs were strongest in species with high interannual variability in seed production and weakened markedly over recent decades in these species. Together, these results show that growth–reproduction trade-offs in woody plants are common but not fixed, and that shifts in reproductive variability under environmental change can alter how trees balance growth and reproduction, with consequences for long-term forest functioning.

## Data description

Main folder:
- **meta_trade-offs.Rproj** 
- **trade-offs_meta_analysis.Rmd** - the Rmd file with the code
- **phyloTreeMetaAN.R** - the R script with the function to plot the reduced phylogenetic tree

DATA_FILES:
- **trade-offs_meta-analysis_data.xlsx** - the spreadsheet with data of all studies screened in the meta-analysis (according to the PRISMA)
- **PRISMA.pdf** - the graph describing the protocol of the meta-analysis
- **Vascular_Plants_rooted.dated** - the phylogenetic tree from Zanne et al. 2014
- **speciesTraits_trade-offs_meta.csv** - data file with the information on species- and genus-level funtional traits
- **worldclim_trade-offs_meta.csv** - data file with climatic variables downloaded from WorldClim for all study locations
- **vpd-locations_trade-offs_meta.csv** - data file with mean values of VapourPressureDeficit derived from TerraClim for all study locations

## References

Zanne, A., Tank, D., Cornwell, W. et al. Three keys to the radiation of angiosperms into freezing environments. 
Nature 506, 89–92 (2014). https://doi.org/10.1038/nature12872

## Contact

For further questions, please contact Maciej Barczyk (maciejkbarczyk@gmail.com) and Michal Bogdziewicz (michalbogdziewicz@gmail.com)
