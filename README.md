"# TraitAssociations" 
To run this analysis - run in this order
1. RowColumnClean-up
2. DataImpute
3. ImputeCheck
4. TraitxTrait_Analysis
5. TraitxIsolationAnalysis
6. Iso_TxT_Analysis
7. ControlForIsolation_NegCheck and/or ControlForIsolation_PosCheck


CarbonMetabolismTrait_Network must be run after TraitxTrait_Analysis (need scripts 1-4 above), but does not need any other the other scripts below that (do not need scripts 5-7).

PhylogeneticSignalAnalysis can be run after ImputeCheck (needs scripts 1-3) does not need any scripts below that (do not need scripts 5-7)

In order to run the VariancePartitioning, you need to run all scripts including the PhylogeneticSignalAnalysis.
