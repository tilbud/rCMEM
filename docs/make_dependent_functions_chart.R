# Make Dependent Functions Chart
#
# This makes a chart showing which functions are called in other functions.

# Make sure Diagrammer is there
require(DiagrammeR, quietly = TRUE)

grViz("
digraph boxes_and_circles {

  graph [overlap = true, fontsize = 10]

  # add function names
  node [shape = box,
        fontname = Helvetica]
  addCohort;animateCohorts; availableSediment; buildHighTideScenario; buildScenarioCurve;   convertProfileAgeToDepth; deliveredSediment3TidalCycle; deliveredSedimentSimple; deliverSedimentFlexibly; calculateDepthOfNonRootVolume; floodTimeFromDatum; massLiveRoots; predictedBiomass; predictLunarNodalCycle; runMemWithCohorts; runToEquilibrium; sedimentInputs; trimCohorts; convertZStarToZ; zToZstar


  # 'edge' statements connect function tree
  calculateDepthOfNonRootVolume->addCohort; 
  massLiveRoots->addCohort; 
  sedimentInputs->addCohort;
  runMemWithCohorts->animateCohorts; 
  predictLunarNodalCycle->buildHighTideScenario;
  deliveredSediment3TidalCycle->deliverSedimentFlexibly; 
  deliveredSedimentSimple->deliverSedimentFlexibly;
  availableSediment->deliveredSediment3TidalCycle;
  floodTimeFromDatum->deliveredSediment3TidalCycle;
  availableSediment->deliveredSedimentSimple;
  massLiveRoots->calculateDepthOfNonRootVolume;
  addCohort->runMemWithCohorts;
  buildHighTideScenario->runMemWithCohorts;
	buildScenarioCurve->runMemWithCohorts; 
	convertProfileAgeToDepth->runMemWithCohorts;
	deliverSedimentFlexibly->runMemWithCohorts; 
	predictedBiomass->runMemWithCohorts; 
	runToEquilibrium->runMemWithCohorts; 
	trimCohorts->runMemWithCohorts; 
	zToZstar->runMemWithCohorts;
  addCohort->runToEquilibrium;
  
}
")

