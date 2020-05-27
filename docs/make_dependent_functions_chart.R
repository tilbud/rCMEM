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
  addCohort;animateCohorts; availableSediment; buildHighTideScenario; buildScenarioCurve;   convertProfile_AgeToDepth; deliveredSediment3TidalCycle; deliveredSedimentSimple; deliverSedimentFlexibly; depthOfNonRootVolume; floodTimeFromDatum; massLiveRoots; predictedBiomass; predictLunarNodalCycle; runMemWithCohorts; runToEquilibrium; sedimentInputs; trimCohorts; zStarToZ; zToZstar


  # 'edge' statements connect function tree
  depthOfNonRootVolume->addCohort; 
  massLiveRoots->addCohort; 
  sedimentInputs->addCohort;
  runMemWithCohorts->animateCohorts; 
  predictLunarNodalCycle->buildHighTideScenario;
  deliveredSediment3TidalCycle->deliverSedimentFlexibly; 
  deliveredSedimentSimple->deliverSedimentFlexibly;
  availableSediment->deliveredSediment3TidalCycle;
  floodTimeFromDatum->deliveredSediment3TidalCycle;
  availableSediment->deliveredSedimentSimple;
  massLiveRoots->depthOfNonRootVolume;
  addCohort->runMemWithCohorts;
  buildHighTideScenario->runMemWithCohorts;
	buildScenarioCurve->runMemWithCohorts; 
	convertProfile_AgeToDepth->runMemWithCohorts;
	deliverSedimentFlexibly->runMemWithCohorts; 
	predictedBiomass->runMemWithCohorts; 
	runToEquilibrium->runMemWithCohorts; 
	trimCohorts->runMemWithCohorts; 
	zToZstar->runMemWithCohorts;
  addCohort->runToEquilibrium;
  
}
")

