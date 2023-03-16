# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 14:25:14 2023

@author: laniq
"""

#all parameters to be tested
# metrics = ["AdvancedKappaStatistic", "AdvancedMattesMutualInformation", 
#            "AdvancedMeanSquares", "AdvancedNormalizedCorrelation", 
#            "CorrespondingPointsEuclideanDistanceMetric", "DisplacementMagnitudePenalty", 
#            "DistancePreservingRigidityPenalty", "KNNGraphAlphaMutualInformation",
#            "MissingStructurePenalty", "NormalizedMutualInformation", 
#            "PolydataDummyPenalty", "StatisticalShapePenalty", "SumSquaredTissueVolumeDifference",
#            "TransformBendingEnergyPenalty", "TransformRigidityPenalty"]

# metrics = ["AdvancedMattesMutualInformation"]

# optimizers = ["AdaGrad", "AdaptiveStochasticGradientDescent", 
#               "AdaptiveStochasticLBFGS", "AdaptiveStochasticVarianceReducedGradient", 
#               "CMAEvolutionStrategy", "ConjugateGradient", "ConjugateGradientFRPR"
#               "FiniteDifferenceGradientDescent", "PreconditionedGradientDescent", 
#               "PreconditionedStochasticGradientDescent", "QuasiNewtonLBFGS",
#               "RegularStepGradientDescent", "RSGDEachParameterApart", 
#               "SimultaneousPerturbation", "StandardGradientDescent"]
# maxNumIterations = range(200, 100, 2100)


finalGridSpacings = range(5,5,45)
numHistogramBins = [16, 32, 64]
numResolutions = range(1,7)
numSpatialSamples = range(1600, 50, 4050)
bSplineInterpOrder = range(1,6)
finalBSplineInterpOrder = range(0,6)


