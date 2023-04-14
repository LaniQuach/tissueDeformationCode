# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 14:25:14 2023

@author: laniq
"""
from numpy import asarray

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

finalGridSpacings = [15,20,25,30,35]
numHistogramBins = [16, 32, 64]
numResolutions = [1,3,5,7,9]
numSpatialSamples = asarray(range(1000, 3000, 500))
bSplineInterpOrder = asarray(range(1,4))
finalBSplineInterpOrder = asarray(range(0,4))

# finalGridSpacings = [28]
# numHistogramBins = [16]
# numResolutions = [6]
# numSpatialSamples = [2048]
# bSplineInterpOrder = [3]
# finalBSplineInterpOrder = [3]
# maxNumIterations = [1000]


