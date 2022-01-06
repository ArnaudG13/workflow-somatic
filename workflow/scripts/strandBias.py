#see https://github.com/tamsen/CallSomaticVariants/blob/master/CallSomaticVariants/StrandBias.cs

import poisson
import math

def getStats(support, coverage, noiseFreq, minDetectableSNP) :
	
	if support == 0 :
		
		#chanceFalsePos = 1
		#chanceVarFreqGreaterThanZero = 0
		#chanceFalseNeg = 0
		
		chanceVarFreqGreaterThanZero = pow((1-minDetectableSNP),coverage)
		chanceFalsePos = 1 - chanceVarFreqGreaterThanZero
		chanceFalseNeg = chanceVarFreqGreaterThanZero
	
	else :
		
		chanceVarFreqGreaterThanZero = poisson.cdf(support-1, coverage * noiseFreq)
		chanceFalsePos = 1 - chanceVarFreqGreaterThanZero
		chanceFalseNeg = poisson.cdf(support, coverage * minDetectableSNP)
		
	return [chanceVarFreqGreaterThanZero, chanceFalsePos, chanceFalseNeg]
	


def strandBiasScore(global_stats, forward_stats, reverse_stats) :
	
	forwardBias = (forward_stats[0] * reverse_stats[1])/global_stats[0]
	reverseBias = (reverse_stats[0] * forward_stats[1])/global_stats[0]
	p = max(forwardBias,reverseBias,1e-10)
	return(10*math.log10(p))
	

def computeStrandBias(alt_p,alt_m,ref_p,ref_m,error_rate=0.01) :
	
	cov_p = alt_p + ref_p
	cov_m = alt_m + ref_m
	
	stats_global = getStats(alt_p+alt_m,cov_p+cov_m,0.01,0.01)
	stats_forward = getStats(alt_p,cov_p,0.01,0.01)
	stats_reverse = getStats(alt_m,cov_m,0.01,0.01)
	
	score = strandBiasScore(stats_global, stats_forward, stats_reverse)
	
	return score

	
