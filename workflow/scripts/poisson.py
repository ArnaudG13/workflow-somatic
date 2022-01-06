#!/usr/bin/env python
# -*- coding: utf-8 -*-

#see https://github.com/tamsen/CallSomaticVariants/blob/master/CallSomaticVariants/StrandBias.cs

import math

global Epsilon 
Epsilon = 1e-20
global Fpmin 
Fpmin = 1e-50
global LanczCutoff 
LanczCutoff = 700
global Itmax 
Itmax = 300

def cdf(numOccurrences, numExpectedOccurrences) :
	
	return IncompleteGammaFunction(numOccurrences + 1.0, numExpectedOccurrences)


def IncompleteGammaFunction(a, x) :

	if x < 0 or a <= 0 :
		return -1.0
		
	if a >= LanczCutoff :
		g = StirlingApproximation(a)
	else :
		g = LanczosApproximation(a)
		
	if x >= a + 1.0 :
		return GammaUsingContinuedFractions(a, x, g)
		
	g = GammaSeries(a, x, g)
	if g < 0 :
		return g
		
	return 1.0-g


def GammaUsingContinuedFractions(a, x, g) :
	
	b = x + 1.0 - a
	c = 1.0 / Fpmin
	d = 1.0 / b
	h = d
	
	for i in range(1,Itmax+1) :
		an = i * (a - i)
		b = b + 2.0
		d = an * d + b
		if abs(d) < Fpmin :
			d = Fpmin
		c = b + an / c
		if abs(c) < Fpmin :
			c = Fpmin
		d = 1.0 / d
		Del = d * c
		h = h * Del
		if abs(Del - 1.0) < Epsilon :
			break
			
	if i > Itmax :
		return -1.0
		
	return (math.exp(a * math.log(x) - x - g) * h)
	
	
def GammaSeries(a, x, g) :
	
	retval = -1.0
	
	if x == 0 :
		return 0.0
	if x < 0 :
		return retval
		
	ap = a
	Sum = 1.0 / a
	Del = Sum
	
	for i in range(1,Itmax+1) :
		ap = ap + 1.0
		Del = Del * (x / ap)
		Sum = Sum + Del
		
		if abs(Del) < (abs(Sum) * Epsilon) :
			retval = Sum * math.exp(a * math.log(x) - x - g)
			break
		
	return retval
	

def LanczosApproximation(p) :
	
	x = p
	tmp = x + 5.5
	tmp = tmp - (x + 0.5) * math.log(tmp)
	
	ser = 1.000000000190015 + 76.18009172947146 / (p + 1.0)
	ser = ser - 86.50532032941678 / (p + 2.0)
	ser = ser + 24.01409824083091 / (p + 3.0)
	ser = ser - 1.231739572450155 / (p + 4.0)
	ser = ser + 0.001208650973866179 / (p + 5.0)
	ser = ser - 5.395239384953E-06 / (p + 6.0)
	
	return (math.log(2.506628274631001 * ser / x) - tmp)
	
	
def StirlingApproximation(n) :
	
	return (0.5 * math.log(2.0 * math.pi) + (0.5 + n) * math.log(n) - n)
	
def ApproximateLNofNFactorial(a) :
	
	if a >= LanczCutoff :
		g = StirlingApproximation(a)
	else :
		g = LanczosApproximation(a)
		
	if a > 100 :
		return g
		
	g = 0;
	for i in range(1,a) :
		g = g + math.log(i)
		
	return g
