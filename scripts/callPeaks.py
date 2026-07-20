#! /python22/Lib

import string,sys,os,math,re,subprocess,random,cmd,time,statistics

def calculate(input, output):

	ratios = []
	f = open(input)
	for thisline in f:
		if (thisline[0] != "p"):
			data = str.split(thisline)
			if (data[2] != "NA"):
				ratios.append(math.log(float(data[2]),2))
			else:
				ratios.append(0)
	f.close()

	currentRatios = ratios
	peakScores = []
	for z in range(10000):
		random.shuffle(currentRatios)
		prior = currentRatios[0]
		for i in range(1,len(currentRatios)):
			if ((prior < 0) and (currentRatios[i] < 0)):
				prior = prior + currentRatios[i]
			elif ((prior > 0) and (currentRatios[i] > 0)):
				prior = prior + currentRatios[i]
			else:
				peakScores.append(prior)
				prior = currentRatios[i]
		peakScores.append(prior)

	ratios = []
	f = open(input)
	for thisline in f:
		if (thisline[0] != "p"):
			data = str.split(thisline)
			if (data[2] != "NA"):
				ratios.append(math.log(float(data[2]),2))
			else:
				ratios.append(0)
	f.close()

	o = open(output, "w")
	o.write("peak_start\tpeak_stop\tpeak_score\tpeak_width\tp_value\n")
	prior = ratios[10]
	peakWidth = 1
	peakStart = 1
	for i in range(10,len(ratios)-10):
		if ((prior < 0) and (ratios[i] < 0)):
			prior = prior + ratios[i]
			peakWidth = peakWidth + 1
		elif ((prior > 0) and (ratios[i] > 0)):
			prior = prior + ratios[i]
			peakWidth = peakWidth + 1
		else:
			found = 0
			notFound = 0
			if (prior < 0):
				for j in range(len(peakScores)):
					if (peakScores[j] <= prior):
						found = found + 1
					else:
						notFound = notFound + 1
			elif (prior > 0):
				for j in range(len(peakScores)):
					if (peakScores[j] >= prior):
						found = found + 1
					else:
						notFound = notFound + 1
			else:
				found = 0
				notFound = 1
			if (prior != 0):
				if (float(found)/float(found+notFound) <= 0.001):
					if (abs(prior) >= 2):
						if (peakWidth >= 5):
							o.write(str(peakStart)+"\t"+str(i)+"\t"+str(prior)+"\t"+str(peakWidth)+"\t"+str(float(found)/float(found+notFound))+"\n")
			peakStart = i+1
			prior = ratios[i]
			peakWidth = 1
	o.close()

calculate(sys.argv[1], sys.argv[2])
sys.exit()

