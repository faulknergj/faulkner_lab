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

calculate("A549_lenti_ratios.txt", "A549_lenti_peaks.txt")
calculate("HCT116_lenti_ratios.txt", "HCT116_lenti_peaks.txt")
calculate("HEK293T_lenti_ratios.txt", "HEK293T_lenti_peaks.txt")
calculate("HEK293T_PB_ratios.txt", "HEK293T_PB_peaks.txt")
calculate("HEK293T_combined_ratios.txt", "HEK293T_combined_peaks.txt")
calculate("HeLa-JVM_lenti_ratios.txt", "HeLa-JVM_lenti_peaks.txt")
calculate("HeLa-JVM_PB_ratios.txt", "HeLa-JVM_PB_peaks.txt")
calculate("HeLa-JVM_combined_ratios.txt", "HeLa-JVM_combined_peaks.txt")
calculate("HTR8_lenti_ratios.txt", "HTR8_lenti_peaks.txt")
calculate("i3Neuron_lenti_ratios.txt", "i3Neuron_lenti_peaks.txt")
calculate("iPSC_lenti_1_ratios.txt", "iPSC_lenti_1_peaks.txt")
calculate("iPSC_lenti_2_ratios.txt", "iPSC_lenti_2_peaks.txt")
calculate("iPSC_lenti_combined_ratios.txt", "iPSC_lenti_combined_peaks.txt")
calculate("MCF7_lenti_ratios.txt", "MCF7_lenti_peaks.txt")
sys.exit()

