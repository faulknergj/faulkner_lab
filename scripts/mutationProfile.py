#! /python22/Lib

import string,sys,os,math,re,subprocess,random,cmd,time,statistics

def generateProfile(input, reference, output):

	L1sequence = ""
	f = open(reference)
	for thisline in f:
		if (thisline[0] != ">"):
			L1sequence = L1sequence + str.strip(thisline)
	f.close()

	maskKey = {}
	maskKey["mismatch"] = {}
	for i in range(1,11):
		maskKey["D"+str(i)] = {}

	##mask start and end of reads where mapping is ambiguous
	maskKey["mismatch"][0] = L1sequence[0]
	maskKey["mismatch"][1] = L1sequence[1]
	maskKey["mismatch"][903] = L1sequence[903]
	maskKey["mismatch"][904] = L1sequence[904]

	for i in range(7):
		for j in range(1,11):
			maskKey["D"+str(j)][i] = L1sequence[i]
		for j in range(1,11):
			maskKey["D"+str(j)][904-i] = L1sequence[904-i]

	##mask positions adjacent to fragment junctions
	for i in range(2,11):
		for j in range(i-1):
			maskKey["D"+str(i)][j] = L1sequence[j]
			maskKey["D"+str(i)][218-j] = L1sequence[218-j]
			maskKey["D"+str(i)][219+j] = L1sequence[219+j]
			maskKey["D"+str(i)][443-j] = L1sequence[443-j]
			maskKey["D"+str(i)][444+j] = L1sequence[444+j]
			maskKey["D"+str(i)][668-j] = L1sequence[668-j]
			maskKey["D"+str(i)][669+j] = L1sequence[669+j]
			maskKey["D"+str(i)][904-j] = L1sequence[904-j]

	##mask positions where homopolymers causing ambiguous mapping
	for i in range(1,11):
		for j in range(i,len(L1sequence)-i):
			matchCount = 0
			for k in range(i):
				if (L1sequence[j-i+k] == L1sequence[j+k]):
					matchCount = matchCount + 1
			if (matchCount == i):
				for k in range(i):
					maskKey["D"+str(i)][j-i+k] = L1sequence[j-i+k]
					maskKey["D"+str(i)][j+k] = L1sequence[j+k]

	result = {}
	for i in range(905):
		result[i] = {}
		result[i]["A"] = 0
		result[i]["C"] = 0
		result[i]["G"] = 0
		result[i]["T"] = 0
		for j in range(1,11):
			result[i]["D"+str(j)] = 0

	readCount = 0
	f = open(input)
	for thisline in f:
		if (thisline[0] != "@"):
			data = str.split(thisline)
			start = int(data[3])
			if (data[1] == "0"):
				proceed = "Y"
				errorCount = 0
				currentResult = []
				markers = {}
				alignment = data[20][5:]+"E"
				for i in range(len(alignment)):
					if (alignment[i] in "=-+*E"):
						markers[i] = alignment[i]
				markerList = list(markers)
				markerList.sort()
				alignmentLength = 0
				start = int(data[3])
				for i in range(len(markerList)):
					if (markers[markerList[i]] == "="):
						for j in range(markerList[i]+1,markerList[i+1]):
							currentResult.append(alignment[j])
							alignmentLength = alignmentLength + 1
					elif (markers[markerList[i]] == "E"):
						pass
					elif (markers[markerList[i]] == "+"):
						errorCount = errorCount + 1
					elif (markers[markerList[i]] == "*"):
						currentResult.append(str.upper(alignment[markerList[i]+2]))
						errorCount = errorCount + 1
						alignmentLength = alignmentLength + 1
					elif (markers[markerList[i]] == "-"):
						errorCount = errorCount + 1
						deletionLength = markerList[i+1]-markerList[i]-1
						if (deletionLength <= 10):
							alignmentLength = alignmentLength
							for j in range(deletionLength):
								currentResult.append("D"+str(deletionLength))
						else:
							proceed = "N"
					else:
						pass
				stop = start + alignmentLength - 1
				if ((errorCount <= 10) and (start <= 20) and (stop >= 885) and (proceed == "Y")):
					readCount = readCount + 1
					for i in range(len(currentResult)):
						result[i+start-1][currentResult[i]] = result[i+start-1][currentResult[i]] + 1
	f.close()

	for i in range(len(result)):
		try:
			maskKey["mismatch"][i]
			result[i]["A"] = "NA"
			result[i]["C"] = "NA"
			result[i]["G"] = "NA"
			result[i]["T"] = "NA"
		except:
			pass
		for j in range(1,11):
			try:
				maskKey["D"+str(j)][i]
				result[i]["D"+str(j)] = "NA"
			except:
				pass

	o = open(output, "w")
	o.write("position\treference\tA\tC\tG\tT\tD1\tD2\tD3\tD4\tD5\tD6\tD7\tD8\tD9\tD10\n")
	for i in range(905):
		o.write(str(i+1)+"\t"+L1sequence[i]+"\t"+str(result[i]["A"])+"\t"+str(result[i]["C"])+"\t"+str(result[i]["G"])+"\t"+str(result[i]["T"]))
		for j in range(1,11):
			o.write("\t"+str(result[i]["D"+str(j)]))
		o.write("\n")
	o.close()

	return result

def combineProfiles(library, reference):

	L1sequence = ""
	f = open(reference)
	for thisline in f:
		if (thisline[0] != ">"):
			L1sequence = L1sequence + str.strip(thisline)
	f.close()

	resultRNA = generateProfile(library+"_RNA.L1_5UTR_aligned.sam", sys.argv[2], library+"_RNA.profile.mask.txt")
	resultDNA = generateProfile(library+"_DNA.L1_5UTR_aligned.sam", sys.argv[2], library+"_DNA.profile.mask.txt")

	sumMismatchRNA = {}
	sumMismatchDNA = {}
	bases = "ACGT"
	sumD1RNA = 0
	sumD1DNA = 0
	sumD5RNA = 0
	sumD5DNA = 0
	sumD10RNA = 0
	sumD10DNA = 0
	sumDXRNA = 0
	sumDXDNA = 0

	##count total times each possible substitution or deletion is observed
	for i in range(len(resultRNA)):
		if (resultRNA[i]["A"] != "NA"):
			for j in range(len(bases)):
				if (bases[j] != L1sequence[i]):
					try:
						sumMismatchRNA[bases[j]+L1sequence[i]] = sumMismatchRNA[bases[j]+L1sequence[i]] + resultRNA[i][bases[j]]
					except:
						sumMismatchRNA[bases[j]+L1sequence[i]] = resultRNA[i][bases[j]]
					try:
						sumMismatchDNA[bases[j]+L1sequence[i]] = sumMismatchDNA[bases[j]+L1sequence[i]] + resultDNA[i][bases[j]]
					except:
						sumMismatchDNA[bases[j]+L1sequence[i]] = resultDNA[i][bases[j]]
		if (resultRNA[i]["D1"] != "NA"):
			sumD1RNA = sumD1RNA + resultRNA[i]["D1"]
			sumD1DNA = sumD1DNA + resultDNA[i]["D1"]

		if (resultRNA[i]["D2"] != "NA"):
			sumDXRNA = sumDXRNA + resultRNA[i]["D2"]
			sumDXDNA = sumDXDNA + resultDNA[i]["D2"]

		if (resultRNA[i]["D3"] != "NA"):
			sumDXRNA = sumDXRNA + resultRNA[i]["D3"]
			sumDXDNA = sumDXDNA + resultDNA[i]["D3"]

		if (resultRNA[i]["D4"] != "NA"):
			sumDXRNA = sumDXRNA + resultRNA[i]["D4"]
			sumDXDNA = sumDXDNA + resultDNA[i]["D4"]

		if (resultRNA[i]["D5"] != "NA"):
			sumD5RNA = sumD5RNA + resultRNA[i]["D5"]
			sumD5DNA = sumD5DNA + resultDNA[i]["D5"]

		if (resultRNA[i]["D6"] != "NA"):
			sumDXRNA = sumDXRNA + resultRNA[i]["D6"]
			sumDXDNA = sumDXDNA + resultDNA[i]["D6"]

		if (resultRNA[i]["D7"] != "NA"):
			sumDXRNA = sumDXRNA + resultRNA[i]["D7"]
			sumDXDNA = sumDXDNA + resultDNA[i]["D7"]

		if (resultRNA[i]["D8"] != "NA"):
			sumDXRNA = sumDXRNA + resultRNA[i]["D8"]
			sumDXDNA = sumDXDNA + resultDNA[i]["D8"]

		if (resultRNA[i]["D9"] != "NA"):
			sumDXRNA = sumDXRNA + resultRNA[i]["D9"]
			sumDXDNA = sumDXDNA + resultDNA[i]["D9"]

		if (resultRNA[i]["D10"] != "NA"):
			sumD10RNA = sumD10RNA + resultRNA[i]["D10"]
			sumD10DNA = sumD10DNA + resultDNA[i]["D10"]

	##calculate the aggregate ratio of normalised ratios between DNA and RNA libraries
	aggregateRatios = []
	for i in range(len(resultRNA)):
		currentRNADX = 0
		currentDNADX = 0
		ratios = []
		if ((resultRNA[i]["A"] != "NA") and (L1sequence[i] != "A")):
			ratios.append((resultRNA[i]["A"]/sumMismatchRNA["A"+L1sequence[i]])/(resultDNA[i]["A"]/sumMismatchDNA["A"+L1sequence[i]]))
		if ((resultRNA[i]["C"] != "NA") and (L1sequence[i] != "C")):
			ratios.append((resultRNA[i]["C"]/sumMismatchRNA["C"+L1sequence[i]])/(resultDNA[i]["C"]/sumMismatchDNA["C"+L1sequence[i]]))
		if ((resultRNA[i]["G"] != "NA") and (L1sequence[i] != "G")):
			ratios.append((resultRNA[i]["G"]/sumMismatchRNA["G"+L1sequence[i]])/(resultDNA[i]["G"]/sumMismatchDNA["G"+L1sequence[i]]))
		if ((resultRNA[i]["T"] != "NA") and (L1sequence[i] != "T")):
			ratios.append((resultRNA[i]["T"]/sumMismatchRNA["T"+L1sequence[i]])/(resultDNA[i]["T"]/sumMismatchDNA["T"+L1sequence[i]]))
		if (resultRNA[i]["D1"] != "NA"):
			ratios.append(float(resultRNA[i]["D1"]/float(sumD1RNA))/float(resultDNA[i]["D1"]/float(sumD1DNA)))
		if (resultRNA[i]["D5"] != "NA"):
			ratios.append(float(resultRNA[i]["D5"]/float(sumD5RNA))/float(resultDNA[i]["D5"]/float(sumD5DNA)))
		if (resultRNA[i]["D10"] != "NA"):
			ratios.append(float(resultRNA[i]["D10"]/float(sumD10RNA))/float(resultDNA[i]["D10"]/float(sumD10DNA)))
		if (resultRNA[i]["D2"] != "NA"):
			currentRNADX = currentRNADX + resultRNA[i]["D2"]
			currentDNADX = currentDNADX + resultDNA[i]["D2"]
		if (resultRNA[i]["D3"] != "NA"):
			currentRNADX = currentRNADX + resultRNA[i]["D3"]
			currentDNADX = currentDNADX + resultDNA[i]["D3"]
		if (resultRNA[i]["D4"] != "NA"):
			currentRNADX = currentRNADX + resultRNA[i]["D4"]
			currentDNADX = currentDNADX + resultDNA[i]["D4"]
		if (resultRNA[i]["D6"] != "NA"):
			currentRNADX = currentRNADX + resultRNA[i]["D6"]
			currentDNADX = currentDNADX + resultDNA[i]["D6"]
		if (resultRNA[i]["D7"] != "NA"):
			currentRNADX = currentRNADX + resultRNA[i]["D7"]
			currentDNADX = currentDNADX + resultDNA[i]["D7"]
		if (resultRNA[i]["D8"] != "NA"):
			currentRNADX = currentRNADX + resultRNA[i]["D8"]
			currentDNADX = currentDNADX + resultDNA[i]["D8"]
		if (resultRNA[i]["D9"] != "NA"):
			currentRNADX = currentRNADX + resultRNA[i]["D9"]
			currentDNADX = currentDNADX + resultDNA[i]["D9"]
		if ((currentRNADX != 0) and (currentDNADX != 0)):
			ratios.append(float(currentRNADX/float(sumDXRNA))/float(currentDNADX/float(sumDXDNA)))
		ratios.sort()
		if (len(ratios) >= 3):
			aggregateRatios.append(statistics.median(ratios))
		else:
			aggregateRatios.append("NA")

	##normalise to fragment median
	fragments = []
	fragments.append([0,218,[],0])
	fragments.append([219,443,[],0])
	fragments.append([444,668,[],0])
	fragments.append([669,904,[],0])

	for i in range(len(aggregateRatios)):
		for j in range(len(fragments)):
			if ((i >= fragments[j][0]) and (i <= fragments[j][1])):
				if (aggregateRatios[i] != "NA"):
					fragments[j][2].append(aggregateRatios[i])
	for i in range(len(fragments)):
		current = fragments[i][2]
		current.sort()
		fragments[i][3] = statistics.median(current)

	allRatios = []
	for i in range(len(aggregateRatios)):
		for j in range(len(fragments)):
			if ((i >= fragments[j][0]) and (i <= fragments[j][1])):
				if (aggregateRatios[i] != "NA"):
					aggregateRatios[i] = float(aggregateRatios[i])/(float(fragments[j][3]))
					allRatios.append(aggregateRatios[i])

	allRatios.sort()
	o = open(library+"_ratios.txt", "w")
	o.write("position\treference\taggregate_ratio\tsliding_5bp_window_aggregate_ratio\n")
	for i in range(len(aggregateRatios)):
		if (aggregateRatios[i] != "NA"):
			proceed = "Y"
			movingAverage = 0
			if ((i >= 2) and (i <= 902)):
				for j in range(i-2,i+3):
					if (aggregateRatios[j] == "NA"):
						proceed = "N"
					else:
						movingAverage = movingAverage + aggregateRatios[j]/statistics.median(allRatios)
				movingAverage = movingAverage / 5
			if (proceed == "Y"):
				o.write(str(i+1)+"\t"+L1sequence[i]+"\t"+str(aggregateRatios[i]/statistics.median(allRatios))+"\t"+str(movingAverage)+"\n")
			else:
				o.write(str(i+1)+"\t"+L1sequence[i]+"\t"+str(aggregateRatios[i]/statistics.median(allRatios))+"\tNA\n")
		else:
			o.write(str(i+1)+"\t"+L1sequence[i]+"\tNA\tNA\n")
	o.close()

combineProfiles(sys.argv[1], sys.argv[2])
sys.exit()

