#Snakemake rules van Tijs van Lieshout en Brecht van den Berg

#Lezen van een bestand met gene id’s.
rule preprocessing:
	input:
		"data/RNA-Seq-counts.txt"
	output:
		"data/geneIDs.txt"
		
	run:
		from itertools import islice
		
		with open(input[0], "r") as f, open(output[0], 'w') as outfile:
				for line in islice(f, 2, None):
					lineLijst = line.split()
					ID = lineLijst[0]
					outfile.write(ID+'\n')

#Herkennen van het type gene id en converteren naar handzame gene identifiers.
rule convertGeneIDs:
	input:
		"data/geneIDs.txt"
	output:
		"data/convertedGeneIDs.txt"
	run:
		from Bio import Entrez
		
		Entrez.email = 'brechtvandenberg20@gmail.com'
		
		with open(input[0], "r") as f, open(output[0], 'w') as outfile:
				for line in f:
					lineLijst = line.split('\n')
					ID = lineLijst[0]
					handle = Entrez.esearch(db="protein", term=ID)
					record = Entrez.read(handle)
					handle.close()
					
					searchID = (record['IdList'])
					searchID = str(searchID[0])
					outfile.write(searchID+'\n')

#De sequentie van het gen
rule getSeq:
	input:
		"data/convertedGeneIDs.txt"
	output:
		"data/seqs.fasta"
	run:
		from Bio import Entrez
		
		Entrez.email = 'brechtvandenberg20@gmail.com'
		
		with open(input[0], "r") as f, open(output[0], 'w') as outfile:
			for line in f:
				lineLijst = line.split('\n')
				ID = lineLijst[0]
				handle = Entrez.efetch(db="nuccore", id=ID, rettype="fasta", retmode="text")
				record = handle.read()
				handle.close()
				
				outfile.write(record+'\n')

#Functionaliteit van het gen
rule getFunction:
	input:
		"data/convertedGeneIDs.txt"
	output:
		"data/functions.txt"
	run:
		from Bio import Entrez
		import re
		
		Entrez.email = 'brechtvandenberg20@gmail.com'
		
		with open(input[0], "r") as f, open(output[0], 'w') as outfile:
			for line in f:
				lineLijst = line.split('\n')
				ID = lineLijst[0]
				handle = Entrez.efetch(db="protein", id=ID, rettype="ft", retmode="text")
				record = handle.read()
				handle.close()
				
				function = re.search('note.*',record).group(0)
				function = function.strip('note\t')
				outfile.write(str(function)+'\n')

#Tussen stap voor het verkrijgen van PubMed ID's
rule getRecords:
	input:
		"data/seqs.fasta",
		"data/convertedGeneIDs.txt"
	output:
		"data/records.txt"
	run:
		from Bio import Entrez
		
		Entrez.email = 'brechtvandenberg20@gmail.com'
		recordList = []
		IDList = []
		counter = 0
		
		with open(input[0], "r") as fasta, open(input[1], "r") as f, open(output[0], 'w') as outfile:
			for line in f:
				lineLijst = line.split('\n')
				ID = lineLijst[0]
				IDList.append(ID)
			for line in fasta:
				if ">" in line:
					if "=" in line:
						lineLijst = line.split("=")
						ID = lineLijst[1]
						handle = Entrez.esearch(db="pubmed", term=ID, rettype="uilist", retmode="text")
						record = handle.read()
						handle.close()
					else:
						record = ""
					recordList.append(record)
					
			for i in IDList:
				outfile.write("Gene: "+IDList[counter]+'\n'+recordList[counter]+'\n'+"-")
				counter+=1

# Haal de PubMed ID's voor alle genen op
rule getPubMed:
	input:
		"data/records.txt"
	output:
		"data/pubMedIDs.txt"
	run:
		with open(input[0], "r") as f, open(output[0], 'w') as outfile:
			for line in f:
				lineLijst = line.split('\n')
				if "<" in line:
					if "<Id>" in line:
						line = line.strip("<Id>")
						line = line[:-6]
						outfile.write(line+'\n')
				else:
					outfile.write(line)

#Zet de genen in volgorde van het aantal PubMed ID's
rule sortPubMed:
	input:
		"data/pubMedIDs.txt"
	output:
		"data/sortedPubMedIDs.txt"
	run:
		counter = 0
		
		with open(input[0], "r") as f, open(output[0], 'w') as outfile:
			content = f.read()
			pmLijst = content.split("-")
			pmLijst.sort(key = len, reverse=True)
			
			for i in pmLijst:
				
				outfile.write(pmLijst[counter])
				counter+=1

#Bepaal per gen in welke pathways deze actief is (KEGG id)
rule getPathway:
	input:
		"data/geneIDs.txt"
	output:
		"data/pathways.txt"
	run:
		from bioservices import KEGG
		
		s = KEGG()
		IDList = []
		counter = 0
		
		with open(input[0], "r") as f, open(output[0], 'w') as outfile:
			for line in f:
				lineLijst = line.split('\n')
				ID = lineLijst[0]
				
				IDList.append(ID)
				pathway = s.get_pathway_by_gene(ID,"lpl")
				if pathway == None:
					pathway = ""
				outfile.write(IDList[counter]+'\t'+str(pathway)+'\n')
				counter+=1

rule results:
	input:
		"data/geneIDs.txt",
		"data/convertedGeneIDs.txt",
		"data/seqs.fasta",
		"data/functions.txt",
		"data/pubMedIDs.txt",
		"data/pathways.txt",
		"data/sortedPubMedIDs.txt"
	output:
		"results.txt"
	run:
		counter = 0
		counter2 = 0
		geneIDLijst = []
		convertedGeneIDLijst = []
		seqLijst = []
		functionLijst = []
		pubMedIDLijst = []
		pathwayLijst = []
		sortedPubMedIDLijst = []
		countPubMedIDLijst = []
		
		with open(input[0], "r") as f1, open(input[1], "r") as f2, open(input[2], "r") as f3, open(input[3], "r") as f4, open(input[4], "r") as f5, open(input [5], "r") as f6, open(input[6], "r") as f7, open(output[0], 'w') as outfile:
			outfile.write("Results of the workflow"+'\n'+'\n')
			
			outfile.write("Genen in volgorde van hoeveelheid voorkomen in artikelen (dalend van veel naar weinig):"+'\n')
			for line in f7:
				if "Gene:" in line:
					countPubMedIDLijst = []
					lineLijst = line.split('\n')
					sortedPubMedID = lineLijst[0]
					sortedPubMedIDLijst.append(sortedPubMedID)
					outfile.write(sortedPubMedIDLijst[counter2]+'\n')
					counter2+=1
			
			outfile.write('\n'+"Rapport per Gen:"+'\n'+'\n')
			for line in f1:
				lineLijst = line.split('\n')
				geneID = lineLijst[0]
				geneIDLijst.append(geneID)
			for line in f2:
				lineLijst = line.split('\n')
				convertedGeneID = lineLijst[0]
				convertedGeneIDLijst.append(convertedGeneID)
			content = f3.read()
			seqLijst = content.split(">")
			for line in f4:
				lineLijst = line.split('\n')
				function = lineLijst[0]
				functionLijst.append(function)
			for line in f5:
				lineLijst = line.split('\n')
				pubMedID = lineLijst[0]
				pubMedIDLijst.append(pubMedID)
			for line in f6:
				lineLijst = line.split('\n')
				pathway = lineLijst[0]
				pathwayLijst.append(pathway)
				
				outfile.write("Gene ID:"+'\t'+geneIDLijst[counter]+'\n'+"Converted Gene ID:"+'\t'+convertedGeneIDLijst[counter]+'\n'+"Sequence:"+'\t'+seqLijst[counter][:-6]+'\n'+"Function:"+'\t'+functionLijst[counter]+'\n'+"PubMed entries:"+'\t'+pubMedIDLijst[counter]+'\n'+"KEGG Pathways:"+'\t'+pathwayLijst[counter]+'\n'+'\n')
				counter+=1

rule visualize:
	output:
		"dag.pdf"
	shell:
		"snakemake --dag results | dot -Tpdf > dag.pdf"