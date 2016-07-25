#################################################################################
#																				#
#	Process multiple react2 datasheets from multiple ReACT runs using R.		#
#	Given a folder with react2 files this will output the total interactions	#
#	found, the unique protein interactions, unique peptide interaction and		#
#	the unique site interactions as well as providing tables representing		#
#	these interactions.															#
#																				#
#	To run: copy all of this script into R, then run command:					#
#	Run React2Processing()														#
#	Then select the appropriate directory.										#
#																				#
#	by Devin Schweppe			last updated: 02/23/2016						#
#																				#
#################################################################################

#open files that will be assimilated into single table
#SET WORKING DIRECTORY WITH THOSE FILES YOU WANT PROCESSESED. 
#02232016: Added FDR determination if react2csv included "-F" option

React2Processing = function(reactdir) {
	
	library(tcltk)
	library(plyr)
#variables
comr2 = NULL

i = NULL
z = NULL
x = NULL
y = NULL

reactdir = NULL
reactdir1 = NULL
r2files = NULL
filex = NULL

upepr2p = NULL
uprotr2p = NULL
usiter2p = NULL

upaccs.site1 = NULL
upaccs.site2 = NULL
siter2p1 = NULL
siter2p2 = NULL

totalft = NULL
totalProt = NULL
protft = NULL
pepft = NULL
siteft = NULL
final.table = NULL


	#Make write.r2p.table true if you want a txt output file of the combined react data (comr2) in the wd.
	
	write.r2p.table = winDialog(type = c("yesnocancel"),"Write Combined React2 Table?")

	reactdir1 = tk_choose.dir(getwd(), "Who? What? Where? When? Give me a ReACT-ion")
	setwd(reactdir1)
	
	#make table containing react2 files from working directory set above
	r2files = list.files(pattern = "\\.react2.xls$")
	r1files = list.files(pattern = "\\.react1.xls$")
	
	comr2 = data.frame()
	
	#combined react2 values from all of these files.
	for (currentFile in r2files){
		filex = read.delim(currentFile,sep="\t",header=T,fill=T)
		filex$File.name = rep(currentFile, nrow(filex))
		if(nrow(comr2)>0){
			#comr2 = merge(comr2,filex,all=T)
			comr2 = rbind.fill(comr2,filex)
		}
		else {
			comr2 = filex
		}
		rm(filex)
	}
	
	#find and remove duplicate peptide-peptide interaction identifications.
	upepr2p = comr2[!duplicated(comr2[c('ms3.pep1','ms3.pep2')]),]

	#Find  Uniprot Accession number between pipes and and site of modification then output the results
	upaccs1=sapply(comr2[c('ms3.prot1')],function(y) {
				substring(
					y,
					first = 4,
					last = gregexpr("\\Q|\\E.",y)[[1]][c(2)]-1
				)
		}
	)
	#bind protein 1 uniport accessions and sites into small table.
	upaccs.site1 = cbind(upaccs1,comr2[c('mod_pos1')])

	#Find  Uniprot Accession number between pipes and and site of modification then output the results
	upaccs2=sapply(comr2[c('ms3.prot2')],function(y) {
			substring(
				y,
				first = 4,
				last = gregexpr("\\Q|\\E.",y)[[1]][c(2)]-1
			)
		}
	)
	
	#find reverse hits
	rev1 = sapply(comr2[c('ms3.prot1')],function(y) {
			ifelse(grepl("rev_.",y),1,0)
		}
	)
	rev2 = sapply(comr2[c('ms3.prot2')],function(y) {
			ifelse(grepl("rev_.",y),1,0)
		}
	)
	
	#determine false discovery rate
	rev.sum = cbind(rev1,rev2,rev1+rev2)
	rev.count = ifelse(rev.sum>0,1,0)
	
	#bind protein 2 uniport accessions and sites into small table.
	upaccs.site2 = cbind(upaccs2,comr2[c('mod_pos2')])

	siter2p1 = do.call(paste,c(upaccs.site1[c(1)],upaccs.site1[c(2)],sep="_"))
	siter2p2 = do.call(paste,c(upaccs.site2[c(1)],upaccs.site2[c(2)],sep="_"))

	prositesr2p = as.data.frame(cbind(siter2p1,siter2p2))

	#find and remove duplicate site-site interaction identifications.
	usiter2p = prositesr2p[!duplicated(prositesr2p[c(1,2)]),]

	#make uniprot accession table for PPI duplicate removal
	protaccsr2p = as.data.frame(cbind(upaccs1,upaccs2))
	
	#find and remove duplicate protein-protein interaction identifications.
	uprotr2p = protaccsr2p[!duplicated(protaccsr2p[c(1,2)]),]
	
	#find total unique proteins
	prottotal = rbind(upaccs1, upaccs2)
	uprottotal = unique(prottotal, incomparables = F)
	
	#Make summary table
	totalft = nrow(comr2)
	totalProt = nrow(uprottotal)
	protft = nrow(uprotr2p)
	pepft = nrow(upepr2p)
	siteft = nrow(usiter2p)
	revft = sum(rev.count)
	fdrft = sum(rev.count)/totalft
	final.table = cbind(
		c(
			"Total # of Identified Interactions",
			"# Unique Protein Accessions",
			"# of Unique PPIs","# of Unique Peptide Interactions",
			"# of Unique Site Interactions",
			"# Reverse Relationships (T_D, D_T, D_D)",
			"False Discovery Rate"
		),
		as.vector(
			c(
				totalft,
				totalProt,
				protft,
				pepft,
				siteft,
				revft,
				fdrft
			)
		)
	)

	#writes comr2 txt file to wd
	if(write.r2p.table == "YES") {
		write.table(comr2,file = paste(reactdir1,"/Combined_React2",format(Sys.time(),"%Y%m%d%H%M%S"),".txt", sep=""), sep="\t",row.names=F)
		}
		
	print(final.table)
	#print(head(upaccs1))
}
