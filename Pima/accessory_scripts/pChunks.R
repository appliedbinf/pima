#!/bin/env Rscript

library(parallel)
library(hash)
library(stringr)
library(grid)
library(gridExtra)
library(optparse)

options(width = 180)

printif = function(string = NULL, condition){
    if (condition) {
        print(string)
    }
}

findPlasmids = function(plasmidPSLFile = NULL, plasmidDatabase = NULL,
    amrPSLFile = NULL, amrDatabase, noAMR = FALSE,
    incPSLFile = NULL, incDatabase, noInc = FALSE,
    outputDirectory = NA, overwrite = TRUE,
    maxTargetLength = 300000,
    minQueryLength = 500,
    makeCircos = FALSE,
    minQueryCoverage = 1/2, minTargetCoverage = 1/2,
    searchDepth = NULL,
    verbosity = 0) {

    ## Verify the arguments
    argumentsGood = TRUE
    if (!file.exists(plasmidPSLFile)) {
        argumentsGood = FALSE
        message(paste('Plasmid PSL file', plasmidPSLFile, 'not found'))
    }
    if (!file.exists(plasmidDatabase)) {
        argumentsGood = FALSE
        message(paste('Plasmid database', plasmidDatabase, 'not found'))
    }
    if (is.na(outputDirectory)) {
        argumentsGood = FALSE
        message('Output directory not given')
    }
    if (file.exists(outputDirectory) && !overwrite) {
        argumentsGood = FALSE
        message(paste('Output directory', outputDirectory, 'already exists.  Add overwrite = TRUE'))
    }
    if (minQueryCoverage < .1 || minQueryCoverage > 1) {
        argumentsGood = FALSE
        message(paste('Minimum query coverage', minQueryCoverage, 'is outside of the range 0.1 <= x <= 1'))
    }
    if (minTargetCoverage < 0.02 || minTargetCoverage > 1) {
        argumentsGood = FALSE
        message(paste('Minimum target coverage', minTargetCoverage, 'is outside of the range 0.1 <= x <= 1'))
    }
    if (!argumentsGood){
        message('There is a problem with the arguments')
        return()
    }
    
    printif(paste('Finding plasmids in', plasmidPSLFile), verbosity > 0)

    ## Keep track of the total score in case we doing a grid search
    totalPlasmidScore = 0
    
    ## Check for the existence of the output directory, remove if it exists
    if (file.exists(outputDirectory)) {
        printif(paste('Removing existing output directory', outputDirectory), verbosity > 1)
        unlink(outputDirectory, recursive = TRUE)
    }
    printif(paste('Making output directory', outputDirectory), verbosity > 1)
    dir.create(outputDirectory)
    outputPrefix = paste0(outputDirectory, "/plasmids")
    
    ## Read in and filter the plasmid hits
    plasmidHits = read.table(plasmidPSLFile, row.names = NULL, header = FALSE, sep = '\t', stringsAsFactors = FALSE, skip = 5)
    colnames(plasmidHits) = c('match', 'mismatch', 'rep_m', 'Ns', 'tgap_c', 'tgap_b',
                'qgap_c', 'qgap_b', 'strand',
                'target', 'tlength', 'tstart', 'tstop',
                'query', 'qlength', 'qstart', 'qstop',
                'blocks', 'block_sizes', 'tstarts', 'qstarts')
    printif(paste("Sequence-plasmid hits:", nrow(plasmidHits)), verbosity > 0)

    plasmidHits = plasmidHits[order(plasmidHits[,'target'], -plasmidHits[,'qlength']), ]

    ## Toss out any hits missing information
    plasmidHits = plasmidHits[complete.cases(plasmidHits),]
    
    ## Toss out very long plasmid sequences -- probably actually genome chunks labeled incorrectly
    veryLongHits = sum(plasmidHits[,'tlength'] >= maxTargetLength)
    printif(paste('Removing', veryLongHits, 'hits greater than', maxTargetLength), verbosity > 0)
    plasmidHits = plasmidHits[plasmidHits[,'tlength'] <= maxTargetLength, ]
    printif(paste("Sequence-plasmid hits after removing very long plasmids:", nrow(plasmidHits)), verbosity > 0)

    ## Toss out very short query sequences -- probably junk or repeats
    veryShortQuery = sum(plasmidHits[,'qlength'] >= minQueryLength)
    printif(paste('Removing', veryShortQuery, 'queries less than', minQueryLength), verbosity > 0)
    plasmidHits = plasmidHits[plasmidHits[,'qlength'] >= minQueryLength, ]
    printif(paste("Sequence-plasmid hits after removing very short queries:", nrow(plasmidHits)), verbosity > 0)

    ## Toss out sequece-plasmid pairs below the coverage cutoff
    sequenceMatches = aggregate(x = plasmidHits[,'match',drop = FALSE],
        by = list(plasmidHits[,'query'], plasmidHits[,'target']), FUN = sum)
    printif(head(sequenceMatches), verbosity > 1)
    printif(paste('Sequence-plasmid pair matches:', paste(dim(sequenceMatches), collapse = 'x')), verbosity > 1)

    sequenceLengths = aggregate(x = plasmidHits[,'qlength', drop = FALSE],
        by = list(plasmidHits[,'query'], plasmidHits[,'target']), FUN = max)
    printif(head(sequenceLengths), verbosity > 1)
    printif(paste('Sequence-plasmid pair lengths:', paste(dim(sequenceLengths), collapse = 'x')), verbosity > 1)
   
    matchingFractions = cbind(sequenceMatches[,c(1,2)], sequenceMatches[,3] / sequenceLengths[,3])
    colnames(matchingFractions) = c('query', 'target', 'fraction')
    printif(head(matchingFractions), verbosity > 1)
    printif(paste('Sequence-plasmid pair fractions:', paste(dim(matchingFractions), collapse = 'x')), verbosity > 1)
    
    matchingFractions = matchingFractions[matchingFractions[,'fraction'] >= minQueryCoverage,]
    printif(head(matchingFractions), verbosity > 1)
    printif(paste('Passing sequence-plasmid pair fractions:', paste(dim(matchingFractions), collapse = 'x')), verbosity > 1)
    
    aboveMinCoverage = apply(matchingFractions, 1, function(i){paste0(i['query'], '|', i['target'])})
    plasmidHits = plasmidHits[apply(plasmidHits, 1, function(i){paste0(i['query'], '|', i['target'])}) %in% aboveMinCoverage, ]
    printif(paste("Sequence-plasmid hits after removing low-coverage hits:", nrow(plasmidHits)), verbosity > 0)

    ## Toss out plasmid sequences below the coverage cutoff
    targetMatches = aggregate(x = plasmidHits[,'match',drop = FALSE],
        by = list(plasmidHits[,'target']), FUN = sum)
    printif(head(targetMatches), verbosity > 1)
    printif(paste('Plasmid matches:', paste(dim(targetMatches), collapse = 'x')), verbosity > 1)

    targetLengths = aggregate(x = plasmidHits[,'tlength', drop = FALSE],
        by = list(plasmidHits[,'target']), FUN = max)
    printif(head(targetLengths), verbosity > 1)
    printif(paste('Plasmid lengths:', paste(dim(targetLengths), collapse = 'x')), verbosity > 1)
   
    matchingFractions = cbind(targetMatches[,1], targetMatches[,2] / targetLengths[,2])
    colnames(matchingFractions) = c('target', 'fraction')
    printif(head(matchingFractions), verbosity > 1)
    printif(paste('Plasmid fractions:', paste(dim(matchingFractions), collapse = 'x')), verbosity > 1)
    
    matchingFractions = matchingFractions[matchingFractions[,'fraction'] >= minTargetCoverage,]
    printif(head(matchingFractions), verbosity > 1)
    printif(paste('Passing plasmid fractions:', paste(dim(matchingFractions), collapse = 'x')), verbosity > 1)
    
    aboveMinCoverage = matchingFractions[, 'target']
    plasmidHits = plasmidHits[plasmidHits[, 'target'] %in% aboveMinCoverage, ]
    printif(paste("Sequence-plasmid hits after removing low-coverage hits:", nrow(plasmidHits)), verbosity > 0)

    
    ## If we're out of sequece-plasmid hits, then stop here
    if (nrow(plasmidHits) == 0) {
        message(paste('Not hits found'))
        return
    }

    ## Find out how much of each query (contig) is covered by each target (plasmid).
    ## Query coverage is constant and does not change as we assign contigs to plasmids
    queryCoverage = hash()
    queryMismatches = hash()
    for (i in 1:nrow(plasmidHits)) {
        if (!(i %% 1000)) {
            printif(paste('Processing hit', i, '/', nrow(plasmidHits)), verbosity > 0)
        }
        
        query = plasmidHits[i,'query']
        target = plasmidHits[i, 'target']

        ## Represent each sequence-plasmid hit as a series of 0/1 vectors that
        if (!has.key(query, queryCoverage)) {
            queryCoverage[[query]] = hash()
            queryMismatches[[query]] = hash()
        }
        if (!has.key(target, queryCoverage[[query]])) {
            queryCoverage[[query]][[target]] = rep(0, times = plasmidHits[i, 'qlength'])
            queryMismatches[[query]][[target]] = 0
        }
        
        blockSizes = as.numeric(unlist(strsplit(x = plasmidHits[i,'block_sizes'], ',')))
        qBlockStarts = as.numeric(unlist(strsplit(x = plasmidHits[i,'qstarts'], ',')))
        
        for (j in 1:length(blockSizes)) {
            queryCoverage[[query]][[target]][qBlockStarts[j]:(qBlockStarts[j]+blockSizes[j])] = 1
        }
        queryMismatches[[query]][[target]] = queryMismatches[[query]][[target]] + plasmidHits[i,'mismatch'] 
    }

    
    ## Pull the full plasmid names from the blast database because BLAT/minimap2 doesn't report them, just the ID's
    targetIDs = plasmidHits[,'target']
    targetIDs = gsub("\\|$", "", targetIDs)
    targetIDs = gsub(".*(\\|.*)$", "\\1", targetIDs)
    noDotIDs = gsub("\\|", "", targetIDs)
    noDotIDs = gsub("(^H[^.]+).[0-9]+$", "\\1", noDotIDs)
    noDotIDs = cbind(noDotIDs)

    targetFile = paste0(outputDirectory, '/targets.tsv', sep = '')
    write.table(file = targetFile, x = noDotIDs, quote = FALSE, row.names = FALSE, col.names = FALSE)
    command = paste('blastdbcmd -db', plasmidDatabase,
        '-entry_batch', targetFile,
        '| grep ">"')
    targetNames = system(command, intern = TRUE)
    printif(paste('Found', length(targetNames), 'target names for', length(targetIDs), 'targets.'), verbosity > 0)
    
    targetNames = gsub('^>.*\\| ', '', targetNames)
    targetNames = gsub('^>[^ ]+', '', targetNames)
    plasmidHits = cbind(plasmidHits, targetIDs, targetNames)
    printif(paste('Named hits:', nrow(plasmidHits)), verbosity > 1)
    
    #Pull just the plasmids out of the larget set of hits, i.e, make sure it has the word 'plasmid' in the description.
    plasmidHits = plasmidHits[grep('plasmid|vector', plasmidHits[,'targetNames'], ignore.case = TRUE), ,drop = FALSE]
    plasmidHits = plasmidHits[!grepl('tig0000|unnamed', plasmidHits[,'targetNames'], ignore.case = TRUE), , drop = FALSE]
    printif(paste("Sequece-plasmid hits after screening by name:", paste(dim(plasmidHits), collapse = 'x')), verbosity > 1)

    ## Stop if there is nothing left
    if (is.null(plasmidHits)) {
        message('Not hits found')
        return()
    }
    if (nrow(plasmidHits) == 0) {
        message('Not hits found')
        return()
    }
                
    ## Clean up the plasmid names -- they look like crap by default.
    plasmidNames = plasmidHits[,'targetNames']
    plasmidNames = gsub(', comp.*', '', plasmidNames)
    plasmidNames = gsub(', contig.*', '', plasmidNames)
    plasmidNames = gsub(', partial.*', '', plasmidNames)
    plasmidNames = gsub('strain ', '', plasmidNames)
    plasmidNames = gsub('^ *', '', plasmidNames)
    plasmidNames = sub('^(cl\\|)(.*?) ', '', plasmidNames)
    plasmidNames = sub('subsp. (.*?) ', '', plasmidNames)
    plasmidNames = sub('serovar (.*?) ', '', plasmidNames)
    plasmidNames = sub('strain ', '', plasmidNames)
    plasmidNames = sub('plasmid$', '', plasmidNames)
    plasmidHits[,'targetNames'] = plasmidNames

    ## Just take the best hit for each plasmid, hence the head, 1 in the agg
    plasmidNames = aggregate(plasmidHits, by  = list(plasmidHits[,'query']), FUN = head, 1)
    plasmidNames = plasmidNames[, 'targetNames', drop = FALSE]

    ## Find the set of plasmid coverage hits for each itteration
    usedContigs = c()
    plasmidMismatches = c()

    ## Order hits by the plasmid ID and the query length
    plasmidHits = plasmidHits[order(plasmidHits[,'target'], -plasmidHits[,'qlength']), ]
    
    ## Iterate, finding plasmids until we run out of usable sequence-plasmid its
    plasmidNumber = 0
    plasmidResults = c()
    while (1) {

        ## Keep track of how many plasmids we have gone over
        plasmidNumber = plasmidNumber + 1
        
        printif(paste('Sequence-plasmid hits left:', nrow(plasmidHits)), verbosity > 1)
        contigToPlasmid = hash()
        plasmidToContig = hash()
        plasmidCoverage = hash()
        plasmidCoverageWithRepeats = hash()
        contigCoverage = hash()
    
        ##Find contig/plasmid plasmid/contig pairs
        if (is.null(plasmidHits)) {
            break
        }
        if (nrow(plasmidHits) == 0) {
            break
        }
        repLengths = c()


        ## Find the coverage of each plasmid in the possible set by the contigs
        for (i in 1:nrow(plasmidHits)) {
        
            query = plasmidHits[i,'query']
            target = plasmidHits[i,'target']
            matches = plasmidHits[i,'match']
            mismatches = plasmidHits[i,'mismatch']
            score = matches - mismatches
            queryLength = plasmidHits[i,'qlength']
            blockSizes = as.numeric(unlist(strsplit(plasmidHits[i, 'block_sizes'], ',')))
            queryStarts = as.numeric(unlist(strsplit(plasmidHits[i, 'qstarts'], ',')))
            targetStarts = as.numeric(unlist(strsplit(plasmidHits[i, 'tstarts'], ',')))

            ## Skip matches which have less than 50% of the bases from the contig on the plasmid -- probably not a good match
            if ((sum(queryCoverage[[query]][[target]]) - queryMismatches[[query]][[target]]) / queryLength <= minQueryCoverage)  {
                next
            }

            targetLength = plasmidHits[i, 'tlength']; targetStart = plasmidHits[i, 'tstart']; targetStop = plasmidHits[i, 'tstop']

            ## Relate this contig to this plasmid
            if (!has.key(query, contigToPlasmid)) {
                contigToPlasmid[[query]] = hash()
                contigCoverage[[query]] = hash()
            }
            if (!has.key(target, contigToPlasmid[[query]])) {
                contigToPlasmid[[query]][[target]] = score
                contigCoverage[[query]][[target]] = rep(0, queryLength)
            } else {
                contigToPlasmid[[query]][[target]] = contigToPlasmid[[query]][[target]] + score
            }
            
            ## Keep track of target(plasmid) coverage by the contigs
            if (!has.key(target, plasmidCoverage)) {
                plasmidCoverage[[target]] = rep(0, targetLength)
                plasmidCoverageWithRepeats[[target]] = rep(0, targetLength)
                plasmidToContig[[target]] = hash()
                plasmidMismatches[target] = 0
            }
            
            penalized = FALSE
            for (j in 1:length(blockSizes)) {

                ## Keep track of all contig alignments to this plasmid, even with repeats
                plasmidCoverageWithRepeats[[target]][targetStarts[j]:(targetStarts[j] + blockSizes[j])] = 1

                ## Skip if this region of the query sequence has already been assigned to this plasmid
                if (sum(contigCoverage[[query]][[target]][queryStarts[j]:(queryStarts[j] + blockSizes[j])] == 0) <= 50) {
                    printif(paste('Sequence', query, 'already used for', target,
                                  '. ', paste0(queryStarts[j], '-', queryStarts[j] + blockSizes[j])), verbosity > 2)
                    next
                }
                if (!penalized) {
                    ## Penalty for every gap, only penalize once per match
                    ## TODO: Penalize for gap length, not just once per gap
                    #plasmidMismatches[target] = plasmidMismatches[target] + (length(blockSizes) - 1) * 100
                    plasmidMismatches[target] = plasmidMismatches[target] + mismatches * 5
                    penalized = TRUE
                }

                plasmidCoverage[[target]][targetStarts[j]:(targetStarts[j] + blockSizes[j])][contigCoverage[[query]][[target]][queryStarts[j]:(queryStarts[j] + blockSizes[j])] == 0] =
                    plasmidCoverage[[target]][targetStarts[j]:(targetStarts[j] + blockSizes[j])][contigCoverage[[query]][[target]][queryStarts[j]:(queryStarts[j] + blockSizes[j])] == 0] + 1
                contigCoverage[[query]][[target]][queryStarts[j]:(queryStarts[j] + blockSizes[j])] = 1
            }

            ## Relate this plasmid to this contig
            if (!has.key(query, plasmidToContig[[target]])) {
                plasmidToContig[[target]][[query]] = score
            } else {
                plasmidToContig[[target]][[query]] = plasmidToContig[[target]][[query]] + score
            }
            if (target == 'NZ_GG692894.1' && query == 'contig_3_0') {
                print('NZ_GG692894.1')
                print(sum(plasmidCoverage[[target]]))
                print(sum(contigCoverage[[query]][[target]]))
            }

        }
        
        ## Get the best set of plasmids out, i.e., the set with the most bases matching between the contig and plasmid
        plasmidScores = c()
        for (thisPlasmid in keys(plasmidCoverage)){
            thisPlasmidScore = sum(plasmidCoverage[[thisPlasmid]])
            plasmidScores = c(plasmidScores, thisPlasmidScore)
            names(plasmidScores)[length(plasmidScores)] = thisPlasmid
        }
        plasmidScores = sort(plasmidScores - plasmidMismatches[names(plasmidScores)], dec = TRUE)

        if (length(plasmidScores) > 0) {
            printif('Highest scoring plasmids', verbosity > 1)
            printif(head(cbind(plasmidScores), 20), verbosity > 1)
        }
        
        ## Stop searching for plasmids if nothing matches well or we're out of hits
        if (length(plasmidScores) == 0) {
            printif('Out of plasmids', verbosity > 0)
            break
        } else if (max(plasmidScores) < 500) {
            printif('Out of min-scoring plasmids', verbosity > 0)
            break
        }
        
        ## For each matching plasmid, ordered by total bases matching the assembly, find the set of corresponding contigs
        plasmidToUse = 1
        if (!is.null(searchDepth) && plasmidNumber <= length(searchDepth)) {
            plasmidToUse = searchDepth[plasmidNumber]
        }
        if (plasmidToUse > length(plasmidScores)) {
            plasmidToUse = length(plasmidScores)
        }
                      
        plasmid = names(plasmidScores)[plasmidToUse]
        print(paste('Plasmid picked', plasmid))
        totalPlasmidScore = totalPlasmidScore + plasmidScores[plasmid]
        printif(paste("Pulling sequences for", plasmid), verbosity > 0)
        
        ## Find contigs what haven't already been given to another plasmid so we can assign them next round
        plasmidContigs = keys(plasmidToContig[[plasmid]])
        unusedContigs = plasmidContigs[!(plasmidContigs %in% usedContigs)]

        ## If no unused contigs that also map to this plasmid then skip it
        if (length(unusedContigs) == 0) {
            printif(paste("No unused sequences for", plasmid), verbosity > 1)
            next
        }

        ## Keep just the rows for this plasmid and which haven't already been used by another plasmid
        plasmidRows = plasmidHits[plasmidHits[,'target'] == plasmid,]
        plasmidRows = plasmidRows[plasmidRows[,'query'] %in% unusedContigs,,drop = FALSE]
        plasmidName = plasmidRows[1, 'targetNames']
        plasmidID = plasmidRows[1, 'targetIDs']
        ## Get the plasmid length
        command = paste('blastdbcmd -db', plasmidDatabase,
            '-entry', plasmidID,
            '-outfmt "%l"')
        plasmidLength = system(command, intern = TRUE)
        plasmidLength = rep(plasmidLength, nrow(plasmidRows))

        ## How many bases from the plasmid are uncovered?
        plasmidMissing = rep(sum(plasmidCoverageWithRepeats[[as.character(plasmidID)]] == 0), nrow(plasmidRows))

        ## Keep track of all of the contigs included
        usedContigs = c(usedContigs, unusedContigs)

        ## How many matching bases for each query sequence?
        thisPlasmidQuerySizes = c()
        thisPlasmidMatches = c()
        for (contig in unusedContigs) {
            thisPlasmidQuerySizes = c(thisPlasmidQuerySizes, length(queryCoverage[[contig]][[plasmid]]))
            thisPlasmidMatches = c(thisPlasmidMatches, sum(queryCoverage[[contig]][[plasmid]]))
        }
        
        ## Add this plasmid's hits onto the growing list of sequence-plasmid hits
        thisPlasmidResults = cbind(unusedContigs, plasmidName, as.character(plasmidID), thisPlasmidQuerySizes, thisPlasmidMatches, plasmidLength, plasmidMissing)
        colnames(thisPlasmidResults) = c('query.name', 'plasmid.name', 'plasmid.accession', 'query.size', 'aligned.bases', 'plasmid.size', 'plasmid.missing')
        thisPlasmidResults = thisPlasmidResults[order(-thisPlasmidMatches), ]
        plasmidResults = rbind(plasmidResults, thisPlasmidResults)
        
        ## Remove the contigs added to this plasmid from the list of plasmid/contig BLAT hits
        plasmidHits = plasmidHits[!(plasmidHits[,'query'] %in% usedContigs),]
        plasmidHits = plasmidHits[!(plasmidHits[,'target'] == plasmid),]
    }

    rownames(plasmidResults) = plasmidResults[,1]
    
    ## Check for the presence of AMR genes in this file
    if (!noAMR) {
        amrBEDFile = paste0(outputDirectory, '/amrMapping.bed')
        command = paste('cat', amrPSLFile,
            '| awk -F \'\\t\' \'($3 >= 80) && ($4 / $14 >= .95){OFS = "\t"; print $2,($9 < $10 ? $9 : $10),($9 < $10 ? $10 : $9),$1,$3/100,($9 < $10 ? "+" : "-")}\'',
            '| sort -k 1,1 -k 2,2n >', amrBEDFile)
        printif(command, verbosity > 1)
        system(command)
        ## Find local overlapping regions
        amrMergedBEDFile = paste0(outputDirectory, '/amrMergedMapping.bed')
        command = paste('bedtools merge -d -30 -i', amrBEDFile, '>', amrMergedBEDFile)
        printif(command, verbosity > 1)
        system(command)
        ## Find the best AMR gene for each region
        amrFinalBEDFile = paste0(outputDirectory, '/amrFinal.bed')
        command = paste('bedtools intersect',
            '-a', amrBEDFile,
            '-b', amrMergedBEDFile,
            '-f .9 -F .9',
            '-wao',
            '| awk \'$7 != "."\'',
            '| awk \'{OFS="\t";locus=$7"\t"$8"\t"$9; if($5 > s[locus]){s[locus]=$5;b[locus] = $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}END{for(i in b){print i,b[i]}}\'',
            '>', amrFinalBEDFile)
        printif(command, verbosity > 1)
        system(command)

        ## Read the AMR results in and add them to the plasmid contigs
        amrResults = read.table(file = amrFinalBEDFile, header = FALSE, row.names = NULL, stringsAsFactors = FALSE, quote = '')
        amrResults[,7] = gsub('(_.*$)|(.*\\|)', '', amrResults[,7])
        amrResults = aggregate(amrResults[ , 7, drop = FALSE], by = list(amrResults[,1]), function(i){paste(i, collapse = ', ')})
        rownames(amrResults) = amrResults[,1]
        amrResults = amrResults[ , 2, drop = FALSE]
        print(amrResults)
        
        plasmidResults = cbind(plasmidResults, rep('', nrow(plasmidResults)))
        colnames(plasmidResults)[ncol(plasmidResults)] = 'amr'
        plasmidResults[rownames(plasmidResults) %in% rownames(amrResults), 'amr'] =
            amrResults[rownames(plasmidResults)[rownames(plasmidResults) %in% rownames(amrResults)], 1]
        
    }
    
    ## Check for the presence of incompatibility groups in this file
    if (!noInc) {
        incBEDFile = paste0(outputDirectory, '/incMapping.bed')
        command = paste('cat', incPSLFile,
            '| awk -F \'\\t\' \'($3 >= 80) && ($4 / $14 >= .95){OFS = "\t"; print $2,($9 < $10 ? $9 : $10),($9 < $10 ? $10 : $9),$1,$3/100,($9 < $10 ? "+" : "-")}\'',
            '| sort -k 1,1 -k 2,2n >', incBEDFile)
        printif(command, verbosity > 1)
        system(command)
        ## Find local overlapping regions
        incMergedBEDFile = paste0(outputDirectory, '/incMergedMapping.bed')
        command = paste('bedtools merge -d -30 -i', incBEDFile, '>', incMergedBEDFile)
        printif(command, verbosity > 1)
        system(command)
        ## Find the best INC group for each region
        incFinalBEDFile = paste0(outputDirectory, '/incFinal.bed')
        command = paste('bedtools intersect',
            '-a', incBEDFile,
            '-b', incMergedBEDFile,
            '-f .9 -F .9',
            '-wao',
            '| awk \'$7 != "."\'',
            '| awk \'{OFS="\t";locus=$7"\t"$8"\t"$9; if($5 > s[locus]){s[locus]=$5;b[locus] = $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}END{for(i in b){print i,b[i]}}\'',
            '>', incFinalBEDFile)
        printif(command, verbosity > 1)
        system(command) 

        ## Read the inc group results in and add them to the plasmid contigs
        incResults = read.table(file = incFinalBEDFile, header = FALSE, row.names = NULL, stringsAsFactors = FALSE, quote = '')
        incResults[,7] = gsub('(_.*$)|(.*\\|)', '', incResults[,7])
        incResults = aggregate(incResults[ , 7, drop = FALSE], by = list(incResults[,1]), function(i){paste(i, collapse = ', ')})
        rownames(incResults) = incResults[,1]
        incResults = incResults[ , 2, drop = FALSE]
        print(incResults)
        
        plasmidResults = cbind(plasmidResults, rep('', nrow(plasmidResults)))
        colnames(plasmidResults)[ncol(plasmidResults)] = 'inc'
        plasmidResults[rownames(plasmidResults) %in% rownames(incResults), 'inc'] =
            incResults[rownames(plasmidResults)[rownames(plasmidResults) %in% rownames(incResults)], 1]
    }
    
    ## Write the plasmid results to file
    plasmidChunksFile = paste0(outputDirectory, '/plasmids.tsv')
    write.table(file = plasmidChunksFile, x = plasmidResults, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
    
    ## Dump a sequence file of potential plasmid contigs
    plasmidSequenceFile = paste0(outputPrefix, '.fna')
    system(paste0('echo "" >', plasmidSequenceFile))
    for (contig in plasmidResults[,'query.name']) {
        command = paste('faidx', paste0('assembly.fasta'), contig, '>>', plasmidSequenceFile)
        print(command)
        #system(command)
    }


    ## Return the total score of this round, in case we are doing a search
    return(totalPlasmidScore)
        
}

pChunks = function(plasmidPSLFile = NULL, plasmidDatabase = NULL,
    amrPSLFile = NULL, amrDatabase, noAMR = FALSE,
    incPSLFile = NULL, incDatabase, noInc = FALSE,
    outputDirectory = NA, overwrite = TRUE,
    maxTargetLength = 300000,
    minQueryLength = 200,
    makeCircos = FALSE,
    minQueryCoverage = 1/2, minTargetCoverage = 1/2,
    searchDepth = c(1),
    threads = 1,
    verbosity = 2) {
    
    print(plasmidDatabase)
    
    ## Verify the arguments
    argumentsGood = TRUE
    if (!file.exists(plasmidPSLFile)) {
        argumentsGood = FALSE
        message(paste('Plasmid PSL file', plasmidPSLFile, 'not found'))
    }
    if (!file.exists(plasmidDatabase)) {
        argumentsGood = FALSE
        message(paste('Plasmid database', plasmidDatabase, 'not found'))
    }
    if (is.na(outputDirectory)) {
        argumentsGood = FALSE
        message('Output directory not given')
    }
    if (file.exists(outputDirectory) && !overwrite) {
        argumentsGood = FALSE
        message(paste('Output directory', outputDirectory, 'already exists.  Add overwrite = TRUE'))
    }
    if (minQueryCoverage < .1 || minQueryCoverage > 1) {
        argumentsGood = FALSE
        message(paste('Minimum query coverage', minQueryCoverage, 'is outside of the range 0.1 <= x <= 1'))
    }
    if (minTargetCoverage < 0.02 || minTargetCoverage > 1) {
        argumentsGood = FALSE
        message(paste('Minimum target coverage', minTargetCoverage, 'is outside of the range 0.1 <= x <= 1'))
    }
    if (!argumentsGood){
        message('There is a problem with the arguments')
        return()
    }

    ## Check for the existence of the output directory, remove if it exists
    if (file.exists(outputDirectory)) {
        printif(paste('Removing existing output directory', outputDirectory), verbosity > 1)
        unlink(outputDirectory, recursive = TRUE)
    }
    printif(paste('Making output directory', outputDirectory), verbosity > 1)
    dir.create(outputDirectory)
    
    ## Default to c(1) for the plasmid search depth
    searchDepths = lapply(searchDepth, function(i){seq(i, 1)})
    searchDepths = as.matrix(expand.grid(searchDepths))
    print(searchDepths)
    
    plasmidScores = mclapply(1:nrow(searchDepths), function(i) {
                               findPlasmids(plasmidPSLFile = plasmidPSLFile, plasmidDatabase = plasmidDatabase,
                                            amrPSLFile = amrPSLFile, amrDatabase, noAMR = noAMR,
                                            incPSLFile = incPSLFile, incDatabase, noInc = noInc,
                                            outputDirectory = paste0(outputDirectory, '/plasmids_', paste(searchDepths[i,], collapse = '_')), overwrite = overwrite,
                                            maxTargetLength = 300000,
                                            minQueryLength = 200,
                                            makeCircos = makeCircos,
                                            minQueryCoverage = 1/2, minTargetCoverage = 1/2,
                                            searchDepth = searchDepths[i,],  ## Search depth i
                                            verbosity = verbosity)
                           },
        mc.cores = threads)
    plasmidScores = unlist(plasmidScores)

    print(cbind(searchDepths, plasmidScores))

    ## Pick out the best set, penalizing for not taking the first
    penalties = (unlist(apply(searchDepths, 1, sum)) - ncol(searchDepths)) * 2500
    print(penalties)
    plasmidScores = plasmidScores - penalties
    print(cbind(searchDepths, plasmidScores))
    bestScoring = which.max(plasmidScores)
    bestScoringDirectory = paste0(outputDirectory, '/plasmids_', paste(searchDepths[bestScoring,], collapse = '_'))

    ## Link to the best scoring files
    files = c('plasmids.tsv', 'amrFinal.bed')
    commands = paste('ln -s ', paste0('"', bestScoringDirectory, '/', files, '"'),
        paste0(outputDirectory, '/', files))
    printif(commands, verbosity >= 2)
    lapply(commands, system)
    
}

    
optionList = list(
    make_option('--plasmid-psl', type = 'character', default = NULL,
                help = 'Plasmid PSL database-v-contig output', metavar = '<PSL_FILE>'),
    make_option('--plasmid-database', type = 'character', default = NULL,
                help = 'Plasmid database', metavar = '<PLASMID_FASTA>'),
    make_option('--amr-database', type = 'character', default = NULL,
                help = 'AMR database', metavar = '<AMR_FASTA>'),
    make_option('--amr-blast', type = 'character', default = NULL,
                help = 'AMR blast output', metavar = '<BLAST_6>'),
    make_option('--output', type = 'character', default = NULL,
                help = 'Output dir', metavar = '<OUTPUT_DIR>'),
    make_option('--threads', type = 'numeric', default = NULL,
                help = 'Output dir', metavar = '<OUTPUT_DIR>'),
    make_option('--no-amr', action = 'store_true', default = FALSE,
                help = 'Don\'t run AMR'),
    make_option('--no-inc', action = 'store_true', default = FALSE,
                help = 'Don\'t run incompatibility groups')
    )

optParser = OptionParser(option_list = optionList)
opt = parse_args(optParser)
print(opt)
pChunks(plasmidPSLFile = opt$'plasmid-psl', plasmidDatabase = opt$'plasmid-database',
        amrPSLFile = opt$'amr-blast', amrDatabase = opt$'amr-database',
        outputDirectory = opt$output,
        threads = opt$threads,
#        searchDepth = c(1,1),
        searchDepth = c(5,5),
        noAMR = TRUE, noInc = TRUE,
        verbosity = 2)



