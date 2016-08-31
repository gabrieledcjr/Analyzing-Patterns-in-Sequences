#!/usr/bin/python -tt
# python seq_identifier.py -p Allels_Trimv2.fa -s MS_1376Translations.fa -d
import os
import sys
import getopt


import GlobalVars

from GetSequences import GetSequences
from DeterminePattern import DeterminePattern
from AssembleFASTAFiles import WriteToFile, WriteToFileNoMatches, MakeDir
from xhtml2pdf import pisa


def usage():
  print "Usage:"
  print "  seq_identifier.py -p <patternFile> -s <sampleFile>"
  print "  seq_identifier.py --min-len=1 --max-gap=3 -p <patternFile> -s <sampleFile>"
  print "  seq_identifier.py [Other Options] [Required Arguments]"
  print
  print "Required Arguments:"
  print "  -p, --pattern=FILE    Pattern file location"
  print "  -s, --sample=FILE     Sample file location"
  print
  print "Optional Arguments:"
  print "  -h, --help            Help"
  print "  -d, --debug           Turn on debug mode"
  print "  --min-len=NUM         Minimum length for pattern match [Default 4]"
  print "  --min-gap=NUM         Minimum gap [Default 1]"
  print "  --max-gap=NUM         Maximum gap [Default 4]"
  print "  --out-pdf=[0 or 1]    Output pdf file [Default 1 (true)]"
  print "  --st-anchor=STRING    Starting anchor [Default KWG]"
  print "  --en-anchor=STRING    Ending anchor [Default GMA]"
  print


def CombineSamePattern(match):
  totalPatterns = len(match)
  if totalPatterns == 1:
    return match
  else:
    startIndex = 0
    while True:
      if match[startIndex][1].num == match[startIndex+1][1].num:
        match[startIndex][1].endIndex = match[startIndex+1][1].endIndex
        match[startIndex][1].size = match[startIndex][1].endIndex - match[startIndex][1].startIndex
        match[startIndex][0].snippets += match[startIndex+1][0].snippets
        match[startIndex][0].endIndex = match[startIndex+1][0].endIndex
        match[startIndex][0].size = match[startIndex][0].endIndex - match[startIndex][0].startIndex
        del match[startIndex+1]
        totalPatterns -= 1
      else:
        startIndex += 1
      if (startIndex+1) == totalPatterns:
        break
  return match


def main(argv):

  # CONFIGURATION VARIABLES
  noMatch = []
  results = {}
  results_ = {}
  minLen = 4        # minimum length of pattern match
  minGap = 2        # minimum gap 1 means don't search for pattern with gap length of 1
  maxGap = 4        # anything over maxGap is considered a no match
  htmlToPdf = True  # set to True if you want pdf output file
  anchor = { "st":"KWG", "en":"GMA" } # set values for st and en if using anchors
  #anchor = { "st":"", "en":"" }       # use this if not using anchors

  try:
    opts, args = getopt.getopt(argv, "hp:s:d",
                               ["help", "debug", "pattern=", "sample=",
                                "min-len=","min-gap=","max-gap=",
                                "out-pdf=","st-anchor=","en-anchor="])

  except getopt.GetoptError:
    usage()
    sys.exit(2)

  for opt, arg in opts:
    if opt in ("-h", "--help"):
      usage()
      sys.exit()
    elif opt in ("-d", "--debug"):
      print "[DEBUG MODE]"
      GlobalVars.DEBUG = True
    elif opt in ("-p", "--pattern"):
      patternFile = arg
    elif opt in ("-s", "--sample"):
      sampleFile = arg
    elif opt == "--min-len":
      minLen = int(arg)
    elif opt == "--min-gap":
      minGap = int(arg)
    elif opt == "--max-gap":
      maxGap = int(arg)
    elif opt == "--out-pdf":
      htmlToPdf = False if arg == '0' else True
    elif opt == "--st-anchor":
      anchor["st"] = arg
    elif opt == "--en-anchor":
      anchor["en"] = arg

  print "Configuration Settings"
  print "----------------------"
  print "Minimum Length: " + str(minLen)
  print "Minimum Gap: " + str(minGap)
  print "Maximum Gap: " + str(maxGap)
  print "PDF Output: " + str(htmlToPdf)
  print "Start Anchor: " + str(anchor["st"])
  print "End Anchor: " + str(anchor["en"])
  print "Sample File: " + str(sampleFile)
  print "\n"

  # Step 1: Retrieve allele patterns
  alleleIds, allelePatterns, numPatterns, _ = GetSequences(patternFile, "fasta", os.getcwd() + "/Library/", True, "Allele Patterns")

  print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
  print "Pattern Sequences:"

  sampFname = sampleFile.split('/')[-1]
  sampFname = sampFname.split('.')[0]
  print sampFname
  cwd = os.getcwd()
  outputPath = cwd + '/Results/' + sampFname
  print outputPath
  MakeDir(outputPath)


  if htmlToPdf:
    html = GlobalVars.HTML_HEAD_TAG
    pdfFile   = open(outputPath + "/" + sampFname + ".pdf", "w+b")
    htmlFile  = open(outputPath + "/" + sampFname + ".htm", "w")
    bodyHTML  = "<span class=\"TEXT\">"
    bodyHTML += "Minimum Length: <strong>" + str(minLen) + "</strong><br />"
    bodyHTML += "Minimum Gap: <strong>" + str(minGap) + "</strong><br />"
    bodyHTML += "Maximum Gap: <strong>" + str(maxGap) + "</strong><br />"
    bodyHTML += "Start Anchor: <strong>" + str(anchor["st"]) + "</strong><br />"
    bodyHTML += "End Anchor: <strong>" + str(anchor["en"]) + "</strong><br />"
    bodyHTML += "Sample File: <strong>" + str(sampleFile) + "</strong><br /><br />"
    bodyHTML += "Pattern Sequences:<br />"
    bodyHTML += "</span>"

  for num in range(0, numPatterns):
    print "%s:\t%s" %(alleleIds[num], allelePatterns[num])

    if htmlToPdf:
      patHTML = "<span class=\"TEXT\">" + str(alleleIds[num]) + ": </span>"
      for letter in allelePatterns[num]:
        patHTML += "<span class=\"" + letter + "\">" + letter + "</span>"
      patHTML += "<br />"
      bodyHTML += patHTML

  if htmlToPdf:
    bodyHTML += "<br />"
  print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

  # Step 2: Load sample sequences and identify duplicates
  sampleIds, sampleSequences, totalSeqs, totalSeqsStCodons =  GetSequences(sampleFile, "fasta", os.getcwd() + "/", False, "Sample Sequences")

  # Step 3: Identify what pattern the sample sequence belongs
  numSamps = len(sampleIds)
  patLen = [len(allelePatterns[i]) for i in range(len(allelePatterns))]

  print "Total number of sequences (with duplicates)    : %d" %(totalSeqs)
  print "Total number of sequences (with no duplicates) : %d" %(numSamps)
  print "Total number of sequences (with stop codons)   : %d" %(totalSeqsStCodons)
  print

  if htmlToPdf:
    bodyHTML += "<span class=\"TEXT\">"
    bodyHTML += "Total number of sequences (with duplicates): " + str(totalSeqs) + "<br />"
    bodyHTML += "Total number of sequences (with no duplicates): " + str(numSamps) + "<br />"
    bodyHTML += "Total number of sequences (with stop codons): " + str(totalSeqsStCodons)
    bodyHTML += "</span><hr /><br />"
    matchHTML = "<div><pdf:nextpage /></div>"

  # Loop through all samples and determine it's allele pattern combination
  for num in range(0, numSamps):
    numberOfSequences = len(sampleIds[num].split(","))
    print sampleIds[num]

    patternMatch, patKey, patHTML = DeterminePattern(sampleSequences[num]
                                                     ,alleleIds
                                                     ,allelePatterns
                                                     ,numPatterns
                                                     ,patLen
                                                     ,minLen
                                                     ,minGap
                                                     ,htmlToPdf
                                                     ,anchor
                                                     )

    if patKey not in results:
      results[patKey] = numberOfSequences
      results_[patKey] = [(numberOfSequences,num)]
    else:
      results[patKey] += numberOfSequences
      results_[patKey].append((numberOfSequences,num))


    if htmlToPdf:
      sampIDHTML = "<a name=\"" + sampleSequences[num] + "\"><span class=\"TEXT\">"
      groupBy = 115
      totalLen = len(sampleIds[num])
      start = 0
      end = groupBy
      while True:
        sampIDHTML += sampleIds[num][start:end]
        if end >= totalLen:
          break
        sampIDHTML += "<br />"
        start = end
        end = end + groupBy
      sampIDHTML += "</a></span><br />"
      matchHTML += sampIDHTML + patHTML + "<br /><br />"


  sortedKeys = sorted(results, key=results.get, reverse=True)
  print "%s\t%s" %("# SEQ", "PATTERN")

  # Printing GOOD MATCH table
  if htmlToPdf:
    resultTable  = "<div><pdf:nextpage /></div>"
    resultTable += "<span class=\"TEXT\">SEQUENCES CONSIDERED A <strong>MATCH</strong></span>"
    resultTable += "<table width=\"100%\" border=\"1px\" cellpadding=\"2px\" class=\"TEXT\">"
    resultTable += "<tr><td width=\"80px\" valign=\"bottom\" height=\"15px\"><strong># OF SEQ</strong></td><td valign=\"bottom\"><strong>PATTERN</strong></td></tr>"

  totalSeq = 0
  noMatchKeys = []
  totalNoMatch = 0
  print "SEQUENCES CONSIDERED A MATCH"
  for key in sortedKeys:

    # Checks that gap lengths are <= maxGap
    tempList = key.split("/")
    noMatch = False
    for t in tempList:
      if (("gap" in t) and (int(t[:-3]) > maxGap)) or ("missing anchor" == t):
        totalNoMatch += results[key]
        noMatchKeys.append(key)
        noMatch = True
        break
    if noMatch:
      continue

    # Print out data for result table for GOOD MATCH
    totalSeq += results[key]
    print "%4d\t%s" %(results[key], key)

    if htmlToPdf:
      sResultTab = "<table cellpadding=\"0px\" width=\"100%\">"
    subResults = sorted(results_[key], reverse=True)
    WriteToFile(sampFname, key, subResults, sampleIds, sampleSequences)
    for sResult in subResults:
      print "\t%4d\t%s" %(sResult[0], sampleSequences[sResult[1]])

      if htmlToPdf:
        sampHTML = ""
        for letter in sampleSequences[sResult[1]]:
          sampHTML += "<span class=\"" + letter + "\">" + letter + "</span>"

        sResultTab += "<tr><td class=\"TEXT\" width=\"50px\"><a href=\"#" + sampleSequences[sResult[1]] + "\">" + str(sResult[0]) + \
                      "</a></td><td>" + sampHTML + "</td></tr>"

    if htmlToPdf:
      sResultTab  += "</table>"
      resultTable += "<tr><td style=\"vertical-align:middle\" align=\"center\">" + str(results[key]) + \
                     "</td><td style=\"vertical-align:middle\">" + str(key) + " " + sResultTab + "</td></tr>"

  if htmlToPdf:
    compute = totalSeqs - totalSeqsStCodons - totalNoMatch
    resultTable += "<tr><td style=\"vertical-align:middle\" align=\"center\"><em>" + str(totalSeq) + \
                   "</em></td><td style=\"vertical-align:middle\"><em>Seq w/ duplicates - Seq w/ stop codons - Seq (no match) = " + str(totalSeqs) + " - " + str(totalSeqsStCodons) + " - " + str(totalNoMatch) + " = " + str(compute) + "</em></td></tr>"
    resultTable += "</table><br />"


  # Printing NO MATCH table
  if htmlToPdf:
    resultTable += "<div><pdf:nextpage /></div>"
    resultTable += "<span class=\"TEXT\">SEQUENCES CONSIDERED A <strong>NO MATCH</strong></span>"
    resultTable += "<table width=\"100%\" border=\"1px\" cellpadding=\"2px\" class=\"TEXT\">"
    resultTable += "<tr><td width=\"80px\" valign=\"bottom\" height=\"15px\"><strong># OF SEQ</strong></td><td valign=\"bottom\"><strong>PATTERN</strong></td></tr>"

  print "SEQUENCES CONSIDERED A NO MATCH"
  for key in noMatchKeys:
    # Print out data for result table for NO MATCH
    print "%4d\t%s" %(results[key], key)

    if htmlToPdf:
      sResultTab = "<table cellpadding=\"0px\" width=\"100%\">"
    subResults = sorted(results_[key], reverse=True)
    WriteToFileNoMatches(sampFname, key, subResults, sampleIds, sampleSequences)
    for sResult in subResults:
      print "\t%4d\t%s" %(sResult[0], sampleSequences[sResult[1]])

      if htmlToPdf:
        sampHTML = ""
        for letter in sampleSequences[sResult[1]]:
          sampHTML += "<span class=\"" + letter + "\">" + letter + "</span>"

        sResultTab += "<tr><td class=\"TEXT\" width=\"50px\"><a href=\"#" + sampleSequences[sResult[1]] + "\">" + str(sResult[0]) + \
                      "</a></td><td>" + sampHTML + "</td></tr>"

    if htmlToPdf:
      sResultTab  += "</table>"
      resultTable += "<tr><td style=\"vertical-align:middle\" align=\"center\">" + str(results[key]) + \
                     "</td><td style=\"vertical-align:middle\">" + str(key) + " " + sResultTab + "</td></tr>"

  if htmlToPdf:
    resultTable += "<tr><td style=\"vertical-align:middle\" align=\"center\"><em>" + str(totalNoMatch) + \
                   "</em></td><td style=\"vertical-align:middle\"><em>Seq w/ gap length > " + str(maxGap) + "</em></td></tr>"
    resultTable += "</table>"


  # Write PDF File
  if htmlToPdf:
    html += bodyHTML + resultTable + "<br /><br />" + matchHTML + "<p style=\"margin-bottom: 15cm;\">&nbsp;</p>"
    pisa.CreatePDF(html, dest=pdfFile)
    htmlFile.write(html)
    htmlFile.close()
    pdfFile.close()

'''
  sortedKeys = sorted(results, key=results.get, reverse=True)
  print "PATTERN   # SEQUENCES"
  for key in sortedKeys:
    print "%-7s   %d" %(key, results[key][0])
    WriteToFile(sampleFile.split('.')[0], key, results[key][1], sampleIds, sampleSequences)
  print

  noMatchLen = len(noMatch)
  if noMatchLen > 0:
    print "Total number of sequences (with no match) : %d" %(totalNoMatches)
    print
    WriteToFileNoMatches(sampleFile.split('.')[0], noMatch, sampleIds, sampleSequences)
    for index in range(0, noMatchLen):
      print "Seq %d: %s" %(index+1, noMatch[index][0])
      print "%s\n" %(noMatch[index][1])
'''


def run(pattern, sample):
  patternFile = pattern
  sampleFile = sample

  #print "[DEBUG MODE]"
  #GlobalVars.DEBUG = True

  # CONFIGURATION VARIABLES
  noMatch = []
  totalNoMatches = 0
  results = {}
  results_ = {}
  minLen = 4        # minimum length of pattern match
  minGap = 0        # minimum gap 1 means don't search for pattern with gap length of 1
  maxGap = 4        # anything over maxGap is considered a no match
  htmlToPdf = True  # set to True if you want pdf output file

  # Step 1: Retrieve allele patterns
  alleleIds, allelePatterns, numPatterns, _ = GetSequences(patternFile, "fasta", os.getcwd() + "/Library/", True, "Allele Patterns")

  print "Sample File: " + str(sampleFile)
  print

  print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
  print "Pattern Sequences:"

  if htmlToPdf:
    pdfFile   = open("Results/" + sampleFile.split('.')[0] + ".pdf", "w+b")
    bodyHTML  = "<span class=\"TEXT\">"
    bodyHTML += "Sample File: <strong>" + str(sampleFile) + "</strong><br /><br />"
    bodyHTML += "Pattern Sequences:<br />"
    bodyHTML += "</span>"

  for num in range(0, numPatterns):
    print "%s:\t%s" %(alleleIds[num], allelePatterns[num])

    if htmlToPdf:
      patHTML = "<span class=\"TEXT\">" + str(alleleIds[num]) + ": </span>"
      for letter in allelePatterns[num]:
        patHTML += "<span class=\"" + letter + "\">" + letter + "</span>"
      patHTML += "<br />"
      bodyHTML += patHTML

  if htmlToPdf:
    bodyHTML += "<br />"
  print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

  # Step 2: Load sample sequences and identify duplicates
  sampleIds, sampleSequences, totalSeqs, totalSeqsStCodons =  GetSequences(sampleFile, "fasta", os.getcwd() + "/", False, "Sample Sequences")

  # Step 3: Identify what pattern the sample sequence belongs
  numSamps = len(sampleIds)
  patLen = [len(allelePatterns[i]) for i in range(len(allelePatterns))]

  print "Total number of sequences (with duplicates)    : %d" %(totalSeqs)
  print "Total number of sequences (with no duplicates) : %d" %(numSamps)
  print "Total number of sequences (with stop codons)   : %d" %(totalSeqsStCodons)
  print

  if htmlToPdf:
    bodyHTML += "<span class=\"TEXT\">"
    bodyHTML += "Total number of sequences (with duplicates): " + str(totalSeqs) + "<br />"
    bodyHTML += "Total number of sequences (with no duplicates): " + str(numSamps) + "<br />"
    bodyHTML += "Total number of sequences (with stop codons): " + str(totalSeqsStCodons)
    bodyHTML += "</span><hr /><br />"

  # Loop through all samples and determine it's allele pattern combination
  for num in range(0, numSamps):
    numberOfSequences = len(sampleIds[num].split(","))
    print sampleIds[num]

    patternMatch, patKey, patHTML = DeterminePattern(sampleSequences[num],
                                                     alleleIds,
                                                     allelePatterns,
                                                     numPatterns,
                                                     patLen,
                                                     minLen,
                                                     minGap,
                                                     htmlToPdf)

    if patKey not in results:
      results[patKey] = numberOfSequences
      results_[patKey] = [(numberOfSequences,num)]
    else:
      results[patKey] += numberOfSequences
      results_[patKey].append((numberOfSequences,num))


    if htmlToPdf:
      sampIDHTML = "<span class=\"TEXT\">"
      groupBy = 115
      totalLen = len(sampleIds[num])
      start = 0
      end = groupBy
      while True:
        sampIDHTML += sampleIds[num][start:end]
        if end >= totalLen:
          break
        sampIDHTML += "<br />"
        start = end
        end = end + groupBy
      sampIDHTML += "</span><br />"
      bodyHTML += sampIDHTML + patHTML + "<br /><br />"


  sortedKeys = sorted(results, key=results.get, reverse=True)
  print "%s\t%s" %("# SEQ", "PATTERN")

  # Printing GOOD MATCH table
  if htmlToPdf:
    resultTable  = "<div><pdf:nextpage /></div>"
    resultTable += "<span class=\"TEXT\"><center>SEQUENCES CONSIDERED A <strong>MATCH</strong></center></span>"
    resultTable += "<table width=\"720px\" border=\"1px\" cellpadding=\"2px\" class=\"TEXT\">"
    resultTable += "<tr><td width=\"80px\" valign=\"bottom\" height=\"15px\"><strong># OF SEQ</strong></td><td valign=\"bottom\"><strong>PATTERN</strong></td></tr>"

  totalSeq = 0
  noMatchKeys = []
  totalNoMatch = 0
  print "SEQUENCES CONSIDERED A MATCH"
  for key in sortedKeys:

    # Checks that gap lengths are <= maxGap
    tempList = key.split("/")
    noMatch = False
    for t in tempList:
      if ("gap" in t) and (int(t[:-3]) > maxGap):
        totalNoMatch += results[key]
        noMatchKeys.append(key)
        noMatch = True
        break
    if noMatch:
      continue

    # Print out data for result table for GOOD MATCH
    totalSeq += results[key]
    print "%4d\t%s" %(results[key], key)

    if htmlToPdf:
      sResultTab = "<table cellpadding=\"0px\" width=\"640px\">"
    subResults = sorted(results_[key], reverse=True)
    for sResult in subResults:
      print "\t%4d\t%s" %(sResult[0], sampleSequences[sResult[1]])

      if htmlToPdf:
        sampHTML = ""
        for letter in sampleSequences[sResult[1]]:
          sampHTML += "<span class=\"" + letter + "\">" + letter + "</span>"

        sResultTab += "<tr><td class=\"TEXT\" width=\"50px\">" + str(sResult[0]) + \
                      "</td><td>" + sampHTML + "</td></tr>"

    if htmlToPdf:
      sResultTab  += "</table>"
      resultTable += "<tr><td style=\"vertical-align:middle\" align=\"center\">" + str(results[key]) + \
                     "</td><td style=\"vertical-align:middle\">" + str(key) + " " + sResultTab + "</td></tr>"

  if htmlToPdf:
    compute = totalSeqs - totalSeqsStCodons - totalNoMatch
    resultTable += "<tr><td style=\"vertical-align:middle\" align=\"center\"><em>" + str(totalSeq) + \
                   "</em></td><td style=\"vertical-align:middle\"><em>Seq w/ duplicates - Seq w/ stop codons - Seq (no match) = " + str(totalSeqs) + " - " + str(totalSeqsStCodons) + " - " + str(totalNoMatch) + " = " + str(compute) + "</em></td></tr>"
    resultTable += "</table>"


  # Printing NO MATCH table
  if htmlToPdf:
    resultTable += "<div><pdf:nextpage /></div>"
    resultTable += "<span class=\"TEXT\"><center>SEQUENCES CONSIDERED A <strong>NO MATCH</strong></center></span>"
    resultTable += "<table width=\"720px\" border=\"1px\" cellpadding=\"2px\" class=\"TEXT\">"
    resultTable += "<tr><td width=\"80px\" valign=\"bottom\" height=\"15px\"><strong># OF SEQ</strong></td><td valign=\"bottom\"><strong>PATTERN</strong></td></tr>"

  print "SEQUENCES CONSIDERED A NO MATCH"
  for key in noMatchKeys:
    # Print out data for result table for NO MATCH
    print "%4d\t%s" %(results[key], key)

    if htmlToPdf:
      sResultTab = "<table cellpadding=\"0px\" width=\"640px\">"
    subResults = sorted(results_[key], reverse=True)
    for sResult in subResults:
      print "\t%4d\t%s" %(sResult[0], sampleSequences[sResult[1]])

      if htmlToPdf:
        sampHTML = ""
        for letter in sampleSequences[sResult[1]]:
          sampHTML += "<span class=\"" + letter + "\">" + letter + "</span>"

        sResultTab += "<tr><td class=\"TEXT\" width=\"50px\">" + str(sResult[0]) + \
                      "</td><td>" + sampHTML + "</td></tr>"

    if htmlToPdf:
      sResultTab  += "</table>"
      resultTable += "<tr><td style=\"vertical-align:middle\" align=\"center\">" + str(results[key]) + \
                     "</td><td style=\"vertical-align:middle\">" + str(key) + " " + sResultTab + "</td></tr>"

  if htmlToPdf:
    resultTable += "<tr><td style=\"vertical-align:middle\" align=\"center\"><em>" + str(totalNoMatch) + \
                   "</em></td><td style=\"vertical-align:middle\"><em>Seq w/ gap length > " + str(maxGap) + "</em></td></tr>"
    resultTable += "</table>"

  # Write PDF File
  if htmlToPdf:
    global html
    html += bodyHTML + resultTable
    pisa.CreatePDF(html, dest=pdfFile)
    pdfFile.close()


if __name__ == '__main__':
  main(sys.argv[1:])
