#!/usr/bin/python -tt
'''
Create on January 5, 2014

@author: Gabriel de la Cruz

'''

import GlobalVars
from DataLibs import *


def DeterminePattern(sample, patIds, patterns, numPatterns, patLen,
                      minLen = 4, minGap = 0, htmlToPdf = False,
                      anchor = {"st":"", "en":""}):
  sampleLength = len(sample)
  matches = {}

  hasAnchor = False
  if anchor["st"] != "" and anchor["en"] != "":
    hasAnchor = True
    anchStIdx = sample.find(anchor["st"])
    anchEnIdx = sample.find(anchor["en"])
    if anchStIdx == -1 or anchEnIdx == -1: # -1 not found
      patternID = "missing anchor"
      print "%s" %(sample)
      print "PATTERN: %s" %(patternID)
      print
      html = None
      if htmlToPdf:
        sampHTML = ""
        for letter in sample:
          sampHTML += "<span class=\"" + letter + "\">" + letter + "</span>"
        sampHTML += "<br />"
        html = sampHTML
        html += "<span class=\"TEXT\">PATTERN: " + str(patternID) + "</span>"
      return ([[(sampleLength, -1, (-1,-1), (0, sampleLength-1))]], patternID, html)

  colSt = 0
  sampLen = sampleLength
  if hasAnchor: # start similarity where the anchor begins
    colSt = anchStIdx
    sampLen = anchEnIdx + 3

  for patNum in range(len(patterns)):
    line = [None] * patLen[patNum]
    for row in range(patLen[patNum]):
      line[row] = [0] * sampleLength
      letP = patterns[patNum][row]
      nxtLetP = -1
      if row+1 != patLen[patNum]:
        nxtLetP = patterns[patNum][row+1]
      col = 0
      for col in range(colSt,sampLen):
        letS = sample[col]
        nxtLetS = -1
        if col+1 != sampLen:
          nxtLetS = sample[col+1]

        # letter match
        if letS == letP:
          if row > 0 and col > colSt and line[row-1][col-1] > 0:
            line[row][col] = line[row-1][col-1] + 1
            if (row+1 == patLen[patNum] or col+1 == sampLen) or \
               (nxtLetS != nxtLetP and line[row][col] >= minLen):
              patStIdx = row-line[row][col] + 1
              patIndices = (patStIdx,row)
              samStIdx = col-line[row][col] + 1
              samIndices = (samStIdx,col)
              length = line[row][col]
              match = (length, patNum, patIndices, samIndices)
              if samStIdx not in matches:
                matches[samStIdx] = []
              matches[samStIdx].append(match)
          else:
              if col + minLen < sampLen and \
                 row + minLen < patLen[patNum]:
                line[row][col] = 1

        # letter don't match
        else:
          # check minimum length lines deep
          if row > 0:
            if line[row-1][col-1] < minLen:
              i = 1
              while row-i >= 0 and col-i >= colSt and i <= minLen:
                if line[row-i][col-i] == 0:
                  break
                line[row-i][col-i] = 0
                i += 1

    # print Similarity Matrix
    if GlobalVars.DEBUG:
      print " %s" %(sample)
      for row in range(patLen[patNum]):
        lineStr = ''.join(str(x if x < 10 else x % 10) if x != 0 else '-' for x in line[row])
        letP = patterns[patNum][row]
        print "%c%s" %(letP, lineStr)

    del line

  #print matches

  print "%s" %(sample)

  # For HTML to PDF
  if htmlToPdf:
    sampHTML = ""
    for letter in sample:
      sampHTML += "<span class=\"" + letter + "\">" + letter + "</span>"
    sampHTML += "<br />"


  if GlobalVars.DEBUG:
    print "FILTER 1:"
    print "-" * sampleLength

  matchKeys = sorted(matches.keys())

  if GlobalVars.DEBUG:
    for key in matchKeys:
      matchList = sorted(matches[key], reverse=True)
      for match in matchList:
        patNum = match[1]
        patStIdx = match[2][0]
        patEnIdx = match[2][1]+1
        samStIdx = match[3][0]
        padding = " " * (samStIdx)
        patStr = patterns[patNum][patStIdx:patEnIdx]
        print padding + patStr,
        print "(%s)" %(patIds[patNum])

  if GlobalVars.DEBUG:
    print "FILTER 2:"
    print "-" * sampleLength

  finMatches = []
  for key in matchKeys:
    matchList = sorted(matches[key], reverse=True)
    maxLength = max(matchList)[0]
    matchList = [match for match in matchList if match[0] == maxLength]

    # check matchList against finMatches for absorption
    match = matchList[0]
    skipList = False
    for m in finMatches:
      fmatch = m[0]
      if match[3][0] >= fmatch[3][0] and match[3][1] <= fmatch[3][1]:
        skipList = True
        break

    if skipList:
      continue

    finMatches.append(matchList[::-1])

    if GlobalVars.DEBUG:
      for match in matchList:
        patNum = match[1]
        patStIdx = match[2][0]
        patEnIdx = match[2][1]+1
        samStIdx = match[3][0]
        padding = " " * (samStIdx)
        patStr = patterns[patNum][patStIdx:patEnIdx]
        print padding + patStr,
        print "(%s)" %(patIds[patNum])


  if GlobalVars.DEBUG:
    print "FILTER 3:"
    print "-" * sampleLength

  sampPad = ['-'] * sampleLength
  tempList = sorted(finMatches, reverse=True)
  for tL in tempList:
    t = tL[0]
    stIdx = t[3][0]
    enIdx = t[3][1]+1
    length = t[0]

    isRemove = True
    count = 0
    for c in sampPad[stIdx:enIdx]:
      if c == '-':
        count += 1
        if count > minGap:
          isRemove = False
          break

    if isRemove:
      finMatches.remove(tL)
      continue

    sampPad[stIdx:enIdx] = ['+'] * (length)

    if GlobalVars.DEBUG:
      print ''.join(c for c in sampPad)

  patternID = ""
  samEnIdx = -1
  finMatchLen = len(finMatches)
  patHTML = ""


  for i in range(finMatchLen):
    matchList = finMatches[i]
    nxSamStIdx = -1
    samStIdx = -1
    samEnIdx = -1
    if i+1 != finMatchLen:
      nxSamStIdx = finMatches[i+1][0][3][0] # get next start index
    subPatID = []
    matchLLen = len(matchList)

    for j in range(matchLLen):
      matchHTML = ""
      match = matchList[j]
      patNum = match[1]
      patStIdx = match[2][0]
      patEnIdx = match[2][1]+1
      samStIdx = match[3][0]
      samEnIdx = match[3][1]+1
      padding = " " * (samStIdx)
      patStr = patterns[patNum][patStIdx:patEnIdx]
      print padding + patStr,
      print "(%s)" %(patIds[patNum])
      subPatID.append(patIds[patNum])

      if htmlToPdf:
        paddingHTML = "<span class=\"SPACE\">" + ("." * (samStIdx)) + "</span>"
        for letter in patStr:
          matchHTML += "<span class=\"" + letter + "\">" + letter + "</span>"
        patNumS = " <span class=\"TEXT\">(" + patIds[patNum] + ")</span>"
        matchHTML_ = paddingHTML + matchHTML + patNumS + "<br />"
        patHTML += matchHTML_


    # check gap at the beginning
    if i == 0 and samStIdx != colSt:
      patternID += str(samStIdx)+"gap" + "/"

    subPatID = ','.join(subPatID)
    if matchLLen > 1:
      patternID += "(" + subPatID + ")/"
    else:
      patternID += subPatID + "/"

    # check for gap/missing
    gap = nxSamStIdx - samEnIdx
    if gap > 0:
      patternID += str(gap)+"gap" + "/"

  # check gap at the end
  gap = sampLen - samEnIdx
  if gap > 0:
    patternID += str(gap)+"gap" + "/"

  print "PATTERN: %s" %(patternID[:-1])
  print

  html = None
  if htmlToPdf:
    html = sampHTML + patHTML
    html += "<span class=\"TEXT\">PATTERN: " + str(patternID[:-1]) + "</span>"

  return (finMatches, patternID[:-1], html)
