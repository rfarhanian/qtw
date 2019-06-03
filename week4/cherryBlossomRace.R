# install.packages("XML")
library(XML)

shons_path <- ''
maryams_path <- ''
ramins_path <- '~/PycharmProjects/qtw/week4/data/Runners/WebPages/'
current_path <- ramins_path
setwd(current_path)

ubase = "http://www.cherryblossom.org/"

menURLs =
  c(
    "results/1999/cb99m.html",
    "results/2000/Cb003m.htm",
    "results/2001/oof_m.html",
    "results/2002/oofm.htm",
    "results/2003/CB03-M.HTM",
    "results/2004/men.htm",
    "results/2005/CB05-M.htm",
    "results/2006/men.htm",
    "results/2007/men.htm",
    "results/2008/men.htm",
    "results/2009/09cucb-M.htm",
    "results/2010/2010cucb10m-m.htm",
    "results/2011/2011cucb10m-m.htm",
    "results/2012/2012cucb10m-m.htm"
  )

womenURLs =
  c(
    "results/1999/cb99f.html",
    "results/2000/Cb003f.htm",
    "results/2001/oof_f.html",
    "results/2002/ooff.htm",
    "results/2003/CB03-F.HTM",
    "results/2004/women.htm",
    "results/2005/CB05-F.htm",
    "results/2006/women.htm",
    "results/2007/women.htm",
    "results/2008/women.htm",
    "results/2009/09cucb-F.htm",
    "results/2010/2010cucb10m-f.htm",
    "results/2011/2011cucb10m-f.htm",
    "results/2012/2012cucb10m-f.htm"
  )

maleUrls = paste(ubase, menURLs, sep = "")
femaleUrls = paste(ubase, womenURLs, sep = "")

extractResTable =
  #
  # Retrieve data from web site,
  # find the preformatted text,
  # and write lines or return as a character vector.
  #
  function(url = "http://www.cherryblossom.org/results/2009/09cucb-F.htm",
           year = 1999,
           sex = "male",
           file = NULL)
  {
    doc = htmlParse(url)
    
    if (year == 1999) {
      pres = getNodeSet(doc, "//pre")
      txt = xmlValue(pres[[1]])
      els = strsplit(txt, "\n")[[1]]
    }
    else if (year == 2000) {
      # Get preformatted text from 4th font element
      # The top file is ill formed so the <pre> search doesn't work.
      ff = getNodeSet(doc, "//font")
      txt = xmlValue(ff[[4]])
      els = strsplit(txt, "\r\n")[[1]]
    }
    else if (year == 2009 & sex == "male") {
      # Get preformatted text from <div class="Section1"> element
      # Each line of results is in a <pre> element
      div1 = getNodeSet(doc, "//div[@class='Section1']")
      pres = getNodeSet(div1[[1]], "//pre")
      els = sapply(pres, xmlValue)
    }
    else {
      # Get preformatted text from <pre> elements
      pres = getNodeSet(doc, "//pre")
      txt = xmlValue(pres[[1]])
      els = strsplit(txt, "\r\n")[[1]]
    }
    if (is.null(file))
      return(els)
    # Write the lines as a text file.
    writeLines(els, con = file)
  }

years = 1999:2012
menTables = mapply(extractResTable,
                   url = maleUrls,
                   year = years,
                   sex = "male")
womenTables = mapply(extractResTable,
                     url = femaleUrls,
                     year = years,
                     sex = "female")
names(menTables) = years
names(womenTables) = years
invisible(sapply(menTables, length))
invisible(sapply(womenTables, length))
save(menTables, file = "CBMenTextTables.rda")
save(womenTables, file = "CBWomenTextTables.rda")

womenTables$'2001'[2:3] <- womenTables$'2002'[2:3]

sapply(c("MenTxt", "WomenTxt"), function(dir)
{
  dir.create(file.path(getwd(), dir))
})


write(x = menTables$'2011', file = "MenTxt/2011.txt")
write(x = menTables$'2012', file = "MenTxt/2012.txt")
write(x = menTables$'2010', file = "MenTxt/2010.txt")
write(x = menTables$'2009', file = "MenTxt/2009.txt")
write(x = menTables$'2008', file = "MenTxt/2008.txt")
write(x = menTables$'2007', file = "MenTxt/2007.txt")
write(x = menTables$'2006', file = "MenTxt/2006.txt")
write(x = menTables$'2005', file = "MenTxt/2005.txt")
write(x = menTables$'2004', file = "MenTxt/2004.txt")
write(x = menTables$'2003', file = "MenTxt/2003.txt")
write(x = menTables$'2002', file = "MenTxt/2002.txt")
write(x = menTables$'2001', file = "MenTxt/2001.txt")
write(x = menTables$'2000', file = "MenTxt/2000.txt")
write(x = menTables$'1999', file = "MenTxt/1999.txt")
write(x = womenTables$'2011', file = "WomenTxt/2011.txt")
write(x = womenTables$'2012', file = "WomenTxt/2012.txt")
write(x = womenTables$'2010', file = "WomenTxt/2010.txt")
write(x = womenTables$'2009', file = "WomenTxt/2009.txt")
write(x = womenTables$'2008', file = "WomenTxt/2008.txt")
write(x = womenTables$'2007', file = "WomenTxt/2007.txt")
write(x = womenTables$'2006', file = "WomenTxt/2006.txt")
write(x = womenTables$'2005', file = "WomenTxt/2005.txt")
write(x = womenTables$'2004', file = "WomenTxt/2004.txt")
write(x = womenTables$'2003', file = "WomenTxt/2003.txt")
write(x = womenTables$'2002', file = "WomenTxt/2002.txt")
write(x = womenTables$'2001', file = "WomenTxt/2001.txt")
write(x = womenTables$'2000', file = "WomenTxt/2000.txt")
write(x = womenTables$'1999', file = "WomenTxt/1999.txt")

# 2012 example race logs

els = readLines("WomenTxt/2012.txt")

# Identify line index for header-data break

eqIndex = grep("^===", els)
eqIndex

first3 = substr(els, 1, 3)
which(first3 == "===")

# Ignore rows above header

spacerRow = els[eqIndex]
headerRow = els[eqIndex - 1]
body = els[-(1:eqIndex)]

# Extract runners' age

headerRow = tolower(headerRow)

ageStart = regexpr("ag", headerRow)
ageStart

age = substr(body, start = ageStart, stop = ageStart + 1)
head(age)

summary(as.numeric(age))

blankLocs = gregexpr(" ", spacerRow)
blankLocs

searchLocs = c(0, blankLocs[[1]])

Values = mapply(substr, list(body),
                start = searchLocs[-length(searchLocs)] + 1,
                stop = searchLocs[-1] - 1)

#Find locations of all blanks in line of '=' characters and extract columns
findColLocs = function(spacerRow) {
  spaceLocs = gregexpr(" ", spacerRow)[[1]]
  rowLength = nchar(spacerRow)
  
  if (substring(spacerRow, rowLength, rowLength) != " ")
    return(c(0, spaceLocs, rowLength + 1))
  else
    return(c(0, spaceLocs))
}

# Extract columns

selectCols =  function(colNames, headerRow, searchLocs)
{
  sapply(colNames,
         function(name, headerRow, searchLocs)
         {
           startPos = regexpr(name, headerRow)[[1]]
           if (startPos == -1)
             return(c(NA, NA))
           
           index = sum(startPos >= searchLocs)
           c(searchLocs[index] + 1, searchLocs[index + 1] - 1)
         },
         headerRow = headerRow, searchLocs = searchLocs)
}

# Testing findColLocs and selectCols functions

searchLocs = findColLocs(spacerRow)
ageLoc = selectCols("ag", headerRow, searchLocs)
ages = mapply(substr, list(body), start = ageLoc[1, ], stop = ageLoc[2,])

summary(as.numeric(ages))

# Create shortened column identifiers and account for when some tables missing columns

shortColNames = c("name", "home", "ag", "gun", "net", "time")

locCols = selectCols(shortColNames, headerRow, searchLocs)

Values = mapply(substr, list(body), start = locCols[1,], stop = locCols[2,])

colnames(Values) = shortColNames
head(Values)

tail(Values)[, 1:3]

# Build wrapper function for column extraction

extractVariables = function(file, varNames = c("name", "home", "ag", "gun", "net", "time"))
{
  eqIndex = grep("^===", file)
  # Extract the two key rows and the real data in the xml
  spacerRow = file[eqIndex]
  headerRow = tolower(file[eqIndex - 1])
  body = file[-(1:eqIndex)]
  
  # Obtain the starting and ending positions of variables
  searchLocs = findColLocs(spacerRow)
  locCols = selectCols(varNames, headerRow, searchLocs)
  
  Values = mapply(substr, list(body), start = locCols[1,], stop = locCols[2,])
  colnames(Values) = varNames
  
  invisible(Values)
}

# Read table lines into R

mfilenames = paste("MenTxt/", 1999:2012, ".txt", sep = "")
menFiles = lapply(mfilenames, readLines)
names(menFiles) = 1999:2012
menFiles[['2009']] <- gsub("�", "", menFiles[['2009']])

# Create list of character matrices containing the column contents for each of the 14 years of data

menResMat = lapply(menFiles, extractVariables)
sapply(menResMat, nrow)

# Read table lines into R

wfilenames = paste("WomenTxt/", 1999:2012, ".txt", sep = "")
womenFiles = lapply(wfilenames, readLines)
names(womenFiles) = 1999:2012
womenFiles[['2009']] <- gsub("�", "", womenFiles[['2009']])

# Create list of character matrices containing the column contents for each of the 14 years of data

womenResMat = lapply(womenFiles, extractVariables)
sapply(womenResMat, nrow)

maleAge = sapply(menResMat, function(x)
  as.numeric(x[, 'ag']))
femaleAge = sapply(womenResMat, function(x)
  as.numeric(x[, 'ag']))

head(menFiles[['2003']])
menFiles[['2006']][2200:2205]

# Updated selectCols takes care of the offset in age column
selectCols = function(shortColNames, headerRow, searchLocs) {
  sapply(shortColNames, function(shortName, headerRow, searchLocs) {
    startPos = regexpr(shortName, headerRow)[[1]]
    if (startPos == -1)
      return(c(NA, NA))
    index = sum(startPos >= searchLocs)
    c(searchLocs[index] + 1, searchLocs[index + 1])
  }, headerRow = headerRow, searchLocs = searchLocs)
}

menResMat = lapply(menFiles, extractVariables)
womenResMat = lapply(womenFiles, extractVariables)

maleAge = sapply(menResMat, function(x)
  as.numeric(x[, 'ag']))
femaleAge = sapply(womenResMat, function(x)
  as.numeric(x[, 'ag']))

sapply(maleAge,  function(x)
  sum(is.na(x)))
sapply(femaleAge,  function(x)
  sum(is.na(x)))

# Example 
maleAge2001 = maleAge[["2001"]]
femaleAge2001 = femaleAge[["2001"]]
grep("^===", menFiles[['2001']])
grep("^===", womenFiles[['2001']])
badAgeIndex = which(is.na(maleAge2001)) + 5
femaleBadAgeIndex = which(is.na(femaleAge2001)) + 5

# new revision of extractVariables handles missing age data

extractVariables = function(file, varNames = c("name", "home", "ag", "gun", "net", "time"))
{
  # Find the index of the row with =s
  eqIndex = grep("^===", file)
  # Extract the two key rows and the data
  spacerRow = file[eqIndex]
  headerRow = tolower(file[eqIndex - 1])
  body = file[-(1:eqIndex)]
  # Remove footnotes and blank rows
  footnotes = grep("^[[:blank:]]*(\\*|\\#)", body)
  if (length(footnotes) > 0)
    body = body[-footnotes]
  blanks = grep("^[[:blank:]]*$", body)
  if (length(blanks) > 0)
    body = body[-blanks]
  
  
  # Obtain the starting and ending positions of variables
  searchLocs = findColLocs(spacerRow)
  locCols = selectCols(varNames, headerRow, searchLocs)
  
  Values = mapply(substr, list(body), start = locCols[1,], stop = locCols[2,])
  colnames(Values) = varNames
  return(Values)
}

menResMat = lapply(menFiles, extractVariables)
womenResMat = lapply(womenFiles, extractVariables)

maleCharTime = menResMat[['2012']][, 'time']
femaleCharTime = womenResMat[['2012']][, 'time']

maleTimePieces = strsplit(maleCharTime, ":")
femaleTimePieces = strsplit(femaleCharTime, ":")

maleTimePieces = sapply(maleTimePieces, as.numeric)
femaleTimePieces = sapply(femaleTimePieces, as.numeric)

maleRunTime = sapply(maleTimePieces,
                  function(x) {
                    if (length(x) == 2)
                      x[1] + x[2] / 60
                    else
                      60 * x[1] + x[2] + x[3] / 60
                  })
femaleRunTime = sapply(femaleTimePieces,
                  function(x) {
                    if (length(x) == 2)
                      x[1] + x[2] / 60
                    else
                      60 * x[1] + x[2] + x[3] / 60
                  })

# Split and process times

convertTime = function(time) {
  maleTimePieces = strsplit(time, ":")
  maleTimePieces = sapply(maleTimePieces, as.numeric)
  sapply(maleTimePieces, function(x) {
    if (length(x) == 2)
      x[1] + x[2] / 60
    else
      60 * x[1] + x[2] + x[3] / 60
  })
}

createDF = function(Res, year, sex)
{
  # Determine which time to use
  useTime = if (!is.na(Res[1, 'net']))
    Res[, 'net']
  else if (!is.na(Res[1, 'gun']))
    Res[, 'gun']
  else
    Res[, 'time']
  
  runTime = convertTime(useTime)
  
  Results = data.frame(
    year = rep(year, nrow(Res)),
    sex = rep(sex, nrow(Res)),
    name = Res[, 'name'],
    home = Res[, 'home'],
    age = as.numeric(Res[, 'ag']),
    runTime = runTime,
    stringsAsFactors = FALSE
  )
  invisible(Results)
}

menDF = mapply(
  createDF,
  menResMat,
  year = 1999:2012,
  sex = rep("M", 14),
  SIMPLIFY = FALSE
)

sapply(menDF, function(x)
  sum(is.na(x$maleRunTime)))

createDF = function(Res, year, sex)
{
  # Determine which time to use
  if (!is.na(Res[1, 'net']))
    useTime = Res[, 'net']
  else if (!is.na(Res[1, 'gun']))
    useTime = Res[, 'gun']
  else
    useTime = Res[, 'time']
  
  
  useTime = gsub("[#\\*[:blank:]]", "", useTime) # Remove footnote symbols from time
  runTime = convertTime(useTime[useTime != ""])
  
  # Drop rows with no time
  Res = Res[useTime != "",]
  
  Results = data.frame(
    year = rep(year, nrow(Res)),
    sex = rep(sex, nrow(Res)),
    name = Res[, 'name'],
    home = Res[, 'home'],
    age = as.numeric(Res[, 'ag']),
    runTime = runTime,
    stringsAsFactors = FALSE
  )
  invisible(Results)
}

menDF = mapply(
  createDF,
  menResMat,
  year = 1999:2012,
  sex = rep("M", 14),
  SIMPLIFY = FALSE
)

sapply(menDF, function(x)
  sum(is.na(x$runTime)))

# Fix missing runTime data issues for 2006

separatorIdx = grep("^===", menFiles[["2006"]])
separatorRow = menFiles[['2006']][separatorIdx]
separatorRowX = paste(substring(separatorRow, 1, 63),
                      " ",
                      substring(separatorRow, 65, nchar(separatorRow)),
                      sep = "")
menFiles[['2006']][separatorIdx] = separatorRowX

menResMat = sapply(menFiles, extractVariables)
menDF = mapply(
  createDF,
  menResMat,
  year = 1999:2012,
  sex = rep("M", 14),
  SIMPLIFY = FALSE
)

separatorIdx = grep("^===", womenFiles[["2006"]])
separatorRow = womenFiles[['2006']][separatorIdx]
separatorRowX = paste(substring(separatorRow, 1, 63),
                      " ",
                      substring(separatorRow, 65, nchar(separatorRow)),
                      sep = "")
womenFiles[['2006']][separatorIdx] = separatorRowX

womenResMat = sapply(womenFiles, extractVariables)
womenDF = mapply(
  createDF,
  womenResMat,
  year = 1999:2012,
  sex = rep("F", length(years)),
  SIMPLIFY = FALSE
)