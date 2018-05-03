library(genbankr)
library(Biostrings)
library(pbapply)

prime_field_re = "^[[:upper:]]+[[:space:]]+"
sec_field_re = "^( {5}|\\t)[[:alnum:]'_]+[[:space:]]+(complement|join|order|[[:digit:]<,])"

strip_fieldname = function(lines) gsub("^[[:space:]]*[[:upper:]]*[[:space:]]+", "", lines)


readGenBank2 <- function(file, text = readLines(file), partial = NA,
                       ret.seq = TRUE, verbose = FALSE) {
  
  if(missing(text)) {
    if(missing(file))
      stop("One of text or file must be specified.")
    if(is(file, "GBAccession"))
      text = .getGBfromNuccore(file)
    else if (!file.exists(file))
      stop("file does not appear to an existing file or a GBAccession object. Valid file or text argument is required.")
  }

  if(is(text, "list"))
    return(pblapply(text, function(txt) readGenBank(text = txt, partial = partial, ret.seq= ret.seq, verbose = verbose)))
  ## we always read in sequence because it is required for variant annotations
  ## we throw it away after if the user set ret.seq=FALSE
  prsed = parseGenBank2(text = text, partial = partial, verbose = verbose,
                       ret.seq = TRUE)
  ret = make_gbrecord(rawgbk = prsed, verbose = verbose)
  if(!ret.seq)
    sequence(ret) = NULL
  ret
}

parseGenBank2 = function(file, text = readLines(file),  partial = NA,
                        verbose = FALSE,
                        ret.anno = TRUE,
                        ret.seq = TRUE) {
  if(!ret.anno && !ret.seq)
    stop("Must return at least one of annotations or sequence.")
  bf = proc.time()["elapsed"]
  if(missing(text) && !file.exists(file))
    stop("No text provided and file does not exist or was not specified. Either an existing file or text to parse must be provided.")
  if(length(text) == 1)
    text = fastwriteread(text)
  
  fldlines = grepl(prime_field_re, text)
  fldfac = cumsum(fldlines)
  fldnames = gsub("^([[:upper:]]+).*", "\\1", text[fldlines])[fldfac]
  
  spl = split(text, fldnames)
  
  resthang = list(LOCUS = genbankr:::readLocus(spl[["LOCUS"]]))
  resthang[["FEATURES"]] = genbankr:::readFeatures(spl[["FEATURES"]],
                                        source.only=!ret.anno,
                                        partial = partial)
  seqtype = genbankr:::.seqTypeFromLocus(resthang$LOCUS)
  resthang$ORIGIN = if(ret.seq)
    genbankr:::readOrigin(spl[["ORIGIN"]],
               seqtype = seqtype)
  else NULL
  
  if(ret.anno) {
    resthang2 = mapply(function(field, lines, verbose) {
      switch(field,
             DEFINITION = genbankr:::readDefinition(lines),
             ACCESSION = genbankr:::readAccession(lines),
             VERSION = genbankr:::readVersions(lines),
             KEYWORDS = genbankr:::readKeywords(lines),
             SOURCE = genbankr:::readSource(lines),
             ## don't read FEATURES, ORIGIN, or LOCUS because they are
             ## already in resthang from above
             NULL)
    }, lines = spl, field = names(spl), SIMPLIFY=FALSE, verbose = verbose)
    resthang2$FEATURES = resthang2$FEATURES[sapply(resthang2$FEATURES,
                                                   function(x) length(x)>0)]
    resthang2 = resthang2[!names(resthang2) %in% names(resthang)]
    resthang = c(resthang, resthang2)
  }
  ##DNAString to DNAStringSet
  origin = resthang$ORIGIN 
  if(ret.seq && length(origin) > 0) {
    typs = sapply(resthang$FEATURES, function(x) x$type[1])
    srcs =  genbankr:::fill_stack_df(resthang$FEATURES[typs == "source"])
    ## dss = DNAStringSet(lapply(GRanges(ranges(srcs), function(x) origin[x])))
    dss = switch(seqtype,
                 bp = DNAStringSet(lapply(ranges(srcs), function(x) origin[x])),
                 aa = AAStringSet(lapply(ranges(srcs), function(x) origin[x])),
                 stop("Unrecognized origin sequence type: ", seqtype)
    )
    names(dss) = sapply(srcs,
                        function(x) as.character(seqnames(x)[1]))
    if(!ret.anno)
      resthang = dss
    else
      resthang$ORIGIN = dss
  } else if (!ret.anno) { ##implies ret.seq is TRUE
    stop("Asked for only sequence (ret.anno=FALSE) from a file with no sequence information")
  }
  af = proc.time()["elapsed"]
  if(verbose)
    message("Done Parsing raw GenBank file text. [ ", af-bf, " seconds ]") 
  resthang
  
  
}



dat <- readGenBank2("./data/UTI89.gbk")
