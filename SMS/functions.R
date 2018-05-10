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

  lineage <- gsub("[ .]", "", prsed[["SOURCE"]][["lineage"]])
  
  all_features <- lapply(prsed[["FEATURES"]], function(ith_feature) {
    if(class(ith_feature) == "data.frame") {
      data.frame(protein_id = ith_feature[["protein_id"]], translation = ith_feature[["translation"]], 
                 stringsAsFactors = FALSE)
    } else {
      NULL
    }
  }) %>% 
    bind_rows() %>% 
    mutate(organism = prsed[["SOURCE"]][["organism"]],
           lineage = paste0(gsub("[ .]", "", prsed[["SOURCE"]][["lineage"]]), collapse = "|"),
           definition = prsed[["DEFINITION"]]) %>% 
    select(definition, organism, lineage, protein_id, translation)
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
  resthang[["FEATURES"]] = readFeatures2(spl[["FEATURES"]],
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
                        function(x) as.character(GenomeInfoDb:::seqnames(x)[1]))
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


readFeatures2 = function(lines, partial = NA, verbose = FALSE,
                         source.only = FALSE) {
  if(substr(lines[1], 1, 8) == "FEATURES")
    lines = lines[-1] ## consume FEATURES line
  fttypelins = grepl(sec_field_re, lines)
  featfactor  = cumsum(fttypelins)
  
  if(source.only) {
    srcfeats = which(substr(lines[fttypelins], 6, 11) == "source")
    keepinds = featfactor %in% srcfeats
    lines = lines[keepinds]
    featfactor = featfactor[keepinds]
  }

  ##scope bullshittery
  chr = "unk"
  
  totsources = length(grep("[[:space:]]+source[[:space:]]+[<[:digit:]]", lines[which(fttypelins)]))
  numsources = 0
  everhadchr = FALSE
  
  do_readfeat = function(lines, partial = NA) {
    
    ## before collapse so the leading space is still there
    type = gsub("[[:space:]]+([[:alnum:]_']+).*", "\\1", lines[1])
    ##feature/location can go across multpiple lines x.x why genbank? whyyyy
    attrstrts = cumsum(grepl("^[[:space:]]+/[^[:space:]]+($|=([[:digit:]]|\"))", lines))
    lines = tapply(lines, attrstrts, function(x) {
      paste(gsub("^[[:space:]]+", "", x), collapse="")
    }, simplify=TRUE)
    
    rawlocstring = lines[1]
    
    rngstr = genbankr:::strip_feat_type(rawlocstring)
    
    ## consume primary feature line
    lines = lines[-1] 
    if(length(lines)) {
      attrs = genbankr:::read_feat_attr(lines)
      
      names(attrs) = gsub("^[[:space:]]*/([^=]+)($|=[^[:space:]].*$)", "\\1", lines)
      if(type == "source") {
        numsources <<- numsources + 1
        if("chromosome" %in% names(attrs)) {
          if(numsources > 1 && !everhadchr)
            stop("This file appears to have some source features which specify chromosome names and others that do not. This is not currently supported. Please contact the maintainer if you need this feature.")    
          everhadchr <<- TRUE
          chr <<- attrs$chromosome
        } else if(everhadchr) {
          stop("This file appears to have some source features which specify chromosome names and others that do not. This is not currently supported. Please contact the maintainer if you need this feature.")
          ## this assumes that if one source has strain, they all will.
          ## Good assumption?
        } else if("strain" %in% names(attrs)) {
          chr <<- if(totsources == 1) attrs$strain else paste(attrs$strain, numsources, sep=":")
        } else {
          chr <<- if(totsources == 1) attrs$organism else paste(attrs$organism, numsources, sep=":")
        }
      }
    } else {
      attrs = list()
    }
    genbankr:::make_feat_gr(str = rngstr, chr = chr, ats = c(type = type, attrs),
                            partial = partial)
    
  }
  
  if(verbose)
    message(Sys.time(), " Starting feature parsing")
  
  resgrs = tapply(lines, featfactor, function(i) {
    do_readfeat(i, partial = partial)
  }, simplify=FALSE)
  if(verbose)
    message(Sys.time(), " Done feature parsing")
  
  resgrs
  
}



my_DT <- function(x)
  datatable(x, escape = FALSE, extensions = 'Buttons',
            filter = "top", rownames = FALSE) %>% 
  formatStyle(1L:ncol(x), color = "black")