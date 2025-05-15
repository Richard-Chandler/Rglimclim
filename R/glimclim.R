##############################################################################
##############################################################################
##############################################################################
######                                                                  ######
######          R wrapper routines for GLIMCLIM                         ######
######                                                                  ######
##############################################################################
##############################################################################
##############################################################################
.onAttach <- function(...) {
 packageStartupMessage("Use 'help(\"Rglimclim-package\")' to get started")
}
##############################################################################
file.OK <- function(filename,status) {
#
#       To check the existence of an input or output file. Arguments:
#
#       filename        The name of the file to be checked
#       status          Either "old" or "new"
#
#       Value: a logical scalar indicating whether or not the 
#       check was passed. If an "old" file doesn't exist, 
#       the routine terminates.
#
 file.ok <- TRUE
 if ( (status=="old") & !file.exists(filename) ) {
  stop(paste("Input file",filename,"doesn't exist."),call.=FALSE)
 }
 if ((status == "new") & file.exists(filename)) {
  cat(paste('\n     ****NOTE**** Output file',filename,
            'already exists.\n     Overwrite it (Y/N, default N)? '))
  option <- readLines(n=1)
  if (!(option %in% c("y","Y"))) file.ok <- FALSE
 }
 file.ok
}
##############################################################################
##############################################################################
##############################################################################
##reformat.scodes <- function(scodes) {
##
##       DEFUNCT FROM RGLIMCLIM V1.4.0
##
##       To reformat a vector of 4-character site codes so that they 
##       can be passed to Fortran. This is necessary because .Fortran()
##       won't pass character arrays(!). However, the Fortran routines
##       need the site codes to identify  observations in data files. 
##       Solution: identify all of the unique characters present in the
##       site codes and concatenate these into a single string; and
##       pass *this* to Fortran along with a numeric array with 4 columns,
##       each picking out the corresponding character in the string. 
##       This is the kind of creativity that comes of desperation ...
##
# chars.used <- unique(unlist(strsplit(scodes,split="")))
# nsites <- length(scodes)
# char.table <- unlist(strsplit(scodes,""))
# char.ids <- t(matrix(match(char.table,chars.used),nrow=4))
# storage.mode(char.ids) <- "integer"
# chars.used <- paste(chars.used,collapse="")
# list(Characters.used=chars.used,IDs=char.ids)
#}
##############################################################################
##############################################################################
##############################################################################
err.msg <- function(filename,line.no,err.no) {
#
#       Input error trapping - called if an error of code err.no 
#       is found on line line.no of file filename. 
#
 if (err.no == 1) {
  msg <- paste("Problem reading file ",filename,". This could be a\n",
               "programming error - please contact Richard ",
               "(r.chandler@ucl.ac.uk).\n",sep="")
 } else if (err.no == 2) {
  msg <- paste("Input error while reading line ",line.no,
               " of file ",filename,".",sep="")
 } else if (err.no == 3) {
  msg <- paste("Line",line.no,"of file",filename,"is out of order")
 } else if (err.no == 4) {
  msg <- paste("There is a problem with input file ",filename,
               " - my guess is that\n  the number of header rows ",
               "is wrongly specified.",sep="")
 } else if (err.no == 5) {
  msg <- paste("Error while reading file ",filename,".",sep="")               
 } else if (err.no == 6) {
  msg <- paste("Mismatch between number of variables in file ",filename,
               " and\n  number of variables named in model definition. Set",
               " nvar.check=FALSE\n  to suppress this warning.",sep="")
 } else if (err.no == 7) {
  msg <- paste("No data found. Check site definitions and data files.",sep="")
 } else if (err.no == 17) {
  msg <- paste("The site attribute requested on line ",line.no,
               " of file ",filename,"\n  has not been defined.",sep="")
 } else if (err.no == 18) {
  msg <- paste("On line ",line.no," of file ",filename,
               " you ask for distance-weighted\n  averages of ",
               "previous days' values, but your siteinfo object",
               "contains\n  fewer than two site attributes so I ",
               "can't calculate distances.",sep="")
 } else if (err.no == 19) {
  msg <- paste("Illegal value of code 1 on line ",line.no,
               " of file ",filename,".\n",sep="")
 } else if (err.no == 20) {
  msg <- paste("Invalid model component on line",line.no,
               "of file",filename)
 } else if (err.no == 21) {
  msg <- paste("Illegal interaction specified on line",line.no,
               "of file",filename)
 } else if (err.no == 22) {
  msg <- paste("On line ",line.no," of file ",filename,
               ", you have asked me\n  to keep track of too ",
               "many previous days' values. I can't\n  remember",
               "more than 10!",sep="")
 } else if (err.no == 23) {
  msg <- paste("Model definition error on line ",line.no," of file ",
               filename,":\n  response variable can be used as a ",
               "covariate only at lags > 0.",sep="")
 } else if (err.no == 25) {
  msg <- paste("In file ",filename,", some covariates are 'trace' ",
               "indicators, but you haven't\n  defined a 'trace' ",
               "threshold in the 'global quantities' section ",
               "of the file.",sep="")
 } else if (err.no == 26) {
  msg <- paste("On line ",line.no," of file ",filename," you define a 'trace' ",
               "threshold. This is not\n  permitted in Gaussian models.",sep="")
 } else if (err.no == 30) {
  msg <- paste("Parameter on line",line.no,"of file",filename,
               "is already defined.")
 } else if (err.no == 31) {
  msg <- paste("On line ",line.no," of file ",filename,
               ", parameters are defined for a function\n  ",
               "  that doesn't exist!",sep="")
 } else if (err.no == 32) {
  msg <- paste("On line ",line.no," of file ",filename,
               ", you haven't specified enough\n  information ",
               "to identify what you want.",sep="")
 } else if (err.no == 33) {
  msg <- paste("The nonlinear parameter definition on line ",
               line.no," of file\n  ",filename," is incorrect. ",
               "It refers either to a covariate that isn't",
               "\n  a main effect, or to a parameter that doesn't ",
               "exist.",sep="")
 } else if (err.no == 34) {
  msg <- paste("Invalid combination of codes on line ",line.no,
               " of file ",filename,sep="")
 } else if (err.no == 35) {
  msg <- paste("On line ",line.no," of file ",filename,", you ask",
               "for estimation of a\n  nonlinear parameter, but ",
               "the corresponding covariate has a\n  zero ",
               "coefficient!",sep="")
 } else if (err.no == 36) {
  msg <- paste("On line ",line.no," of file ",filename,", you've",
               "asked me to estimate\n  a parameter that should be ",
               "regarded as known.",sep="")
 } else if (err.no == 37) {
  msg <- paste("Nonlinear transformations have been requested in file ",
               filename,", but the parameters\n  of the transformations ",
               "have not all been defined.",sep="")
 } else if (err.no == 38) {
  msg <- paste("You have asked for orthogonal series representation ",
               "of a function over an interval\n  whose upper limit is ",
               "not greater than its lower limit (line ",line.no,
               " of file ",filename,")",sep="")
 } else if (err.no == 39) {
  msg <- paste("Site number ",line.no," has an attribute outside ",
               "the range you have specified\n  for orthogonal series ",    
               "representation in file ",filename,".",sep="")
 } else if (err.no == 40) {
#
#       This was an "inadequate storage" error - not needed now. 
#
 } else if (err.no == 41) {
  msg <- "Problem with specification of Box-Cox transforms - see message above"
 } else if (err.no == 42) {
  msg <- paste("Attempt to log a non-positive number when calculating",
              "covariates -\n  check model specification")
 } else if (err.no == 45) {
  msg <- paste("On line ",line.no," of file ",filename,", you specify",
               "a spatial structure\n  different from the one you ",
               "specified previously in the same file!",sep="")
 } else if (err.no == 46) {
  msg <- paste("In file ",filename,", you've asked for spatial dependence",
               "\nmodel ",line.no," - I don't have one of those!",sep="")
 } else if (err.no == 47) {
  msg <- paste("On line ",line.no," of file ",filename,", you specify ",
               "a spatial structure that is\n  inconsistent with the ", 
               "type of model you're using.",sep="")
 } else if (err.no == 48) {
  msg <- paste("You may not specify dispersion parameters in ",
               "normal-heteroscedastic models\n  (offence committed on ",
               "line ",line.no," of file ",filename,").",sep="")
 } else if (err.no == 49) {
  msg <- paste("You may not specify global values or spatial structures in ",
               "the\n  dispersion part of a model (offence committed on ",
               "line ",line.no," of file\n  ",filename,").",sep="")
 } else if (err.no == 50) {
  msg <- paste("On line ",line.no," of file ",filename,
               " you've specified a second\n  dispersion parameter",
               " - I can only deal with one.",sep="")
 } else if (err.no == 51) {
  msg <- paste("In file ",filename,", you haven't defined ",
               "a dispersion parameter.",sep="")
 } else if (err.no == 52) {
  msg <- paste("On line ",line.no," of file ",filename,
               ", you've specified a dispersion\n   parameter for a ",
               "model that doesn't have one.",sep="")
 } else if (err.no == 60) {
  msg <- paste("You've asked me to compute weighted averages of ",
               "previous days' values,\n  but you haven't given me any ",
               "parameters for the weighting scheme!",sep="")
# } else if (err.no == 65) {
#  msg <- paste("In file ",filename,", you haven't defined all the parameters",
#               "\n  needed for your chosen spatial dependence structure",sep="")
# } else if (err.no == 66) {
#  msg <- paste("In file ",filename,", you have asked for a constant inter-site",
#               "\n  correlation with magnitude greater than 1!",sep="")
 } else if (err.no == 67) {
  msg <- paste("Estimation of spatial correlation model failed. Ignore any",
               "\n  Fortran error above - perhaps there are no pairs of sites",
               "\n  with simultaneous observations? Examine file containing",
               "\n  (latent) correlations to see what's going on.",sep="")
 } else if (err.no == 68) {
  msg <- "Problem with spatial correlation model - see message above."
 } else if (err.no == 70) {
  msg <- "Missing external predictor data for part of requested period"
 } else if (err.no == 90) {
  msg <- paste("Data file contains non-binary values of response variable ",
                "for a logistic model.\n  Did you mean to set response.check=1?",
                sep="")
 } else if (err.no == 91) {
  msg <- paste("Data file contains non-positive values of response variable ",
                "for a gamma model.\n  Did you mean to set response.check=1?",
                sep="")
 } else if (err.no == 100) {
  msg <- paste("Current parameter estimates yield a fitted probability ",
                "of either 0 or 1 (to\n  machine precision). This isn't ",
                "going to converge - I'm bailing out.\n  Try different ",
                "starting values.",sep="")
 } else if (err.no == 101) {
  msg <- paste("Model seems singular (see above - row of zero diagonal is",
               "number of","\nparameter causing the problem). This",
               "could be due to a covariate that is\na linear combination",
               "of other covariates in the model, or to an\ninappropriately", 
               "chosen nonlinear parameter (e.g. an exponential decay",
               "rate\nthat is too large). Please check and try again.")
 } else if (err.no == 200) {
  msg <- paste("Problem with correlations in file ",filename,
               ".\n  See message above for details.",sep="")
 } else if (err.no == 210) {
  msg <- "Problem with external data - see message above for details."
 } else if (err.no == 300) {
  msg <- paste("Problem with model specification in file ",filename,
               ".\n  See message above for details.",sep="")
 } else if (err.no == 301) {
  msg <- paste("Problem with model specification in file ",filename,
               ". This relates\n  to external predictors - most likely",
               " the data files don't contain all\n  the predictors",
               " requested.",sep="")
 } else if (err.no == 400) {
  msg <- paste("\n  Simulation has generated a value outside the range ",
               "(-99.99,999.99):\n  this won't fit in the daily output ",
               "file. Possible causes include:",
               "\n\n  1. Dynamical instability caused by lagged values ",
               "of the response\n  variable. Consider different ",
               "transformations of the lagged responses,\n  and check ",
               "for interactions with other covariates that show ",
               "strong\n  time trends.\n  2. The model may be ",
               "overfitted, with some large coefficients\n  causing ",
               "instability. In this case, simplify the model and ",  
               "try again.","\n  3. If the problem ",
               "occurs while carrying out imputations, it\n  may ",
               "relate to an ill-conditioned residual correlation ",
               "matrix\n  associated with high correlations between ",
               "pairs of sites that are\n  close together. Consider ",
               "either removing some of these sites, or\n  using a ",
               "different correlation structure (including a nugget ",
               "term\n  may help).",
               "\n  4. Finally, it may be that the response ",
               "variable can naturally\n  take such values. If so, ",
               "change the measurement units (e.g.\n  divide ",
               "all values of the response by 10), refit the model\n  ",
               "and try again.",sep="")
 } else if (err.no == 997) {
  msg <- paste("Fortran storage inadequate - please contact Richard",
               "(r.chandler@ucl.ac.uk)")
 } else if (err.no == 998) {
  msg <- paste("\r        ")
 } else if (err.no == 999) {
  msg <- paste("**\n  ** Coding error found in package source - please", 
               "contact Richard\n  ** (r.chandler@ucl.ac.uk)\n  **")
 }
 msg
}
##############################################################################
##############################################################################
##############################################################################
model.num <- function(model.text,which.part) {
#
#       Returns the numeric code associated with a particular
#       model type (and, if the type is "normal-heteroscedastic",
#       whether it's the mean or dispersion component)
#
 if (model.text == "logistic") {
  z <- as.integer(1)
 } else if (model.text == "normal") {
  z <- as.integer(10)
 } else if (model.text == "gamma") {
  z <- as.integer(11)
 } else if (model.text == "normal-heteroscedastic") {
  if (which.part == "mean") {
   z <- as.integer(20)
  } else if (which.part == "dispersion") {
   z <- as.integer(21)
  } else {
   stop("which.part must be either 'mean' or 'dispersion'")
  }
 } else {
  stop(paste("Model <",model.text,"> not implemented yet",sep=""))
 }
 z
}
##############################################################################
##############################################################################
##############################################################################
OpenTempFile <- function() {
#
#       To open a temporary file in R, that can be accessed from Fortran
#       by passing an integer to the Fortran code. This is a bit clunky,
#       but it's the only way I can think of to get round the fact that
#       the passing of character strings (e.g. filenames) between R and 
#       Fortran was deprecated in R v3.6.1.
#
#       Arguments: none
#       Value: a list containing two components:
#         Name    - the name of the file that has been opened. This has
#                   the form RGLCtmp_####, where #### is .. 
#         Number: an integer between 1000 and 9999
#
 TmpFileAllocated <- FALSE
 TmpFlNo <- 1000
 while(!TmpFileAllocated) {
  TmpFlName <- paste("RGLCtmp_",TmpFlNo,sep="")
  if (!file.exists(TmpFlName)) {
   cat("### This file can safely be deleted unless Rglimclim is running ###\n", 
       file=TmpFlName)
   TmpFileAllocated <- TRUE
  } else {
   TmpFlNo <- TmpFlNo+1
   if (TmpFlNo > 9999) stop("Something's wrong: too many temporary files open")
  }
 }
 list(Name=TmpFlName, Number=as.integer(TmpFlNo))
}
##############################################################################
##############################################################################
##############################################################################
OpenFile <- function(filename,scratch=FALSE,Nfields,direct=TRUE,
                                                 formatted=FALSE) {
#
#       Finds an available Fortran I/O unit, connects it
#       and returns the unit number. Arguments:
#
#       filename        Name of an ASCII file to connect.
#       scratch         If TRUE, opens a Fortran scratch file
#                       (and in this case the argument filename, 
#                       if present, is ignored)
#       Nfields         The number of fields per record of 
#                       an unformatted scratch file
#       direct          Open scratch file for direct access
#       formatted       Open formatted scratch file
#
 if (scratch) { 
  acc.code <- ifelse(direct,1,2) # should be 1 for direct access, 2 for sequential
  form.code <- ifelse(formatted,1,2) # 1 for formatted, 2 for unformatted
  z <- .Fortran("SCROPN",ACODE=as.integer(acc.code),FCODE=as.integer(form.code),
                NFIELDS=as.integer(Nfields),FILNO=as.integer(0))$FILNO
 } else {
  LenFlNm <- as.integer(nchar(filename))
  TempFile <- OpenTempFile()
  cat(paste(filename,"\n",sep=""), file=TempFile$Name, append=TRUE)
  z <- .Fortran("RGLCFileConnect",TmpFlNo=as.integer(TempFile$Number),
                VarLen=LenFlNm,UnitNo=as.integer(0),Ifail=as.integer(0))
  unlink(TempFile$Name)
  if (z$Ifail == 0) z <- z$UnitNo else z <- -1
 }
 z
}
##############################################################################
##############################################################################
##############################################################################
 CloseFiles <- function(UnitNos,error=FALSE) {
#
#       Closes all Fortran I/O units
#
 .Fortran("CloseFiles",UnitNos=as.integer(UnitNos),
          Nunits=as.integer(length(UnitNos)),PACKAGE="Rglimclim")
 if (error) {
  stop(paste("Fortran got confused about which files were open, probably",
             "due to\n  an earlier error. All open connections have been",
             "closed - please try\n  your last command again."))
 }
 NULL
}
##############################################################################
##############################################################################
##############################################################################
set.maxlag <- function(max.lag,covcode,which.response,
                                     allow.incomplete.averages) {
#
#       Sets up the number of previous days' values required for each variable
#       in order to include a case when fitting models. Arguments:
#
#       max.lag         User-requested values - calculated from 
#                       nprev.required in GLCfit
#       covcode         The elements of the "covcode" part of a GLC.modeldef
#                       object that relate to lagged values
#       which.response  Which is the response variable? 
#       allow.incomplete.averages       See header for GLCfit 
#
 weighting.scheme <- floor((covcode / 1e4)) %% 10
 varnum <- floor(covcode / 1e6)
 lag.wanted <- covcode %% 1000
 tmp <- tapply(lag.wanted,INDEX=list(varnum),FUN=max)
#
#       If we're allowing incomplete averages then the only variables
#       for which we (possibly) need to set additional constraints are 
#       the response, and those for which we're doing some unweighted
#       estimation
#
 if (allow.incomplete.averages) {
  vars.to.constrain <- (names(tmp) == as.character(which.response)) |
                       tapply(weighting.scheme,INDEX=list(varnum),
                              FUN=function(x) {any(x==0)})
  tmp <- tmp[vars.to.constrain]
 }
 max.lag[as.numeric(names(tmp))] <- pmax(tmp,max.lag[as.numeric(names(tmp))])
 max.lag
}
##############################################################################
##############################################################################
##############################################################################
GLCfit <- function(model.def,dispersion.def,response.check=1,max.iter=Inf,
                   data.file="gaugvals.dat",missval=-99.99,which.response=1,
                   nprev.required,nvar.check=TRUE,
                   allow.incomplete.averages=FALSE,siteinfo,model.type,
                   external.files=c("yr_preds.dat","mn_preds.dat",
                   "dy_preds.dat"),diagnostics=2,verbosity=2,
                   resid.file,cor.file) {
#
#       Close any open Fortran file handles in case something went wrong
#       previously
#
 CloseFiles(10:1000)
#
#       Check site information
#
 z <- check.siteinfo(siteinfo,model.def$max.pars)
 nsites <- z$nsites
 nattr <- z$nattr
 sitinf <- z$sitinf
#
#       Check that modeldef is a valid model definition object
#
 if (class(model.def)[1] != "GLC.modeldef") {
  stop("model.def must be an object of class GLC.modeldef")
 }
 if(missing(model.type)) model.type <- model.def$model.type
#
#       Check that missval is numeric (cause of trouble for one user)
#
 if(!(is.numeric(missval))) stop(paste("missval must be numeric (currently ",
                                                class(missval),")",sep=""))
#
#       For heteroscedastic normal model, check first that a dispersion
#       model has been defined (either via the dispersion component of
#       model.def, or by explicitly supplying a dispersion.def argument); 
#       then cross-check the mean and dispersion components and ensure
#       that the components of model.def contain everything necessary
#       for both mean and dispersion models. 
#
 if (model.type == "normal-heteroscedastic") {
  if ( !inherits(model.def$dispersion,"GLC.modeldef") & 
       (missing(dispersion.def)) ) {
   stop("No dispersion structure defined in heteroscedastic normal model")
  } 
  if ( inherits(model.def$dispersion,"GLC.modeldef") & 
       !missing(dispersion.def) ) {
   stop(paste("Dispersion structure defined twice in heteroscedastic normal ",
              "model - \nonce in ",substitute(dispersion.def)," and once in ",
              substitute(model.def),"$dispersion.",sep=""))
  } 
  if (!missing(dispersion.def)) {
   if (!inherits(dispersion.def,"GLC.modeldef")) {
    stop("dispersion.def is not an object of class GLC.modeldef")
   }
   model.def$dispersion <- dispersion.def
  }
  model.def <- CheckJointModel(model.def)
 }
#
#       Now check files and allocate Fortran unit numbers. CheckFiles 
#       returns NULL if there's a problem
#
 UnitNos <- CheckFiles(data.file,resid.file,external.files,cor.file,
                      model.def,diagnostics,msg=(verbosity>0))
 if (is.null(UnitNos)) return(NULL)
 phi <- model.def$dispersion
 if (isTRUE(is.na(phi)) | !is.numeric(phi)) phi <- -1
 izero <- as.integer(0)
#
#       Figure out how many variables are in the data file, and check that
#       this matches the number of variable names in the model definition
#
 nvars <- (nchar(scan(data.file,what="",sep="\n",n=1,quiet=TRUE))-12)/6
 if (abs(nvars-round(nvars)) > 1e-6) stop(err.msg(data.file,err.no=5))
 if (nvar.check & abs(nvars-length(model.def$var.names)) > 1e-6) {
  warning(err.msg(data.file,err.no=6))
 }
#
#       Set up the mxlag argument to the Fortran code. Start with 
#       nprev.required, and then go through the max.lag component of
#       the model definition to see if there were any terms there
#       with higher lags. This also gives the opportunity to check 
#       that the model definition doesn't ask for more variables than
#       are present in the data file. 
#
 if (missing(nprev.required)) {
  CloseFiles(UnitNos)
  stop("You must specify a value for nprev.required") 
 } else if (length(nprev.required) == 1) {
  max.lag <- rep(nprev.required,nvars)
 } else if (length(nprev.required) == nvars) {
  max.lag <- nprev.required
 } else {
  CloseFiles(UnitNos)
  stop(paste("nprev.required is the wrong length (it should be either a ",
             "scalar, or\n  a vector with an entry for each variable in ",
             data.file,")",sep=""))
 }
 if (length(model.def$max.lag) > nvars) {
  CloseFiles(UnitNos)
  stop(paste("Model definition refers to more variables than are present\n",
             " in",data.file))
 }
 codes <- model.def$covcode[(model.def$Np[3]+1):model.def$Np[4]]
 codes <- codes[floor(codes / 1e6) > 0]
 max.lag <- set.maxlag(max.lag=max.lag,covcode=codes,
            which.response=which.response,
            allow.incomplete.averages=allow.incomplete.averages)
 max.lag <- as.integer(max.lag)
#
#     Need to know for what years are data available, to pass 
#     dimensions through to the Fortran code. A quick scan through
#     the data file is all that's needed, and is quick
#
 tmp <- as.numeric(substr(scan(data.file,what="",flush=TRUE,quiet=TRUE),1,4))
 year.range <- range(tmp)
#
#       Now open Fortran scratch files: one to store response and 
#       covariate information, the other to accumulate information
#       for residual analysis. Each record in the first needs a field
#       for each covariate; then one each for constant, response, case
#       weight, year, month and day. For the second we need Site, year, 
#       month, day, response, fitted and linear predictor. 
#
 UnitNos[90] <- OpenFile(scratch=TRUE,Nfields=6+max(model.def$Np[6:7]),
                         direct=TRUE,formatted=FALSE)
 if (UnitNos[90] < 0) CloseFiles(UnitNos[-90],error=TRUE)
 UnitNos[91] <- OpenFile(scratch=TRUE,Nfields=10,direct=TRUE,formatted=FALSE)
 if (UnitNos[91] < 0) CloseFiles(UnitNos[-91],error=TRUE)
#
#       Call the fitting routine. NB it is embedded within tryCatch() so
#       that file units can be disconnected if the user does Ctrl-C 
#       in the middle (this ensures that the "finally" argument is
#       executed in case this occurs). 
#
 if (model.type == "normal-heteroscedastic") {
  UnitNos.disp <- UnitNos
  UnitNos.disp[90] <- 
          OpenFile(scratch=TRUE,Nfields=6+max(model.def$dispersion$Np[6:7]),
                                                direct=TRUE,formatted=FALSE)
  if (UnitNos.disp[90] < 0) CloseFiles(UnitNos.disp[-90],error=TRUE)
  UnitNos.disp[91] <- 
          OpenFile(scratch=TRUE,Nfields=10,direct=TRUE,formatted=FALSE)
  if (UnitNos.disp[91] < 0) CloseFiles(UnitNos.disp[-91],error=TRUE)
  OpenUnits <- unique(c(UnitNos,UnitNos.disp))
  z <- tryCatch(JointFit.internal(model.type,max.iter,response.check,nsites,
                       model.def,nattr,sitinf, siteinfo$Site.codes, 
                       phi,year.range,diagnostics,UnitNos,UnitNos.disp,
                       data.file,nvars,missval,which.response,siteinfo,
                       max.lag,allow.incomplete.averages,verbosity),
                finally=CloseFiles(OpenUnits))
  CloseFiles(OpenUnits)
#
#       NB covcode potentially got modified by the Fortran code - reset
#       to avoid problems if the model is used again
#
  z$dispersion$covcode <- dispersion.def$covcode
 } else {
  z <- tryCatch(GLCfit.internal(model.type,max.iter,response.check,nsites,
                       model.def,nattr,sitinf, siteinfo$Site.codes, 
                       phi,year.range,diagnostics,UnitNos,TRUE,TRUE,
                       data.file,nvars,missval,which.response,siteinfo,
                       max.lag,allow.incomplete.averages,verbosity),
                finally=CloseFiles(UnitNos))
  CloseFiles(UnitNos)
 }
 z$covcode <- model.def$covcode
 z$var.names <- model.def$var.names
#
#       Also, ensure that parameter indices explicitly accommodate 
#       dispersion parameters (not required when reading model 
#       definitions)
 if (is.numeric(z$dispersion)) {
  if(isTRUE(z$dispersion > 0)) {
   z$Np[9] <- 1
   z$beta[z$df] <- z$dispersion
  }
 }
#
#       Finally, add file information, and indication of response
#       variable, to the object
#
 z$which.response <- which.response
 files.list <- list(Data=data.file,External=external.files)
 if ( (model.def$spmod >= 1) & (model.def$spmod <= 20) ) {
  files.list$Correlation <- cor.file
 } 
 if (diagnostics == 2) {
  files.list$Residuals <- resid.file
 }
 z$filenames <- files.list
 z
}
##############################################################################
##############################################################################
##############################################################################
GLCfit.internal <- function(model,max.iter,response.check,nsites,
                      model.def,nattr,sitinf, SiteCodes, 
                      phi,year.range,diagnostics,UnitNos,Read.needed,
                      Do.fit,data.file,nvars,missval,which.response,siteinfo,
                      max.lag,allow.incomplete.averages,verbosity,
                      DateRange=c(-1,-1)) {
##############################################################################
#       Workhorse for GLCfit(). Not for user consumption. Arguments:
#
#       model           Code for model being fitted. 
#       max.iter        Maximum number of IWLS iterations
#       response.check  Passed through from GLCfit()
#       nsites          Number of sites
#       model.def       GLC model definition object
#       nattr           Number of site attributes defined
#       sitinf          Array containing site information
#       SiteCodes       Vector of 4-character site codes
#       phi             Dispersion parameter
#       year.range      Range of years present in data file
#       diagnostics     Passed through from GLCfit()
#       UnitNos         Array of Fortran unit numbers
#       Read.needed     Indicates whether the Fortran routine needs to 
#                       read the data from data.file (TRUE), or whether
#                       the required scratch file has already been written
#                       (FALSE)
#       Do.fit          Indicates whether to compute parameter estimates
#                       (TRUE) or not (FALSE) - the latter being relevant
#                       if all that is required is to read the data and
#                       produce a scratch file, or to perform residual
#                       analysis for an existing model.
#       data.file       Name of input data file
#       which.response  Which of the variables in the data file should be 
#                       taken as the response
#       siteinfo        Passed through from GLCfit()
#       max.lag         Only fit to observations with this many previous
#                       days' observations
#       verbosity       Controls verbosity of output (larger values: more 
#                       verbose, 0=work silently)
#       DateRange       A 2-element vector with elements that, if
#                       non-negative, should be of the form of YYYYMMDD;
#                       in this case, only observations between these
#                       two dates inclusive will be used for fitting
#                       (this is for use when fitting joint mean-variance
#                       models where the subset of observations for
#                       which covariate information is available might
#                       differ between the mean and dispersion models)
##############################################################################
 izero <- as.integer(0)
 nyears <- as.integer(diff(year.range)+1)
 model.code <- model.num(model)
 nobs <- as.integer(ifelse(is.null(model.def$nobs),0,model.def$nobs))
#
#       Ensure dimensioning of pwtidx matches that in Fortran 
#       declarations
#
 pwtidx <- expand.pwtidx(model.def,nvars)
#
#       Pass everything over to Fortran. For versions of R from 3.6.1 onwards, 
#       the only way to get character info from R into Fortran is to write to 
#       a temporary file and then read it back in from the Fortran code - VERY 
#       annoying. 
#
 TempFile <- OpenTempFile()
 write(SiteCodes, file=TempFile$Name, append=TRUE)
 z <- .Fortran("glmfit",model=model.code,
          maxit=as.integer(ifelse(is.infinite(max.iter),-1,max.iter)),
          ycheck=as.integer(response.check),nsites=as.integer(nsites),
          nvars=as.integer(nvars),missval=missval,
          RespIdx=as.integer(which.response),
          AllowIncAvge=as.integer(allow.incomplete.averages),
          mxp=model.def$max.pars,nattr=as.integer(nattr),sitinf=sitinf,
          TmpFlNo=TempFile$Number, NP=as.integer(model.def$Np), 
          beta=model.def$beta,covcode=model.def$covcode,two=model.def$twoway,
          three=model.def$threeway,theta=model.def$theta,
          CovNaive=diag(rep(0,model.def$max.pars+1)),
          CovRobust=diag(rep(0,model.def$max.pars+1)),LogLik=0,Deviance=0,
          nobs=nobs,p=izero,df=izero,sitxfm=model.def$sitxfm,
          fouidx=model.def$fourier.idx,legidx=model.def$legendre.idx,
          mxlag=as.integer(max.lag),pwtidx=model.def$prevwt.idx,
          phi=phi,ybar=0,sy=0,errbar=0,mse=0,Rsq=0,prbar=0,spr=0,sepr=0,
          arbaro=0,saro=0,arbare=0,sare=0,bypred=matrix(0,nrow=3,ncol=10),
          prmnth=matrix(0,nrow=12,ncol=4),pryear=matrix(0,nrow=nyears,ncol=4),
          prsite=matrix(0,nrow=nsites,ncol=4),nyears=nyears,
          firstyr=as.integer(year.range[1]),spmod=model.def$spmod,
          rho=model.def$rho,glbcod=model.def$global.codes,
          glbval=model.def$global.vals,DoFit=as.integer(Do.fit),
          DoResid=as.integer(diagnostics),UnitNos=as.integer(UnitNos),
          NeedRead=as.integer(Read.needed),DateRange=as.integer(DateRange),
          verbosity=as.integer(verbosity),ifail=izero,PACKAGE="Rglimclim",
          NAOK=TRUE)
 unlink(TempFile$Name)
#
#       Check and assemble the results. NB the degrees of freedom 
#	count *all* parameters estimated, whereas the Fortran routines
#	don't count the constant terms - hence add 1 to df in the next
#	command
#
 if (z$ifail != 0) {
  stop(err.msg(data.file,err.no=z$ifail),call.=FALSE)
 }
 zz <- list(model.type=model,Np=z$NP,iext=model.def$iext,
            model.title=model.def$model.title,max.pars=z$mxp,beta=z$beta,
            covcode=z$covcode,twoway=z$two,threeway=z$three,theta=z$theta,
            cov.naive=z$CovNaive,cov.robust=z$CovRobust,nobs=z$nobs,
            df=z$p+1,df.resid=z$df,LogLik=z$LogLik,Deviance=z$Deviance,
            nlstat=model.def$nlstat,sitxfm=z$sitxfm,fourier.idx=z$fouidx,
            legendre.idx=z$legidx,max.lag=z$mxlag,prevwt.idx=z$pwtidx,
            weighting.scheme=model.def$weighting.scheme,
            dispersion=z$phi,spmod=z$spmod,rho=z$rho,global.codes=z$glbcod,
            global.vals=z$glbval,model.labels=model.def$model.labels,
            siteinfo=siteinfo,DateRange=z$DateRange)
 if (isTRUE(zz$dispersion >= 0)) zz$df <- zz$df+1
 zz$df.resid <- zz$nobs-zz$df
 if (!(diagnostics == 0)) {
  zz$Residuals <- list(ybar=z$ybar,sy=z$sy,errbar=z$errbar,MSE=z$mse,Rsq=z$Rsq,
                       Pearson=list(Mean=z$prbar,SD=z$spr,SEmean=z$sepr,
                                    Month.table=as.data.frame(z$prmnth),
                                    Site.table=as.data.frame(z$prsite),
                                    Year.tables=as.data.frame(z$pryear)),
                       Anscombe=list(Obs.mean=z$arbaro,Obs.sd=z$saro,
                                     Exp.mean=z$arbare,Exp.sd=z$sare),
                       Prob.table=z$bypred)
  tmp <- zz$Residuals$Pearson
#
#       NB the Fortran code skips over groups with no observations, but
#       initialises all its arrays to zero - hence standard errors come
#       back incorrectly as zero for these groups. This needs to be 
#       corrected here.  
#       
  tmp$Month.table[tmp$Month.table[,1] == 0,2:4] <- NA
  tmp$Year.table[tmp$Year.table[,1] == 0,2:4] <- NA
  tmp$Site.table[tmp$Site.table[,1] == 0,2:4] <- NA
  rownames(tmp$Month.table) <- month.abb
  names(tmp$Month.table) <- names(tmp$Year.table) <- 
                            c("# of days","Mean","Std Dev","S.E. mean")
  tmp$Site.table <- cbind(tmp$Site.table[,1],siteinfo$Site.names,
                          tmp$Site.table[,-1])
  rownames(tmp$Site.table) <- siteinfo$Site.codes
  colnames(tmp$Site.table) <- c("# of days","Name","Mean","Std Dev","S.E. mean")
  rownames(tmp$Year.table) <- year.range[1]:year.range[2]
  zz$Residuals$Pearson <- tmp
  if (!is.null(zz$Residuals$Prob.table)) {
   rownames(zz$Residuals$Prob.table) <- 
      c("Observed prop.","Expected prop.","Number of cases")
   colnames(zz$Residuals$Prob.table) <- 1:10
  }
 }
 class(zz) <- "GLC.modeldef"
 zz
}
##############################################################################
##############################################################################
##############################################################################
JointFit.internal <- function(model.type,max.iter,response.check,nsites,
                       model.def,nattr,sitinf, SiteCodes, 
                       phi,year.range,diagnostics,UnitNos,UnitNos.disp,
                       data.file,nvars,missval,which.response,siteinfo,
                       max.lag,allow.incomplete.averages,verbosity) {
 mean.def <- model.def; mean.def$model.type <- "normal"
 disp.def <- model.def$dispersion; disp.def$model.type <- "gamma"
#
#       Set up scratch file for dispersion model (at this point it 
#       will have the wrong values of the response). This will also
#       set up things like degrees of freedom, which are required below.
#       Need to check for the possibility that the mean and dispersion
#       models have different external predictors in and hence need
#       to ensure that they're fitted to the same set of cases.
#       Also do *not* check that response values are strictly positive
#       here (they're just dummy values at this stage)
#
 izero <- as.integer(0)
 disp.def <- GLCfit.internal(disp.def$model.type,max.iter=0,response.check=-1,
                      nsites,disp.def,nattr,sitinf, SiteCodes, 
                      phi,year.range,diagnostics=0,UnitNos.disp,
                      Read.needed=TRUE,Do.fit=FALSE,data.file,
                      nvars,missval,which.response,siteinfo,max.lag,
                      allow.incomplete.averages,verbosity)
#
#       Preliminary fit of mean model with equal weights
# 
 if (verbosity > 0) cat("### Preliminary fit of mean model ###\n")
 mean.def <- GLCfit.internal(mean.def$model.type,max.iter,response.check,nsites,
                      mean.def,nattr,sitinf, SiteCodes, 
                      phi,year.range,diagnostics=1,UnitNos,
                      Read.needed=TRUE,Do.fit=TRUE,data.file,
                      nvars,missval,which.response,siteinfo,max.lag,
                      allow.incomplete.averages,verbosity,
                      DateRange=disp.def$DateRange)
 if (!isTRUE(all.equal(mean.def$DateRange,disp.def$DateRange))) {
  disp.def <- GLCfit.internal(disp.def$model.type,max.iter,response.check,nsites,
                       disp.def,nattr,sitinf, SiteCodes, 
                       phi,year.range,diagnostics=0,UnitNos.disp,
                       Read.needed=TRUE,Do.fit=FALSE,data.file,
                       nvars,missval,which.response,siteinfo,max.lag,
                       allow.incomplete.averages,verbosity,
                       DateRange=mean.def$DateRange) 
 }
 if (mean.def$nobs != disp.def$nobs) {
#
#       If the configuration of external covariates means that there are
#       different numbers of cases available for the mean and dispersion
#       models, sort it out ... NB the "rbind" in the Fortran call is 
#       correct since the rows and columns in R and Fortran seem to get
#       interchanged. Also, note that the "residual" file for the 
#       dispersion model is passed as Unit 0 because this doesn't exist. 
#
  z <- .Fortran("ScratchSync",UnitNos=rbind(as.integer(UnitNos[90:91]),
                 as.integer(c(UnitNos.disp[90],0))),
                 Nobs=as.integer(c(mean.def$nobs,disp.def$nobs)),
                 P=as.integer(c(mean.def$df,disp.def$df)),Ifail=izero,
                 PACKAGE="Rglimclim")
  if (z$Ifail != 0) stop(err.msg(err.no=z),call.=FALSE)
  mean.def$nobs <- z$Nobs[1]
  disp.def$nobs <- z$Nobs[2]
 }
 OldLogL <- mean.def$LogLik
 NewLogL <- -Inf
 while (abs(NewLogL-OldLogL) > 1e-4) {
  OldLogL <- NewLogL
#
#       Define the squared residuals from the mean model as the responses
#       for the dispersion, and fit
#
  z <- .Fortran("ScratchUpdate",FromUnit=as.integer(UnitNos[91]),
                 ToUnit=as.integer(UnitNos.disp[90]),
                 UpdateType=as.integer(1),Nobs=as.integer(mean.def$nobs),
                 P=as.integer(disp.def$df),Ifail=izero,
                 PACKAGE="Rglimclim")$Ifail
  if (z != 0) stop(err.msg(err.no=z),call.=FALSE)
  if (verbosity > 0) cat("### Fitting / updating dispersion model ###\n")
  disp.def <- GLCfit.internal(disp.def$model.type,max.iter,response.check=0,
                       nsites,disp.def,nattr,sitinf, SiteCodes, 
                       phi,year.range,diagnostics=0,
                       UnitNos.disp,Read.needed=FALSE,Do.fit=TRUE,data.file,
                       nvars,missval,which.response,siteinfo,max.lag,
                       allow.incomplete.averages,verbosity)
#
#       Write the inverse of the fitted variances as case weights in the
#       mean model, and refit
#
  z <- .Fortran("ScratchUpdate",FromUnit=as.integer(UnitNos.disp[91]),
                 ToUnit=as.integer(UnitNos[90]),
                 UpdateType=as.integer(2),Nobs=as.integer(disp.def$nobs),
                 P=as.integer(mean.def$df),Ifail=izero,
                 PACKAGE="Rglimclim")$Ifail
  if (z != 0) stop(err.msg(err.no=z),call.=FALSE)
  if (verbosity > 0) cat("### Updating mean model ###\n")
  mean.def <- GLCfit.internal(mean.def$model.type,max.iter,response.check=0,
                       nsites,mean.def,nattr,sitinf, SiteCodes, 
                       phi,year.range,diagnostics,UnitNos,
                       Read.needed=FALSE,Do.fit=TRUE,data.file,
                       nvars,missval,which.response,siteinfo,max.lag,
                       allow.incomplete.averages,verbosity)
  NewLogL <- mean.def$LogLik
 }
#
#	Degrees of freedom: put them all in both models (and don't count the
#	dispersion parameter estimates, which *are* counted by GLCfit.internal)
#
 mean.def$df <- mean.def$df - 1
 disp.def$df <- disp.def$df - 1
 mean.def$df <- mean.def$df + disp.def$df
 mean.def$df.resid <- mean.def$nobs - mean.def$df
 disp.def$df.resid <- mean.def$df.resid
#
#       Get rid of components from disp.def that are likely to be 
#       confusing or redundant; put it in the dispersion component of
#       mean.def; and return
#
 els.to.drop <- names(disp.def) %in% c("LogLik","Deviance","dispersion",
                                       "spmod","rho","global.codes",
                                       "global.vals","siteinfo")
 disp.def <- disp.def[!els.to.drop]
 class(disp.def) <- "GLC.modeldef"
 mean.def$dispersion <- disp.def
 mean.def$model.type <- model.type
 mean.def
}
##############################################################################
##############################################################################
##############################################################################
CheckFiles <- function(data.file,resid.file,external.files,cor.file,
                       model.def,diagnostics,fitting=TRUE,which.model=1,
                       msg=TRUE) {
##############################################################################
#       To check files and allocate Fortran unit numbers. All arguments
#       except for fitting, msg and which.model are passed through from 
#       GLCfit() - see header there for details. fitting is used to 
#       figure out whether we need to check that correlation files can 
#       be overwritten (for fitting routines) or must be present (for
#       simulation). which.model is always 1 except when defining the 
#       second model in a simulation setup that requires a pair of 
#       models with different correlation structures( e.g. simulating
#       precipitation with occurrence and intensity models). msg just 
#       indicates whether or not to print "Checking files" when the 
#       routine is called. 
#
#       The function returns a vector of length 100, containing the unit
#       numbers that have been allocated. The elements of this vector are
#       as follows:
#
#       1     Input file containing daily time series data
#       2     Output data file for residual analysis (fitting routines)
#       3     Input file containing external annual predictors
#       4     Input file containing external monthly predictors
#       5     Input file containing external daily predictors
#       6     File containing spatial correlations. This is output by 
#             the fitting routines, and taken as input when simulating
#       7     File containing a second set of spatial correlations for
#             use in simulation settings where two models are used. This 
#             gets allocated if which.model=2. 
#      20     Daily output file for simulated data.
#      21     Monthly output file for simulated data.
#      90     Scratch file for covariate and response information (fitting)
#      91     Scratch file for accumulating residual information (fitting)
#      95     Scratch file for storing all daily values (simulation)
#
#       NB the CloseFiles() in the return() options are needed to close
#       any units that have been opened in the interim. The CloseFiles()
#       function returns NULL. 
##############################################################################
 if (msg) cat("\nChecking files ...\n")
 UnitNos <- rep(0,100)
 if (!file.OK(data.file,"old")) return(CloseFiles(UnitNos))
 UnitNos[1] <- OpenFile(data.file) 
 if (diagnostics > 1) {
  if (missing("resid.file")) {
   CloseFiles(UnitNos)
   stop("In GLCfit() you must specify a value for resid.file if diagnostics==2",
        call.=FALSE)
  }
  if (!file.OK(resid.file,"new")) return(CloseFiles(UnitNos))
  UnitNos[2] <- OpenFile(resid.file) 
 } else {
  if (!missing("resid.file")) {
   warning("In GLCfit(), arguments resid.file and diagnostics are incompatible",
           call.=FALSE)
  }
 }
 if (!missing(model.def)) {
  for (i in 1:3) {
   if (model.def$iext[i] != 0) { 
    if (!file.OK(external.files[i],"old")) return(CloseFiles(UnitNos))
    UnitNos[2+i] <- OpenFile(external.files[i]) 
   }
  }
  if ( (model.def$spmod >= 1) & (model.def$spmod <= 20) ) {
   if (fitting) {
    if (!file.OK(cor.file,"new")) return(CloseFiles(UnitNos))
   } else {
#
#       For simulation, only need the correlation file if empirical 
#       inter-site correlations are being used
#
    if (model.def$spmod == 1) {
     if(!file.OK(cor.file,"old")) return(CloseFiles(UnitNos))
    }
   }
   if (fitting | model.def$spmod ==1) UnitNos[5+which.model] <- OpenFile(cor.file) 
  }
 }
 UnitNos
}
##############################################################################
##############################################################################
###s###########################################################################
CheckJointModel <- function(model.def) {
#
#       To check the mean and dispersion parts of a joint 
#       mean-dispersion model (the mean part being held in model.def
#       and the dispersion part in model.def$dispersion); and update
#       the mean part so that all necessary files etc. get opened
#       while fitting it. 
#
 mean.def <- model.def
 disp.def <- model.def$dispersion
#
#       Check that the two models have identical siteinfo components
#
 if (!identical(mean.def$siteinfo,disp.def$siteinfo)) {
  stop(paste("In GLCfit(), mean and dispersion models have",
             "different siteinfo components"),.Call=FALSE)
 }
#
#       Check that limits for orthogonal series representation of
#       site effects are the same in the two models. The logic 
#       below is that these limits correspond to elements of theta
#       that are less than 1e9 in magnitude (this being the missing value
#       code for theta) and correspond to site effects for which there
#       is either a Legendre or Fourier index entry ...
#
 mean.ortab <- cbind(mean.def$theta,mean.def$covcode,mean.def$legendre.idx,
                     mean.def$fourier.idx)[1:mean.def$Np[1],,drop=FALSE]
 mean.ortab <- 
   mean.ortab[rowMeans(abs(mean.ortab[,1:3,drop=FALSE]) < 1e9)>0,,drop=FALSE]
 mean.ortab <- 
   mean.ortab[rowMeans(mean.ortab[,5:6,drop=FALSE]) > 0,,drop=FALSE]
 disp.ortab <- cbind(disp.def$theta,disp.def$covcode,disp.def$legendre.idx,
                     disp.def$fourier.idx)[1:disp.def$Np[1],,drop=FALSE]
 disp.ortab <- 
   disp.ortab[rowMeans(abs(disp.ortab[,1:3,drop=FALSE]) < 1e9)>0,,drop=FALSE]
 disp.ortab <- 
   disp.ortab[rowMeans(disp.ortab[,5:6,drop=FALSE]) > 0,,drop=FALSE]
#
#       ... and for which the same covariate is present in both
#       models ...
# 
 shared.seffs <- intersect(mean.ortab[,4],disp.ortab[,4])
 mean.ortab <- mean.ortab[mean.ortab[,4] %in% shared.seffs,,drop=FALSE]
 disp.ortab <- disp.ortab[disp.ortab[,4] %in% shared.seffs,,drop=FALSE]
 
 for (covcode in shared.seffs) {
  mean.row <- (mean.ortab[,4] == covcode) & (mean.ortab[,5] > 0) # Legendre
  disp.row <- (disp.ortab[,4] == covcode) & (disp.ortab[,5] > 0) 
  if ( (sum(mean.row) > 0) & (sum(disp.row) > 0) ) {
   if (!identical(mean.ortab[mean.row,1:3],disp.ortab[disp.row,1:3])) {
    stop(paste("In GLCfit(). mean and dispersion models have different limits", 
               "for Legendre\npolynomial representation of the same site",
               "effects."),.Call=FALSE)
   }
  }
  mean.row <- (mean.ortab[,4] == covcode) & (mean.ortab[,6] > 0) # Fourier
  disp.row <- (disp.ortab[,4] == covcode) & (disp.ortab[,6] > 0) 
  if ( (sum(mean.row) > 0) & (sum(disp.row) > 0) ) {
   if (!identical(mean.ortab[mean.row,1:3],disp.ortab[disp.row,1:3])) {
    stop(paste("In GLCfit(), mean and dispersion models have different limits",
               "for Fourier\nrepresentation of the same site effects."),
               .Call=FALSE)
   }
  }
 }
#
#       Ensure that all external files are flagged in main model 
#       definition
#
 mean.def$iext <- pmax(mean.def$iext,disp.def$iext)
#
#       That max.lag is the same for both models
#
 mean.def$max.lag <- pmax(mean.def$max.lag,disp.def$max.lag)
 disp.def$max.lag <- mean.def$max.lag
 mean.def$dispersion <- disp.def
 mean.def
}
##############################################################################
##############################################################################
##############################################################################
expand.pwtidx <- function(model,nvars) {
#
#       To expand an array of indices relating to weighted averages of
#       lagged values, from the size currently held in an R object
#       to the dimensions required by the underlying Fortran code
#
 pwtdim <- dim(model$prevwt.idx)
 z <- array(dim=pmax(pwtdim,c(0,0,nvars)))
 z[1:pwtdim[1],1:pwtdim[2],1:pwtdim[3]] <- model$prevwt.idx
 storage.mode(z) <- "integer"
 z
}
##############################################################################
##############################################################################
##############################################################################
read.siteinfo <- function(site.file,nhead=38,
                                    which.coords=NULL,regions=NULL) {
#
#       To read a site information file and produce a siteinfo
#       object suitable for passing to GLCfit
#
 if (!(site.file %in% list.files())) {
  stop(paste("Input file",site.file,"doesn't exist."))
 }
 if (!is.null(which.coords) & length(which.coords) > 2) {
  stop(paste("At most two attributes can be used to define ",
             "inter-site distances\n(supplied argument which.coords ",
             "has length",length(which.coords),")",sep=""))
 }
#
#       Read all lines except header into a single character vector
#
 file.input <- scan(site.file,what=character(),sep="\n",
                                            quiet=TRUE,skip=nhead)
#
#       Read attribute descriptions, remove these from the 
#       input vector and order so that the spatial co-ordinates
#       are the first two attributes
#
 nattr <- as.numeric(file.input[1])
 if (is.null(which.coords) & nattr > 0) {
  which.coords <- 1:min(nattr,2)
 } else if (!is.null(which.coords) & 
            length(which.coords) > length(nattr)) {
  warning(paste("which.coords has more elements than # of site",
                "attributes defined"))
 }
 file.input <- file.input[-1]
 if (length(file.input) < nattr) {
  stop(paste("End of input file ",site.file," before ",
             "specified # (",nhead,")\nof site attributes ",
             "was found"),sep="")
 }
 attrtxt <- file.input[1:nattr]
 attrtxt <- c(attrtxt[which.coords],attrtxt[-which.coords])
 file.input <- file.input[-(1:(nattr+1))] 
#
#       Read site names and attributes - there are two rows in the
#       file per site.
#
 nsites <- length(file.input) / 2
 if (nsites < 1) stop(paste("No sites found in file",site.file))
 if (nsites != round(nsites)) {
  stop(paste("Incorrect number of rows in site definition section of",
             "file",site.file))
 }
 sitetxt <- file.input[2*(1:nsites)-1]
 file.input <- file.input[2*(1:nsites)]
 scodes <- substr(file.input,1,4)
 regcodes <- as.numeric(substr(file.input,5,9))
 regcodes[is.na(regcodes)] <- 0
 sitinf <- matrix(nrow=nsites,ncol=nattr)
 rownames(sitinf) <- scodes
 colnames(sitinf) <- substr(attrtxt,1,16)
 if (any(nchar(attrtxt) > 16) & any(duplicated(colnames(sitinf)))) {
  warning(paste("Some attribute text was truncated when making column names",
            "\n  for 'Attribute.values' component, and some column names are",
            "\n  duplicated as a result. To make this warning go away, change",
            "\n  the attribute text in file",site.file)) 
 }
 for (i in 1:nattr) {
  sitinf[,i] <- as.numeric(substr(file.input,10*i,(10*i)+9))
 }
#
#       Set the co-ordinates as the first two attributes
#
 siteinf <- cbind(sitinf[,which.coords],sitinf[,-which.coords]) 
#
#       And finally, make some region names
#
 regnames <- make.regnames(max(regcodes),regions)
 
 if (regnames$warn) warning(paste("More regions found in file",
                                  site.file,"than in 'regions' argument"))
 z <- list(Nsites=nsites,Site.names=sitetxt,Site.codes=scodes,
      Regions=regcodes,Region.names=regnames$names,Nattr=nattr,
      Attribute.names=attrtxt,Attribute.values=sitinf)
 class(z) <- c("siteinfo","list")
 z
}
##############################################################################
##############################################################################
##############################################################################
make.siteinfo <- function(site.data,coord.cols=NULL,site.codes=NULL,
                          site.names=NULL,region.col=NULL,
                          attr.cols=NULL,attr.names=NULL,regions=NULL) {
#
#       To extract site information from a data frame and produce 
#       a siteinfo object suitable for passing to GLCfit
#
 if (!is.data.frame(site.data)) stop("site.data must be a data frame")
 nsites <- nrow(site.data)
#
#       Extract region codes
#
 if (is.null(region.col)) {
  regcodes <- rep(0,nsites)
  name.cols <- numeric(0)
 } else {
  if (length(region.col) > 1) stop("region.col must contain a a single column number")
  regcodes <- site.data[,region.col]
  if (max(regcodes) > nsites) {
   stop(paste("Highest region number (",max(regcodes),") exceeds number of sites (",
              nsites,"):\n  this will cause Fortran memory allocation problems.",
              sep=""))
  }
  name.cols <- region.col
 }
#
#       Define 4-character site codes. Default is the row names
#       of site.data. Could also be a column of site.data, or 
#       a character vector provided separately. Apart from this, 
#       don't make any automatic decisions here - they need to be 
#       right so that we can identify cases in the data file later 
#       on. 
#
 if (is.null(site.codes)) {
  scodes <- rownames(site.data)
 } else if (is.numeric(site.codes)) {
  scodes <- as.character(site.data[,site.codes[1]])
  if (length(site.codes) > 1) {
   warning(">1 column of site codes specified - only the first was used")
  }
  name.cols <- c(name.cols,site.codes)
 } else if (is.character(site.codes)) {
  if (length(site.codes) != nsites) {
   stop("Length of site.codes doesn't match number of sites")
  }
  scodes <- site.codes
 } else {
  stop("Invalid value of site.codes")
 }
 if (!all(nchar(scodes) == 4)) {
  errmsg <- 
    "Site codes must contain exactly 4 characters including spaces"
  if (is.null(site.codes)) errmsg <- paste(errmsg,
     "\n(default value of rownames(site.data) was requested)")
  stop(errmsg)
 }
#
#       Now site names
#
 if (is.null(site.names)) {
  stop("site.names must be specified")
 } else if (is.numeric(site.names)) {
  sitetxt <- as.character(site.data[,site.names[1]])
  if (length(site.names) > 1) {
   warning(">1 column of site names specified - only the first was used")
  }
  name.cols <- c(name.cols,site.names)
 } else if (is.character(site.names)) {
  if (length(site.names) != nsites) {
   stop("Length of site.names doesn't match number of sites")
  }
  sitetxt <- site.names
 } else {
  stop("Invalid value of site.names")
 }
#
#       Finally, the attributes and associated text. Default
#       is to take all columns of site.data that are *not* 
#       site names or codes - but start by taking *all* 
#       columns so as to assist with locating columns containing
#       co-ordinates
#
 if (!is.null(attr.cols) & any(name.cols %in% attr.cols)) {
  stop(paste("Some columns in attr.cols have been specified as",
             "containing\n  site information"))
 }
 if (is.null(coord.cols)) {
  which.coords <- numeric(0)
 } else if (length(coord.cols) > 2) {
  stop(paste("At most two attributes can be used to define ",
             "inter-site distances\n(supplied argument coord.cols ",
             "has length",length(coord.cols),")",sep=""))
 } else {
  which.coords <- coord.cols
 }
 if (any(name.cols %in% which.coords)) {
  stop(paste("Some columns in coord.cols have been specified as",
             "containing\n  site information"))
 }
 if (is.null(attr.cols)) {
  which.attr <- 1:ncol(site.data)
 } else {
  which.attr <- attr.cols
 }
 cols.to.drop <- match(name.cols,which.attr)
 cols.to.drop <- cols.to.drop[!is.na(cols.to.drop)]
 if (length(cols.to.drop) > 0) {
  which.attr <- which.attr[-cols.to.drop]
 }
 if (!all(which.coords %in% which.attr)) {
  stop(paste("coord.cols contains columns that are not specified",
             "in attr.cols"))
 }
 sitinf <- as.matrix(site.data[,which.attr])
 mode(sitinf) <- "numeric" # In case there are integer columns which screw up the Fortran
 if (is.null(attr.names)) {
  attrtxt <- colnames(site.data)[which.attr]
 } else {
  attrtxt <- attr.names
  if (length(attrtxt) != ncol(sitinf)) {
   stop("Length of attr.names doesn't match number of site attributes")
  }
 }
#
#       Put the coordinate columns first - if the user didn't
#       define any coordinate columns, use the given attributes
#       if there are at most two of them, but otherwise stop
#       with an error
#
 nattr <- ncol(sitinf)
 if (is.null(coord.cols)) {
  if (nattr <= 2) {
   which.coords <- which.attr
  } else {
   stop(paste("If more than two site attributes have been defined, you",
              "must specify\n  which ones are spatial coordinates"))
  }
 }
 if (length(which.coords) > 0) {
  first.cols <- match(which.coords,which.attr) 
 } else {
  first.cols <- 1
 }
 sitinf <- cbind(sitinf[,first.cols,drop=FALSE],sitinf[,-first.cols,drop=FALSE])
 attrtxt <- c(attrtxt[first.cols],attrtxt[-first.cols])
 rownames(sitinf) <- scodes
 colnames(sitinf) <- substr(attrtxt,1,16)
 if (any(nchar(attrtxt) > 16) & any(duplicated(colnames(sitinf)))) {
  if (is.null(attr.names)) {
   warning(paste("Some attribute text was truncated when making column names",
             "\n  for 'Attribute.values' component, and some column names are",
             "\n  duplicated as a result. To make this warning go away, change",
             "\n  the column names in site.data")) 
  } else {
   warning(paste("Some attribute text was truncated when making column names",
             "\n  for 'Attribute.values' component, and some column names are",
             "\n  duplicated as a result. To make this warning go away, change",
             "\n  the attribute text in attr.names")) 
  }
 }
#
#       And finally, make some region names
#
 regnames <- make.regnames(max(regcodes),regions)
 if (regnames$warn) warning(paste("Site information and 'regions' argument",
                                  " suggest different numbers of regions"))
 z <- list(Nsites=nsites,Site.names=as.character(sitetxt),
      Site.codes=as.character(scodes),Regions=regcodes,
      Region.names=regnames$names,Nattr=nattr,
      Attribute.names=as.character(attrtxt),Attribute.values=sitinf)
 class(z) <- c("siteinfo","list")
 z
}
##############################################################################
##############################################################################
##############################################################################
print.siteinfo <- function(x,...) {
#
#       To display the attributes that have been defined in a siteinfo
#       object
#
 if (x$Nsites == 1) {
  cat("1 site defined, with the following attributes:\n\n")
 } else {
  cat(paste(x$Nsites,"sites defined, each with the following attributes:\n\n"))
 }
 cat(paste(paste(1:x$Nattr,". ",x$Attribute.names,sep=""),collapse="\n"))
 cat("\n")
 invisible(NULL)
}
##############################################################################
##############################################################################
##############################################################################
make.regnames <- function(nregs,regions) {
#
#       To make some region names to go into a siteinfo object. If the
#       regions argument is non-null, it is taken to be a "regions"
#       object; otherwise some default names are created. nregs is 
#       the minimum number of regions required. Names are created
#	up to the larger of nregs or the largest number defined in 
#	the regions argument. 
#
 warn <- FALSE
 maxreg <- nregs
 if (!is.null(regions)) {
  if (is.numeric(regions$Region)) maxreg <- max(nregs,regions$Region)
 }
#
#	Default names
#
 if (maxreg > 0) {
  region.names <- c("Whole area",paste("Region",1:maxreg))
 } else {
  region.names <- "Whole area"
 }
 if (!is.null(regions)) {
  warn <- !(maxreg==nregs & maxreg==length(unique(regions$Region))-1 )
  region.names[1+regions$Region] <- as.character(regions$Name)
 }
 invisible(list(names=region.names,warn=warn))
}
##############################################################################
##############################################################################
##############################################################################
read.regiondef <- function(region.file,nhead=19) {
#
#       To read a region definition file and produce a reginfo
#       object
#
 if (!(region.file %in% list.files())) {
  stop(paste("Input file",region.file,"doesn't exist."))
 }
#
#       Read all lines except header into a single character vector.
#       NB we're expecting some text for region 0 - the whole area -
#       so the number of subregions is 1 less than the number of 
#       lines read.
#
 file.input <- scan(region.file,what=character(),sep="\n",
                                            quiet=TRUE,skip=nhead)
 nregs <- length(file.input) - 1
 if (nregs < 0) stop(paste("No regions defined in file",region.file))
 regcodes <- as.numeric(substr(file.input,1,5))
 regnames <- substr(file.input,6,nchar(file.input))
 reg.frame <- data.frame(Region=regcodes,Name=regnames)
 reg.frame <- reg.frame[order(reg.frame$Region),]
 if (!all(sort(reg.frame$Region) == 0:nregs)) {
  stop(paste("Regions incorrectly specified:",region.file,
             "should contain definitions\n",
             " for regions 0 (for whole area) through",nregs))
 }
 reg.frame
}
##############################################################################
##############################################################################
##############################################################################
define.regions <- function(names,codes=NULL) {
#
#       To produce a reginfo object from a vector of region names
#
 nregs <- length(names)
 if (is.null(codes)) codes <- 0:(nregs-1)
 if (length(codes) != nregs) stop("codes and names have different lengths")
 reg.frame <- data.frame(Region=codes,Name=names)
 reg.frame <- reg.frame[order(reg.frame$Region),]
 if (!all(reg.frame$Region == 0:(nregs-1))) {
  stop(paste("Regions incorrectly specified: you must provide",
             "names\n for regions 0 (for whole area) through",nregs))
 }
 reg.frame
}
##############################################################################
##############################################################################
##############################################################################
read.modeldef <- 
 function(model.file,nhead=46,model.type,which.part="mean",siteinfo,
          var.names,which.response=1,
          external.files=c("yr_preds.dat","mn_preds.dat","dy_preds.dat"),
          sim=FALSE,oldGlimClim.warning=TRUE) {
#
#       To read a model definition from file. Arguments:
#
#       model.file      File from which to read
#       nhead           Number of header lines in file. The 
#                       default value of 38 is the number of 
#                       header rows in 'siteinfo.def' files for
#                       previous versions of GLIMCLIM.
#       model.type      One of "normal", "gamma", "logistic" or
#                       "normal-heteroscedastic". 
#       which.part      If model.type="normal-heteroscedastic",
#                       indicates whether we're defining the 
#                       mean (which.part="mean") or dispersion 
#                       (which.part="dispersion") components. 
#       siteinfo        An object of class siteinfo, resulting
#                       from a call to read.siteinfo() or 
#                       make.siteinfo(). Used for checking that
#                       any required site attributes have 
#                       been defined
#       var.names       Vector of variable names (used to construct 
#                       labels, in multivariate models in particular) 
#       which.response  Index number of response variable in data file 
#       external.files  Names of files containing "external" covariate
#                       information (not used unless external 
#                       covariates are explicitly specified in 
#                       the model definition)
#       sim             Indicates whether we're defining a model
#                       for simulation or fitting purposes. 
#                       For simulating some model types, a
#                       dispersion parameter should be specified 
#                       in the model definition file (this is not
#                       necessary for model fitting purposes). 
#                       Conversely, if we're fitting models then
#                       we'll need to calculate derivatives of
#                       log-likelihoods with respect to any 
#                       parameters that enter nonlinearly into
#                       the "linear predictor". Defaults to FALSE.
#
 if (!file.exists(model.file)) stop(paste(model.file,"not found"))
 if (!inherits(siteinfo,"siteinfo")) {
  stop("Wrongly structured siteinfo argument - use make.siteinfo()")
 }
 if (missing(model.type)) stop("model.type must be supplied")
#
#       Minimum number of codes required for each component type
#
 Ncodes <- c(1,1,1,1,2,3,2,1,0,1)
#
#       Initialise flags and dimensions
#
 max.compnt <- length(Ncodes)
 Np <- rep(0,max.compnt)
 trace.needed <- spmod <- 0
 nvars <- which.response	# Know there must be at least which.response variables!
 dispersion <- NA
 iext <- rep(0,3)
 old.compnt <- -1
#
#       Read input from file
#
 file.input <- scan(model.file,what=character(), sep="\n",
                                        quiet=TRUE, skip=nhead+1)
#
#       Model title is the first line in what remains
#
 model.title <- file.input[1]; file.input <- file.input[-1]
 model.title <- sub(' +$', '', model.title)
#
#       And we know what is the maximum number of parameters, so
#       can dimension things accordingly. NB when passing to 
#       Fortran, the coefficient vector is indexed from 0:MXP,
#       so max.pars here corresponds to MXP. Also, Fortran crashes
#       if MXP is zero or if it is less than (NATTR + number of 
#       nonlinear transformations of site attributes). Since the 
#       number of such nonlinear transformations can't be greater
#       than the number of lines in the definition file, that 
#       gives us an upper bound.
#
 max.pars <- max(length(file.input)-1+siteinfo$Nattr,1)
 beta <- rep(NA,max.pars+1)
 global.codes <- rep(-1,max.pars)
 fourier.idx <- legendre.idx <- rep(0,max.pars)
 global.vals <- rho <- rep(1e9,max.pars)
 threeway <- matrix(0,nrow=max.pars,ncol=3)
 nlstat <- matrix(-1,nrow=max.pars,ncol=3)
 theta <- matrix(1e9,nrow=max.pars,ncol=3)
 sitxfm <- twoway <- matrix(0,nrow=max.pars,ncol=2)
 covcode <- weighting.scheme <- rep(NA,max.pars)
 prevwt.idx <- array(0,dim=c(max.pars,4,1))
#
#       Now loop over the remaining lines of input 
#
 warned <- FALSE
 for (N in 0:min(max.pars,length(file.input)-1)) {
  i <- N+1
  code <- as.numeric(substring(file.input[i],c(1,6,16,21,26),
                                                 c(5,15,20,25,30)))
  code[is.na(code)] <- 0
  compnt <- code[1]; coeff <- code[2]; code <- code[3:5]
#
#       Update according to the value of COMPNT - if it's zero then we've
#       got the constant and there's nothing more to do, otherwise need
#       to check input, evaluate code and update Np. NB indexing of the
#       beta argument in FORTRAN starts at 0, whereas it starts at
#       1 here (zero array indices not allowed in R)
#
  if (model.type == "normal-heteroscedastic") {
   if (compnt == 9) stop(err.msg(model.file,N+3+nhead,48))
   if (compnt %in% c(8,10) & which.part == "dispersion") {
    stop(err.msg(model.file,N+3+nhead,49))
   }
  } 
  if (compnt <= 6) beta[N+1] <- coeff

  if (compnt == 0) {
   if (N != 0) stop(err.msg(model.file,N+3+nhead,2))
   next
  }
#       
#       Don't allow users to get their input file in the wrong order
#       (it will almost inevitably screw up the interactions). Also
#	check that the value of COMPNT is valid.
#
  if (compnt < old.compnt) stop(err.msg(model.file,N+3+nhead,3))
  if (compnt > max.compnt) stop(err.msg(model.file,N+3+nhead,20))
#
#	Now check that at least the minimum required amount of 
#	information has been given. The only situations in which a 
#	CODE of zero is allowed in a field that would normally be
#	required is if we have the default spatial structure, if
#	we're defining the dispersion parameter or if we're
#       doing contemporaneous regression on other variables. The 
#       "if"s here are a bit inelegant, but they work!
#
  if ( !(compnt %in% c(4,9)) | ((code[1]!=0) & (compnt != 9)) )  {
   if (code[Ncodes[compnt]] == 0) {
    stop(err.msg(model.file,N+3+nhead,32))
   }
  }
  old.compnt <- compnt
#
#	At this stage we update everything except NP(7) - NP(10)
#	(NP(7)=0 means something later on, and NP(9) (dispersion
#	parameter), NP(8) (global quantities) and NP(10) (spatial 
#       structure) are unrelated to the rest of the NPs.
#
  if (compnt <= 6) Np[compnt:6] <- Np[compnt:6] + 1
  
  if (compnt <= 4) {
#
#	Possible error trap!
#
   if (N <= 0) {
    stop(err.msg(model.file,0,4))
   }
   covcode[N] <- code[1]
#
#	For site effects, check that the requested attribute has been 
#       defined. Then record if any nonlinear transformations
#	have been defined, and remember what site attributes they
#	correspond to. Also, allow estimation of nonlinear parameters
#	here.
#
   if (compnt == 1) {
    if (code[1] > siteinfo$Nattr) {
     stop(err.msg(model.file,N+3+nhead,17))    
    }
    if (code[2] > 0) {
     sitxfm[N,] <- code[1:2]
     nlstat[N,] <- c(0,0,NA)
    }
   }
#
#	Also allow nonlinear parameters to be estimated for trend effects,
#       and note whether there are any external predictors on an annual
#       timescale; if so, check that the required input file is present.
#       The sign() thing is to allow for the fact that the user may 
#       ask for negatively lagged external predictors (i.e. future
#       values) - it was ISIGN(COVCODE(N),CODE(2)) in Fortran, which 
#       transfers the sign of the second argument to the first (with
#       convention that a + sign is transferred if the second arg is
#       zero - R transfers a zero in this case, hence the addition of 
#       1e-6)
#
   if (compnt == 2) {
    if (code[1] %in% 4:49) stop(err.msg(model.file,N+3+nhead,19)) 
    nlstat[N,] <- c(0,0,NA)
    if (code[1] > 50) {
     iext[1] <- 1
     if (!file.OK(external.files[1],"old")) return()
     covcode[N] <- (1000*code[2]) + (abs(covcode[N])*sign(code[2]+1e-6))
    }
   }
#
#     Monthly effects - note whether any external ones have been 
#     requested.
#
   if (compnt == 3) {
    if (code[1] %in% c(9:10,23:49)) stop(err.msg(model.file,N+3+nhead,19)) 
    if (code[1] > 50) {
     iext[2] <- 1
     if (!file.OK(external.files[2],"old")) return()
     covcode[N] <- (1000*code[2]) + (abs(covcode[N])*sign(code[2]+1e-6))
    }
   }
   if (compnt == 4) {
#
#       For 'previous days' effects, keep track of the maximum lag required
#       for each variable. We also need to keep track of (a) which 
#       variable is being considered (b) whether any nonlinear transformations
#       are being defined. This is done by coding in COVCODE e.g. a value
#       of 2005003 means variable 2, transformation 5, 3 days ago. This
#       definition is different from that in the old (non-R) version of
#       Glimclim where code[3] was taken to indicate the number of previous
#       days' values that were required to be present: therefore, if a
#       value of code[3] is present a warning message will be issued to 
#       alert the user (this can be suppressed using the oldGlimClim.warning
#       argument to this function). Finally, for compatibility with old
#       versions, if code[3] is zero or blank then previous values of 
#       the *response* variable will be used. 
#
    if (code[1] %in% c(11:20,29:30,43:49)) stop(err.msg(model.file,N+3+nhead,19)) 
    if (code[1] <= 10) {
     if (code[1] == 0 & code[3] %in% c(0,which.response)) {
      stop(err.msg(model.file,N+3+nhead,23))
     }
     if (code[3] > 0) {
      if (oldGlimClim.warning & !warned) {
       warning(paste("When defining contributions of previous days' values,", 
                     "meaning of code3\n  is different in Rglimclim",
                     "to that in previous versions of GlimClim\n ",
                     "(see manual for details). Set",
                     "oldGlimClim.warning=FALSE to suppress\n  this",
                     "warning."
              ),immediate.=TRUE)
       warned <- TRUE
      }
      nvars <- max(nvars,code[3])
#
#       Ensure we have enough space to record information on weighting
#       schemes for all variables if necessary
#
      prevwt.idx <- array(0,dim=c(max.pars,4,nvars))
      covcode[N] <- covcode[N] + (1000*code[2]) + (1e6*code[3])
     } else {
      covcode[N] <- covcode[N] + (1000*code[2]) + (1e6*which.response)
     }
     if (covcode[N] >= 2^31) stop(err.msg(model.file,N+3+nhead,2))  
#
#     Mark any request for a trace value - so we can check that the
#     threshold gets defined.
#
     if ( (code[2] %% 10) == 4) trace.needed <- 1
#
#     Are we doing spatial averaging of previous days' values, that
#     involves nonlinear functions of inter-site distance? If so,
#     make sure we've got at least 2 site attributes defined (otherwise
#     we can't calculate distances).
#
     if ( (floor(code[2]/10) %% 10) > 1) {
      if (siteinfo$Nattr < 2) stop(err.msg(model.file,N+3+nhead,18))
      weighting.scheme[N] <- floor(code[2]/10) %% 10
      nlstat[N,] <- rep(0,3)
     }
    } else {
#
#     Check no transformations requested for anything other than 
#     previous days' effects (may be 'daily' effects other than 
#     previous days)
#
     if (code[2] > 0) stop(err.msg(model.file,N+3+nhead,31))
     if (code[1] > 50) {
      iext[3] <- 1
      if (!file.OK(external.files[3],"old")) return()
      covcode[N] <- (1000*code[2]) + (abs(covcode[N])*sign(code[2]+1e-6))
     }
    }
   }
  } else if (compnt == 5) {
#
#       2-way interactions - we check that the terms requested 
#       are in fact main effects
#
   if ( ( !((code[1] >= 1) & (code[1] <= Np[4])) ) | 
        ( !((code[2] >= 1) & (code[2] <= Np[4])) ) ) {
    stop(err.msg(model.file,N+3+nhead,21))
   }
   twoway[N,] <- code[1:2]
  } else if (compnt == 6) {
#
#       And 3-way interactions
#
   if ( ( !((code[1] >= 1) & (code[1] <= Np[4])) ) |
        ( !((code[2] >= 1) & (code[2] <= Np[4])) ) |
        ( !((code[3] >= 1) & (code[3] <= Np[4])) ) ) {
    stop(err.msg(model.file,N+3+nhead,21))
   }
   threeway[N,] <- code[1:3]
  } else if (compnt == 7) {
#
#	Now parameters in nonlinear functions. code[1] will be the
#	predictor where the nonlinear function was requested
#	(NB we know what the function is, because it was defined
#	by the request); code[2] - either 1 or 2 - the number
#	of the parameter we're currently defining, and code[3], if
#	present and non-zero, indicates that we're supposed to estimate
#	the parameter rather than take it as known (the default). 
#       The check that theta[code[1],code[2]] is not NA is to
#       ensure that we haven't already defined it. If the parameter is 
#       to be estimated rather than assumed known, we duplicate it 
#       in beta, as then the IWLS procedure is easier to implement.
#       Useful to have in theta as well, for when we come to 
#       recalculate nonlinear functions.
#
   if ( ( !((code[1] >= 1) & (code[1] <= Np[4])) ) |
        ( !((code[2] >= 1) & (code[2] <= 3)) ) ) {
    stop(err.msg(model.file,N+3+nhead,33))
   } else if (abs(theta[code[1],code[2]]) < 1e8) {
    stop(err.msg(model.file,N+3+nhead,30))
   } else if (nlstat[code[1],code[2]] != 0) {
    stop(err.msg(model.file,N+3+nhead,31))
   }
#
#	We've found a nonlinear parameter, so Np[7] is no longer 
#	zero. For nonlinear parameters, covcode is used to 
#	identify (a) the covariate number with which it is
#	associated (b) the parameter number for that covariate.
#	The formula is covcode = 1000a + b. E.g. covcode = 3002
#	means covariate 3, parameter 2. We copy it both to the 
#	parameter vector (for updating estimates) and to the
#	theta array (easier bookkeeping when calculating covariates)
#	Note that if we're estimating a nonlinear parameter, the 
#	associated beta can't be zero. Also, parameters for 
#	orthogonal series representations of site attributes are 
#	always treated as known, so check that we haven't been 
#	asked to estimate these. There are additional complications
#	for orthogonal series and weighted averages of previous days'
#       values, like I'm only going to allow limits to be specified 
#       once.
#
   if (Np[7] == 0) Np[7] <- Np[6]
#
#	Orthogonal series options are fiddly and relegated to 
#	a separate function.
#
   if ( (sitxfm[code[1],2] >= 11) &
        (sitxfm[code[1],2] <= 30) ) {
    tmp <- orth.set(sitxfm,code,fourier.idx,theta,
                    siteinfo,model.file,nhead,N,coeff)
    theta <- tmp$theta; fourier.idx <- tmp$ortidx
   } else if ( (sitxfm[code[1],2] >= 31) &
               (sitxfm[code[1],2] <= 40) ) {
    tmp <- orth.set(sitxfm,code,legendre.idx,theta,
                    siteinfo,model.file,nhead,N,coeff)
    theta <- tmp$theta; legendre.idx <- tmp$ortidx
#
#     Ditto weighted averages of previous days' values
#
   } else if (isTRUE(weighting.scheme[code[1]] > 1)) { 
    varnum <- floor(covcode[code[1]] / 1e6)
    if (varnum > nvars) stop("Programming error: varnum > nvars - tell Richard")
    tmp <- lagweights.def(weighting.scheme,code,varnum,prevwt.idx,
                          theta,model.file,nhead,N,coeff) 
    theta <- tmp$theta; prevwt.idx <- tmp$prevwt.idx
#
#	No such complications with any other nonlinear functions
#
   } else {
    theta[code[1],code[2]] <- coeff
   }
#
#	Now see if the parameter is to be estimated & deal accordingly
#       (weighting schemes: copy covariate code into index for 
#       derivatives as well). Remember that indexing for beta starts at
#       0 in the Fortran code - hence need to add 1 to all the indices
#       here
#
   if ( (code[3] != 0) & !sim ) {
    nlstat[code[1],code[2]] <- code[3]
    if (abs(beta[code[1]+1]) <= 0) stop(err.msg(model.file,N+3+nhead,35))
    Np[7] <- Np[7] + 1
    covcode[Np[7]] <- (1000*code[1]) + code[2]
    beta[Np[7]+1] <- coeff  
    if (isTRUE(weighting.scheme[code[1]] > 1)) {
#
#       NB in the Fortran code, the columns of prevwt.idx
#       are numbered from 0 to 3 rather than from 1 to 4
#       as here - hence code[2]+1 in next line. 
#
     prevwt.idx[weighting.scheme[code[1]],code[2]+1,varnum] <- Np[7]
    }
   }
  } else if (compnt == 8) {
#
#     `Global' quantities. Each can be defined at most once. We've
#     already checked that code[1] has been defined. We also record
#     whether or not a trace threshold has been defined, by setting
#     trace.needed  to 0 if it has (translates as `no longer require 
#     the threshold to be defined'). Note that trace thresholds 
#     cannot be defined for Gaussian models
#
   Np[8] <- Np[8] + 1
   global.codes[Np[8]] <- (1000*code[1]) + code[2]
   if (abs(global.vals[code[1]]) <= 1e8) {
    stop(err.msg(model.file,N+3+nhead,30))
   }
   if (code[1] == 1) {
    if (model.type %in% c("normal","normal-heteroscedastic")) {
     stop(err.msg(model.file,N+3+nhead,26))
    }
    if (code[2] == 0) stop(err.msg(model.file,N+3+nhead,32))
    if ((code[2] < 0) | (code[2] > 3)) {
                        stop(err.msg(model.file,N+3+nhead,34))
    }
    if (code[2] == 1) trace.needed <- 0
   } else stop(err.msg(model.file,N+3+nhead,34))
   global.vals[code[1]] <- coeff
  } else if (compnt == 9) {
#
#	Dispersion parameter. No more than 1! (and none for a logistic
#       model)
#
   if (Np[9] == 1) stop(err.msg(model.file,N+3+nhead,50))
   if (model.type == "logistic") stop(err.msg(model.file,N+3+nhead,52))
   Np[9] = Np[9] + 1
   dispersion <- coeff
  } else if (compnt == 10) {
#
#       Here's specification of spatial structure. We check (a) that all
#	entries relating to spatial structure attempt to define the same
#	model (b) that no parameter is defined more than once (c) that the
#	spatial structure being defined is appropriate for the type of
#	model being used (variable model.type). We use a similar coding 
#       for COVCODE e.g. COVCODE = 3002 means covariance model 3, parameter 
#       2. Distinguish between `default' and `requested' selection of model 
#       0 (`requested') gets set to -1 for now).
#
   if ( !(code[1] %in% c(0:8,21,22)) ) stop(err.msg(model.file,code[1],46))
   if (spmod != 0) {
    if ((spmod != code[1]) & !((spmod == -1) & (code[1] == 0)) ) {
     stop(err.msg(model.file,N+3+nhead,45))
    } else if ( (code[2] == 0) & (spmod > 1) ) {
     stop(err.msg(model.file,N+3+nhead,32))
    }
   }
   if (code[2] != 0) {
    if (rho[code[2]] < 1e8) stop(err.msg(model.file,N+3+nhead,30))
   }
#
#	Spatial structures 1 through 20 are generic, but 21 through 
#       40 are only valid for logistic model (model.type=1).
#
   if ( (model.type != "logistic") & (floor((code[1]-1)/20) == 1) ) {
    stop(err.msg(model.file,N+3+nhead,47))
   }
   if (code[1] < 1) {
    spmod <- -1
   } else if (code[1] != 1) {
    if (code[2] == 0) stop(err.msg(model.file,N+3+nhead,32))
    rho[code[2]] <- coeff
    Np[10] <- Np[10]+ 1
    spmod <- code[1]
   } else spmod <- code[1]
  } else {
#
#	Trap invalid values of COMPNT
#
   stop(err.msg(model.file,N+3+nhead,20))
  }
 }
#
#	Reset 'requested' default spatial structure to model 0, and 
#	check that a dispersion parameter has been defined if we're in
#	'simulation' mode.
#
 if (spmod == -1) spmod <- 0
 if ( (Np[9] == 0) & sim ) stop(err.msg(model.file,0,51))
#
#     Check that trace threshold has been defined if it's required
#
 if (trace.needed == 1) stop(err.msg(model.file,0,25))
#
#       Calculate the maximum lag referenced in any of the "previous
#       day" codes.
#
  max.lag <- rep(0,nvars)
  prevcodes <- covcode[(Np[3] + 1):Np[4]]
  prevcodes <- prevcodes[floor(prevcodes/1e+06) > 0]
  max.lag <- set.maxlag(max.lag,prevcodes,which.response,TRUE)
  max.lag <- as.integer(max.lag)
#
#     Check that all required parameters for spatial structures
#     have been defined, and that correlations are < 1
#
 rho <- rho.check(rho,spmod,model.file)
#
#       Coerce all pointers and indices to integer
#
 storage.mode(covcode) <- "integer"
 storage.mode(twoway) <- "integer"
 storage.mode(threeway) <- "integer"
 storage.mode(nlstat) <- "integer"
 storage.mode(sitxfm) <- "integer"
 storage.mode(fourier.idx) <- "integer"
 storage.mode(legendre.idx) <- "integer"
 storage.mode(prevwt.idx) <- "integer"
 storage.mode(weighting.scheme) <- "integer"
 storage.mode(global.codes) <- "integer"

#
#       Extract the necessary descriptive text from the Fortran code
#       (this requires creativity!)
#
 nsites <- as.integer(siteinfo$Nsites)
 nattr <- as.integer(siteinfo$Nattr)
 sitinf <- array(dim = c(nsites, max.pars, 4))
 sitinf[,1:nattr, 1] <- siteinfo$Attribute.values

 UnitNos <- rep(0,10)
 for (i in 1:3) {
  if (iext[i] != 0) { 
   if (!file.OK(external.files[i],"old")) return()
   UnitNos[2+i] <- OpenFile(external.files[i])
  }
 }
 storage.mode(UnitNos) <-  "integer"

#
#      For versions of R from 3.6.1 onwards, the only way to get character
#      info from R into Fortran is to write to a temporary file and then
#      read it back in from the Fortran code - VERY annoying. 
#
 TempFile <- OpenTempFile()
 write(siteinfo$Site.codes, file=TempFile$Name, append=TRUE)
 text.strings <- tryCatch(capture.output(invisible(.Fortran("WriteLabels",
          model=model.num(model.type,which.part),nsites=nsites,
          mxp=as.integer(max.pars),nattr=nattr,sitinf=sitinf,
          TmpFlNo=TempFile$Number,np=as.integer(Np),
          covcode=covcode,nvars=as.integer(dim(prevwt.idx)[3]),
          theta=theta,sitxfm=sitxfm,fouidx=fourier.idx,
          legidx=legendre.idx,pwtidx=prevwt.idx,wtschm=weighting.scheme,
          spmod=as.integer(spmod),glbcod=global.codes,UnitNos=UnitNos,
          PACKAGE="Rglimclim",NAOK=TRUE)
 )),finally=CloseFiles(UnitNos))
 unlink(TempFile$Name)
 
 ifail <- grep("Error code ",text.strings)
 if (length(ifail) > 0) { 
  ifail <- as.numeric(substr(text.strings[ifail[1]],12,
                             nchar(text.strings[ifail[1]])))
  stop(err.msg(model.file,0,ifail))
 }
#
#       And create arrays of labels for use when printing output
#
 if (missing(var.names)) {
  var.names <- paste("Y",1:nvars,sep="")
  if (length(var.names) == 1) var.names <- "Y"
 }
 model.labels <- make.labels(text.strings,siteinfo$Attribute.names,
                             var.names,covcode,Np)
 
 z <- list(model.type=model.type,Np=as.integer(Np),iext=as.integer(iext),
      model.title=model.title,max.pars=as.integer(max.pars),
      max.lag=max.lag,beta=beta,covcode=covcode,twoway=twoway,
      threeway=threeway,theta=theta,nlstat=nlstat,sitxfm=sitxfm,
      fourier.idx=fourier.idx,legendre.idx=legendre.idx,
      prevwt.idx=prevwt.idx,weighting.scheme=weighting.scheme,
      dispersion=dispersion,spmod=as.integer(spmod),rho=rho,
      global.codes=global.codes,global.vals=global.vals,
      model.labels=model.labels,siteinfo=siteinfo,var.names=var.names,
      which.response=which.response)
 class(z) <- "GLC.modeldef"
 z
}
##############################################################################
##############################################################################
##############################################################################
rho.check <- function(rho,spatial.model,model.file) {
#
#     To check that the correct number of parameters has been defined 
#     for a spatial dependence structure, and that they satisfy any 
#     necessary constraints. Arguments:
#
#     rho               Array of parameters
#     spatial.model     Selected dependence model
#     model.file        Name of file from which we're reading 
#                       model definition
#

#
#       Start by setting up arrays giving the numbers of parameters
#       required by different model specifications. NB the specifications
#       themselves are defined in the FORTRAN code, in routine CORSET
#       for correlation-based structures (currently spatial.model %in% 1:8);
#       and in SPEST / WD21 / WD22 (for models 21 and 22). 
#
 npar <- rep(-1,40)
 npar[c(1:8,21,22)] <- c(0,1,1,2,2,3,2,3,1,3)
 ll <- rep(-1e6,40); ul <- -ll
 if (spatial.model %in% c(0,1,21,22)) {
# do nothing! 
 } else if (spatial.model == 2) {
  ll[1] <- -1 ; ul[1] <- 1
 } else if (spatial.model %in% 3:4) {
  ll[1] <- 0
  if (spatial.model ==4) {
   ll[2] <- -1 ; ul[2] <- 1   
  } 
 } else if (spatial.model %in% 5:6) {
  ll[1:2] <- 0
  if (spatial.model == 6) {
   ll[3] <- -1 ; ul[3] <- 1      
  } 
 } else if (spatial.model == 7) {
  ll[1:2] <- 0
  ul[2] <- 1
 } else if (spatial.model == 8) {
  ll[1:3] <- 0
  ul[3] <- 1
 } else {
  if (missing(model.file)) {
   stop("Unknown spatial dependence structure in model definition")
  } else {
   stop(paste("Unknown spatial dependence structure requested in file",
              model.file))
  }
 }
#
#       Check that the requested spatial model exists in the software
#
 if (spatial.model > 0) {
  if (npar[spatial.model] < 0) {
   if (missing(model.file)) {
    stop("Unknown spatial dependence structure in model definition")
   } else {
    stop(paste("Unknown spatial dependence structure requested in file",
               model.file))
   } 
  }
 }
#
#       Now check that all necessary parameters have been defined
#
 needed.pars <- numeric(0)
 undefined.pars <- 1:length(rho)
 if (spatial.model > 0) {
  if (npar[spatial.model] > 0) {
   needed.pars <- 1:npar[spatial.model]
   undefined.pars <- (npar[spatial.model]+1):length(rho)
  }
  if (spatial.model==22) { # Parameters 2 & 3 are optional for this one: set to defaults if undefined
   needed.pars <- 1 
   if (rho[2] > 1e8) rho[2] <- 0.01
   if (rho[3] > 1e8) rho[3] <- 1
  }
 }
 if (any(rho[needed.pars] > 1e8) | any(rho[undefined.pars] < 1e8)) {
  if (missing(model.file)) {
   stop(paste(" Wrong number of parameters for spatial dependence",
              "structure\nin model definition"))
  } else {
   stop(paste(" Wrong number of parameters supplied in file",model.file,
              "\nfor spatial dependence structure"))
  }
 }
#
#       And that they all satisfy constraints 
#
 if (any(rho[needed.pars] < ll[needed.pars] | 
         rho[needed.pars] > ul[needed.pars])) {
  if (missing(model.file)) {
   stop(paste("One or more spatial dependence structure parameters is",
              "out of bounds."))
  } else {
   stop(paste("In file",model.file,",one or more spatial dependence",
              " structure\n parameters is out of bounds.",sep=""))
  }
 }
 rho
}
##############################################################################
##############################################################################
##############################################################################
orth.set <- function(sitxfm,code,ortidx1,theta1,
                     siteinfo,filename,nhead,N,coeff) {
##############################################################################
#
#	To set limits for orthogonal series representations of site
#	effects, and check that these have been defined correctly.
#	This stuff really belongs in the MDLSET routine, but it's so
#	fiddly that it's easier to read if relegated to its own
#	personal space. Arguments:
#
#	sitxfm	- indicates which nonlinear functions of site properties
#		  have been selected.
#       code    - Identifiers of attribute being defined.
#	ortidx1	- Index array to track attributes for which orthogonal
#		  series representations are being defined, on entry.
#	theta	- parameters in nonlinear functions on entry.
#	siteinfo- Array of site information, used to check validity of
#		  limits. 
#       filename- Input file name. Used for identifying where errors are
#		  found.
#       nhead   - No. of header lines in file.
#	N	- Used to identify whereabouts in the input file we are,
#		  for error tracking.
#       coeff   - value we're currently trying to allocate. 
#
#       Value: a list containing components ortidx and theta, which
#       are updated versions of ortidx1 and theta1
#
##############################################################################
#
#	First check: all parameters here should be regarded as known
#
 if (code[3] != 0) stop(err.msg(filename,N+3+nhead,36))
 ortidx <- ortidx1; theta <- theta1
#
#	Next: has this limit been defined somewhere else?
#
 if (ortidx[sitxfm[code[1],1]] != 0) {
  if ( (theta[ortidx[sitxfm[code[1],1]],1] < 1e8) & 
       (theta[ortidx[sitxfm[code[1],1]],2] < 1e8) ) {
   stop(err.msg(filename,N+3+nhead,30))
  }
 }
#
#	If we've got here, things are looking basically OK so set the 
#	parameter: if it's the first one being defined for this attribute,
#	record it.
#
 if ( ortidx[sitxfm[code[1],1]] == 0) {
  ortidx[sitxfm[code[1],1]] <- code[1]
  oset <- 1
 } else oset <- 2
 theta[ortidx[sitxfm[code[1],1]],code[2]] <- coeff
#
#	If both parameters have been defined, some further checks ...
#
 if (oset == 2 ) {
#
#	Are they in the right order? (don't allow UL < LL as it may 
#	indicate some other inadvertent cock-up). 
#
  LL <- theta[ortidx[sitxfm[code[1],1]],1]
  UL <- theta[ortidx[sitxfm[code[1],1]],2]
  if ( LL >= UL ) stop(err.msg(filename,N+3+nhead,38))
#
#	Are all sites within the range? If not, representation is 
#	invalid. 
#
   outside <- (siteinfo$Attribute.values[,sitxfm[code[1],1]] <= LL) |
              (siteinfo$Attribute.values[,sitxfm[code[1],1]] >= UL)
   if (any(outside)) {
    stop(err.msg(filename,((1:siteinfo$Nsites)[outside])[1],39))
   }
 }
 list(theta=theta,ortidx=ortidx)
}
##############################################################################
##############################################################################
##############################################################################
lagweights.def <- function(weight.scheme,code,varnum,prevwt.idx1,
                           theta1,filename,nhead,N,coeff) {
##############################################################################
#
# To define parametrisation of weighted averages of previous days'
# values, and check that these have been defined correctly. Same
# idea as orth.set routine above. Arguments:
#
# weight.scheme - indicates which weighting scheme has been selected.
# code          - Identifiers of attribute being defined, read from input
#		  file.
# varnum        - Number of variable for which weights are being defined.
# prevwt.idx1   - Index array to track predictors for which weighting
#		  parameters are being defined, on entry. NB the 
#                 columns of this are indexed from 1 to 4 in R, 
#                 but from 0 to 3 in the underlying Fortran code.
# theta1        - parameters in nonlinear weighting functions, on entry
# filename      - Input file name. Used for identifying where errors are
#		  found.
# nhead         - No. of header lines in file. 
# N             - Used to identify whereabouts in the input file we are,
#		  for error tracking.
# coeff         - value we're currently trying to allocate. 
##############################################################################
#
#	Has this parameter been defined somewhere else? (NB
#       refs to column 1 of prevwt.idx here correspond to column
#       0 in the Fortran code). code[1] is the number of the 
#       covariate to which the weighting scheme is attached.
#       We might need to "grow" the prevwt.idx array if the 
#       current variable hasn't been seen before.
#
 cursize <- dim(prevwt.idx1)
 prevwt.idx <- array(dim=pmax(cursize,c(0,0,varnum)))
 prevwt.idx[1:cursize[1],1:cursize[2],1:cursize[3]] <- prevwt.idx1
 theta <- theta1
 if (prevwt.idx[weight.scheme[code[1]],1,varnum] != 0) {
  if ( (theta[prevwt.idx[weight.scheme[code[1]],1,varnum],code[2]] < 1e8) ) {
   stop(err.msg(filename,N+3+nhead,30))
  }
 }
#
#	If we've got here, things are looking basically OK so set the 
#	parameter: if it's the first one being defined for this variable,
#	record the covariate to which it's attached. Otherwise make
#       a note of the model parameter number.
#
 if ( prevwt.idx[weight.scheme[code[1]],1,varnum] == 0) {
  prevwt.idx[weight.scheme[code[1]],1,varnum] <- code[1]
 }
 theta[prevwt.idx[weight.scheme[code[1]],1,varnum],code[2]] <- coeff

 list(prevwt.idx=prevwt.idx,theta=theta)
}
##############################################################################
##############################################################################
##############################################################################
make.labels <- function(text.strings,attr.names,var.names,covcode,Np) {
##############################################################################
#
#   To generate arrays of labels for model output. Arguments:
#
#   text.strings A character vector arising from a call to 
#                capture.output(invisible(.Fortran("WriteLabels",...))).
#                Don't ask.
#   attr.names   A vector of site attribute names
#   var.names    A vector of names for the daily weather variables
#   which.response      Index number of response variable
#   covcode      Coding of covariates in the model
#   Np           A vector indicating the size of each component of
#                an Rglimclim model.
#
#   Value: a list of five character vectors, containing labels for
#   different components of the model
#
##############################################################################
 text.type <- substr(text.strings,1,10)
#
#       Covariates (replace site attribute text as appropriate)
#
 wanted <- (text.type == "Covariate ")
 tmparr <- strsplit(text.strings[wanted],"|| ",fixed=TRUE)
 tmparr <- matrix(unlist(tmparr),ncol=3,byrow=TRUE)[,-1,drop=FALSE]
 cov.nos <- as.integer(tmparr[,1])
 if (!identical(cov.nos,1:(Np[4]+1))) {
  stop("Programming error 1 - tell Richard (richard@stats.ucl.ac.uk)")
 }
 xlabels <- tmparr[,-1]
 wanted <- grep("~~~",xlabels)
 if (length(wanted) > 0) {
  site.text <- xlabels[wanted]
  cutout <- gregexpr("~~~",site.text)
  if (!all(sapply(cutout,FUN=length) == 2)) {
   stop("Programming error 2 - tell Richard (richard@stats.ucl.ac.uk)")
  } 
  cutout <- matrix(unlist(cutout),ncol=2,byrow=TRUE)
  cutout[,2] <- cutout[,2] + 3   # Cut out the final ~~ (yes, I know 
                                 # it should be cutout[,2] + 2, but
                                 # wait to see what happens later)
  attr.initxt <- substr(site.text,cutout[,1],cutout[,2])
  attr.sec <- gregexpr("#",attr.initxt)
  if (!all(sapply(attr.sec,FUN=length) == 2)) {
   stop("Programming error 3 - tell Richard (richard@stats.ucl.ac.uk)")
  } 
  attr.sec <- matrix(unlist(attr.sec),ncol=2,byrow=TRUE)
  attr.nos <- as.numeric(substr(attr.initxt,attr.sec[,1]+1,attr.sec[,2]-1))
  site.text <- paste(substr(site.text,1,cutout[,1]-1),attr.names[attr.nos],
                     substr(site.text,cutout[,2],nchar(site.text)),sep="")
  xlabels[wanted] <- site.text
 }
#
#       Lagged covariate values (NB xlabels starts with the intercept,
#       but covcode omits this, hence the -1)
#
 wanted <- grep("###",xlabels)
 if (length(wanted) > 0) {
  varnums <- floor(covcode[wanted-1] / 1e6)
  for (i in 1:length(wanted)) {
   xlabels[wanted[i]] <- gsub("###",var.names[varnums[i]],
                              xlabels[wanted[i]])
  }
 }
#
#       Nonlinear transformations
#
 wanted <- (text.type == "Nonlinear ")
 tmparr <- strsplit(text.strings[wanted],"|| ",fixed=TRUE)
 if (length(tmparr) > 0) {
  tmparr <- matrix(unlist(tmparr),ncol=4,byrow=TRUE)[,-1,drop=FALSE]
  cov.nos <- as.integer(tmparr[,1])
  par.nos <- as.integer(tmparr[,2])
  nonlin.txt <- matrix(tmparr[,3],ncol=3,byrow=TRUE)
  if (nrow(nonlin.txt) != Np[4]) {
   stop("Programming error 4 - tell Richard (richard@stats.ucl.ac.uk)")
  }
  nonlin.txt[grep("NOT USED",nonlin.txt)] <- NA
  nonlin.txt[grep("Undefined",nonlin.txt)] <- NA
 } else {
  nonlin.txt <- NULL
 }
#
#       Global quantities
#
 wanted <- (text.type == "Global    ")
 tmparr <- strsplit(text.strings[wanted],"|| ",fixed=TRUE)
 if (length(tmparr) > 0) {
  tmparr <- matrix(unlist(tmparr),ncol=3,byrow=TRUE)[,-1,drop=FALSE]
  global.txt <- tmparr[,2]
 } else {
  global.txt <- NULL
 }
#
#       Spatial dependence structure
#
 wanted <- (text.type == "Corstruct ")
 if (sum(wanted) != 1) {
  stop("Programming error 5 - tell Richard (richard@stats.ucl.ac.uk)")
 }
 tmparr <- unlist(strsplit(text.strings[wanted],"|| ",fixed=TRUE))
 if (length(tmparr) != 2) {
  stop("Programming error 6 - tell Richard (richard@stats.ucl.ac.uk)")
 }
 corstruct.txt <- tmparr[2]
#
#       Parameters in spatial dependence structure
#
 wanted <- (text.type == "Corparams ")
 tmparr <- strsplit(text.strings[wanted],"|| ",fixed=TRUE)
 if (length(tmparr) > 0) {
  tmparr <- matrix(unlist(tmparr),ncol=3,byrow=TRUE)[,-1,drop=FALSE]
  corpars.txt <- tmparr[,2]
 } else {
  corpars.txt <- NULL
 }
 list(xlabels=xlabels,nonlin.txt=nonlin.txt,global.txt=global.txt,
      corstruct.txt=corstruct.txt,corpars.txt=corpars.txt)
}
##############################################################################
##############################################################################
##############################################################################
write.modeldef <- function(x,file,check.file=TRUE,mean.only=FALSE) {
##############################################################################
#
#       To write model definition files that can be read subsequently
#       using read.modeldef(). Arguments:
#
#       x               An object of class GLC.modeldef, for example 
#                       resulting from a call to GLCfit(). 
#       file            Name(s) of file(s) to write to. If x is a joint
#                       mean-variance model then this should be a vector of
#                       length two (one for the mean and one for the 
#                       dispersion component); otherwise a single 
#                       character string
#       check.file      If TRUE, a check will be made that file does not
#                       already exist; if it does, user will be prompted
#                       to overwrite it. 
#       mean.only       If TRUE, only those parts of the definition 
#                       corresponding to the mean component of the model
#                       will be written.
#
##############################################################################
 if (x$model.type=="normal-heteroscedastic") {
  if (length(file) != 2) {
   stop("'file' must contain 2 filenames for normal-heteroscedastic models")
  }
 }
 if (check.file) {
  for (curfil in file) {if (!file.OK(curfil,"new")) return()}
 }
 OutText <- capture.output(invisible(.Fortran("WriteHeader",
                                              HeadType=as.integer(1),
                                              PACKAGE="Rglimclim")))
#
#       Obsolete line, followed by model description
#
 CurLine <- paste(paste(rep("%",20),collapse=""),
                  "LINE NOT CURRENTLY USED",
                  paste(rep("%",20),collapse=""),sep="")
 OutText <- c(OutText,CurLine,x$model.title)
#
#       Model specification. First the constant:
#
 Compnt <- 0
 CurLine <- paste(formatC(Compnt,digits=0,width=5,format="d"),
                  formatC(x$beta[1],digits=4,width=10,format="f"),
                  paste(rep(" ",15),collapse=""),
                  char.trunc(x$model.labels$xlabels[1],40),sep="")
 OutText <- c(OutText,CurLine)
#
#       Now site effects. There may be 2 codes here - if not, the second
#       one is a zero, it's not a nonlinear xformation of site
#       attributes and x$covcode contains the correct value. If
#       x$sitxfm[i,2] isn't zero, we're working with something nonlinear,
#       x$covcode *doesn't* contain the correct value and we have
#       to use x$sitxfm[i,1] in its place. A bit messy ...
#
 Compnt <- 1
 if (x$Np[1] > 0) {
  for (i in 1:x$Np[1]) {
   if (x$sitxfm[i,2] == 0) {
    CurLine <- paste(formatC(Compnt,digits=0,width=5,format="d"),
                     formatC(x$beta[i+1],digits=4,width=10,format="f"),
                     formatC(x$covcode[i],digits=0,width=5,format="d"),
                     paste(rep(" ",10),collapse=""),
                     char.trunc(x$model.labels$xlabels[i+1],40),
                     formatC(i,digits=0,width=3,format="d"),sep="")
   } else {
    CurLine <- paste(formatC(Compnt,digits=0,width=5,format="d"),
                     formatC(x$beta[i+1],digits=4,width=10,format="f"),
                     formatC(x$sitxfm[i,1],digits=0,width=5,format="d"),
                     formatC(x$sitxfm[i,2],digits=0,width=5,format="d"),
                     paste(rep(" ",5),collapse=""),
                     char.trunc(x$model.labels$xlabels[i+1],40),
                     formatC(i,digits=0,width=3,format="d"),sep="")    
   } 
   OutText <- c(OutText,CurLine)  
  } 
 }
#
#       Year and month effects have only 1 argument, unless there 
#       are external factors involved
#
 for (Compnt in 2:3) {
  if (x$Np[Compnt] > x$Np[Compnt-1]) {
   for (i in (x$Np[Compnt-1]+1):x$Np[Compnt]) {
    tmp <- abs(x$covcode[i]) %% 1000 
    if (tmp <= 50) {
     CurLine <- paste(formatC(Compnt,digits=0,width=5,format="d"),
                      formatC(x$beta[i+1],digits=4,width=10,format="f"),
                      formatC(x$covcode[i],digits=0,width=5,format="d"),
                      paste(rep(" ",10),collapse=""),
                      char.trunc(x$model.labels$xlabels[i+1],40),
                      formatC(i,digits=0,width=3,format="d"),sep="")
    } else {
     CurLine <- paste(formatC(Compnt,digits=0,width=5,format="d"),
                 formatC(x$beta[i+1],digits=4,width=10,format="f"),
                 formatC(tmp,digits=0,width=5,format="d"),
                 formatC(floor(x$covcode[i]/1000),digits=0,width=5,format="d"),
                 paste(rep(" ",5),collapse=""),
                 char.trunc(x$model.labels$xlabels[i+1],40),
                 formatC(i,digits=0,width=3,format="d"),sep="")
    }
   OutText <- c(OutText,CurLine)  
   }
  }
 }
#
#       Day effects are a bit complicated - they're split into
#       `Previous days' and `external':
#
 Compnt <- 4
 code <- rep(NA,3)
 if (x$Np[4] > x$Np[3]) {
  for (i in (x$Np[3]+1):x$Np[4]) {
   code[1] <- abs(x$covcode[i]) %% 1000
#
#     Here are the previous days - write max.lag out for the first one.
#
   if (code[1] <= 10) {
    code[2] <- floor(x$covcode[i] / 1000) %% 1000
    code[3] <- x$covcode[i] %/% 1e6
    CurLine <- paste(formatC(Compnt,digits=0,width=5,format="d"),
                formatC(x$beta[i+1],digits=4,width=10,format="f"),
                formatC(code[1],digits=0,width=5,format="d"),
                formatC(code[2],digits=0,width=5,format="d"),
                formatC(code[3],digits=0,width=5,format="d"),
                char.trunc(x$model.labels$xlabels[i+1],40),
                formatC(i,digits=0,width=3,format="d"),sep="")
   } else {
#
#     Here are external factors
#
    if (code[1] <= 50) {
     CurLine <- paste(formatC(Compnt,digits=0,width=5,format="d"),
                 formatC(x$beta[i+1],digits=4,width=10,format="f"),
                 formatC(code[1],digits=0,width=5,format="d"),
                 paste(rep(" ",10),collapse=""),
                 char.trunc(x$model.labels$xlabels[i+1],40),
                 formatC(i,digits=0,width=3,format="d"),sep="")
    } else {
     code[2] <- floor(x$covcode[i] / 1000)
     CurLine <- paste(formatC(Compnt,digits=0,width=5,format="d"),
                 formatC(x$beta[i+1],digits=4,width=10,format="f"),
                 formatC(code[1],digits=0,width=5,format="d"),
                 formatC(code[2],digits=0,width=5,format="d"),
                 paste(rep(" ",5),collapse=""),
                 char.trunc(x$model.labels$xlabels[i+1],40),
                 formatC(i,digits=0,width=3,format="d"),sep="")
    }
   }
   OutText <- c(OutText,CurLine)  
  }
 } 
#
#       For interactions, output appropriate codes; abandon model.labels
#
 Compnt <- 5
 if (x$Np[5] > x$Np[4]) {
  for (i in (x$Np[4]+1):x$Np[5]) {
   CurLine <- paste(formatC(Compnt,digits=0,width=5,format="d"),
                   formatC(x$beta[i+1],digits=4,width=10,format="f"),
                   formatC(x$twoway[i,1],digits=0,width=5,format="d"),
                   formatC(x$twoway[i,2],digits=0,width=5,format="d"),
                   paste(rep(" ",5),collapse=""),
                   "2-way interaction: covariates ",
                   formatC(x$twoway[i,1],digits=0,width=2,format="d"),
                   " and ",
                   formatC(x$twoway[i,2],digits=0,width=2,format="d"),sep="")
   OutText <- c(OutText,CurLine)    
  }
 }
 Compnt <- 6
 if (x$Np[6] > x$Np[5]) {
  for (i in (x$Np[5]+1):x$Np[6]) {
   CurLine <- paste(formatC(Compnt,digits=0,width=5,format="d"),
                   formatC(x$beta[i+1],digits=4,width=10,format="f"),
                   formatC(x$threeway[i,1],digits=0,width=5,format="d"),
                   formatC(x$threeway[i,2],digits=0,width=5,format="d"),
                   formatC(x$threeway[i,3],digits=0,width=5,format="d"),
                   "3-way interaction: covariates ",
                   formatC(x$threeway[i,1],digits=0,width=2,format="d"),", ",
                   formatC(x$threeway[i,2],digits=0,width=2,format="d"),
                   " and ",
                   formatC(x$threeway[i,3],digits=0,width=2,format="d"),sep="")
   OutText <- c(OutText,CurLine)    
  }
 }
#
#       Parameters in nonlinear transformations. Rather than look at 
#       definitions in positions Np[6] to Np[7], it's easier to go 
#       through the thetas themselves, directly (if you don't believe me,
#       try it and see what happens!). One slight complication:
#       for weighted averages of previous days' values, some
#       elements of theta have been used as as temporary storage -
#       so we have to check that the current row genuinely requires
#       something to be output.
#
 Compnt <- 7
 if (x$Np[4] > 0) {
  for (i in 1:x$Np[4]) {
   wt.scheme <- floor(x$covcode[i] / 10000) %% 10
   varnum <- floor(x$covcode[i] / 1e6)
#
#       Next two lines: clumsy but the second isn't defined if the 
#       first is FALSE
#
   nl.notneeded <- (i > x$Np[3]) & (wt.scheme > 1)
   if (nl.notneeded) nl.notneeded <- (i != x$prevwt.idx[wt.scheme,1,varnum])
   if (!nl.notneeded) {
    for (j in 1:3) {
     if (x$theta[i,j] < 1e8) {
      CurLine <- paste(formatC(Compnt,digits=0,width=5,format="d"),
                       formatC(x$theta[i,j],digits=4,width=10,format="f"),
                       formatC(i,digits=0,width=5,format="d"),
                       formatC(j,digits=0,width=5,format="d"),
                       formatC(x$nlstat[i,j],digits=0,width=5,format="d"),
                       "Parameter ",j," in transformation of covariate ",
                       i,sep="")
      OutText <- c(OutText,CurLine)
     }
    }
   }
  }
 }
 if (mean.only) {
  write(OutText,file=file[1])
  return()
 }
#
#     Global quantities, if present
#
 Compnt <- 8
 if (x$Np[Compnt] > 0) {
  for (i in 1:x$Np[Compnt]) {
   tmp <- floor(x$global.codes[i]/1000)
   tmp2 <- x$global.codes[i] %% 1000
   if (tmp == 1) {
    CurLine <- paste(formatC(Compnt,digits=0,width=5,format="d"),
                     formatC(x$global.vals[i],digits=4,width=10,format="f"),
                     formatC(tmp,digits=0,width=5,format="d"),
                     formatC(tmp2,digits=0,width=5,format="d"),
                     paste(rep(" ",5),collapse=""),
                     char.trunc(x$model.labels$global.txt[i],40),sep="")
    OutText <- c(OutText,CurLine)
   }
  }
 }
#
#	Dispersion parameter, if required
#
 Compnt <- 9
 if (is.numeric(x$dispersion)) {
  if (x$dispersion >= 0) {
   CurLine <- paste(formatC(Compnt,digits=0,width=5,format="d"),
                    formatC(x$dispersion,digits=4,width=10,format="f"),
                    paste(rep(" ",15),collapse=""),"Dispersion parameter",sep="")
   OutText <- c(OutText,CurLine) 
  }
 }
#
#	And spatial structure
#
 Compnt <- 10
 if (x$spmod == 1) {
  CurLine <- paste(formatC(Compnt,digits=0,width=5,format="d"),
                   formatC(0,digits=4,width=10,format="f"),
                   formatC(x$spmod,digits=0,width=5,format="d"),
                   paste(rep(" ",10),collapse=""),
                   "Observed residual correlation structure",sep="")
  OutText <- c(OutText,CurLine) 
 } else if (x$Np[Compnt] > 0) {
  for (i in 1:x$Np[Compnt]) {
   CurLine <- paste(formatC(Compnt,digits=0,width=5,format="d"),
                    formatC(x$rho[i],digits=4,width=10,format="f"),
                    formatC(x$spmod,digits=0,width=5,format="d"),
                    formatC(i,digits=0,width=5,format="d"),
                    paste(rep(" ",5),collapse=""),"Parameter ",
                    i," in spatial dependence model",sep="")
   OutText <- c(OutText,CurLine) 
  }
 }   
 write(OutText,file=file[1])
#
#       Write dispersion info for joint mean-variance models
#
 if (x$model.type == "normal-heteroscedastic") {
  write.modeldef(x$dispersion,file=file[2],check.file=FALSE,mean.only=TRUE)
 }
 invisible(NULL)
} 
##############################################################################
##############################################################################
##############################################################################
write.GLCexternal <- function(x,var.names,file,check.file=TRUE) {
##############################################################################
#
#       To write a data file of "external predictors" that can subsequently
#       be read by the model fitting and simulation routines. Arguments:
#
#       x               A data frame containing predictor information at a
#                       daily, monthly or yearly time scale. If the time
#                       scale is daily then x must have columns named
#                       "Year", "Month" and "Day" (case-sensitive) in 
#                       addition to the predictor data. If monthly, it must
#                       have columns "Year" and "Month"; and if yearly 
#                       it must have a column "Year". All other columns
#                       will be taken to contain data values.
#       var.names       A character vector of names for the non-date 
#                       variables. This defaults to the corresponding
#                       variable names in x.
#       file            Name of output file.
#       check.file      Check whether output file exists?
#
##############################################################################
 if (!is.data.frame(x)) stop("x must be a data frame")
 tscale.labels <- c("Year","Month","Day")
 dates.present <- tscale.labels %in% names(x)
 if (all(dates.present)) {
  tscale <- "daily"
  HeadType <- 4
 } else if (all(dates.present[1:2])) {
  tscale <- "monthly"
  HeadType <- 3
 } else if (dates.present[1]) {
  tscale <- "yearly"
  HeadType <- 2
 } else {
  stop("I can't figure out the required time scale from the variable names in x")
 }
 if (check.file) {
  if (!file.OK(file,"new")) return()
 }
#
#       Write header section
#
 OutText <- capture.output(invisible(.Fortran("WriteHeader",
                                              HeadType=as.integer(HeadType),
                                              PACKAGE="Rglimclim")))
 pred.cols <- (1:ncol(x))[!(names(x) %in% tscale.labels)]
 Npreds <- length(pred.cols)
 if (missing(var.names)) {
  vnames <- names(x)[pred.cols]
 } else {
  if (length(var.names) == Npreds) {
   vnames <- var.names
  } else {
   stop("Length of var.names doesn't match number of predictors in x")
  }
 }
 OutText <- c(OutText,Npreds,vnames)
 CurLine <- paste(paste(rep("*",25),collapse=""),
                  "END OF PREDICTOR DEFINITION",
                  paste(rep("*",25),collapse=""),sep="")
 OutText <- c(OutText,CurLine)
 write(OutText,file=file)
#
#       Now data section (replace missing values with -9999.9)
#
 if (tscale=="daily") {
  OutText <- formatC(10000*x$Year+100*x$Month+x$Day,
                           digits=0,width=8,format="d")
  OutText <- paste(OutText," ",sep="")
 } else if (tscale=="monthly") {
  OutText <- paste(formatC(x$Year,digits=0,width=4,format="d"),
                   formatC(x$Month,digits=2,width=2,format="d"))
  OutText <- paste(OutText,"  ",sep="")
 } else if (tscale=="yearly") {
  OutText <- formatC(x$Year,digits=0,width=4,format="d")
  OutText <- paste(OutText,"     ",sep="")
 }
 xx <- x
 xx[is.na(xx)] <- -9999.9
 for (i in 1:Npreds) {
  OutText <- paste(OutText,formatC(x[,pred.cols[i]],width=10,format="f"),sep="")
 }
 write(OutText,file=file,append=TRUE)
 invisible(NULL)
} 
##############################################################################
##############################################################################
##############################################################################
write.GLCdata <- function(x,date.cols=1:3,site.col=4,data.cols,
                          missval=-99.99,file,check.file=TRUE) {
##############################################################################
#
#       To write a daily data file of response variables that can subsequently
#       be read by the model fitting and simulation routines. Arguments:
#
#       x               A data frame containing the data to be written.
#       date.cols       A numeric vector of length 3, containing the 
#                       numbers of the columns containing respectively 
#                       the values of year, month and day in x. 
#       site.col        Number of column containing 4-character site codes
#       data.cols       Vector of column numbers containing data. Defaults
#                       to all columns except those containing date and site
#                       information. 
#       missval         Code to use for missing values (this will be used to 
#                       replace NA values in x)
#       file            Name of output file.
#       check.file      Check whether output file exists?
#
##############################################################################
 if (!is.data.frame(x)) stop("x must be a data frame")
#
#       Check output file if requested
#
 if (check.file) {
  if (!file.OK(file,"new")) return()
 }
 if (any(nchar(as.character(x[,site.col])) > 4)) {
  stop("x must contain 4-character site codes in column site.col")
 }
#
#       Identify date and site columns, and replace NA values with 
#       missing value code
#
 Dat.Id <- x[,c(date.cols,site.col)]
 if (missing(data.cols)) {
  datcols <- (1:ncol(x))[-c(date.cols,site.col)]
 } else {
  datcols <- data.cols
 }
 z <- x[,datcols,drop=FALSE]
 all.missing <- apply(z,MARGIN=1,FUN=function(x) {all(is.na(x))})
 z[is.na(z)] <- missval
 z <- z[!all.missing,,drop=FALSE]
 Dat.Id <- Dat.Id[!all.missing,,drop=FALSE]
#
#       Write data
#
 OutText <- paste(formatC(Dat.Id[,1],digits=0,width=4,format="d"),
                  formatC(Dat.Id[,2],digits=0,width=2,format="d"),
                  formatC(Dat.Id[,3],digits=0,width=2,format="d"),
                  format(Dat.Id[,4],width=4,justify="right"),sep="")
 for (i in 1:ncol(z)) {
  OutText <- paste(OutText,formatC(z[,i],width=6,digits=2,format="f"),sep="")
 }
 write(OutText,file=file)
 invisible(NULL)
} 
##############################################################################
##############################################################################
##############################################################################
read.GLCdata <- function(file,missval=-99.99) {
##############################################################################
#
#       To read a daily data file of response variables (e.g. for
#       postprocessing of simulation results produced by GLCsim). 
#       This is considerably quicker than using read.fwf().
#       Arguments:
#
#       file            The filename
#       missval         Code to use for missing values (these will be 
#                       replaced by NA)
#
##############################################################################
 if (!file.OK(file,"old")) return()
 z <- scan(file,what="",sep="\n")
 nvars <- nchar(z)
 if (length(table(nvars)) > 1) {
  stop(paste("Rows of",file,"have different lengths"))
 }
 nvars <- (nvars[1]-12)/6
 field.lengths <- c(4,2,2,4,rep(6,nvars))
#
# Need to replicate the start and end positions of each field to 
# match the length of z; to save on memory, store them explicitly
# as integers
#
 ends <- rep(as.integer(cumsum(field.lengths)),each=length(z))
 starts <- ends+as.integer(1)-rep(as.integer(field.lengths),each=length(z))
 z <- matrix(substring(z,first=starts,last=ends),ncol=nvars+4)
 z <- as.data.frame(z,stringsAsFactors=FALSE)
 for (i in (1:ncol(z))[-4]) z[,i] <- as.numeric(z[,i])
 z[,-4][z[,-4]==missval] <- NA
 names(z) <- c("Year","Month","Day","Site",paste("Var",1:nvars,sep=""))
 z
} 
##############################################################################
##############################################################################
##############################################################################
char.trunc <- function(x,width) {
##############################################################################
#
#       Truncates the string x to an exact width
#
##############################################################################
 formatC(substr(x,1,width),width=width,flag="-")
}
##############################################################################
##############################################################################
##############################################################################
print.GLC.modeldef <- function(x,scr.width=NULL,global.warn=TRUE,
                               mean.only=FALSE,which.se="robust",...) {
##############################################################################
#
#       print() function for objects of class GLC.modeldef
#
##############################################################################
 if (is.null(scr.width)) scr.width <- as.numeric(options("width"))-2 
 title.print(x$model.title)
 if (!is.null(x$which.response)) {
  cat(paste("\nResponse variable:",x$var.names[x$which.response],"\n\n"))
 }
#
#       Set up column widths for printing, etc. 
#
 colnames.ct <- c(" ","Coefficient","Std Err","T-stat","Pr(|T|>t)")  
 test.type <- ifelse(x$model.type %in% c("normal","gamma"),"t","z")
 if (test.type=="z") colnames.ct[4:5] <- c("Z-stat","Pr(|Z|>z)")
 if (is.null(x$cov.robust)) {
  nc.ct <- 2
 } else {
  nc.ct <- 5
  if (which.se == "robust") {
   covmat <- x$cov.robust 
  } else if (which.se == "naive") {
   covmat <- x$cov.naive
  } else {
   stop("which.se must be either 'naive' or 'robust'")
  }
 }
 label.widths <- scr.width - ((nc.ct-1)*8) - 4
 prevwt.flag <- FALSE
#
#       Np in the x object does *not* include the constant
#       in the # of parameters
#
 np <- x$Np
#
#       Main effects
#
 idx <- 1:(1+np[4])
 tmp.table <- data.frame(matrix(nrow=length(idx),ncol=nc.ct))
 names(tmp.table) <- colnames.ct[1:nc.ct]
 if (np[4] > 0) {
  rownames(tmp.table) <- c("",as.character(1:np[4]))
 } else {
  rownames(tmp.table) <- ""
 } 
 num.widths <- ceiling(log(max(1,np[4]),base=10))
 tmp.table[,1] <- substr(x$model.labels$xlabels,1,
                                        label.widths-num.widths-nc.ct)
 tmp.table[,1] <- format(tmp.table[,1],width=label.widths-num.widths-nc.ct,
                                                           justify="right")
 tmp.table[,2] <- x$beta[idx]
 if (nc.ct > 2) {
  tmp.table[,3] <- sqrt(diag(covmat)[idx])
  tmp.table[,4] <- tmp.table[,2] / tmp.table[,3]
  if (test.type=="z") {
   tmp.table[,5] <- 2*pnorm(abs(tmp.table[,4]),lower.tail=FALSE)
  } else {
   tmp.table[,5] <- 2*pt(abs(tmp.table[,4]),
                         df=x$df.resid,lower.tail=FALSE)  
  }
  tmp.table[,3:4] <- formatC(as.matrix(tmp.table[,3:4]),
                             format="f",width=8,digits=4)
  tmp.table[,5] <- format.pval(tmp.table[,5],digits=4)
 }
 tmp.table[,2] <- formatC(tmp.table[,2],format="f",width=11,digits=4)
 cat("Main effects:\n-------------\n")
 print(tmp.table,row.names=TRUE)
#
#       Two-way interactions
#
 if (np[5] > np[4]) {
  idx <- (2+np[4]):(1+np[5])
  tmp.table <- data.frame(matrix(nrow=2*length(idx),ncol=nc.ct))
  names(tmp.table) <- colnames.ct[1:nc.ct]
  coef.rows <- 2*(idx-np[4]-1)-1
  tmp.text <- x$model.labels$xlabels[x$twoway[idx-1,1]+1]
  tmp.table[coef.rows,1] <- 
                  paste("    ",substr(tmp.text,1,label.widths-5-nc.ct))
  tmp.text <- x$model.labels$xlabels[x$twoway[idx-1,2]+1]
  tmp.table[coef.rows+1,1] <- 
                  paste("with",substr(tmp.text,1,label.widths-5-nc.ct))
  tmp.table[,1] <- format(tmp.table[,1],
                          width=label.widths-nc.ct,justify="right")
  tmp.table[coef.rows,2] <- x$beta[idx]
  NA.pos <- is.na(tmp.table[,-1])
  if (nc.ct > 2) {
   tmp.table[coef.rows,3] <- sqrt(diag(covmat)[idx])
   tmp.table[coef.rows,4] <- tmp.table[coef.rows,2] / 
                             tmp.table[coef.rows,3]
   if (test.type=="z") {
    tmp.table[coef.rows,5] <- 
              2*pnorm(abs(tmp.table[coef.rows,4]),lower.tail=FALSE)
   } else {
    tmp.table[coef.rows,5] <- 2*pt(abs(tmp.table[coef.rows,4]),
                                   df=x$df.resid,lower.tail=FALSE)  
   }
   NA.pos <- is.na(tmp.table[,-1])
   tmp.table[,3:4] <- formatC(as.matrix(tmp.table[,3:4]),
                              format="f",width=8,digits=4)
   tmp.table[,5] <- format.pval(tmp.table[,5],digits=4)
  }
  tmp.table[,2] <- formatC(tmp.table[,2],format="f",width=11,digits=4)
  tmp.table[,-1][NA.pos] <- ""
  cat("\nTwo-way interactions:\n---------------------\n")  
  print(tmp.table,row.names=FALSE)  
 }
#
#       Three-way interactions
#
 if (np[6] > np[5]) {
  idx <- (2+np[5]):(1+np[6])
  tmp.table <- data.frame(matrix(nrow=3*length(idx),ncol=nc.ct))
  names(tmp.table) <- colnames.ct[1:nc.ct]
  coef.rows <- 3*(idx-np[5]-1)-2
  tmp.text <- x$model.labels$xlabels[x$threeway[idx-1,1]+1]
  tmp.table[coef.rows,1] <- 
                  paste("    ",substr(tmp.text,1,label.widths-5-nc.ct))
  tmp.text <- x$model.labels$xlabels[x$threeway[idx-1,2]+1]
  tmp.table[coef.rows+1,1] <- 
                  paste("with",substr(tmp.text,1,label.widths-5-nc.ct))
  tmp.text <- x$model.labels$xlabels[x$threeway[idx-1,3]+1]
  tmp.table[coef.rows+2,1] <- 
                  paste("and",substr(tmp.text,1,label.widths-4-nc.ct))
  tmp.table[,1] <- format(tmp.table[,1],
                          width=label.widths-nc.ct,justify="right")
  tmp.table[coef.rows,2] <- x$beta[idx]
  NA.pos <- is.na(tmp.table[,-1])
  if (nc.ct > 2) {
   tmp.table[coef.rows,3] <- sqrt(diag(covmat)[idx])
   tmp.table[coef.rows,4] <- tmp.table[coef.rows,2] / 
                             tmp.table[coef.rows,3]
   if (test.type=="z") {
    tmp.table[coef.rows,5] <- 
              2*pnorm(abs(tmp.table[coef.rows,4]),lower.tail=FALSE)
   } else {
    tmp.table[coef.rows,5] <- 2*pt(abs(tmp.table[coef.rows,4]),
                                   df=x$df.resid,lower.tail=FALSE)  
   }
   NA.pos <- is.na(tmp.table[,-1])
   tmp.table[,3:4] <- formatC(as.matrix(tmp.table[,3:4]),
                                      format="f",width=8,digits=4)
   tmp.table[,5] <- format.pval(tmp.table[,5],digits=4)
  }
  tmp.table[,2] <- formatC(tmp.table[,2],format="f",width=11,digits=4)
  tmp.table[,-1][NA.pos] <- ""
  cat("\nThree-way interactions:\n-----------------------\n")  
  print(tmp.table,row.names=FALSE)  
 }
#
#	If there are *no* nonlinear parameters at all in the model, 
#       np[7] is zero; if there *are* nonlinear parameters but 
#       they're treated as known, np[7]=np[6], but we still need to say
#       something about them!
#
 if ( (np[7]-np[6] >= 0) & (np[6] != 0) ) {
  cat(paste("\nParameters in nonlinear transformations:",
            "\n----------------------------------------\n",sep=""))
#
#	For nonlinear parameters, it's easier to go through the covariates
#	themselves looking for nonlinear parameters - otherwise we end up
#	doubling the output if there are 2 parameters if we're not
#	careful! The way to find out is if a value for the first THETA has
#	been specified. The additional IF is in case there isn't a second
#	THETA. We also flag whether or not there's a parameter associated 
#       with a distance-weighted average of previous days' values - if
#       so we'll add a footnote to tell the user where the distances are 
#       coming from. Here, the formatting is quite fiddly so there is
#       not much to be gained from avoiding loops.
#
  idx <- 1:np[4]
  nl.covs <- (x$theta[idx,1] < 1e8)
  wtschm <- xfrm <- rep(NA,np[4])
  wtschm[nl.covs] <- floor(x$covcode[idx][nl.covs]/10000) %% 10
  xfrm[nl.covs] <- floor(x$covcode[idx][nl.covs]/1000) %% 10
  if (np[4] > np[3]) idx <- (1+np[3]):np[4] else idx <- NULL
  if (!is.null(idx)) {
   prevwt.flag <- any(wtschm[idx] > 1, na.rm=TRUE)
#
#       The next line returns a vector containing TRUE for any 
#       covariate to which a set of weighting parameters is 
#       attached (for any of the variables). The next line 
#       works by recycling the elements of idx columnwise; the 
#       outcome of this calculation is a 2-d array containing at
#       most one TRUE in each row (because the columns 
#       correspond to different variables, and a weighting
#       parameter can't simultaneously belong to more than one variable).
#
   tmp <- (idx == x$prevwt.idx[wtschm[idx],1,,drop=FALSE])
   if (any(tmp,na.rm=TRUE)) { 
    wanted.rows <- apply(tmp,MARGIN=1,FUN=any)
   } else {
    wanted.rows <- rep(FALSE,nrow(tmp))
   }
   nl.covs[idx][!wanted.rows] <- FALSE
  }
  for (i in (1:np[4])[nl.covs]) {
#
#       i+1 in next line is because the first label is for the 
#       constant
#
   cat(paste(substr(x$model.labels$xlabels[i+1],1,scr.width-5),
             ":\n",sep=""))
   pars.wanted <- (x$theta[i,] < 1e8)
   tmp.text <- format(paste("    ",
               substr(x$model.labels$nonlin.txt[i,pars.wanted],
                          1,scr.width-39),sep=""),width=scr.width-34,
                      justify="right")
   nl.known <- x$nlstat[i,pars.wanted] == 0
   tmp.text[nl.known] <- paste(tmp.text,": ",
                               formatC(x$theta[i,pars.wanted],
                                      width=11,digits=4,format="f"),
                               " (prespecified)",sep="")
   if (is.null(x$cov.naive)) {
    tmp.text[!nl.known] <- paste(tmp.text[!nl.known],": ",
                                 formatC(x$theta[i,pars.wanted],
                                   width=11,digits=4,format="f"),sep="")
   } else {
    nlpars.wanted <- !is.na(x$covcode) & (x$covcode %/% 1000 == i) &
                     (1:length(x$covcode)) %in% (np[6]+1):np[7]
    nlpars.wanted <- (1:length(nlpars.wanted))[nlpars.wanted]
#
#       In the next line, have to remember that covmat is indexed to 
#       align with beta (which starts from 0 in the Fortran code) but
#       theta is not - hence adding 1 to the subscript for covmat. Why
#       do I do this to myself?!
#   
    tmp.text[!nl.known] <- paste(tmp.text[!nl.known],": ",
                                 formatC(x$theta[i,pars.wanted],
                                      width=11,digits=4,format="f"),
                                 " (Std Err: ",
                                 formatC(sqrt(diag(covmat)[1+nlpars.wanted]),
                                      width=9,digits=4,format="f"),
                                 ")",sep="")   
   }
   cat(tmp.text,sep="\n")
  }
 }
 if (mean.only) return(invisible(x))
#
#       Global quantities
#
 if (np[8] > 0) {
  cat("\nGlobal quantities:\n------------------\n")
  idx <- 1:np[8]
  glob.id <- floor(x$global.codes[idx] / 1000)
  tmp.text <- paste("     ",
      substr(x$model.labels$global.txt[glob.id],1,scr.width-16),": ",
      formatC(x$global.vals[glob.id],
              format="f",width=11,digits=4),sep="")
  cat(tmp.text,sep="\n")
 }
#
#       Dispersion parameter
#
 if (x$model.type == "normal-heteroscedastic") {
  if (inherits(x$dispersion,"GLC.modeldef")) {
   cat(paste("\n",paste(rep("<",scr.width),collapse = ""),"\n",sep = ""))
   cat("Dispersion structure:\n---------------------\n")
   print(x$dispersion,global.warn=FALSE,scr.width=scr.width,mean.only=TRUE)
   cat(paste(paste(rep(">",scr.width),collapse = ""),"\n", sep = ""))
  }
 } else if (is.numeric(x$dispersion)) {
  if (isTRUE(x$dispersion > 0)) {
   cat(paste("\nDispersion parameter: ",
             formatC(x$dispersion,format="f",width=11,digits=4),
             "\n",sep=""))
  } else {
   cat("\nNo dispersion parameters defined\n") 
  }
 } else {
  cat("\nNo dispersion parameters defined\n")
 }
#
#       And spatial structure
#
 cat("\nSpatial dependence structure:\n-----------------------------\n")
 cat(substr(x$model.labels$corstruct.txt,1,scr.width-5),sep="\n")
 if (x$spmod != 0 & x$Np[10] > 0) {
  idx <- 1:x$Np[10]
  tmp.text <- paste("     ",
      substr(x$model.labels$corpars.txt[idx],1,scr.width-16),
      ": ",formatC(x$rho[idx],format="f",width=11,digits=4),sep="")
  cat(tmp.text,sep="\n") } 
#
#       Finally, notes if relevant
#
 if (prevwt.flag | (x$spmod > 2 & x$spmod < 20)) {
  cat(paste("\n**** NOTE: distance-dependent weights (for averages ",
            "of previous days'\nvalues and / or correlations) are ",
            "based on distances calculated from\n     ",
            substr(x$siteinfo$Attribute.names[1],1,scr.width-10),
            "\nand ",
            substr(x$siteinfo$Attribute.names[2],1,scr.width-10),
            "\n",sep=""))
 } 
 if (np[8] < 1 & x$model.type %in% c("logistic","gamma") & global.warn) {
  warning("\nNo global quantities (trace thresholds etc.) defined")
 }
 invisible(x)
}
##############################################################################
##############################################################################
##############################################################################
logLik.GLC.modeldef <- function(object,...) {
#
#       Log-likelihood extraction for Rglimclim model fits
#
 if (!inherits(object,"GLC.modeldef")) {
  stop("argument must be the result of a call to GLCfit")
 }
 if (is.null(object$LogLik)) return(NULL)
 z <- object$LogLik
 attr(z,which="df") <- object$df
 attr(z,which="nobs") <- object$nobs
 class(z) <- "logLik"
 z
}
##############################################################################
##############################################################################
##############################################################################
title.print <- function(model.title) {
#
#       Utility function to print model titles
#
 if (is.null(model.title)) {
  print.tit <- "UNTITLED MODEL"
 } else {
  print.tit <- model.title
 }
 cat(paste(print.tit, "\n", paste(rep("=", nchar(print.tit)), 
     collapse = ""), "\n", sep = ""))
 NULL
}
##############################################################################
##############################################################################
##############################################################################
summary.GLC.modeldef <- function(object, tables=c("month","site","year"), ...) {
#
#       summary function for Rglimclim model fits. This just
#       extracts the relevant components from the fit, and 
#       returns a list that is suitable for passing to 
#       the corresponding print method
#
 if (!inherits(object,"GLC.modeldef")) {
  stop("argument must be the result of a call to GLCfit")
 }
 if (is.null(object$nobs)) {
  warning(paste("No summary available for object",substitute(object),
                "- run it through GLCfit first"))
  return(NULL)
 }
#
#       Originally had a separate print method, but that was silly. 
#       Code below is cribbed from the old print method so, to
#       avoid errors in incorporating this into the summary 
#       method, the next line creates a list of the form that
#       was originally passed to the print method
#
 x <- list(title=object$model.title,type=object$model.type,
           N=object$nobs,df=object$df,df.resid=object$df.resid,
           LogLik=object$LogLik,deviance=object$Deviance,
           dispersion=object$dispersion,Residuals=object$Residuals)
 title.print(x$title)
 if (!is.null(object$which.response)) {
  cat(paste("\nResponse variable:",object$var.names[object$which.response],"\n\n"))
 }
 cat(paste("Model of type '",x$type,"', fitted to ",
           x$N," observations\n\n",sep=""))
 cat(paste(formatC(paste("# of parameters estimated:",x$df),
                   width=39,flag="-"),
           formatC(paste("Independence log-likelihood:",round(x$LogLik,2)),
                   width=30,flag="-"),
           "\n",sep=""))
 cat(paste(formatC(paste("Residual degrees of freedom:",x$df.resid),
                   width=39,flag="-"),
           formatC(paste("Deviance:",round(x$deviance,2)),width=30,flag="-"),
           "\n",sep=""))
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<              Model-specific section here             <<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 if (!is.null(x$Residuals)) {
  if (x$type == "logistic") {
   cat(paste("Mean squared error (mean Brier score): ",
             round(x$Residuals$MSE,4),"\n\n",sep=""))
  } else {
   cat(paste(formatC(paste("Mean of observations:",round(x$Residuals$ybar,3)),
                     width=39,flag="-"),
             formatC(paste("Std dev of observations:",round(x$Residuals$sy,3)),
                     width=30,flag="-"),
             "\n",sep="")) 
   if (x$type == "normal") {
    cat(paste(formatC(paste("Mean error (observed-fitted):",
                            format(round(x$Residuals$errbar,4),nsmall=4)),
                            width=39,flag="-"),
              formatC(paste("Residual standard deviation:",
                            round(sqrt(x$dispersion),4)),width=30,flag="-"),
              "\n",sep="")) 
   } else {
    cat(paste(formatC(paste("Mean error (observed-fitted):",
                            format(round(x$Residuals$errbar,4),nsmall=4)),
                            width=39,flag="-"),
              formatC(paste("Root mean squared error:",
                            round(sqrt(x$Residuals$MSE),4)),width=30,flag="-"),
              "\n",sep=""))   
   }
   cat(paste("Percentage of variance explained: ",
             format(round(100*x$Residuals$Rsq,1),nsmall=1),"\n",sep=""))
  }
 }
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 if (x$type == "gamma") {
  cat(paste("\nCommon shape parameter for gamma distributions: ",
            round(1/x$dispersion,3),"\n",sep=""))
 } else if (x$type == "logistic") {
  cat("\nNo dispersion parameters estimated for this model\n")
 }
 if (!is.null(x$Residuals)) {
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#<<<<<              Model-specific section here             <<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  if (x$type == "logistic") {
   cat(paste("\nPearson residuals: mean ",
             format(round(x$Residuals$Pearson$Mean,4),nsmall=4)," (std err ",
             format(round(x$Residuals$Pearson$SEmean,4),nsmall=4),
             "), standard deviation ",
             format(round(x$Residuals$Pearson$SD,4),nsmall=4),"\n",sep=""))
   tmp <- formatC(x$Residuals$Prob.table,format="f",digits=3)
   tmp[3,] <- formatC(x$Residuals$Prob.table[3,],format="f",digits=0)
   cat(paste("\nOccurrence frequencies vs forecasts:\n",
             "------------------------------------\n",
             "                                      ",
             "Forecast decile\n",sep=""))
   print(as.data.frame(tmp))
  } else if (x$type == "gamma") {
   cat(paste("\nPearson residuals: mean ",
             format(round(x$Residuals$Pearson$Mean,3),nsmall=3)," (std err ",
             format(round(x$Residuals$Pearson$SEmean,3),nsmall=3),
             "), std dev ",
             format(round(x$Residuals$Pearson$SD,3),nsmall=3),
             " (expected ",format(round(sqrt(x$dispersion),3),nsmall=3),
             ")\n",sep=""))  
   cat(paste("Anscombe residuals: mean = ",
             format(round(x$Residuals$Anscombe$Obs.mean,3),nsmall=3),
             " (expected ",
             format(round(x$Residuals$Anscombe$Exp.mean,3),nsmall=3),
             "), sd = ",
             format(round(x$Residuals$Anscombe$Obs.sd,3),nsmall=3),
             " (expected ",
             format(round(x$Residuals$Anscombe$Exp.sd,3),nsmall=3),
             ")\n",sep=""))
  }
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ntabs <- length(tables)
  if (ntabs > 0) {
   for (i in 1:ntabs) {
    if (tables[i] == "month") {
     cat(paste("\nPearson residual summaries by month\n",
               "-----------------------------------\n",sep=""))
     wanted.rows <- x$Residuals$Pearson$Month.table[,1] > 0
     print(round(x$Residuals$Pearson$Month.table[wanted.rows,],digits=3))
    } else if (tables[i] == "site") {
     cat(paste("\nPearson residual summaries by site\n",
               "----------------------------------\n",sep=""))
     tmp <- x$Residuals$Pearson$Site.table
     tmp$Name <- substr(as.character(tmp$Name),1,30)
     tmp[,-2] <- round(tmp[,-2],3)
     wanted.rows <- tmp[,1] > 0
     print(tmp[wanted.rows,])
    } else if (tables[i] == "year") {
     cat(paste("\nPearson residual summaries by year\n",
               "----------------------------------\n",sep=""))
     wanted.rows <- x$Residuals$Pearson$Year.table[,1] > 0
     print(round(x$Residuals$Pearson$Year.table[wanted.rows,],3))  
    }
   }
   if (any(!is.na(match(c("month","year"),tables)))) {
    cat(paste("\n***NOTE*** Standard errors for monthly and",
              "annual tables above are\n           adjusted",
              "for possible inter-site dependence.\n"))
   }
  }
 }
 invisible(NULL)
}
##############################################################################
##############################################################################
##############################################################################
plot.GLC.modeldef <- function(x,which.plots=1:2,sd.plots=TRUE,
              site.options=list(add.to.map=FALSE,scale=NULL,axis.labels=NULL,
                                coord.cols=1:2,site.labels="all"),
              titles=TRUE,distance.units,plot.cols=gray(c(0.4,0)),...) {
 if (!inherits(x,"GLC.modeldef")) {
  stop("argument must be the result of a call to GLCfit")
 }
 if (is.null(x$Residuals)) {
  warning(paste("Object",substitute(x),"contains no residuals",
                "- nothing to plot"))
  return(NULL)
 }
#
#       Next if {} block is a way to allow the user to specify just
#       part of the site.options() list
#
 if (!missing(site.options)) {
  if (!is.list(site.options)) stop("site.options must be a list")
  site.opt.default <- eval(formals(plot.GLC.modeldef)$site.options)
  els.to.replace <- names(site.opt.default) %in% names(site.options)
  site.options <- c(site.opt.default[!els.to.replace],site.options)
 }
 monthly.data <- x$Residuals$Pearson$Month.table
 yearly.data <- x$Residuals$Pearson$Year.table
 names(monthly.data) <- names(yearly.data) <- c("Ndays","Mean","SD","SEmean")
 monthly.data$Month <- 1:12
 yearly.data$Year <- as.numeric(rownames(yearly.data))
 site.res <- x$Residuals$Pearson$Site.table
 names(site.res) <- c("Ndays","Name","Mean","SD","SEmean")
 site.res$SITE <- rownames(site.res)
 sites <- as.data.frame(x$siteinfo$Attribute.values)
 sites$SITE <- rownames(sites)
#
#       Theoretical std dev of Pearson residuals is model-dependent
#
 pr.sd <- 1
 if (x$model.type %in% c("normal","gamma")) {
  pr.sd <- sqrt(x$dispersion)
 } 
 if (1 %in% which.plots) {
  monthly.ul95 <- 1.96*monthly.data$SEmean
  monthly.ll95 <- -monthly.ul95
  plotlims <- 1.05*max(monthly.ul95,abs(monthly.data$Mean),na.rm=TRUE)
  plotlims <- c(-plotlims,plotlims)
  plot(monthly.data$Month,monthly.data$Mean,type="l",xlab="Month",
	ylab="Mean",ylim=plotlims,xlim=c(0.5,12.5),las=1,tck=-0.05,
	cex=0.9,lwd=1)
  lines(monthly.data$Month,monthly.ul95,lty=2,lwd=1)
  lines(monthly.data$Month,monthly.ll95,lty=2,lwd=1)
  axis(1,at=1:12,labels=F,pos=0,tck=0.0)
  if (titles) title("Monthly residual means",cex=0.6)
  
  if (sd.plots) {
   plotlims <- c(0,1.1*max(pr.sd,monthly.data$SD,na.rm=TRUE))
   plot(monthly.data$Month,monthly.data$SD,type="l",xlab="Month",
	ylab="Std Dev",ylim=plotlims,xlim=c(0.5,12.5),las=1,tck=-0.05,
	cex=0.9,lwd=1)
   axis(1,at=1:12,labels=F,pos=pr.sd,tck=0.0)
   if (titles) title("Monthly residual standard deviations",cex=0.6)
  }
 }

 if (2 %in% which.plots) {
  yearly.ul95 <- 1.96*yearly.data$SEmean
  yearly.ll95 <- -yearly.ul95
  plotlims <- 1.05*max(yearly.ul95,abs(yearly.data$Mean),na.rm=TRUE)
  plotlims <- c(-plotlims,plotlims)
  plot(yearly.data$Year,yearly.data$Mean,type="l",xlab="Year",
	ylab="Mean",ylim=plotlims,las=1,tck=-0.05,cex=0.9,lwd=1)
  lines(yearly.data$Year,yearly.ul95,lty=2,lwd=1)
  lines(yearly.data$Year,yearly.ll95,lty=2,lwd=1)
  lines(range(yearly.data$Year),c(0,0))
  if (titles) title("Annual residual means",cex=0.6)

  if (sd.plots) {
   plotlims <- c(0,1.1*max(pr.sd,yearly.data$SD,na.rm=TRUE))
   plot(yearly.data$Year,yearly.data$SD,type="l",xlab="Year",
	ylab="Std Dev",ylim=plotlims,las=1,tck=-0.05,cex=0.9,lwd=1)
   lines(range(yearly.data$Year),c(pr.sd,pr.sd))
   if (titles) title("Annual residual standard deviations",cex=0.6)
  }
 }

 if (3 %in% which.plots) {
  site.res <- merge(sites, site.res, by=c("SITE"))
  site.rec <-  site.res[site.res$Nday > 0,]
  xcol <- 1 + site.options$coord.cols[1] # Add 1 to column numbers because
  ycol <- 1 + site.options$coord.cols[2] # merge inserted SITE column

  if (!site.options$add.to.map) {
   if (is.null(site.options$axis.labels)) {
    ax.labs <- x$siteinfo$Attribute.names[site.options$coord.cols]
   } else {
    ax.labs <- site.options$axis.labels
    if (length(ax.labs) !=2) {
     stop("site.options$axis.labels should be a vector of length 2")
    }
   }
   plot(site.res[,xcol],site.res[,ycol],type="n",
                                      xlab=ax.labs[1],ylab=ax.labs[2])
   box(lwd=2)
  }

  std.res <- site.rec$Mean / site.rec$SEmean

#
#       Default symbol size aims to give a reasonable visual 
#       impression - this depends upon the co-ordinate
#       units (numeratore in expression below), the number of 
#       sites and the magnitude of the residuals themselves
#       (denominator)
#
  symscl <- mean(diff(range(site.res[,xcol])),
                  diff(range(site.res[,ycol]))) / 
             (max(20,sqrt(nrow(site.res))) * mean(abs(std.res)))
  if (!is.null(site.options$scale)) symscl <- symscl * site.options$scale

#
#       Positive residuals ...
#
  pos.res <- (site.rec$Mean > 0)
  x.wanted <- site.rec[pos.res,xcol]
  y.wanted <- site.rec[pos.res,ycol]
  z.wanted <- std.res[pos.res]
  points(x.wanted,y.wanted,pch=20,col="black",cex=0.8)

  sigres <- (z.wanted > 1.96)

  if (any(!sigres)) {
   circles(x.wanted[!sigres],y.wanted[!sigres],z=symscl*sqrt(z.wanted[!sigres]),lwd=1)
  }
  if (any(sigres)) {
   circles(x.wanted[sigres],y.wanted[sigres],z=symscl*sqrt(z.wanted[sigres]),lwd=3)
   if (site.options$site.labels == "significant") {
    text(x.wanted[sigres],y.wanted[sigres],(site.rec$SITE[pos.res])[sigres])
   }
  }
  if (site.options$site.labels == "all") {
   text(x.wanted,y.wanted,site.rec$SITE[pos.res])
  }

#
#       ... and negative.
#
  x.wanted <- site.rec[!pos.res,xcol]
  y.wanted <- site.rec[!pos.res,ycol]
  z.wanted <- std.res[!pos.res]
  points(x.wanted,y.wanted,pch=20,col="black",cex=0.8)

  sigres <- (z.wanted < -1.96)

  if (any(!sigres)) {
   circles(x.wanted[!sigres],y.wanted[!sigres],z=symscl*sqrt(-z.wanted[!sigres]),
           lwd=1,lty=2)
  }
  if (any(sigres)) {
   circles(x.wanted[sigres],y.wanted[sigres],z=symscl*sqrt(-z.wanted[sigres]),
           lwd=3,lty=2)
   if (site.options$site.labels == "significant") {
    text(x.wanted[sigres],y.wanted[sigres],(site.rec$SITE[!pos.res])[sigres])
   }
  }
  if (site.options$site.labels == "all") {
   text(x.wanted,y.wanted,site.rec$SITE[!pos.res])
  }
  if (titles) title("Mean Pearson residuals by site",cex=0.6)
 }
 if (4 %in% which.plots) {
#
#       For quantile-quantile plots, plot largest and smallest 1000 points
#       along with the middle 1000 points (with large datasets, if you plot
#       all of the data then everything gets overplotted and the only effect 
#       is to increase the file size and slow down printing.
#
  if (x$model.type == "logistic") {
   warning(paste("Q-Q plot (plot 4) cannot be produced for logistic regression",
                 "model -\n  use table of observed and expected proportions in",
                 "summary() of model\n  to check probability structure."))
  } else if (is.null(x$filenames$Residuals)) {
    warning(paste("Q-Q plot (plot 4) cannot be produced for this model because",
                  "residuals\n  were not saved - rerun GLCfit() with diagnostics=2",
                  "to correct this."))
  } else if (!file.exists(x$filenames$Residuals)) {
    warning(paste("Q-Q plot cannot be produced because residual file",
                  x$filenames$Residuals,"\n  is not present - rerun GLCfit() to",
                  "generate it again."))
  } else {
   resid.data <- scan(x$filenames$Residuals,
                      what=list(Site=character(),Year=numeric(),Month=numeric(),
                                Day=numeric(),Obs=numeric(),Mu=numeric(),
                                Sigma=numeric()),
                      skip=1,nlines=x$nobs,quiet=TRUE)
   resid.data <- data.frame(Y=resid.data$Obs,Mu=resid.data$Mu,
                            Sigma=resid.data$Sigma)
   prob.grid <- ppoints(x$nobs)
   if (x$nobs > 3000) {
    points.to.plot <- c(1:1000,round(seq(1001,x$nobs-1001,length.out=1000),
                        (x$nobs-1000):x$nobs))
    prob.grid <- prob.grid[points.to.plot]
   } else {
    points.to.plot <- 1:x$nobs
   }
   if (x$model.type %in% c("normal","normal-heteroscedastic")) {
    e <- (resid.data$Y-resid.data$Mu) / resid.data$Sigma
    qq <- qnorm(prob.grid)
    xlab = "Quantiles of standard normal"
   } else if (x$model.type == "gamma") {
    e <- resid.data$Y / resid.data$Mu
    qq <- qgamma(prob.grid,shape=1/x$dispersion, rate=1/x$dispersion)
    xlab <- expression(paste("Quantiles of ",Gamma(alpha,alpha)))
    xlab <- substitute(paste("Quantiles of ",Gamma,"(",list(a),",",list(a),")"),
                       list(a=round(1/x$dispersion,2)))
   }
   e <- sort(e)[points.to.plot]
   plot(qq,e,pch=20,xlab=xlab,ylab="Sample quantiles",col=plot.cols[1])
   abline(0,1,lty=2,col=plot.cols[2])
   if (titles) title("Q-Q plot of standardised residuals")
  }
 }
 if (5 %in% which.plots) {
  if (is.null(x$filenames$Correlation)) {
   warning("No correlation file in model definition - plot 5 not produced") 
  } else {
   obs.corr <- read.table(x$filenames$Correlation,header=TRUE)
   obs.corr$Corr[abs(obs.corr$Corr) > 1] <- NA  # Missing values
   d <- sqrt(obs.corr$Xsep^2 + obs.corr$Ysep^2)
   xlab <- "Distance"
   if (!missing(distance.units)) {
    xlab <- paste(xlab," (",distance.units,")",sep="")
   }
#
#       For the plot of inter-site correlations, use opacity on devices
#       that support it so that points corresponding to pairs of sites
#       with very few observations are very faint.
#
   point.cols <- col2rgb(rep(plot.cols[1],length(obs.corr$Corr)),alpha=TRUE)
   point.cols[4,] <- 255 * obs.corr$N / max(obs.corr$N)
   point.cols <- rgb(point.cols[1,],point.cols[2,],point.cols[3,],
                     point.cols[4,],maxColorValue=255)
   plot(d,obs.corr$Corr,pch=20,xlab=xlab,ylab="Correlation",col=point.cols)
   d.grid <- seq(0,max(d),length.out=100)
   if (x$spmod == 2) {
    rho <- rep(x$rho[1],length(d.grid))
   } else if (x$spmod == 3) {
    phi <- x$rho[1]
    rho <- exp(-phi*d.grid)
   } else if (x$spmod == 4) {
    phi <- x$rho[1]
    alpha <- x$rho[2]
    rho <- alpha + (1-alpha) * exp(-phi*d.grid)
   } else if (x$spmod == 5) {
    phi <- x$rho[1]
    kappa <- x$rho[2]
    rho <- exp(-phi*(d.grid^kappa))
   } else if (x$spmod == 6) {
    phi <- x$rho[1]
    kappa <- x$rho[2]
    alpha <- x$rho[3]
    rho <- alpha + (1-alpha) * exp(-phi*(d.grid^kappa))
   } else if (x$spmod == 7) {
    phi <- x$rho[1]
    lambda <- x$rho[2]
    rho <- lambda * exp(-phi*d.grid)
   } else if (x$spmod == 8) {
    phi <- x$rho[1]
    kappa <- x$rho[2]
    lambda <- x$rho[3]
    rho <- lambda * exp(-phi*(d.grid^kappa))
   }
   lines(d.grid,rho,col=plot.cols[2],lwd=2)
   if (titles) title("Inter-site correlations",cex=0.6)
  }
 }
 invisible(NULL)
}
##############################################################################
##############################################################################
##############################################################################
circles <- function(x,y,z,...) {
#
#       To add circles to a plot. This is a replacement for the corresponding 
#       feature in the symbols() command, for which control over line 
#       types / thicknesses etc. doesn't seem to work under Windows. Bleargh.
#       Arguments:
#
#       x       Vector of x-co-ordinates of circle centres
#       y       Vector of y co-ordinates
#       z       Controls sizes of circles to plot. If the x- and y-dimensions
#               are on the same scale, this will correspond to radii in 
#               user co-ordinates. Otherwise, the circles are scaled so that
#               they appear as circles on the plot and hence z doesn't 
#               have a direct interpretation in user co-ordinates
#
 plot.lims <- par()$usr
 plot.size <- par()$pin
 plot.ratio <- ( (plot.lims[2]-plot.lims[1]) * (plot.size[2]) ) /
               ( (plot.lims[4]-plot.lims[3]) * (plot.size[1]) )
 n.points <- 101
 theta <- c(seq(0,2*pi,length.out=n.points-1),NA)
 tmp <- data.frame(x=rep(x,each=n.points),y=rep(y,each=n.points),
                   z=rep(z,each=n.points),theta=rep(theta,length(x)),
                   x.circ=rep(NA,n.points*length(x)),
                   y.circ=rep(NA,n.points*length(x)))
 tmp$x.circ <- tmp$x + (tmp$z * cos(tmp$theta))
 tmp$y.circ <- tmp$y + (tmp$z * sin(tmp$theta) / plot.ratio)
 lines(tmp$x.circ,tmp$y.circ,...)
 NULL
}
##############################################################################
##############################################################################
##############################################################################
anova.GLC.modeldef <- function(object,...,include.constraints=FALSE,
                               nonlin.warning=TRUE) {
 model.list <- list(object,...)
 nmodels <- length(model.list)
 if (nmodels < 2) stop("Need at least two models to make comparisons")
 LogLiks <- sapply(model.list,function(x) ifelse(is.null(x$LogLik),NA,x$LogLik))
 if (any(is.na(LogLiks))) {
  stop(paste("Log-likelihoods haven't been calculated for all models -",
             "fit them first\n  using GLCfit() before comparing them"))
 }
 classes <- sapply(model.list,function(x) x$model.type)
 if (!all(classes == classes[1])) {
  stop("All models must have the same value of model.type")
 }
 nobs <- sapply(model.list,function(x) x$nobs)
 if (!all(nobs == nobs[1])) {
  stop("Models were fitted to different numbers of observations")
 }
 modelnames <- sapply(model.list,function(x) x$model.title)
 res.df <- sapply(model.list,function(x) x$df.resid)
#
#       Work through the models from largest to smallest, and check
#       that they are consecutively nested. For each successive pair, 
#       we define a matrix indicating which parameters in the larger 
#	model have been omitted or fixed in the smaller one. We then 
#	create a corresponding list in which everything is related 
#	back to the first model.
#
 p <- sapply(model.list,function(x) x$df)
 if (any(duplicated(p))) stop("Models are not all nested")
 model.order <- order(p,decreasing=TRUE)
 
 Deltalist <- vector(length=nmodels-1,mode="list")
 Deltalist.full <- Deltalist
 theta.twid <- Deltalist
 p1 <- p[model.order[1]]
 theta.hat <- model.list[[model.order[1]]]$beta[1:p1] 
 if (classes[1] == "normal-heteroscedastic") {
  p1.disp <- model.list[[model.order[1]]]$dispersion$df
  p1.mean <- p1 - p1.disp
  theta.hat[(p1.mean+1):p1] <- 
        model.list[[model.order[1]]]$dispersion$beta[1:p1.disp] 
 }
 for (i in 2:nmodels) {
  Deltalist[[i-1]] <- get.nesting(model.list[[model.order[i-1]]],
                                  model.list[[model.order[i]]],
                                  nonlin.warning)
  Deltalist.full[[i-1]] <- get.nesting(model.list[[model.order[1]]],
                                  model.list[[model.order[i]]],
                                  nonlin.warning=FALSE)
  pcur <- p[model.order[i]]
  theta.cur <- model.list[[model.order[i]]]$beta[1:pcur] 
  if (classes[1] == "normal-heteroscedastic") {
   pcur.disp <- model.list[[model.order[i]]]$dispersion$df
   pcur.mean <- pcur - pcur.disp
   theta.cur[(pcur.mean+1):pcur] <- 
        model.list[[model.order[i]]]$dispersion$beta[1:pcur.disp] 
  }
  tmp <- rep(NA,p1)
  cols.fixed <- colSums(Deltalist.full[[i-1]]$Delta) > 0
  tmp[cols.fixed] <- Deltalist.full[[i-1]]$del0
#
#	The next lines take the estimates from the smaller model, and 
#	put them in the correct positions of the theta.twid vector
#	(a vector corresponding to the constrained estimates for the
#	larger model, with all coefficients in the same order). 
#
#	Note that this is the only place where this permutation is
#	needed, because all of the other calculations below are done
#	with quantities computed from the largest model.
#
  tmp[Deltalist.full[[i-1]]$TermIndices] <- theta.cur
  theta.twid[[i-1]] <- tmp
 }

#
#	Now create a data frame with test statistics etc.
# 
 z <- as.data.frame(matrix(nrow=nmodels,ncol=7))
 names(z) <- c("Resid DF","DF2-DF1","LogL","LLR",
               "p","Robust LLR","Robust p")
 rownames(z) <- c("M1",paste("M",1:(nmodels-1)," vs M",2:nmodels,sep=""))
 z[,1] <- res.df[model.order]
 z[,2] <- c(NA,diff(z[,1]))
 z[,3] <- LogLiks[model.order]
 z[,4] <- c(NA,-diff(z[,3]))
 z[,5] <- pchisq(2*z[,4],df=z[,2],lower.tail=FALSE)
#
#	Here are the adjusted likelihood ratios. Everything is related
#	back to the largest model in the first instance, then take 
#	differences to get the successive nesting
#
 z[1,6] <- 0
 if (classes[1] != "normal-heteroscedastic") {
  Hadj.inv <- model.list[[model.order[1]]]$cov.robust[1:p1,1:p1]
  Hind <- solve(model.list[[model.order[1]]]$cov.naive[1:p1,1:p1])
 } else {
  Hadj.inv <- matrix(0,nrow=p1,ncol=p1); Hind <- Hadj.inv
  Hadj.inv[1:p1.mean,1:p1.mean] <- 
       model.list[[model.order[1]]]$cov.robust[1:p1.mean,1:p1.mean]
  Hadj.inv[(p1.mean+1):p1,(p1.mean+1):p1] <- 
       model.list[[model.order[1]]]$dispersion$cov.robust[1:p1.disp,1:p1.disp]
  Hind[1:p1.mean,1:p1.mean] <- 
       model.list[[model.order[1]]]$cov.naive[1:p1.mean,1:p1.mean]
  Hind[(p1.mean+1):p1,(p1.mean+1):p1] <- 
       model.list[[model.order[1]]]$dispersion$cov.naive[1:p1.disp,1:p1.disp]
  Hind <- solve(Hind)
 }
 Hadj <- solve(Hadj.inv)
 for (i in 2:nmodels) {
  tmp1 <- (Deltalist.full[[i-1]]$Delta %*% theta.hat) - 
                         (Deltalist.full[[i-1]]$del0)
  tmp2 <- Deltalist.full[[i-1]]$Delta %*% Hadj.inv %*%
                        t(Deltalist.full[[i-1]]$Delta)
  numer <- t(tmp1) %*% solve(tmp2) %*% tmp1
  tmp1 <- theta.hat - theta.twid[[i-1]]
  denom <- t(tmp1) %*% Hind %*% tmp1
  z[i,6] <- (z[1,3]-z[i,3]) * numer / denom
 }
 z[,6] <- c(NA,diff(z[,6]))
 z[,7] <- pchisq(2*z[,6],df=z[,2],lower.tail=FALSE)
 attr(z,which="titles") <- modelnames[model.order]
 if (include.constraints) {
  constraints <- vector(nmodels-1,mode="list")
  names(constraints) <- paste("Constraints imposed in moving from model",
                              (1:(nmodels-1)),"to model",(2:nmodels))
  for (i in 1:(nmodels-1)) {
   tmp <- data.frame(matrix(nrow=1,ncol=length(Deltalist[[i]]$del0)))
   names(tmp) <- paste("Parameter",
        (1:ncol(Deltalist[[i]]$Delta))[colSums(Deltalist[[i]]$Delta)==1])
   rownames(tmp) <- c("Value")
   tmp[1,] <- Deltalist[[i]]$del0
   constraints[[i]] <- tmp
  }
  attr(z,which="constraints") <- constraints
 }
 class(z) <- c("anova.GLC.modeldef","data.frame")
 z
}
##############################################################################
##############################################################################
##############################################################################
get.nesting <- function(model1,model2,nonlin.warning=TRUE) {
#
#       Given two GLC.modeldef objects model1 and model2 (the first
#       being the more complex), figures out what linear restrictions on
#       the parameters have been imposed in model1 to get to model2, as 
#	well as the parameter numbers in model1 corresponding to the 
#       parameters in model2. The linear restrictions are of the form 
#
#       Delta theta = del0
#
#       where Delta is a matrix and del0 is a vector
#
 p <- model1$df
 q <- model1$df-model2$df
 if (q < 0) stop("model1 has fewer parameters than model2")
#
#       Set up Delta and del0; then we'll delete rows corresponding
#       to terms that *are* in model2
#
 Delta <- diag(rep(1,p))
 del0 <- rep(0,p)
#
#       Go through the specification for model2, for each parameter 
#       identifying the corresponding parameter in model1. Start with
#       main effects. NB the constant term in every model is included
#       in all the df calculations, but not in covcode, Np etc. NB also
#       the meaning of covcode is different for different model 
#       components - this means we need to take account of both the 
#       component and the value of covcode in each case. Site effects
#       are particularly fiddly because the nonlinear transformation
#       info is stored in sitxfm rather than in covcode. :-( 
#
 terms.idx <- rep(NA,model2$df-1)
 compnt.type1 <- pmax(0, # pmax because Np[7] is 0 if no nonlinear pars
                      c(model1$Np[1],diff(model1$Np[1:7]),model1$Np[8:10]))
 compnt.type1 <- rep(1:10,compnt.type1)
 unique.code1 <- 100*model1$covcode[1:length(compnt.type1)] + compnt.type1
 tmp <- which(compnt.type1 == 1)  # Numeric vector of indices
 unique.code1[tmp] <- 10000*model1$sitxfm[tmp,2] + unique.code1[tmp]
 compnt.type2 <- pmax(0,
                      c(model2$Np[1],diff(model2$Np[1:7]),model2$Np[8:10]))
 compnt.type2 <- rep(1:10,compnt.type2)
 unique.code2 <- 100*model2$covcode[1:length(compnt.type2)] + compnt.type2
 tmp <- which(compnt.type2 == 1)
 unique.code2[tmp] <- 10000*model2$sitxfm[tmp,2] + unique.code2[tmp]
#
#	Next lines were replaced in version 1.3-7: they don't work when
#	the user changes the order of main effects between models.
#
# tmp <- (unique.code1[1:model1$Np[4]] %in% unique.code2[1:model2$Np[4]])
# if (sum(tmp) < model2$Np[4]) stop("smaller model is not nested in larger one")
# maineff.idx <- (if (model2$Np[4]>=1) 1:model2$Np[4] else numeric(0))
# terms.idx[maineff.idx] <- (1:model1$Np[4])[tmp]
tmp <- numeric(0)
 if (model2$Np[4]>=1) {
  tmp <- match(unique.code2[1:model2$Np[4]], unique.code1[1:model1$Np[4]])
  if (any(is.na(tmp))) {
   stop("Models are not nested: main effects in smaller model do not all appear in larger one")
  }
 }
 maineff.idx <- (if (model2$Np[4]>=1) 1:model2$Np[4] else numeric(0))
 terms.idx[maineff.idx] <- tmp
#
#       Now 2-way interactions. Care needed here because the numbering of 
#       main effects might be different in the two models. Next two lines
#       convert the indices of interacting main effects in model2 into 
#       the corresponding indices for model1
#
 if (model2$Np[5] > model2$Np[4]) {
  tmp <- model2$twoway[(model2$Np[4]+1):model2$Np[5],,drop=FALSE]
  tmp <- matrix((terms.idx[maineff.idx])[tmp],nrow=nrow(tmp))
  for (i in 1:nrow(tmp)) {
   tmp2 <- model1$twoway[,1] %in% tmp[i,] &
           model1$twoway[,2] %in% tmp[i,]
   if (!any(tmp2)) {
     stop(paste("Models are not nested: larger model contains no interaction",
                "\n  corresponding to that between covariates",
                paste(model2$twoway[model2$Np[4]+i,], collapse=" and "), 
                "in smaller model."))
   }
   terms.idx[model2$Np[4]+i] <- (1:length(tmp2))[tmp2]
  }
 }
#
#       3-way interactions, similarly.
#
 if (model2$Np[6] > model2$Np[5]) {
  tmp <- model2$threeway[(model2$Np[5]+1):model2$Np[6],,drop=FALSE]
  tmp <- matrix((terms.idx[maineff.idx])[tmp],nrow=nrow(tmp))
  for (i in 1:nrow(tmp)) {
   tmp2 <- model1$threeway[,1] %in% tmp[i,] &
           model1$threeway[,2] %in% tmp[i,] &
           model1$threeway[,3] %in% tmp[i,]
   if (!any(tmp2)) {
     stop(paste("Models are not nested: larger model contains no interaction",
                "\n  corresponding to that between covariates (",
                paste(model2$threeway[model2$Np[5]+i,], collapse=", "), 
                ") in smaller model.", sep=""))
   }
   terms.idx[model2$Np[5]+i] <- (1:length(tmp2))[tmp2]
  }
 }
#
#       And parameters in nonlinear transformations. Code here is a bit
#       messy due to the way that information on nonlinear transformations
#       is dispersed around the model specification. Need to know the 
#       following:
#
#       - Only unknown parameters are accounted for in (e.g.) covcode and Np
#       - Corresponding elements of covcode are (1000*a)+b, where a
#         is the covariate to which the nonlinear parameter belongs
#         and b is the number of the parameter being defined
#       - nlstat(a,b) gives the status of this nonlinear parameter
#
#       First, identify parameters that are estimated in both models ...
#
 if (model2$Np[7] > model2$Np[6]) {
  nlpars2 <- (model2$Np[6]+1):model2$Np[7]
  nlcovs2 <- model2$covcode[nlpars2] %/% 1000
  which.nl2 <- model2$covcode[nlpars2] %% 1000
  tmp <- (1000*terms.idx[nlcovs2]) + which.nl2
  tmp2 <- (tmp %in% model1$covcode)
  if (!all(tmp2)) stop(paste("smaller model is not nested in larger ",
                             "one (problem is with\nspecification of",
                             "nonlinear transformations"))
  nlpars1 <- (model1$Np[6]+1):model1$Np[7]
  terms.idx[nlpars2][tmp2] <- nlpars1[tmp2]
 }  
#
#       Now figure out what are the values of del0 corresponding to 
#       any parameters that were estimated in model1 but fixed in 
#       model2. This is tricky because they're not recorded in covcode 
#       for model2. It's done by making tables corresponding to nlstat, 
#       showing where the fixed and estimated parameters are for both
#       models
#
 if (model1$Np[7] > model1$Np[6]) {
  main1 <- 1:model1$Np[4]
  nltab1 <- (model1$nlstat[main1,] >= 0) &        # All nonlinear pars in model1
            (abs(model1$theta[main1,]) < 1e9) & 
            !duplicated(model1$weighting.scheme[main1],incomparables=NA) &
            !duplicated(model1$legendre.idx[main1],incomparables=0) &
            !duplicated(model1$fourier.idx[main1],incomparables=0) 
  fixed1 <- nltab1 & (model1$nlstat[main1,] == 0) # And those that are fixed
  nltab2 <- (model2$nlstat[maineff.idx,] >= 0) &  # Ditto, model2
            (abs(model2$theta[maineff.idx,]) < 1e9) &   
            !duplicated(model2$weighting.scheme[maineff.idx],
                                             incomparables=NA) &
            !duplicated(model2$legendre.idx[maineff.idx],
                                             incomparables=0) &
            !duplicated(model2$fourier.idx[maineff.idx],
                                             incomparables=0) 
  fixed2 <- nltab2 & (model2$nlstat[maineff.idx,drop=FALSE] == 0)
#
#       Next line gives a table, the same size as fixed2, that is TRUE
#       whenever one of those parameters was estimated in model1
#
  which.refixed <- nltab1[terms.idx[maineff.idx],,drop=FALSE] & 
                   !fixed1[terms.idx[maineff.idx],,drop=FALSE] & fixed2
#
#       And the *next *line is the value of covcode that we should find
#       somewhere in model1 for those parameters
#
  covcode.refixed <- 1000*terms.idx[row(which.refixed)[which.refixed]] +
                     col(which.refixed)[which.refixed] 
#
#       Now go through each parameter, find the corresponding entry in 
#       model1 and set the value of del0 accordingly (NB at this point
#       the first element of del0 corresponds to the constant term in the
#       model, so need to shift everything up a position - hence insert
#       entry in appropriate position of del0[2:pp]). pp here is equal
#       to p, except for normal-heteroscedastic models where it is the 
#       number of parameters in the *mean* part of the model (for the
#       dispersion part, the value of model.type is stored as "gamma").
#
  pp <- p
  if (model1$model.type=="normal-heteroscedastic") pp <- p-model1$dispersion$df
  if (length(covcode.refixed) > 0) {
   for (i in 1:length(covcode.refixed)) {
    tmp <- (model1$covcode[1:(pp-1)] == covcode.refixed[i])
    tmp[is.na(tmp)] <- FALSE    # covcode contains NAs
    if (sum(tmp) == 0) {
     stop(paste("smaller model has a fixed nonlinear parameter that",
                "doesn't appear in\n  larger model"))
    }
    del0[(2:pp)][tmp] <- model2$theta[which.refixed][i]
   }
  }
#
#       Finally in this section, deal with any nonlinear parameters 
#	defined in model1 that don't appear anywhere in model2. Since we
#	need to insert a value in del0 for such parameters, default 
#	behaviour is to keep the value from model1 and to issue a 
#	warning. The warning can be suppressed by setting nonlin.warning to
#	FALSE.
#
  nl1.unknown <- nltab1
  nl1.unknown[terms.idx[maineff.idx],][nltab2] <- FALSE
  if (any(nl1.unknown)) {
   if (nonlin.warning) {
    warning(paste("In anova.GLC.modeldef, larger model has a nonlinear",
                  "parameter\n  attached to a covariate that doesn't appear",
                  "in smaller model. For\n  calculation purposes, I'm taking",
                  "this parameter in the smaller model\n  to be fixed at its",
                  "value from the larger. This is probably unimportant.\n\n",
                  " Set nonlin.warning=FALSE to suppress this warning."),
                  call.=FALSE,immediate.=TRUE)
   }
   covcode.dropped <- 1000*row(nl1.unknown)[nl1.unknown] +
                           col(nl1.unknown)[nl1.unknown]  
   for (i in 1:length(covcode.dropped)) {
    tmp <- (model1$covcode[1:(pp-1)] == covcode.dropped[i])
    tmp[is.na(tmp)] <- FALSE    # covcode contains NAs
#
#	Easiest to get the value from beta, which has the constant at
#	the beginning - hence remove first element before selecting elements
#
    del0[(2:pp)][tmp] <- model1$beta[(2:pp)][tmp]
   }
  }
 } else if (model2$Np[7] > model2$Np[6]) {
  stop("smaller model is not nested in larger one")
 }
#
#	Dispersion parameters, for models that have them
#
 if (model1$Np[9] > 0 & model1$model.type != "normal-heteroscedastic") {
  terms.idx[max(model2$Np[6:7])+1] <- max(model1$Np[6:7])+1
 }
#
#       At this stage, we're sorted for all model types except
#       "normal-heteroscedastic" - for this one, we need to go
#       through the dispersion component as well. This is slightly
#       clumsy - get the Delta and del0 components from that
#       nesting, then the columns of Delta with a "1" in them are the
#       parameters that were discarded in the simpler model
#
 if (model1$model.type == "normal-heteroscedastic") {
  disp.nesting <- get.nesting(model1$dispersion,model2$dispersion,
                                                      nonlin.warning)
  tmp <- max(model1$Np[6:7]) + 
         (1:model1$dispersion$df)[colSums(disp.nesting$Delta) == 0]
  terms.idx[(max(model2$Np[6:7])+1):(model2$df-1)] <- tmp
  tmp <- max(model1$Np[6:7]) + 1 +  
         (1:model1$dispersion$df)[colSums(disp.nesting$Delta) == 1]
  del0[tmp] <- disp.nesting$del0
 }
#
#       Now remove the appropriate rows from Delta and del0. NB adding
#       1 to everything (and removing the first row in all cases) is
#       necessary to account for the constant term in parameter vectors 
#       and covariance matrices
#
 terms.idx <- c(1,terms.idx+1)
 Delta <- Delta[-terms.idx,,drop=FALSE]
 del0 <- del0[-terms.idx]
#
#	Return Delta and del0, as well as the indices of the parameters in
#	model1 corresponding to each term in model2
#
 list(Delta=Delta,del0=del0,TermIndices=terms.idx)
}
##############################################################################
##############################################################################
##############################################################################
print.anova.GLC.modeldef <- function(x,...) {
#
#	To print an anova.GLC.modeldef object (the result of a call to
#	anova.GLC.modeldef).
#
 cat("Comparison of nested models\n---------------------------\n")
 nmodels <- nrow(x)
 cat(paste("\nModel ",1:nmodels,": ",attr(x,which="titles"),sep=""))
 cat("\n\n")
 tmp <- as.data.frame(x)
 tmp[,c(5,7)] <- format.pval(tmp[,c(5,7)])
 tmp[,c(4,6)] <- round(tmp[,c(4,6)],3)
 tmp[is.na(tmp) | tmp == "NA"] <- ""
 print(tmp)
 constraints <- attr(x,which="constraints") 
 if(!is.null(constraints)) {
  cat("\nConstraints:\n")
  for (i in 1:(nmodels-1)) {
   tmp <- constraints[[i]]
   rownames(tmp) <- paste("M",i," vs M",i+1,sep="")
   tmp[1,] <- formatC(as.numeric(tmp[1,]),digits=5,width=12,format="f")
   print(tmp)
  }
 }
}
##############################################################################
##############################################################################
##############################################################################
GLCsim <- function(modeldefs,siteinfo,start,end,nsims,impute.until=end,
                   output=c("daily","monthly"),which.regions=0,
                   which.daily=1:nsims,daily.start=start,daily.end=end,
                   data.file,external.files,simdir,file.prefix,missval=-99.99) {
##############################################################################
#
#       To simulate from a fitted model or set of models. Arguments:
#
#       modeldefs       - Either an object defining a model for a single
#                         variable (for univariate simulation), or a list
#                         of such objects, each representing a different 
#                         variable (multivariate simulation). For most
#                         variables, these individual objects will be 
#                         of class GLC.modeldef; the exception is precipitation
#                         for which separate GLC.modeldef objects are
#                         required for occurrence and intensity. A 
#                         precipitation model must itself be specified as a 
#                         list with named components "Occurrence" and
#                         "Intensity", which are GLC.modeldef objects 
#                         for "logistic" and "gamma" GLMs respectively. 
#       siteinfo        - A siteinfo object containing information about
#                         the sites to be used in the simulation. This
#                         should give all of the attribute values for each
#                         site that are required by the models in modeldefs,
#                         and should usually be generated by a call to 
#                         make.siteinfo(). 
#       start           - Start date for simulation, in format YYYYMM. 
#                         The first day simulated is the first of the month.
#       end             - End date for simulation, similarly. The last
#                         day simulated is the last of the month.
#       nsims           - Number of simulations
#       impute.until    - A date, in the form YYYYMM. If present then 
#                         simulations up to and including the last day
#                         of this month will be conditioned on all 
#                         available observations
#       output          - Chooses whether to produce daily output files,
#                         monthly output files or both
#       which.regions   - A vector containing the numbers of regions for 
#                         which monthly summaries should be produced (if
#                         monthly outputs have been requested). NB 0 is the
#                         entire area.
#       which.daily     - Vector of simulation numbers for which to produce
#                         daily output (if daily output files have been 
#                         requested). Defaults to all simulations
#       daily.start     - Start date for daily output, in form YYYYMM.
#                         Default is start of simulation
#       daily.end       - End date for daily output. Default is end of 
#                         simulation
#       data.file       - File from which to take input data. If not supplied,
#                         the routine will take the filename from the model 
#                         objects in modeldefs.
#       external.files  - Character vector of length 3, giving names of 
#                         files from which to take "external" predictor
#                         data (yearly, monthly and daily) to drive the 
#                         simulations. If not supplied, the routine will
#                         take the file names from the model objects in 
#                         modeldefs.
#       simdir          - Name of directory in which to store the output
#                         files. This will be interpreted as a pathname 
#                         relative to the current working directory. 
#                         The routine removes any trailing directory
#                         separators so, for example, simdir="TestSim",
#                         simdir="TestSim/" and simdir="TestSim\" all have 
#                         the same effect
#       file.prefix     - Output files are named in a structured way as,
#                         for example, AshdownSim_Daily_Sim0001.dat. Here,
#                         "AshdownSim" is a user-specified character string 
#                         defined via the file.prefix argument, and
#                         the remainder of the filename is generated
#                         automatically. By default, file.prefix is the
#                         same as simdir with any leading dots and slashes
#                         removed.
#       missval         - The missing value code in the input file. 
#
##############################################################################
#
#       Store the state of the random number generator
#
 if (exists(".Random.seed")) {
  RNGstate <- list(RNGkind=RNGkind(),seed=.Random.seed)
 } else {
  RNGstate <- list(RNGkind=RNGkind(),seed=NULL)
 }
#
#       Close any open Fortran file handles in case something went wrong
#       previously
#
 CloseFiles(10:1000)
#
#       Figure out whether we're doing univariate or multivariate
#       simulation, and the number(s) (if any) of the model(s) corresponding
#       to precipitation. NB allow more than one - there may be other 
#       "precipitation-like" variables, I guess. 
#
 multivariate <- TRUE
 precip.model <- -1
 if (inherits(modeldefs,"GLC.modeldef")) {
  multivariate <- FALSE
 } else if (length(modeldefs)==2 & !is.null(names(modeldefs)) &
            all(names(modeldefs) %in% c("Occurrence","Intensity"))) {
  multivariate <- FALSE
  precip.model <- 1
 } else {
  GLCmodels <- sapply(modeldefs,FUN=function(x) {inherits(x,"GLC.modeldef")})
  precip.model <- (1:length(modeldefs))[!GLCmodels]
  for (i in precip.model) {
   if (!all(names(modeldefs[[i]]) %in% c("Occurrence","Intensity"))) {
    stop(paste("Element",i,"of modeldefs has incorrect structure or is",
               "wrongly named"))
   }
  }
 }
#
#       Check that all required input files are present. Precipitation
#       needs to be handled separately. For univariate (non-precip)
#       simulation, need to wrap the GLC.modeldef list into an "outer" 
#       list for the loop to work - this is what the first two lines are
#       for below. For precip models, also check that the values of any
#       "global" quantities are the same for both occurrence and amounts. 
#
 cat("\nChecking files ...\n")
 mlist <- modeldefs
 if (!multivariate) mlist <- list(mlist) 
 for (model in mlist) {
  if (inherits(model,"GLC.modeldef")) {
   datflnm <- model$filenames$Data
   if (!missing(data.file)) datflnm <- data.file
   exflnms <- if(missing(external.files)) {
               model$filenames$External } else external.files
   tmp <- CheckFiles(data.file=datflnm,external.files=exflnms,
                     cor.file=model$filenames$Correlation,
                     model.def=model,diagnostics=0,fitting=FALSE,msg=FALSE)
   CloseFiles(tmp) 
  } else {
   datflnm <- model$Occurrence$filenames$Data
   if (!missing(data.file)) datflnm <- data.file
   exflnms <- if(missing(external.files)) {
               model$Occurrence$filenames$External } else external.files
   tmp <- CheckFiles(data.file=datflnm,external.files=exflnms,
                     cor.file=model$Occurrence$filenames$Correlation,
                     model.def=model$Occurrence,diagnostics=0,
                     fitting=FALSE,msg=FALSE)
   CloseFiles(tmp) 
   datflnm <- model$Intensity$filenames$Data
   if (!missing(data.file)) datflnm <- data.file
   exflnms <- if(missing(external.files)) {
               model$Intensity$filenames$External } else external.files
   tmp <- CheckFiles(data.file=datflnm,external.files=exflnms,
                     cor.file=model$Intensity$filenames$Correlation,
                     model.def=model$Intensity,diagnostics=0,
                     fitting=FALSE,msg=FALSE)
   CloseFiles(tmp)
#
#       Also check that the same global values have been defined for
#       both occurrence and amounts in a precipitation model (NB
#       the use of model$Occurrence in both subsetting commands is
#       deliberate and *correct*).
#
   if (any(model$Occurrence$global.codes[1:model$Occurrence$Np[8]] !=
           model$Intensity$global.codes[1:model$Occurrence$Np[8]]) ) {
    stop(paste("Precipitation occurrence and intensity models have",
               "different definitions of global quantities"))
   }
   tmp <- floor(model$Occurrence$global.codes[1:model$Occurrence$Np[8]] / 1000)
   tmp <- model$Occurrence$global.vals[tmp] - model$Intensity$global.vals[tmp]
   if (any(abs(tmp) > 1e-6)) {
    stop(paste("Precipitation occurrence and intensity models have",
               "different definitions of global quantities"))
   }
#
#       And that both models reference the same variable
#
   if (model$Occurrence$which.response != model$Intensity$which.response) {
    stop(paste("Precipitation occurrence and intensity models reference",
               "different\n  response variables"))
   }
  }
 }
#
#       Worth checking that the data file is the same for all models 
#       (all other files might conceivably be different because models
#       might use different external predictors, correlation structures
#       etc.)
#
 if (missing(data.file)) {
  tmp <- unlist(mlist,recursive=TRUE)
  tmp <- tmp[grep("filenames.Data",names(tmp))]
  if (length(unique(tmp)) > 1) {
   stop("models in modeldefs were fitted using different data files")
  }
  datflnm <- unique(tmp)
 } else {
  datflnm <- data.file
 }
#
#       Figure out how many variables are in the data file, and what their
#       names are (take these just from the first model in the list)
#
 nvars <- (nchar(scan(datflnm,what="",sep="\n",n=1,quiet=TRUE))-12)/6
 if (abs(nvars-round(nvars)) > 1e-6) stop(err.msg(datflnm,err.no=5))
#
#       Check that any required covariates are being simulated. Also,
#       If we're doing multivariate simulation, check that modeldefs 
#       is self-consistent by going through each definition in turn to 
#       check that there are no circular dependencies (also check 
#       that the user hasn't supplied two models for the same variable -
#       there are some things you just wouldn't believe possible unless
#       you had students available to demonstrate your lack of 
#       imagination). And extract the observed mean of each variable, 
#       for use in initialising simulations where there are missing data 
#       values. Also variable names, from the first model in the list. 
#       And model types (to be stored in the final object).
#
 dependencies <- diag(0,nrow=nvars)
 mean.table <- rep(0,nvars)
 which.vars <- numeric(0)
 var.names <- NULL
 model.types <- character(nvars)
 for (model in mlist) {
  if (inherits(model,"GLC.modeldef")) {
   which.vars <- c(which.vars,model$which.response)
   covs.needed <- floor(model$covcode[(model$Np[3]+1):model$Np[4]] / 1e6)
   covs.needed <- covs.needed[covs.needed > 0]
   dependencies[model$which.response,covs.needed] <- 1
   mean.table[model$which.response] <- model$Residuals$ybar
   if (is.null(var.names)) {
     var.names <- model$var.names
   } else if (!identical(var.names,model$var.names)) {
     stop("models have different var.names components - check original definitions")
   }
   model.types[model$which.response] <- model$model.type
  } else {
   tmp <- model$Occurrence
   which.vars <- c(which.vars,tmp$which.response)
   covs.needed <- floor(tmp$covcode[(tmp$Np[3]+1):tmp$Np[4]] / 1e6)
   covs.needed <- covs.needed[covs.needed > 0]
   dependencies[tmp$which.response,covs.needed] <- 1
   tmp <- model$Intensity
   covs.needed <- floor(tmp$covcode[(tmp$Np[3]+1):tmp$Np[4]] / 1e6)
   covs.needed <- covs.needed[covs.needed > 0]
   dependencies[tmp$which.response,covs.needed] <- 1
#
#       For precipitation-like variables, set the default initial
#       value to zero if this is the median of the distribution, and
#       to the overall mean otherwise. NB if the model object 
#       doesn't 
#
   if (is.null(tmp$Residuals$ybar)) {
    warning(paste("Model object doesn't contain means of variables - simulations",
                  "will be initialised\n  with zero unless file",data.file,
                  "contains appropriate values\n"),immediate.=TRUE)
    mean.table[tmp$which.response] <- 0
   } else {
    mean.table[tmp$which.response] <- tmp$Residuals$ybar
    if (tmp$Residuals$ybar < 0.5) {
     mean.table[tmp$which.response] <- mean.table[tmp$which.response] *
                                       model$Intensity$Residuals$ybar
    }
   }
   if (is.null(var.names)) {
     var.names <- model$Occurrence$var.names
   } else if (!identical(var.names,model$Occurrence$var.names)) {
     stop("models have different var.names components - check original definitions")
   }
   model.types[model$Occurrence$which.response] <- "logistic-gamma"
  }
 }
 if(any(duplicated(which.vars))) {
  stop(paste("You have specified multiple models for variable",
             var.names[unique(which.vars[duplicated(which.vars)])[1]]))
 }
 vars.needed <- colSums(dependencies) > 0
 vars.needed <- (1:nvars)[vars.needed]
 vars.needed <- vars.needed[!(vars.needed %in% which.vars)]
 if (length(vars.needed) > 0) {
  stop(paste("This simulation requires models for the following additional ",
             "variables\n  in file ",datflnm,": variable number(s) ",
             paste(vars.needed,collapse=", "),sep=""))
 }
 diag(dependencies) <- 0 # Direct dependencies are allowed
 dependencies <- dep.calc(dependencies)
 if (any(diag(dependencies) > 0)) {
  tmp <- (1:ncol(dependencies))[diag(dependencies) > 0]
  if (length(tmp) == 2) {
   tmp <- paste(tmp,collapse=" and ") 
  } else {
   tmp <- paste(paste(tmp[-length(tmp)],collapse=", "),"and",tmp[length(tmp)])
   stop(paste("Circular dependencies",
              "(e.g. X depends on Y, Y depends on Z, Z depends on X)\n",
              " found between elements",tmp,"of modeldefs"))
  }
 }
#
#       Now check the siteinfo argument: if it's missing then set
#       it to the siteinfo element of the model objects (after 
#       checking that they're all the same)
#
 if (missing(siteinfo)) {
  if (inherits(mlist[[1]],"GLC.modeldef")) {
   simsites <- mlist[[1]]$siteinfo
  } else {
   simsites <- mlist[[1]]$Occurrence$siteinfo
  }
  for (model in mlist) {
   if (inherits(model,"GLC.modeldef")) {
    tmp <- identical(model$siteinfo,simsites)
   } else {
    tmp <- identical(model$Occurrence$siteinfo,simsites)
   }
   if (!tmp) stop("No siteinfo argument provided and component models disagree")
  }
 } else {
  simsites <- siteinfo
 }
#
#       Simulate variables in order starting with those that have the
#       fewest dependencies; and re-order the list of models correspondingly.
#       The logic here is that the which.vars vector is currently in the 
#       order obtained by reading the models in mlist one at a time;
#       applying the sam re-ordering to which.vars and to mlist should 
#       therefore do the trick. 
#
 sim.order <- order(rowSums(dependencies[which.vars,,drop=FALSE]))
 mlist <- mlist[sim.order]
 sim.order <- which.vars[sim.order]
#
#       Now check site and region information
#
 z <- check.siteinfo(simsites,gen.sitinf=FALSE)
 nsites <- z$nsites
 nattr <- z$nattr
 nregions <- max(length(simsites$Region.names),1)-1
 if (any(which.regions > nregions | which.regions < 0)) {
  stop("Output requested for undefined region - check value of which.regions")
 }
 if (any(simsites$Regions > nregions)) {
  stop(paste("One or more sites has a region number greater than the number ", 
             "of\n  named regions: check the site information for your simulation.",
             sep=""))
 }
 regs.wanted <- rep(0,nsites+1) # add an element at the start for the whole area
 regs.wanted[1:(nregions+1)] <- (0:nregions) %in% which.regions
 region.codes <- simsites$Regions
#
#       And that start and end dates seem valid
#
 if (missing(start) | missing(end)) stop("Start and end dates must be supplied")
 if (start < 0 | start > 999999) {
  stop("Invalid start date supplied (format should be YYYYMM)")
 }
 if (end < 0 | end > 999999) {
  stop("Invalid end date supplied (format should be YYYYMM)")
 }
 if (start >= end) stop("Start date must precede end date")
#
#       Check requested output options
#
 if (!all(output %in% c("daily","monthly"))) {
  stop("invalid value of 'output' argument")
 }
 if ("daily" %in% output) {
  if (daily.start < 0 | daily.start > 99999999) {
   stop("Invalid start date supplied for daily output")
  }
  if (daily.end < 0 | daily.end > 99999999) {
   stop("Invalid end date supplied for daily output")
  }
  if (daily.start >= daily.end) {
   stop("For daily output, start date must precede end date")
  }
 }
#
#       Create simdir if necessary and check that there are no files with 
#       names beginning file.prefix. If there *are*, and if the user
#       is happy to proceed, delete them all before proceeding (this is 
#       to avoid confusion if you do one run with 10 simulations and
#       another with 100 simulations - in this case the original files
#       won't be overwritten because the simulation numbers will appear
#       in the filenames as 01, 02, ..., 10 rather than 001, 002, ..., 100.
#
 tmp <- nchar(simdir)
 while (substr(simdir,tmp,tmp) %in% c("/","\\")) {
  simdir <- substr(simdir,1,tmp-1)
  tmp <- tmp - 1
 }
 if (!file.exists(simdir)) dir.create(simdir)
 if (missing(file.prefix)) {
  file.prefix <- simdir
  while (substr(file.prefix,1,1) %in% c(".","/","\\")) {
   file.prefix <- substr(file.prefix,2,nchar(file.prefix))
  }
 }
 tmp <- list.files(path=simdir)
 tmp2 <- grep(paste("^",file.prefix,sep=""),tmp)
 if (length(tmp2) > 0) {
  cat(paste('\n  ****NOTE**** Directory ',simdir,
            ' already contains files with\n  names beginning "',
            file.prefix,'". If you proceed, these will be deleted.\n',
            '  OK to continue (Y/N, default N)? ',sep=""))
  if (!(readLines(n=1) %in% c("y","Y"))) return(NULL)
  file.remove(paste(simdir,"/",tmp[tmp2],sep=""))
  }
#
#       Initialise Fortran random number seed (this is done using an 
#       R-generated random number, so the Fortran random number generator 
#       can be set to a repeatable state via a call to set.seed() in R
#       before calling this routine.
# 
 Fortran.seed <- as.integer(floor(runif(1)*(2^31-1)))
 Fortran.seed <- .Fortran("zbqlini",seed=Fortran.seed)$seed
#
#       Here are the simulations: one set of output files per simulation.
#       Before starting, read any values to be conditioned upon in 
#       imputation and write these to a direct-access scratch file; 
#       and at the start of each simulation, open another scratch
#       file to hold the complete daily data for that simulation (NB 
#       we have already checked that all of the component models access 
#       the same data file). Next two lines find the maximum lag required 
#       by any of the models
#
 tmp <- unlist(mlist,recursive=TRUE)
 max.lag <- as.integer(max(as.numeric(tmp[grep("max.lag",names(tmp))])))
 OpenUnits <- CheckFiles(data.file=datflnm,diagnostics=0,
                                             fitting=FALSE,msg=FALSE)
 OpenUnits[95] <- OpenFile(scratch=TRUE,Nfields=nvars+4)
 if (OpenUnits[95] < 0) CloseFiles(OpenUnits[95],error=TRUE)
 storage.mode(OpenUnits) <- "integer"
 DateLims <- c(start,end); storage.mode(DateLims) <- "integer"
#
#  Next lines: write site codes to a temporary file, as the only way to 
#  pass character information from R to Fortran (Rv3.6.1 onwards).
#
 TempFile <- OpenTempFile()
 write(simsites$Site.codes, file=TempFile$Name, append=TRUE)
 tmp <- .Fortran("ReadData2",UnitNos=OpenUnits,MissVal=missval,
                  NSITES=as.integer(nsites),NVARS=as.integer(nvars),
                  TmpFlNo=TempFile$Number,DateRange=DateLims,
                  ImputeTill=as.integer(impute.until),MaxLag=max.lag,
                  IFAIL=as.integer(0))
 if (tmp$IFAIL != 0) {
  stop("Error while reading data - see message above for details")
 }
#
#       Output a list of variable names
# 
 if (multivariate) {
  cat(paste("\nSimulating the following variables from file ",datflnm,":\n\n",sep=""))
  for (i in sort(which.vars)) {
   cat(paste("\t",i,".\t",var.names[i],"\n",sep="")) 
  }
 }
 cat("\n")
 for (sim in 1:nsims) {
  cat(paste("Simulation ",sim," of ",nsims,": opening files ...\r",sep=""))
  in.unit <- as.integer(95) # Read original data when simulating first variable
  OpenUnits[96] <- OpenFile(scratch=TRUE,Nfields=nvars+4)
  if (OpenUnits[96] < 0) CloseFiles(c(OpenUnits),error=TRUE)
  for (i in 1:length(mlist)) {
   if (multivariate) {
    out.text <- paste("Simulation ",sim," of ",nsims,", variable ",sim.order[i],
              " (simulation order: ",paste(sim.order,collapse="-"),
              ") ...",sep="")
   } else {
    out.text <- paste("Simulation ",sim," of ",nsims,": simulating ...",sep="")
   }
   cat(paste(out.text,"\r",sep=""))
   model <- mlist[[i]]
#
#       Set output option for this simulation
#
   OutOpt <- (2*as.numeric("daily" %in% output & sim %in% which.daily)) + 
                                         as.numeric("monthly" %in% output)
#
#       Open the required data files etc.for this model, and sort out the 
#       arguments for the Fortran call. This is very tedious ... 
#
   if (inherits(model,"GLC.modeldef")) {
    exflnms <- if(missing(external.files)) {
                model$filenames$External } else external.files
    UnitNos <- CheckFiles(data.file=datflnm,external.files=exflnms,
                      cor.file=model$filenames$Correlation,
                      model.def=model,diagnostics=0,fitting=FALSE,msg=FALSE)
    CorFileName <- model$filenames$Correlation
    Np1 <- as.integer(model$Np)
    Mxp1 <- as.integer(model$max.pars)
    sitinf1 <- check.siteinfo(simsites,Mxp1)$sitinf
    Sitxfm1 <- model$sitxfm
    Covcode1 <- model$covcode; Covcode1[is.na(Covcode1)] <- -1
    Two1 <- model$twoway
    Three1 <- model$threeway
    Fouid1 <- model$fourier.idx
    Legid1 <- model$legendre.idx
    Pwtidx1 <- expand.pwtidx(model,nvars)
    Spmod1 <- model$spmod 
    Sitxfm1 <- model$sitxfm
    Beta1 <- model$beta; Beta1[is.na(Beta1)] <- 1e9
    Theta1 <- model$theta
    Rho1 <- model$rho
    global.codes <- model$global.codes
    global.vals <- model$global.vals
    if (model$model.type == "normal-heteroscedastic") {
     model.arg <- c(model.num(model$model.type,"mean"),
                    model.num(model$model.type,"dispersion"))
     Phi1 <- Phi2 <- 1
     Np2 <- as.integer(model$dispersion$Np)
     Mxp2 <- as.integer(model$dispersion$max.pars)
     sitinf2 <- check.siteinfo(simsites,Mxp2)$sitinf
     Covcode2 <- model$dispersion$covcode; Covcode2[is.na(Covcode2)] <- -1
     Two2 <- model$dispersion$twoway
     Three2 <- model$dispersion$threeway
     Fouid2 <- model$dispersion$fourier.idx
     Legid2 <- model$dispersion$legendre.idx
     Pwtidx2 <- expand.pwtidx(model$dispersion,nvars)
     Sitxfm2 <- model$dispersion$sitxfm
     Beta2 <- model$dispersion$beta; Beta2[is.na(Beta2)] <- 1e9
     Theta2 <- model$dispersion$theta
    } else {
     model.arg <- c(model.num(model$model.type),0)
     Phi1 <- model$dispersion; Phi2 <- 1
     Mxp2 <- as.integer(1)
     Np2 <- as.integer(rep(0,10))
     sitinf2 <- array(0,dim=c(nsites,Mxp2,4))
     Covcode2 <- -1
     Two2 <- Sitxfm2 <- matrix(rep(0,2),nrow=Mxp2)
     Three2 <- matrix(rep(0,3),nrow=Mxp2)
     Fouid2 <- Legid2 <- 0
     Pwtidx2 <- array(0,dim=c(Mxp2,4,nvars))
     Beta2 <- rep(1e9,Mxp2+1)
     Theta2 <- matrix(rep(1e9,3),nrow=Mxp2)
    }
    Spmod2 <- as.integer(0)
    Rho2 <- rep(1e9,Mxp2)
    is.precip <- as.integer(0)
   } else {
    exflnms <- if(missing(external.files)) {
                model$Occurrence$filenames$External } else external.files
    UnitNos <- CheckFiles(data.file=datflnm,external.files=exflnms,
                      cor.file=model$Occurrence$filenames$Correlation,
                      model.def=model$Occurrence,diagnostics=0,
                      fitting=FALSE,msg=FALSE)
    CorFileName <- model$Occurrence$filenames$Correlation
#
#       Any files that have already been connected will not be reconnected
#       in the following command; all that happens is that a -1 gets
#       returned in the corresponding position of the tmp vector
#
    exflnms <- if(missing(external.files)) {
                model$Intensity$filenames$External } else external.files
    tmp <- CheckFiles(data.file=datflnm,external.files=exflnms,
                      cor.file=model$Intensity$filenames$Correlation,
                      model.def=model$Intensity,diagnostics=0,
                      fitting=FALSE,which.model=2,msg=FALSE)   
    if (any(tmp != UnitNos & tmp > 0 & UnitNos > 0)) {
     stop("Input files for occurrence and intensity models don't match")
    }
    UnitNos <- pmax(UnitNos,tmp)
    if (!is.null(model$Intensity$filenames$Correlation)) {
     if (!is.null(CorFileName)) {
      CorFileName <- paste(CorFileName,"or",model$Intensity$filenames$Correlation)
     } else {
      CorFileName <- model$Intensity$filenames$Correlation
     }
    }
    model.arg <- c(model.num(model$Occurrence$model.type),
                   model.num(model$Intensity$model.type))
    is.precip <- as.integer(1)
    Np1 <- as.integer(model$Occurrence$Np)
    Mxp1 <- as.integer(model$Occurrence$max.pars)
    sitinf1 <- check.siteinfo(simsites,Mxp1)$sitinf
    Sitxfm1 <- model$Occurrence$sitxfm
    Covcode1 <- model$Occurrence$covcode; Covcode1[is.na(Covcode1)] <- -1
    Two1 <- model$Occurrence$twoway
    Three1 <- model$Occurrence$threeway
    Fouid1 <- model$Occurrence$fourier.idx
    Legid1 <- model$Occurrence$legendre.idx
    Pwtidx1 <- expand.pwtidx(model$Occurrence,nvars)
    Spmod1 <- model$Occurrence$spmod 
    Sitxfm1 <- model$Occurrence$sitxfm
    Beta1 <- model$Occurrence$beta; Beta1[is.na(Beta1)] <- 1e9
    Theta1 <- model$Occurrence$theta
    Rho1 <-model$Occurrence$rho
    Phi1 <- model$Occurrence$dispersion
    Np2 <- as.integer(model$Intensity$Np)
    Mxp2 <- as.integer(model$Intensity$max.pars)
    sitinf2 <- check.siteinfo(simsites,Mxp2)$sitinf
    Sitxfm2 <- model$Intensity$sitxfm
    Covcode2 <- model$Intensity$covcode; Covcode2[is.na(Covcode2)] <- -1
    Two2 <- model$Intensity$twoway
    Three2 <- model$Intensity$threeway
    Fouid2 <- model$Intensity$fourier.idx
    Legid2 <- model$Intensity$legendre.idx
    Pwtidx2 <- expand.pwtidx(model$Intensity,nvars)
    Spmod2 <- model$Intensity$spmod 
    Sitxfm2 <- model$Intensity$sitxfm
    Beta2 <- model$Intensity$beta; Beta2[is.na(Beta2)] <- 1e9
    Theta2 <- model$Intensity$theta
    Rho2 <-model$Intensity$rho
    Phi2 <- model$Intensity$dispersion
    global.codes <- model$Occurrence$global.codes
    global.vals <- model$Occurrence$global.vals
   }
   OpenUnits <- pmax(OpenUnits,UnitNos)
   Rho1 <- rho.check(Rho1,Spmod1)
   Rho2 <- rho.check(Rho2,Spmod2)

   storage.mode(model.arg) <- "integer"
   DateLims <- c(start,end); OutLims <- c(daily.start,daily.end)
   storage.mode(DateLims) <- storage.mode(OutLims) <- "integer"
   storage.mode(Covcode1) <- storage.mode(Covcode2) <- "integer"
   storage.mode(Two1) <- storage.mode(Two2) <- "integer"
   storage.mode(Three1) <- storage.mode(Three2) <- "integer"
   storage.mode(Fouid1) <- storage.mode(Fouid2) <- "integer"
   storage.mode(Legid1) <- storage.mode(Legid2) <- "integer"
   storage.mode(Sitxfm1) <- storage.mode(Sitxfm2) <- "integer"
   storage.mode(global.codes) <- "integer"
   storage.mode(OpenUnits) <- "integer"
#
#       Do the simulation. 
#
   z <- tryCatch(.Fortran("simulate",model=model.arg,IsPrecip=is.precip,
                  Nsites=as.integer(nsites),nattr=as.integer(nattr),
                  Nvars=as.integer(nvars),CurVar=as.integer(sim.order[i]),
                  MXP1=as.integer(Mxp1),MXP2=as.integer(Mxp2),
                  Sitinf1=sitinf1,Sitinf2=sitinf2,TmpFlNo=TempFile$Number,
                  NP1=Np1,NP2=Np2,
                  COVCODE1=Covcode1,COVCODE2=Covcode2,TWO1=Two1,TWO2=Two2,
                  THREE1=Three1,THREE2=Three2,FOUID1=Fouid1,FOUID2=Fouid2,
                  LEGID1=Legid1,LEGID2=Legid2,PWTID1=Pwtidx1,PWTID2=Pwtidx2,
                  SPMOD1=as.integer(Spmod1),SPMOD2=as.integer(Spmod2),
                  SITXFM1=Sitxfm1,SITXFM2=Sitxfm2,GLBCOD=global.codes,
                  MaxLag=max.lag,BETA1=Beta1,BETA2=Beta2,THETA1=Theta1,
                  THETA2=Theta2,PHI1=Phi1,PHI2=Phi2,RHO1=Rho1,RHO2=Rho2,
                  GLBVAL=global.vals,InitVals=mean.table,DateLims=DateLims,
                  UnitNos=OpenUnits,ReadFrom=in.unit,IFAIL=as.integer(0),
                  PACKAGE="Rglimclim",NAOK=TRUE)$IFAIL,
                 finally=CloseFiles(OpenUnits[-(95:96)]))
   CloseFiles(OpenUnits[-(95:96)])
   if (z!=0) {
    CloseFiles(OpenUnits)
    if (z == 200) {
     stop(err.msg(err.no=z,filename=CorFileName))
    } else if (z == 400) {
     stop(err.msg(err.no=z))
    } else {
     stop(err.msg(err.no=z))
    }
   }
   in.unit <- as.integer(96) # Read data generated so far when simulating 
                             # remaining variables
   OpenUnits[-(95:96)] <- 0
  }
#
#       Output the results for this simulation ... 
#
  if (OutOpt >= 2) {
    fnm <- make.simfilename(simdir,file.prefix,"Daily",sim,nsims)
    OpenUnits[20] <- OpenFile(filename=fnm, direct=FALSE)
    if (OpenUnits[20] < 0) {
      stop(paste("Failed to connect file ", fnm, ". Does the folder exist?", sep=""))
    }
  }
  if (OutOpt %in% c(1,3)) {
    fnm <- make.simfilename(simdir,file.prefix,"Monthly",sim,nsims)
    OpenUnits[21] <- OpenFile(filename=fnm, direct=FALSE)
    if (OpenUnits[21] < 0) {
      stop(paste("Failed to connect file ", fnm, ". Does the folder exist?", sep=""))
    }
  }
  storage.mode(OpenUnits) <- "integer"
  storage.mode(regs.wanted) <- storage.mode(region.codes) <- "integer"
  cat(paste(c(rep(" ",nchar(out.text)),"\r"),collapse=""))
  cat(paste("Simulation ",sim," of ",nsims,": writing output ...\r",sep=""))
  z <- tryCatch(.Fortran("WriteData",UnitNos=OpenUnits,MissVal=missval,
                 NSITES=as.integer(nsites),NVARS=as.integer(nvars),
                 TmpFlNo=TempFile$Number,OUTOPT=as.integer(OutOpt),
                 DateLims=DateLims,OutLims=OutLims,WHICHREGS=regs.wanted,
                 REGCODE=region.codes,NREGS=as.integer(nregions),
                 IFAIL=as.integer(0),PACKAGE="Rglimclim")$IFAIL,
                 finally=CloseFiles(OpenUnits[-95]))

  if (z!=0) {
   CloseFiles(OpenUnits)
   stop(err.msg(err.no=z))
  }
#
#       Close all files except the scratch file containing the data
#
  CloseFiles(OpenUnits[-95])
  OpenUnits[-95] <- 0
 }
 unlink(TempFile$Name) # This was opened earlier
 CloseFiles(OpenUnits)
 cat("\n")
#
#       And create a GLCsim object containing summary information about the
#       simulation
#
 z <- list(DataFile=datflnm,var.numbers=sort(which.vars),var.names=var.names,
           modeldefs=modeldefs,model.types=model.types,SimDir=simdir,
           FilePrefix=file.prefix,RNGstate=RNGstate,
           siteinfo=simsites,start=start,end=end,nsims=nsims,
           impute.until=impute.until,output=output,which.regions=which.regions,
           which.daily=which.daily,daily.start=daily.start,daily.end=daily.end,
           missval=missval)
 class(z) <- "GLCsim"
 z
}
##############################################################################
##############################################################################
##############################################################################
dep.calc <- function(dependencies) {
#
#       To convert a matrix of direct dependencies between variables
#       into a matrix of indirect dependencies (e.g. if V1 depends on 
#       V2 and V2 depends on V3, then V1 depends indirectly on V3)
#
 n <- nrow(dependencies)
 d <- dependencies
 z <- matrix(0,nrow=n,ncol=n)
 while (!identical(d,z)) {
  z <- d
  for (i in 1:n) {
   tmp <- d[i,] > 0; tmp[i] <- TRUE
   d[i,] <- as.numeric(colSums(d[tmp,,drop=FALSE]) > 0)
  }
 }
 z
}
##############################################################################
##############################################################################
##############################################################################
check.siteinfo <- function(siteinfo,max.pars,gen.sitinf=TRUE) {
#
#       To check that site information in object siteinfo is defined
#       correctly; and also to set up data structures required by
#       the site effects defined in a model with max.pars parameters
#       (if gen.sitinf=TRUE)
#
 if (!is.list(siteinfo)) {
  stop("siteinfo must be a list (use make.siteinfo() to create it)")
 }
 sitinf.needed <- c("Nsites","Site.names","Site.codes","Nattr",
                    "Attribute.names","Attribute.values")
 present <- sitinf.needed %in% names(siteinfo)
 if (!all(present)) {
  stop(c("The following components are missing from siteinfo:\n",
          paste(sitinf.needed[!present]," ")))
 }
 nsites <- as.integer(siteinfo$Nsites)
 nattr <- as.integer(siteinfo$Nattr)
#
#       First slice of sitinf contains site attributes; remaining
#       slices are for potential derivatives in nonlinear 
#       transformations
#
 sitinf=NULL
 if (gen.sitinf) {
  sitinf <- array(dim=c(nsites,max.pars,4))
  sitinf[,1:nattr,1] <- siteinfo$Attribute.values
 }
 z <- list(nsites=nsites,nattr=nattr,sitinf=sitinf)
 z
}
##############################################################################
##############################################################################
##############################################################################
make.simfilename <- function(simdir,file.prefix,tscale="Daily",simno,nsims) {
##############################################################################
#
#       To generate the name of an output file for simulated data.
#       file.prefix identifies the simulation; tscale is either "Daily"
#       or "Monthly", simno is the number of the current simulation and
#       nsims is the total number of simulations (used to ensure that
#       simno is padded with an appropriate number of zeroes throughout).
#
##############################################################################
 simno.format <- formatC(simno,digits=0,format="d",
                         width=floor(log(nsims,base=10))+1,flag="0")
 paste(simdir,"/",file.prefix,"_",tscale,"_Sim",simno.format,".dat",sep="")
}
##############################################################################
##############################################################################
##############################################################################
print.GLCsim <- function(x, ...) {
##############################################################################
#
#       Print method for objects of class GLCsim
#
##############################################################################
 title.text <- "Object of class GLCsim:"
 cat(paste(title.text,"\n",
           paste(rep("=",nchar(title.text)),collapse=""),"\n",sep=""))
 cat(paste("\nVariables taken from data file",x$DataFile,"\n"))
 cat("Variables simulated:\n")
  for (i in x$var.numbers) {
   cat(paste("\t",i,".\t",x$var.names[i],
             " (model type: ",x$model.types[i],")\n",sep="")) 
  }
 cat(paste("\nSimulation period: ",x$start %% 100,"/",x$start %/% 100,
           " to ",x$end %% 100,"/",x$end %/% 100,sep=""))
 if (x$impute.until >= x$start) {
 cat(paste(" (imputations from ",x$start %% 100,"/",
           x$start %/% 100," to ",x$impute.until %% 100,"/",
           x$impute.until %/% 100,")\n",sep=""))
 } else {
  cat("\nNo imputation performed\n") 
 }
 cat(paste(x$nsims," realisations generated\n",sep=""))
 outopt.text <- 
 cat(paste("\nOutput files generated: ",paste(x$output,collapse=" and "),
           "\n",sep=""))
 if ("daily" %in% x$output) {
  cat(paste("Daily output written from ",x$daily.start %% 100,"/",
            x$daily.start %/% 100," to ",x$daily.end %% 100,"/",
            x$daily.end %/% 100,"\n",sep=""))
 }
 if ("monthly" %in% x$output) {
  cat("Monthly summaries written for the following regions:\n")
  for (j in sort(x$which.regions)) {
   cat(paste("\t",j,"\t",x$siteinfo$Region.names[j+1],"\n",sep=""))
  }
 }
 cat(paste("Output directory: ",getwd(),"/",
            gsub("^.\\/","",x$SimDir),"\n",sep=""))
 cat(paste("Prefix for output filenames: ",x$FilePrefix,"\n",sep=""))
 invisible(NULL)
}
##############################################################################
##############################################################################
##############################################################################
summary.GLCsim <- function(object,which.variables,which.sites,which.regions,
                           which.timescales,thresholds,season.defs,...) {
##############################################################################
#
#       Summary method for objects of class GLCsim. This processes the 
#       output files generated by the GLCsim() routine and calculates 
#       tables of summary statistics for each simulation. These tables are
#       not necessarily designed to be read by the user; but the routine
#       is called by the plot method to produce its plots. Arguments:
#
#       which.variables Defaults to all of the variables produced in the 
#                       simulation
#       which.sites     Defaults to all of the sites for which data were
#                       simulated; otherwise a character vector of 
#                       4-character site codes. Summaries will be 
#                       produced for all of these sites.
#       which.regions   Defaults to all available regions. If non-null, summaries
#                       will be calculated for time series of daily averages
#                       over the selected regions
#       which.timescalesEither "daily", "monthly" or c("daily","monthly"). 
#                       Default is to process all of the output files
#                       produced in the simulation
#       thresholds      A vector the same length as which.variables. If
#                       present and if the element corresponding to 
#                       a particular variable is non-NA, the proportion of 
#                       exceedances of the corresponding threshold will
#                       be calculated for that variable.
#       season.defs     A list defining months to be grouped together to
#                       form "seasons" for monthly summaries. If not 
#                       specified, separate summaries are produced for each 
#                       month of the year. 
#
#       Value: a list with two components entitled "Daily" and "Monthly". 
#       The "Daily" component is itself a list with components "Sites" and 
#       "Regions". Each of these is a 5-dimensional array with dimensions 
#       simulation number, site, variable, month, statistic; the 
#       statistics produced are mean, standard deviation, maximum, minimum,
#       autocorrelation at lags 1-3, proportion of threshold exceedances 
#       (if a threshold is set), conditional means and standard deviations 
#       (again, if a threshold is set) and inter-variable correlations. 
#       The "Monthly" component is a 5-dimensional arrays with dimensions 
#       simulation number, region, variable, year and season. The values 
#       are the annual time series of means for each variable within the 
#       corresponding seasons. 
#
##############################################################################
#
#       Define variables for summary
#
 obj.name <- match.call()$object 
 vars.wanted <- object$var.numbers
 if (!missing(which.variables)) {
  if (!all(which.variables %in% object$var.numbers)) {
   stop(paste("Summary requested for variable that was not simulated in",
              obj.name))
  }
  vars.wanted <- which.variables
 }
 var.names <- object$var.names[vars.wanted]
#
#       Now sites ...
#
 sites.wanted <- object$siteinfo$Site.codes
 if (!missing(which.sites)) {
  if (!all(which.sites %in% object$siteinfo$Site.codes)) {
   stop(paste("Summary requested for site that is not present in '",
              obj.name,"'",sep=""))
  }
  sites.wanted <- which.sites
 }
#
#       ... regions (NB numbering starts at 0 for whole area)
#
 regions.wanted <- NULL
 if (!missing(which.regions)) {
  regions.wanted <- which.regions
  if (!all(regions.wanted %in% c(0,object$siteinfo$Regions))) {
   stop(paste("Summary requested for region that is not defined in '",
              obj.name,"'",sep=""))
  }
  if (!all(regions.wanted %in% object$which.regions)) {
   stop(paste("Not all requested regions have monthly data stored in '",
              obj.name,"'",sep=""))
  }
 } else {
  regions.wanted <- object$which.regions
 }
 region.names <- object$siteinfo$Region.names[1+regions.wanted]
 Nregions <- length(regions.wanted)
#
#       ... timescales ...
#
 timescales.wanted <- object$output
 if (!missing(which.timescales)) {
  timescales.wanted <- which.timescales
  if (!all(timescales.wanted %in% c(0,object$output))) {
   stop(paste("Summary requested for time scale that is not present in",
              obj.name))
  }
 }
#
#       ... thresholds ...
#
 wanted.thresholds <- rep(-1e100,length(vars.wanted))
 if (!missing(thresholds)) {
  if (length(thresholds) != length(vars.wanted)) {
   stop("Length of thresholds vector doesn't match number of variables requested")
  }
  wanted.thresholds[!is.na(thresholds)] <- thresholds[!is.na(thresholds)]
 }
#
#       ... and seasons. For the subsequent processing, when a season spans
#       the year end it will be convenient to label the months in the 
#       new year starting 13 (the code below also allows for "seasons" that span 
#       more than one year end)
#
 wanted.seasons <- as.list(1:12)
 if (!missing(season.defs)) {
  wanted.seasons <- season.defs
 }
 season.names <- sapply(wanted.seasons,
                        FUN=function(x) {
                         paste(month.abb[x],collapse=", ")
                        })
 wanted.seasons <- lapply(wanted.seasons,
                          FUN=function(x) {
                           x + (12 * cumsum(c(0,x[-1] < x[-length(x)])))
                          })
#
#       Now create the structure of the result
#
 z <- list(Daily=NULL,Monthly=NULL)
 thresh.text <- as.character(wanted.thresholds)
 thresh.text[wanted.thresholds <= -1e100] <- ""
 var.table <- data.frame(Variable=var.names,Threshold=thresh.text,
                         row.names=vars.wanted,stringsAsFactors=FALSE)
 reg.table <- NULL
 years <- NULL
 if (Nregions > 0) {
  reg.table <- data.frame(Region=region.names,row.names=regions.wanted,
                          stringsAsFactors=FALSE)
 }
 if ("daily" %in% timescales.wanted) {
  stat.names <- c("Mean","Std Dev","Max","Min","ACF1","ACF2","ACF3")
  if (any(wanted.thresholds > -1e100)) {
   stat.names <- c(stat.names,c("P(exceed threshold)","Conditional mean",
                              "Conditional std dev"))
  }
  if (length(var.names) > 1) {
   stat.names <- c(stat.names,paste("Correlation with",var.names))
  }
  z$Daily <- list(Sites=NULL,Regions=NULL)
   if (!is.null(sites.wanted)) {
   z$Daily$Sites <- array(dim=c(object$nsims,12,length(sites.wanted),
                                length(var.names),length(stat.names)))
   dimnames(z$Daily$Sites) <- list(Simulation=1:object$nsims,
                                   Month=1:12,Site=sites.wanted,
                                   Variable=var.names,
                                   Statistic=stat.names)
  }
  if (!is.null(regions.wanted)) {
   z$Daily$Regions <- array(dim=c(object$nsims,12,Nregions,
                                  length(var.names),length(stat.names)))
   dimnames(z$Daily$Regions) <- list(Simulation=1:object$nsims,
                                     Month=1:12,Region=region.names,
                                     Variable=var.names,
                                     Statistic=stat.names)
  }
 } else {
 }
 if ("monthly" %in% timescales.wanted) {
  years <- (object$start %/% 100):(object$end %/% 100)
  nyears <- length(years)
  if (!is.null(regions.wanted)) {
   z$Monthly <- array(dim=c(object$nsims,Nregions,
                          length(var.names),nyears,length(season.names)))
   dimnames(z$Monthly) <- list(Simulation=1:object$nsims,
                                       Region=region.names,
                                       Variable=var.names,
                                       Year=years,
                                       "Season (months)" = season.names)
  }  
 }
#
#       Prepare variables needed for Fortran call
#
 Nsites <- as.integer(object$siteinfo$Nsites)
 site.indicators <- rep(0,Nsites)
 site.indicators[match(sites.wanted,object$siteinfo$Site.codes)] <- 1
 Nvars <- as.integer(length(object$var.names))
 var.indicators <- rep(0,Nvars)
 var.indicators[vars.wanted] <- 1
 regions.matrix <- matrix(0,nrow=Nsites,ncol=Nregions)
 if (Nregions > 0) {
  for (i in 1:Nregions) {
   if (regions.wanted[i] == 0) {
    regions.matrix[,i] <- 1
   } else {
    regions.matrix[object$siteinfo$Regions==regions.wanted[i],i] <- 1
   }
  }
 }
 storage.mode(regions.matrix) <- "integer"
 stat.indicators <- c(rep(TRUE,10),rep(FALSE,Nvars))
 if (length(vars.wanted) > 1) {                 # only correlations
  stat.indicators[10+vars.wanted] <- TRUE       # with wanted vars
 }
 if (all(wanted.thresholds <= -1e100)) stat.indicators[8:10] <- FALSE
#
#       Now loop through the simulations and calculate the summaries (as usual,
#       need a temporary file to pass site codes to Fortran)
#
 cat("\n")
 TempFile <- OpenTempFile()
 write(object$siteinfo$Site.codes, file=TempFile$Name, append=TRUE)
 for (sim in 1:object$nsims) {
  cat(paste("\rWorking on simulation",sim,"of",object$nsims,"...\r"))
  if ("daily" %in% timescales.wanted) {
   simfilename <- make.simfilename(object$SimDir,object$FilePrefix,
                                   "Daily",sim,object$nsims)
   if (!file.exists(simfilename)) stop(paste("File",simfilename,"doesn't exist"))
   simfilno <- OpenFile(filename=simfilename,direct=FALSE)
   current.stats <- array(0,dim=c(12,Nsites+Nregions,Nvars,10+Nvars))
   stat.table <- tryCatch(.Fortran("DailyStats",FilNo=as.integer(simfilno),
                          Nsites=Nsites,Nvars=Nvars,Nreg=Nregions,
                          WantedSites=as.integer(site.indicators),
                          WantedVars=as.integer(var.indicators),
                          Thresholds=wanted.thresholds,
                          RegionDefs=regions.matrix,
                          TmpFlNo=TempFile$Number,
                          MoStats=current.stats,
                          Ifail=as.integer(0)),
                          finally=CloseFiles(simfilno))
   CloseFiles(simfilno)
   if (!stat.table$Ifail == 0) stop(err.msg(stat.table$Ifail))
   stat.table$MoStats[abs(stat.table$MoStats) > 9.999e99] <- NA
   z$Daily$Sites[sim,,,,] <- 
     stat.table$MoStats[,(1:Nsites)[object$siteinfo$Site.codes %in% sites.wanted],
                       vars.wanted,stat.indicators,drop=FALSE]
   if (Nregions > 0) {
#    wanted.col <- Nsites+as.numeric(0%in%regions.wanted)+regions.wanted # Corrected in v1.3-7
    wanted.col <- Nsites+1:Nregions
    z$Daily$Regions[sim,,,,] <- stat.table$MoStats[,wanted.col,vars.wanted,
                                                 stat.indicators,drop=FALSE]   
   }
  }
  if ("monthly" %in% timescales.wanted) {
   if (Nregions > 0) {
    simfilename <- make.simfilename(object$SimDir,object$FilePrefix,
                                    "Monthly",sim,object$nsims)
#
#       Read data and process each variable in turn, each time putting
#       the data for the current variable into an array (Region * Year * Month)
#       and looping over the seasons (can't think of a neat way of avoiding
#       loops here ...). Note that the for (i in 1:min(max.yrspan,nyears)-1) loop 
#       appends extra "columns" to the temporary array in case of seasons
#       that span year-ends. The min() here is in case the simulation doesn't
#       span enough years to get started ...
#
    current.data <- read.fwf(simfilename,widths=c(4,4,rep(7,max(vars.wanted)*13)))
    current.data <- current.data[current.data[,2] %in% regions.wanted,]
    current.data[current.data == object$missval] <- NA
    max.yrspan <- ceiling(max(unlist(wanted.seasons))/12)
    for (i in 1:length(vars.wanted)) {
     v <- vars.wanted[i]
     tmp <- array(dim=c(Nregions,nyears,12*max.yrspan))
     tmp[,,1:12] <- array(unlist(current.data[,2 + (13*(v-1)+(1:12))]))
     if ((max.yrspan > 1) & (nyears > 1)) {
      for (j in 1:(min(max.yrspan,nyears)-1)) {
       tmp[,1:(nyears-j),(j*12)+(1:12)] <- tmp[,-(1:j),1:12]
      }
     }
     for (j in 1:length(wanted.seasons)) {
      cur.seas <- wanted.seasons[[j]]
#
#       indexing is sim, region, variable, year, season
#
      z$Monthly[sim,,i,,j] <- apply(tmp[,,cur.seas,drop=FALSE],
                                                     MARGIN=1:2,FUN=mean)
     }
    }
   }
  }
 }
 unlink(TempFile$Name)
 cat("\rAll done.                                                      \n")
 attr(z,which="GLCsim.object") <- obj.name
 attr(z,which="var.table") <- var.table
 attr(z,which="thresholds") <- wanted.thresholds
 attr(z,which="region.table") <- reg.table
 attr(z,which="Years") <- years
 class(z) <- "summary.GLCsim"         
 z
}
##############################################################################
##############################################################################
##############################################################################
print.summary.GLCsim <- function(x,...) {
##############################################################################
#
#       print method for objects of class summary.GLCsim. This just produces
#       a concise printout of the summaries that are present in the object, 
#       and can be useful for (e.g.) finding the relevant names that should 
#       be used in a call to plot.summary.GLCsim.
#
##############################################################################
#
#       Now create the structure of the result and write a summary of this
#       structure to screen
#
 title.text <- paste("Simulation '",attr(x,which="GLCsim.object"),
                     "' - summaries calculated:",sep="")
 cat(paste("\n",paste(rep("#",nchar(title.text)),collapse=""),"\n",
           title.text,"\n",
           paste(rep("#",nchar(title.text)),collapse=""),"\n",sep=""))
 var.table <- attr(x,which="var.table")
 wanted.thresholds <- attr(x,which="thresholds")
 cat("\nVariables:\n")
 cat("----------\n")
 if (is.null(x$Daily) | all(wanted.thresholds <= -1e100)) {
  print(var.table[,1,drop=FALSE])
 } else {
  print(var.table)
  if (any(wanted.thresholds <= -1e100)) {
   cat("-----\nNote: thresholds are displayed only when defined\n")
  }
 }
 if (!is.null(x$Daily)) {
  cat("\nDaily summaries, by month:\n==========================\n")
  if (!is.null(x$Daily$Sites)) {
   stat.names <- dimnames(x$Daily$Sites)$Statistic
  } else if (!is.null(x$Daily$Regions)) {
   stat.names <- dimnames(x$Daily$Regions)$Statistic
  }  
  if (!is.null(x$Daily$Sites)) {
   cat("\nSites:\n------\n")
   write.CharVecAsTable(dimnames(x$Daily$Sites)$Site,options()$width-5)
  } else {
   cat("\n##### Individual site summaries not requested\n")
  }
  if (!is.null(x$Daily$Regions)) {
   cat("\nRegions:\n--------\n")
   print(attr(x,which="region.table"))
  } else {
   cat("\n##### No regional summaries requested for daily data\n")
  }
  cat("\nSummary statistics calculated:\n------------------------------\n")
  write.CharVecAsTable(stat.names,options()$width-5)
  if (any(wanted.thresholds > -1e100)) {
   cat("-----\n'Conditional' statistics are conditional on threshold exceedance\n")
  }
 } else {
  cat("\n#####Daily summaries not requested or unavailable\n")
 }
 if (!is.null(x$Monthly)) {
  cat("\nMonthly / seasonal means:\n=========================\n")
  write.CharVecAsTable(dimnames(x$Monthly)$`Season (months)`,options()$width-5)
 } else {
  cat("\n##### Monthly summaries not requested or unavailable\n")
 }
 NULL
}
##############################################################################
##############################################################################
##############################################################################
plot.summary.GLCsim <- function(x,imputation,quantiles,
                                which.variables,which.sites,which.regions,
                                which.timescales,which.stats,which.seasons,
                                plot.titles,ylabs,colours.sim="greyscale",
                                colour.obs="black",...) {
##############################################################################
#
#       To plot a summary of objects of class GLCsim. Arguments:
#       
#       x               - A summary for an object of class GLCsim
#       imputation      - An optional summary for another object of class
#                         GLCsim, that is treated as a set of imputations
#                         conditioned upon all available data
#       quantiles       - A set of quantiles to plot for the simulated
#                         distributions. Default is 
#                         c(0,0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99,1)
#       which.variables - Character vector selecting the variables for 
#                         which summaries are to be plotted. These should 
#                         be specified as the variable *names* - note that
#                         this differs from the usage in both the summary 
#                         and plot methods for objects of class GLCsim, 
#                         where the variables should be specified by number
#       which.sites     - Character vector specifying the sites for which
#                         summaries are to be plotted (4-character site
#                         codes)
#       which.regions   - Character vector specifying the names of regions
#                         for which summaries are to be plotted
#       which.timescales- Either "daily", "monthly" or c("daily","monthly").
#       which.stats     - For summaries of daily data, a character vector
#                         specifying the summary statistics to be plotted
#       which.seasons   - For summaries of monthly data, a character vector 
#                         specifying the names of the seasons for which 
#                         summaries are to be plotted.
#       plot.titles     - A vector of titles for the plots. It is the 
#                         user's responsibility to ensure that the 
#                         length of this vector matches the number of
#                         plots actually produced. If omitted, the routine
#                         will construct plot titles automatically.
#       ylabs           - A vector of y axis labels for the plots. Ditto.
#       colours.sim     - Either "greyscale" (the default), "colour"
#                         or a vector of valid colour specifiers, of
#                         length length(quantiles)-1. These colours
#                         will be used to shade the simulated distributions.
#                         "greyscale" produces plots that are suitable for
#                         inclusion in printed material; "colour" 
#                         uses a default colour scale generated using the 
#                         rainbow() command. 
#       colour.obs      - The colour to use for the imputation envelope
#                         if an imputation object is supplied. 
#       ...             - Additional arguments to plot functions.
#
#       The default for the selection arguments which.xxxxx is to include 
#       all the values for which information is present in the summary 
#       object. The routine plots the daily summaries first if requested, 
#       and then the monthly summaries. For each set of summaries, it loops
#       over variables, then sites / regions, then individual statistics 
#       (for daily summaries) or seasons (for monthly summaries). No attempt
#       is made to control the layout of plots; the user should therefore
#       set up the plot layout before calling the function.
#
##############################################################################
 obj.name <- match.call()$x
 if (!missing(imputation)) {
  if (class(imputation)[1] != "summary.GLCsim") {
   stop("imputation must be an object of class summary.GLCsim")
  }
 }
 imputation.name <- match.call()$imputation
 if (missing(quantiles)) {
  q.vec <- c(0,0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99,1)
 } else {
  if (any(quantiles < 0 | quantiles > 1)) {
   stop("Quantiles must be between 0 and 1")
  }
  q.vec <- quantiles
 }
#
#       Set up colour scale for plotting
#
 npolys <- length(q.vec)-1
 if (isTRUE(colours.sim == "greyscale")) {
  colscl <- grey(1 - (1:npolys)/(2*npolys))
 } else if (isTRUE(colours.sim == "colour")) {
  colscl <- rainbow(1+1.4*npolys)[-1]
 } else {
  colscl <- colours.sim
  if (length(colscl) != npolys) {
   stop("Number of colours doesn't match number of requested quantiles")
  }
 }
#
#       Find candidate values for the choice of variables, sites, regions etc.
#       There's a bit of redundancy in the code below, but it's quick and 
#       the end result is that everything is set correctly.
#
 sites.wanted <- regions.wanted <- stats.wanted <- seasons.wanted <- NULL
 if (!is.null(x$Daily$Sites)) {
  vars.wanted <- dimnames(x$Daily$Sites)$Variable
  sites.wanted <- dimnames(x$Daily$Sites)$Site
  stats.wanted <- dimnames(x$Daily$Sites)$Statistic
 } 
 if (!is.null(x$Daily$Regions)) {
  vars.wanted <- dimnames(x$Daily$Regions)$Variable
  regions.wanted <- dimnames(x$Daily$Regions)$Region
  stats.wanted <- dimnames(x$Daily$Regions)$Statistic
 } 
 if (!is.null(x$Monthly)) {
  vars.wanted <- dimnames(x$Monthly)$Variable
  regions.wanted <- dimnames(x$Monthly)$Region
  seasons.wanted <- dimnames(x$Monthly)$Season
 }
 timescales.wanted <- 
            c("daily","monthly")[c(!is.null(x$Daily),!is.null(x$Monthly))]
#
#       Choose variables ...
#
 if (!missing(which.variables)) {
  tmp <- (which.variables %in% vars.wanted)
  if (!all(tmp)) {
   warning(paste("Plots have been requested for variables that are not ",
                 "present\n  in object '",obj.name,"' - plots will be ",
                 "produced only for\n  variables that are actually present\n",
                 sep=""))
  }
  vars.wanted <- which.variables[tmp]
 }
#
#       ... locations ...
#
 if (!missing(which.sites)) {
  tmp <- (which.sites %in% sites.wanted)
  if (!all(tmp)) {
   warning(paste("Plots have been requested for sites that are not ",
                 "present\n  in object '",obj.name,"' - plots will be ",
                 "produced only for\n  sites that are actually present\n",
                 sep=""))
  }
  sites.wanted <- which.sites[tmp]
 }
 if (!missing(which.regions)) {
  tmp <- (which.regions %in% regions.wanted)
  if (!all(tmp)) {
   warning(paste("Plots have been requested for regions that are not ",
                 "present\n  in object '",obj.name,"' - plots will be ",
                 "produced only for\n  regions that are actually present\n",
                 sep=""))
  }
  regions.wanted <- which.regions[tmp]
 }
#
#       ... timescales ...
#
 if (!missing(which.timescales)) {
  tmp <- (which.timescales %in% timescales.wanted)
  if (!all(tmp)) {
   warning(paste("Plots have been requested a time scale that is not ",
                 "present\n  in object '",obj.name,"' - plots will be ",
                 "produced only for\n  time scales that are actually present\n",
                 sep=""))
  }
  timescales.wanted <- which.timescales[tmp]
 }
#
#       ... and statistics / seasons ...
#
 if (!missing(which.stats)) {
  tmp <- (which.stats %in% stats.wanted)
  if (!all(tmp)) {
   warning(paste("Plots have been requested for statistics that are not ",
                 "present\n  in object '",obj.name,"' - plots will be ",
                 "produced only for\n  statistics that are actually present\n",
                 sep=""))
  }
  stats.wanted <- which.stats[tmp]
 }
 if (!missing(which.seasons)) {
  tmp <- (which.seasons %in% seasons.wanted)
  if (!all(tmp)) {
   warning(paste("Plots have been requested for seasons that are not ",
                 "present\n  in object '",obj.name,"' - plots will be ",
                 "produced only for\n  seasons that are actually present\n",
                 sep=""))
  }
  seasons.wanted <- which.seasons[tmp]
 }
#
#       Now loop through and produce plots. NB at the point where 
#       cur.table is defined (which should be an array of dimension
#       nsims*12) need to explicitly set the array size in case 
#       someone has only one month or one simulation (can't use
#       drop=FALSE, because this retains *all* dimensions AAARGH)
#
 plot.num <- 0
 for (cur.var in vars.wanted) {
  if ("daily" %in% timescales.wanted) {
   for (site in sites.wanted) {
    for (stat in stats.wanted) {
     plot.num <- plot.num + 1
     if (missing(plot.titles)) {
      cur.title <- paste("Site ",site,", variable ",cur.var,
                         ":\n",stat,sep="")
      } else {
      cur.title <- plot.titles[plot.num]
     }
     if (missing(ylabs)) {
      cur.label <- "Quantiles"
     } else {
      cur.label <- ylabs[plot.num]
     }
     cur.table <- array(x$Daily$Sites[,,site,cur.var,stat],
                        dim=dim(x$Daily$Sites)[1:2])
     q.tab <- apply(cur.table,MARGIN=2,FUN=quantile,probs=q.vec,na.rm=TRUE)
     q.tab.imputed <- NULL
     if (!missing(imputation)) {
      cur.table <- try(array(imputation$Daily$Sites[,,site,cur.var,stat],
                             dim=dim(imputation$Daily$Sites)[1:2]),
                       silent=TRUE)
      if (class(cur.table)[1] == "try-error") {
       warning(paste("Object '",imputation.name,"' doesn't contain information ",
                     "for\nvariable",cur.var,", site ",site,", statistic ",
                     "stat",sep=""))
      } else {
       q.tab.imputed <- apply(cur.table,MARGIN=2,FUN=quantile,
                                                 probs=c(0,1),na.rm=TRUE)
      }
     }
     quantileplot(q.tab,q.tab.imputed,"monthly",cur.label,cur.title,colscl,
                                                            colour.obs,...)
    }
   }
   for (region in regions.wanted) {
    for (stat in stats.wanted) {
     plot.num <- plot.num + 1
     if (missing(plot.titles)) {
      cur.title <- paste(region,", mean of variable ",cur.var,
                         ":\n",stat,sep="")
     } else {
      cur.title <- plot.titles[plot.num]
     }
     if (missing(ylabs)) {
      cur.label <- "Quantiles"
     } else {
      cur.label <- ylabs[plot.num]
     }
     cur.table <- array(x$Daily$Regions[,,region,cur.var,stat],
                        dim=dim(x$Daily$Regions)[1:2])
     q.tab <- apply(cur.table,MARGIN=2,FUN=quantile,probs=q.vec,na.rm=TRUE)
     q.tab.imputed <- NULL
     if (!missing(imputation)) {
      cur.table <- try(array(imputation$Daily$Regions[,,region,cur.var,stat],
                             dim=dim(imputation$Daily$Regions)[1:2]),silent=TRUE)
      if (class(cur.table)[1] == "try-error") {
       warning(paste("Object '",imputation.name,"' doesn't contain information ",
                     "for\nvariable",cur.var,", region ",region,", statistic ",
                     stat,sep=""))
      } else {
       q.tab.imputed <- apply(cur.table,MARGIN=2,FUN=quantile,
                                                 probs=c(0,1),na.rm=TRUE)
      }
     }
     quantileplot(q.tab,q.tab.imputed,"monthly",cur.label,cur.title,colscl,
                                                            colour.obs,...)
    }
   }
  }
  if ("monthly" %in% timescales.wanted) {
   if (length(attr(x,which="Years")) < 2) {
    stop("Annual plots can't be produced with < 2 years' data")
   }
   for (region in regions.wanted) {
    for (season in seasons.wanted) {
     plot.num <- plot.num + 1
     if (missing(plot.titles)) {
      cur.title <- paste(region,", mean of ",cur.var,
                         ",\n",season,sep="")
     } else {
      cur.title <- plot.titles[plot.num]
     }
     if (missing(ylabs)) {
      cur.label <- "Quantiles"
     } else {
     cur.label <- ylabs[plot.num]
     }
     cur.table <- array(x$Monthly[,region,cur.var,,season],
                        dim=dim(x$Monthly)[c(1,4)])
     q.tab <- apply(cur.table,MARGIN=2,FUN=quantile,probs=q.vec,na.rm=TRUE)
     q.tab.imputed <- NULL
     if (!missing(imputation)) {
      if (!identical(dimnames(imputation$Monthly)[4], dimnames(x$Monthly)[4])) {
       stop("objects 'x' and 'imputation' contain summaries for different years")
      }
      cur.table <- try(array(imputation$Monthly[,region,cur.var,,season],
                             dim=dim(imputation$Monthly)[c(1,4)]),silent=TRUE)
      if (class(cur.table)[1] == "try-error") {
       warning(paste("Object '",imputation.name,"' doesn't contain seasonal ",
                     "information\nfor variable",cur.var,", region ",region,
                     ",\nseason (",season,")",sep=""))
      } else {
       q.tab.imputed <- apply(cur.table,MARGIN=2,FUN=quantile,
                                                 probs=c(0,1),na.rm=TRUE)
      }
     }
     quantileplot(q.tab,q.tab.imputed,"yearly",cur.label,cur.title,colscl,
                                    colour.obs,years=attr(x,which="Years"),...)
    }
   }
  }
 }
 invisible(NULL)
}
##############################################################################
##############################################################################
##############################################################################
write.CharVecAsTable <- function(x,width) {
#
#       To write the elements of a character vector as a table
#
 cols.needed <- floor(width/(7+max(nchar(x))))
 x.table <- matrix("",nrow=cols.needed,                         # nrow and ncol
                      ncol=ceiling(length(x)/cols.needed))      # are *right* -
                                                                # will transpose
 x.table[1:length(x)] <- x
 x.table <- t(x.table)
 write.table(format(x.table,width=max(nchar(x.table))+1),
              quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
}
##############################################################################
##############################################################################
##############################################################################
quantileplot <- function(sim,obs,timescale,ylabel,plot.title,colscl,
                                                  colour.obs,years,...) {
#
#       To produce "caterpillar plots" of monthly summary statistics or
#       annual seasonal means, optionally overlaying an imputation envelope.
#
 ylim <- extendrange(c(as.numeric(sim),as.numeric(obs)))
 if (timescale == "monthly") {
  x.grid <- 1:12
  x.poly <- c(x.grid,rev(x.grid))
  xlab="Month"
 } else if (timescale == "yearly") {
  x.grid <- years
  x.poly <- c(x.grid,rev(x.grid))
  xlab="Year"
 }
 plot(x.grid,sim[1,],type="n",xlab=xlab,ylab=ylabel,main=plot.title,
      ylim=ylim,...)
 for (j in 1:(nrow(sim)-1)) {
  y.poly <- c(sim[j,],rev(sim[j+1,]))
  polygon(x.poly[!is.na(y.poly)],y.poly[!is.na(y.poly)],col=colscl[j])
 }
 if (!is.null(obs)) {
  y.poly <- c(obs[1,],rev(obs[2,]))
  polygon(x.poly[!is.na(y.poly)],y.poly[!is.na(y.poly)],col=colour.obs,
          border=colour.obs,lwd=2)  
 }
 invisible(NULL)
}

