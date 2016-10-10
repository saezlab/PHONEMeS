#Classes:
#background K/P-S interactions class
setClass("KPSbg",
         representation=representation(
           interactions="data.frame",
           species="character"))
#GMMres class: res is the results from the GMM (list with a slot per peptide)
#IDmap is a data frame that maps the ppetide IDs to the S.cc (sites IDs from the network)
#resFC is the same as res but with a FC column
setClass("GMMres",
         representation=representation(
           res="list",
           IDmap="data.frame",
           resFC="list"))
#GMMbyCond
setClass("GMMbyCond",
         representation=representation(
           dataBC="list",
           speciesP="list"))
#PKNlist
setClass("PKNlist",
         representation=representation(
           interactions="data.frame",
           interactionsD="data.frame",
           species="character",
           sinks="character",
           integrators="character",
           intermediates="character"))
#OptParam
setClass("OptParam",
         representation=representation(
           sizeP="numeric",
           nG1="numeric",
           cstart="logical",
           intgAsintm="logical",
           nScripts="numeric",
           absTol="logical",
           tol="numeric",
           cap="numeric",
           resN="numeric",
           nG="numeric"))
#nG1 is the number of models per script, nG is the number of generations           
#ModelList
setClass("ModelList",
         representation=representation(
           model="data.frame",
           scores="numeric",
           scoresList="list",
           intAnd="logical",
           intgAnd="logical"))
##Checks: 
#KPSbg: I need the data frame to have colums K.ID, S.ID, S.cc, SID
#everything else is up to the user
setValidity("KPSbg",
            function(object){
              isValid<-TRUE
              if(length(colnames(object@interactions)) > 0){
                isValid<-all(c("K.ID", "S.ID", "S.cc", "SID") %in% colnames(object@interactions))
              }
              if (!isValid)
                cat("You interactions data frame must contain columns named K.ID, S.ID, S.cc and SID")
              if(!is.null(object@species) && isValid){
                isValid<-all(object@species %in% c(object@interactions$K.ID, object@interactions$S.cc))
              }
              if (!isValid)
                cat("species must match the names in the K.ID and S.cc columns")
              return(isValid)
            })
#GMMres: res field needs to have named slots (a peptide per ID), 
#and each slot has colums "clus" and a column "status"
#IDmap field needs to have a "dataID" column and a "S.cc" column
setValidity("GMMres",
            function(object){
              isValid<-TRUE
              if(!length(object@res) == 0 && is.null(names(object@res))){
                cat("The res field must be a named list")
                isValid<-FALSE
              }
              if(!length(object@res) == 0 && !is.null(colnames(object@res[1]))){
                colTable<-table(unlist(lapply(object@res, function(x){return(colnames(x))})))
                isValid<-all(c("clus", "status") %in% names(colTable)) && all(colTable == dim(object@res)[1])
                if(!isValid)
                  cat("The res field must contain matrices/data frames with identical column names, which must include clus and status")
              }
              if(!dim(object@IDmap)[1] == 0 && isValid){
                isValid<-all(c("dataID", "S.cc") %in% colnames(object@IDmap))
                if(!isValid) cat("IDmap field must contain columns dataID and S.cc")
              }
              return(isValid)
            })
#KPSbg getter and setter
setGeneric("interactions", function(object) standardGeneric("interactions"))
setMethod("interactions", "KPSbg", function(object) object@interactions)
setGeneric("interactions<-", function(object,value) standardGeneric("interactions<-"))
setReplaceMethod("interactions", signature(object="KPSbg",
                                 value="data.frame"),
                 function(object, value) {
                   object@interactions <- value
                   if (validObject(object))
                     return(object)
                 })
setGeneric("species", function(object) standardGeneric("species"))
setMethod("species", "KPSbg", function(object) object@species)
setGeneric("species<-", function(object,value) standardGeneric("species<-"))
setReplaceMethod("species", signature(object="KPSbg",
                                           value="character"),
                 function(object, value) {
                   object@species <- value
                   if (validObject(object))
                     return(object)
                 })
#KPSbg show method
setMethod("show",
          "KPSbg",
          function(object) {
            cat("Object of class",class(object),"\n")
            cat("Interactions:",dim(interactions(object))[1],"\n")
            cat("species:",length(species(object)),"\n")
            cat("K/P:",length(unique(interactions(object)[,"K.ID"])) ,"\n")
            cat("Sites:",length(unique(interactions(object)[,"S.cc"])) ," on ", length(unique(interactions(object)[,"S.ID"])), "substrates \n")
          })
#GMMres getter and setter
setGeneric("res", function(object) standardGeneric("res"))
setMethod("res", "GMMres", function(object) object@res)
setGeneric("res<-", function(object,value) standardGeneric("res<-"))
setReplaceMethod("res", signature(object="GMMres",
                                           value="list"),
                 function(object, value) {
                   object@res <- value
                   if (validObject(object))
                     return(object)
                 })
setGeneric("IDmap", function(object) standardGeneric("IDmap"))
setMethod("IDmap", "GMMres", function(object) object@IDmap)
setGeneric("IDmap<-", function(object,value) standardGeneric("IDmap<-"))
setReplaceMethod("IDmap", signature(object="GMMres",
                                    value="data.frame"),
                 function(object, value) {
                   object@IDmap <- value
                   if (validObject(object))
                     return(object)
                 })
setGeneric("resFC", function(object) standardGeneric("resFC"))
setMethod("resFC", "GMMres", function(object) object@resFC)
setGeneric("resFC<-", function(object,value) standardGeneric("resFC<-"))
setReplaceMethod("resFC", signature(object="GMMres",
                                           value="list"),
                 function(object, value) {
                   object@res <- value
                   if (validObject(object))
                     return(object)
                 })                 
#GMMres show method
setMethod("show",
          "GMMres",
          function(object) {
            cat("Object of class",class(object),"\n")
            cat("GMM results list:",length(res(object)),"peptides \n")
            cat("Conditions:",length(rownames(res(object)[[1]]))-1,"\n")
            cat("ID map:",!(dim(IDmap(object))[1] == 0) ,"\n")
          })
#GMMbyCond getter and setter
setGeneric("dataBC", function(object) standardGeneric("dataBC"))
setMethod("dataBC", "GMMbyCond", function(object) object@dataBC)
setGeneric("dataBC<-", function(object,value) standardGeneric("dataBC<-"))
setReplaceMethod("dataBC", signature(object="GMMbyCond",
                                  value="list"),
                 function(object, value) {
                   object@dataBC <- value
                   return(object)
                 })
setGeneric("speciesP", function(object) standardGeneric("speciesP"))
setMethod("speciesP", "GMMbyCond", function(object) object@speciesP)
setGeneric("speciesP<-", function(object,value) standardGeneric("speciesP<-"))
setReplaceMethod("speciesP", signature(object="GMMbyCond",
                                    value="list"),
                 function(object, value) {
                   object@speciesP <- value
                   return(object)
                 })
#####
#GMMbyCond show method
setMethod("show",
          "GMMbyCond",
          function(object) {
            cat("Object of class",class(object),"\n")
            cat("GMM results list, by condition:",dim(dataBC(object)[[1]][[1]])[1],"peptides \n")
            cat("Conditions:",length(dataBC(object)),"\n")
            sP<-unlist(lapply(speciesP(object), function(x){lapply(x,length)}))
            cat("Species perturbed:" ,"\n")
            for(i in 1:length(sP)){
              cat("\t",names(sP)[i],":", sP[i] ,"\n")
            }
          })
####
#PKNlist getter and setter
setMethod("interactions", "PKNlist", function(object) object@interactions)
setGeneric("interactions<-", function(object,value) standardGeneric("interactions<-"))
setReplaceMethod("interactions", signature(object="PKNlist",
                                           value="data.frame"),
                 function(object, value) {
                   object@interactions <- value
                   return(object)
                 })
setGeneric("interactionsD", function(object) standardGeneric("interactionsD"))
setMethod("interactionsD", "PKNlist", function(object) object@interactionsD)
setGeneric("interactionsD<-", function(object,value) standardGeneric("interactionsD<-"))
setReplaceMethod("interactionsD", signature(object="PKNlist",
                                  value="data.frame"),
                 function(object, value) {
                   object@interactionsD <- value
                   return(object)
                 })
setMethod("species", "PKNlist", function(object) object@species)
setReplaceMethod("species", signature(object="PKNlist",
                                      value="character"),
                 function(object, value) {
                   object@species <- value
                 })
setGeneric("sinks", function(object) standardGeneric("sinks"))
setMethod("sinks", "PKNlist", function(object) object@sinks)
setGeneric("sinks<-", function(object,value) standardGeneric("sinks<-"))
setReplaceMethod("sinks", signature(object="PKNlist",
                                  value="character"),
                 function(object, value) {
                   object@sinks <- value
                   return(object)
                 })
setGeneric("integrators", function(object) standardGeneric("integrators"))
setMethod("integrators", "PKNlist", function(object) object@integrators)
setGeneric("integrators<-", function(object,value) standardGeneric("integrators<-"))
setReplaceMethod("integrators", signature(object="PKNlist",
                                    value="character"),
                 function(object, value) {
                   object@integrators <- value
                   return(object)
                 })
setGeneric("intermediates", function(object) standardGeneric("intermediates"))
setMethod("intermediates", "PKNlist", function(object) object@intermediates)
setGeneric("intermediates<-", function(object,value) standardGeneric("intermediates<-"))
setReplaceMethod("intermediates", signature(object="PKNlist",
                                          value="character"),
                 function(object, value) {
                   object@intermediates <- value
                   return(object)
                 })

#PKNlist show method
setMethod("show",
          "PKNlist",
          function(object) {
            cat("Object of class",class(object),"\n")
            cat("Interactions:",dim(interactions(object))[1],"\n")
            cat("Direct interactions with targets:",dim(interactionsD(object))[1],"\n")
            cat("Species:",length(species(object)),"\n")
            cat("Sinks:",length(sinks(object)) ,"\n")
            cat("Integrators:",length(integrators(object)), "\n")
            cat("Intermediates:",length(intermediates(object)), "\n")
          })
#OptParam getter and setter
setGeneric("sizeP", function(object) standardGeneric("sizeP"))
setMethod("sizeP", "OptParam", function(object) object@sizeP)
setGeneric("sizeP<-", function(object,value) standardGeneric("sizeP<-"))
setReplaceMethod("sizeP", signature(object="OptParam",
                                            value="numeric"),
                 function(object, value) {
                   object@sizeP <- value
                   return(object)
                 })
setGeneric("nG1", function(object) standardGeneric("nG1"))
setMethod("nG1", "OptParam", function(object) object@nG1)
setGeneric("nG1<-", function(object,value) standardGeneric("nG1<-"))
setReplaceMethod("nG1", signature(object="OptParam",
                                    value="numeric"),
                 function(object, value) {
                   object@nG1 <- value
                   return(object)
                 })
setGeneric("cstart", function(object) standardGeneric("cstart"))
setMethod("cstart", "OptParam", function(object) object@cstart)
setGeneric("cstart<-", function(object,value) standardGeneric("cstart<-"))
setReplaceMethod("cstart", signature(object="OptParam",
                                  value="logical"),
                 function(object, value) {
                   object@cstart <- value
                   return(object)
                 })
setGeneric("intgAsintm", function(object) standardGeneric("intgAsintm"))
setMethod("intgAsintm", "OptParam", function(object) object@intgAsintm)
setGeneric("intgAsintm<-", function(object,value) standardGeneric("intgAsintm<-"))
setReplaceMethod("intgAsintm", signature(object="OptParam",
                                  value="logical"),
                 function(object, value) {
                   object@intgAsintm <- value
                   return(object)
                 })
setGeneric("nScripts", function(object) standardGeneric("nScripts"))
setMethod("nScripts", "OptParam", function(object) object@nScripts)
setGeneric("nScripts<-", function(object,value) standardGeneric("nScripts<-"))
setReplaceMethod("nScripts", signature(object="OptParam",
                                         value="numeric"),
                 function(object, value) {
                   object@nScripts <- value
                   return(object)
                 })
setGeneric("absTol", function(object) standardGeneric("absTol"))
setMethod("absTol", "OptParam", function(object) object@absTol)
setGeneric("absTol<-", function(object,value) standardGeneric("absTol<-"))
setReplaceMethod("absTol", signature(object="OptParam",
                                       value="logical"),
                 function(object, value) {
                   object@absTol <- value
                   return(object)
                 })
setGeneric("tol", function(object) standardGeneric("tol"))
setMethod("tol", "OptParam", function(object) object@tol)
setGeneric("tol<-", function(object,value) standardGeneric("tol<-"))
setReplaceMethod("tol", signature(object="OptParam",
                                     value="numeric"),
                 function(object, value) {
                   object@tol <- value
                   return(object)
                 })
setGeneric("cap", function(object) standardGeneric("cap"))
setMethod("cap", "OptParam", function(object) object@cap)
setGeneric("cap<-", function(object,value) standardGeneric("cap<-"))
setReplaceMethod("cap", signature(object="OptParam",
                                  value="numeric"),
                 function(object, value) {
                   object@cap <- value
                   return(object)
                 })
setGeneric("resN", function(object) standardGeneric("resN"))
setMethod("resN", "OptParam", function(object) object@resN)
setGeneric("resN<-", function(object,value) standardGeneric("resN<-"))
setReplaceMethod("resN", signature(object="OptParam",
                                  value="numeric"),
                 function(object, value) {
                   object@resN <- value
                   return(object)
                 })   
setGeneric("nG", function(object) standardGeneric("nG"))
setMethod("nG", "OptParam", function(object) object@nG)
setGeneric("nG<-", function(object,value) standardGeneric("nG<-"))
setReplaceMethod("nG", signature(object="OptParam",
                                  value="numeric"),
                 function(object, value) {
                   object@nG <- value
                   return(object)
                 })                 
#OptParam show method
setMethod("show",
          "OptParam",
          function(object) {
            cat("Object of class",class(object),"\n")
            cat("Size penalty:",sizeP(object),"\n")
            cat("Number of models per job:",nG1(object),"\n")
            cat("Number of parallel jobs:", nScripts(object),"\n")
            cat("Number of generations:", nG(object),"\n")
            cat("Force in direct substrates of drug targets:",cstart(object),"\n")
            cat("Tolerance:",tol(object), ifelse(absTol(object), "(absolute)", "(relative"),"\n")
            cat("Cap for weights adjustment:",cap(object),"\n")
            cat("Treat integrators as intermediates:",intgAsintm(object) ,"\n")
            cat("Index of the optimisation:",resN(object) ,"\n")
          })

#ModelList getter and setter
setGeneric("model", function(object) standardGeneric("model"))
setMethod("model", "ModelList", function(object) object@model)
setGeneric("model<-", function(object,value) standardGeneric("model<-"))
setReplaceMethod("model", signature(object="ModelList",
                                    value="data.frame"),
                 function(object, value) {
                   object@model <- value
                   return(object)
                 })
setGeneric("scores", function(object) standardGeneric("scores"))
setMethod("scores", "ModelList", function(object) object@scores)
setGeneric("scores<-", function(object,value) standardGeneric("scores<-"))
setReplaceMethod("scores", signature(object="ModelList",
                                    value="numeric"),
                 function(object, value) {
                   object@scores <- value
                   return(object)
                 })
setGeneric("scoresList", function(object) standardGeneric("scoresList"))
setMethod("scoresList", "ModelList", function(object) object@scoresList)
setGeneric("scoresList<-", function(object,value) standardGeneric("scoresList<-"))
setReplaceMethod("scoresList", signature(object="ModelList",
                                     value="list"),
                 function(object, value) {
                   object@scoresList <- value
                   return(object)
                 })
setGeneric("intAnd", function(object) standardGeneric("intAnd"))
setMethod("intAnd", "ModelList", function(object) object@intAnd)
setGeneric("intAnd<-", function(object,value) standardGeneric("intAnd<-"))
setReplaceMethod("intAnd", signature(object="ModelList",
                                         value="logical"),
                 function(object, value) {
                   object@intAnd <- value
                   return(object)
                 })
setGeneric("intgAnd", function(object) standardGeneric("intgAnd"))
setMethod("intgAnd", "ModelList", function(object) object@intgAnd)
setGeneric("intgAnd<-", function(object,value) standardGeneric("intgAnd<-"))
setReplaceMethod("intgAnd", signature(object="ModelList",
                                     value="logical"),
                 function(object, value) {
                   object@intgAnd <- value
                   return(object)
                 })