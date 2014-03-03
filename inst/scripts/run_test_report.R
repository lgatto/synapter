library("synapter")
library("synapterdata")

fasfile <- synapterdata::getFasta()
mstrfile <- synapterdata::getMaster()
pepfile <- synapterdata::getMSeFinalPeptide()[2]
pep3dfile <- synapterdata::getMSePep3D()[2]

inlist <- list(quantpeptide = pepfile,
               quantpep3d = pep3dfile,
               identpeptide = mstrfile,
               fasta = fasfile)    

syndat <- Synapter(inlist, master = TRUE)


(v <- packageVersion("synapter"))
outdir <- paste0("Report_", v)
if (!file.exists(outdir)) dir.create(outdir)

set.seed(1L)
synres <-
    synergise(object = syndat, master = TRUE,
              outputdir = outdir, 
              span = 0.05,
              grid.subset = 0.5)

htmlReport <- paste0("file:///",
                     file.path(getwd(), outdir, "index.html"))
browseURL(htmlReport)


load("Report_1.5.2/SynapterObject.rda", verbose = TRUE)
obj2 <- obj$copy()
load("Report_1.5.3/SynapterObject.rda", verbose = TRUE)
obj3 <- obj$copy()

xx <- ls(obj2@.xData)

eql <- sapply(xx, function(x)
              all.equal(get(x, obj2@.xData), get(x, obj3@.xData)))

eql[eql != "TRUE"]
##         DbFastaFile         SynapterLog             Version 
## "1 string mismatch" "1 string mismatch" "1 string mismatch" 
##         OK                                          OK

all.equal(obj3@.xData$SynapterLog[-1], obj2@.xData$SynapterLog[-1])
