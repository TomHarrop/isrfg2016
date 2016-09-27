## Title slide

* resource for studying the basic molecular mechanisms involved in panicle branching

## Inflorescence architecture

* rice inflorescence has a main axis called the rachis
* the primary branches emerge from nodes on the rachis
* secondary branches similarly attached to primary branches
* spikelets are short branches on the primary and secondary branches
* each spikelet bears one grain and spikelets are the individual floral units of the rice inflorescence

## Branching component of yield

* phenotyping data from around 20 accessions of african and asian wild and domesticated rice
* primary and secondary branches correlate to the number of spikelets produced
* illustrates that branching is a component of grain yield

## Inflorescence development

* Branching complexity is specified within the first 8–10 days of panicle development
* After induction of flowering the shoot apical meristem (SAM) is converted to a rachis meristem (RM)
* Next primary branch meristems (PBM) are established peripheral to the RM
* After formation of primary branches primary branches elongate (ePBM)
* During elongation, the apex retains meristem activity and can give rise to axillary meristems (AM)
* AM may differentiate into secondary and higher-order branches or be converted directly to spikelet meristem
* Spikelet meristem is the determinate phase of development
* After spikelet differentiation complexity is fixed, and the branches elongate flowers develop
* Two stages where phase transition of meristems conditions branching:
	* Timing of rachis meristem abortion determines the number of primary branches.
    * Transition of indeterminate branch meristems to determinate spikelet meristems specifies the complexity of branching
	* delays in transitions may lead to reiteration of branching (more yield potential)
* APO1/APO2, UFO/LFY homologs, maintain BM identity and suppress SM identity. Increase in APO1 expression leads to prolonged indeterminate meristem fate and delays in SM specification
* TAW1 promotes BM and suppresses SM. in o/e mutants BM formation is reiterated. in KD spikelets form precociously.
* PAP2 (MADS32, SEP) positively regulates SM identity

## Molecular control of meristem identity

* There are several examples of genes and interactions that have been studied individually and been found to control meristem identity transitions
* the landscape of gene expression is different between meristem types
* no mechanistic understanding of dynamic gene–gene interactions that control these patterns
* do the members of this transcriptional network control branching and yield, and how were they affected by domestication?
* to answer this broad question we started by describing gene expression in different meristem types

## Laser microdissection (LMD)

* Used O. sativa japonica cv. Nipponbare
* Prepared and dissected RM, PBM ePBM/AM and SM samples

## LMD results

* Dissected around 1 square mm of 8 um sections
* Used three true biological replicates at each stage: multiple panicles from different individuals were used for each sample
* Recovered an average of less than 20 ng of RNA per sample
* RNA integrity number higher than 7 for each sample. Good results for LMD.
* Used Ovation RNA-Seq System V2 to amplify RNA and produce libraries
* Sequenced around 70 M reads per library, resulting in an average of more than 20 M reads within annotated genes

## LMD Dataset

* Resource is publicly available
* Look up your favourite genes if you are interested in panicle development
* Find out more information about genes e.g. from GWAS, maybe narrow down your list of candidates
* More complicated analysis: 
	* Correlate coexpression patterns with TF binding sites
	* Try to model a gene regulatory network for meristem identity!
* Easiest way to find it is to go to EvoRepRice github page:
* Link to raw data
* Link to paper with some supplementary data
* Can also email me or the corresponding authors

## Coexpressed genes

* As an overview of the data:
	* fuzzy c-means clustering of transformed read count of detected genes
	* switch between apical and axillary meristems, followed by gradual changes during transitions between PBM and SM (corresponding to transition of axillary meristem from indeterminate to determinate state)

## O. sativa G1-Like genes

* As an example of how the data can be used:
* Looked more closely in the clusters and found three G1-Like genes
* G1-Like is a small family of transcriptional activators
* G1L5 is overexpressed in the TAWAWA1 mutant, causing increased secondary branch and spikelet production

## G1L5 targes

* G1L5 activates two SVP-like MADS transcription factors, which may promote branch meristem identity (delay SM determination)
* these two genes were also coexpressed with G1L1, G1L2 and G1L5
* interesting: TAW1 / SVP genes are supposed to promote BM identity but all had an expression peak in the RM
* potential additional role in apical meristem

## Oryza G1-Like genes

* also measured expression of G1-Like family in **whole panicles** from African and Asian wild and domesticated rice
* Along with G1L1, G1L2 and G1L5, we detected G1L3 and G1L4 at PBM and SM stage
* These genes were not detected in the LMD dataset, suggesting expression outside the meristem

## More about G1-Like genes

* Israr's poster

## Information for GWAS candidates

* to finish: using the dataset to investigate candidate gene lists
* Primary branch number QTL from a large GWAS on 242 tropical sativa accessions
* Out of 25 genes inside the region, 10 were detected in the LMD dataset
* mostly uncharacterised genes
* six were differentially expressed with respect to stage, and several were expressed most strongly in the RM, which could be consistent with a role in primary branch formation

## LMD dataset

* Given examples of using the LMD dataset for analysing genes of interest and working with candidate gene lists
* Lots of other things you could do with this data
* These landscapes of expression are basically the final outputs of complex interactions between many TFs that coordinate the expression of each gene
* What we've done is measured the "steady states" of gene expression in different meristem types
* Would be v. interesting to model the interactions in the regulatory components that produce these steady states
* Obviously with a view to testing the network e.g. with mutants






