# BeTAE
A local clone of MG-RAST 

First and foremost, this is not a tool that outperforms everything else. 
These few hundred lines of code have been produced because we had literally thousand of metgaenomic samples to process.
I know that it sounds stupid to "reinvente the wheel", but this time, it was more than necessary. MG-RAST has been a gold standard of metagenomic annotation and is extremly usefull to annotate a few samples. Unfortunatly, it does not scale up well with thousand of samples and you'll probably get stucked forever in the queue.

So, BeTAE is "just" a local clone of MG-RAST. It has been designed to run on high computing services and can process a single sample of metagenomic raw reads ~ 30GB in ~4 hours (including cleaning, processing, annotation and assembly).

Here are the different steps (copy of material and methods used for publication):

" Paired-end raw reads of 150 base pairs (bp) were further analyzed using a home-made bioinformatics pipeline. Firstly, quality trimming of raw reads was performed by the SolexaQA v3.1.7.1 program with default settings [66]. Trimmed reads shorter than 75 nt were removed for further analysis. Artificial duplicate removal was performed using an in-house script based on the screening of identical leading 20 bp. From the trimmed high-quality reads, gene fragments were predicted using FragGeneScan-Plus v3.0 [67]. Cd-hit v4.8.1 was applied to cluster predicted protein fragments at a 90% similarity [68]. One representative of each cluster was used for a similarity search on the M5nr database (https://github.com/MG-RAST/myM5NR) using the Diamond engine [69]. For assessment of taxonomic affiliation of gene fragments encoding proteins, we took into account best hits (minimal e-value of 1 × 10−5) combined with a last common ancestor approach."

If you read the last word in details, you know that you need to install third party programs :

  SolexaQA
  FragGeneScan-Plus
  Cd-Hit
  Diamond
  M5nr database
  
  and assembly tools like MegaHit, SPAdes, ...
