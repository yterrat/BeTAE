# BeTAE
A local clone of MG-RAST 

First and foremost, this is not a tool that outperforms everything else. 
These few hundred lines of code have been produced because we had literally thousand of metgaenomic samples to process.
I know that it sounds stupid to "reinvente the wheel", but this time, it was more than necessary. MG-RAST has been a gold standard of metagenomic annotation and is extremly usefull to annotate a few samples. Unfortunatly, it does not scale up well with thousand of samples and you'll probably get stucked forever in the queue.
So, BeTAE is "just" a local clone of MG-RAST. It has been designed to run on high computing services and can process a single sample of metagenomic raw reads ~ 30GB in ~4 hours (including cleaning, processing, annotation and assembly).
