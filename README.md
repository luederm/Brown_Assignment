To run pipeline run the following command:
luigi --module coverage RunAll --r1 <read 1> --r2 <read 2> --x <Reference fasta> --t <number of threads> --pandoc <Path to pandoc> [--Q <Quality cutoff>] [--local-scheduler]

--r1 and --r2: the paths to the paired end reads (in fastq format w/ phred33 encoding) from the MiSeq sequencer
Note: This pipeline does not perform adaptor or barcode trimming. This should be done ahead of time.

--x: Path to the reference genome (should be in fasta format)

--t: Number of threads (only certain portions of pipeline multi-threaded)

--pandoc: Path to pandoc. You can find this by opening up Rstudio and typing: Sys.getenv("RSTUDIO_PANDOC")

--Q: [optional] Filter out reads with mean quality score less than this argument. Should be integer.

--local-scheduler: [optional] Use this argument unless you plan to use a luigi central scheduler

