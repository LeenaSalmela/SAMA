# MAGA

Misassembly Avoidance Guaranteed Assembler

## Prerequisites
In order to build MAGA you need to the following software or library packages:

  * CMake version 2.6.3 or higher
  * GCC version 4.7 or higher
  * Google's [sparsehash](https://github.com/sparsehash/sparsehash)
  * ZLIB (optional)
  * recent [Boost](https://www.boost.org/) libraries (for approximate inference)
  * GMP library (for approximate inference)

For approximate inference computations the Detox method included in MAGA relies on an adaptation of the Loopy belief propagation implementation of the [libDAI](https://bitbucket.org/jorism/libdai/src/master/) C++ libary. 
The classes from libDAI that are necessary for Loopy Belief Propagation with our own additions of more recent message passing schemes are provided in this repository.

In order to run MAGA you also need to install and compile [BCALM 2](https://github.com/GATB/bcalm).

## Compilation

Clone the MAGA Github repository:

```bash
git clone https://github.com/LeenaSalmela/MAGA
```

Next, compile the C++ code as follows:

```bash
cd MAGA
mkdir build
cd build
cmake ..
make
```

After compilation, the `maga` binary can be found in:

```bash
MAGA/build/src/maga
```

## Usage

We use BCALM 2 to build a compacted de Bruijn graph from sequencing data. After that MAGA is used to produce contigs with a guaranteed correctness.

### 1. Generating the de Bruijn graph using BCALM 2
The first step is to prepare a manifest file, e.g. named `reads.mf`, in which the input FASTQ files are listed. For example, suppose the reads are stored in two FASTQ files named `reads_1.fastq` and `reads_2.fastq`, then the `reads.mf` file would simply list these input files, one file per line:

```
reads_1.fastq
reads_2.fastq
```
The set of reads must be in FASTQ (`.fq`/`.fastq`) format.

Using BCALM 2, one can then create a de Bruijn graph as follows:

```bash
bcalm -in reads.mf -kmer-size 31 -abundance-min 2
```

For performance reasons we recommend to remove some sequencing errors already here by only including k-mers with abundance at least 2 (the -abundance-min 2 parameter).

When BCALM 2 is finished, it produces the file `reads.unitigs.fa`, in which the de Bruijn graph's unitig node sequences as well as their connectivity are encoded in a FASTA-format. This file serves as input for MAGA.

### 2. Generating contigs with correctness guarantee

MAGA uses the Detox method to estimate k-mer and k+1-mer distribution in the genome. Then it computes the abundancy thresholds and generates contigs.

The pipeline takes as input a de Bruijn graph in the FASTA-format of BCALM 2 (`reads.unitigs.fa`) as well as the original read set (`reads.mf`). MAGA can then be run as follows:

```bash
maga -abundance-min 2 -misassembly-likelihood 1e-03 reads.unitigs.fa reads.mf
```

The two mandatory input arguments are `reads.unitigs.fa` and `reads.mf`, in that order. All other flags or command line options should be listed prior to these two input arguments.

When BCALM 2 was run using the `-abundance-min`  option, it is best to provide this information also to MAGA. When MAGA uses the Detox method to fit the error model to the data, it then knows that all k-mers that occur fewer times than the value provided by `abundance-min` are missing. It takes this into account when fitting the model to the data.

Upon completion, the contigs are stored in file `output.fa`.

### 3. Advanced user information 
The MAGA pipeline consists of 3 stages. Each stage outputs a number of intermediate files. When re-running MAGA (e.g. using different parameter settings), the existence of these intermediate files is checked. If they exist, stages 1 or 2 may be skipped. In order to force re-running stages 1, 2 or 3, simply remove the appropriate files: `rm *.st1` for stage 1 or  `rm *.st2` for stage 2. We will now explain the stages in more detail below.

#### Stage 1 (part of the Detox method)
The input file `reads.unitigs.fa`, produced by BCALM 2 is read. Even though BCALM 2 provides average k-mer counts for each node (unitig), it does not provide (k+1)-mer counts for the arcs.

So after reading the de Bruijn graph into memory, Detox needs to stream through all reads in order to establish both the appropriate node and arc counts. The results are written to `reads.unitigs.fa.st1`, a binary file containing the de Bruijn graph annotated with k-mer coverage. This file is a prerequisite for stage 2.

In principle, one should only re-run stage 1 when switching between different values of k.

#### Stage 2 (part of the Detox method)
In stage 2, a k-mer and (k+1)-mer mixture model is fitted to the data. By default, Detox selects a number of random nodes and random arcs from the de Bruijn graph on which the model is then trained. The number of nodes and arcs can be specified using the `-em-train-size` (default = 10.000) option. Larger values lead to larger runtimes, but model parameters will be estimated more accurately.

Both models are fitted using the expectation-maximization (EM) algorithm. During the E-step, the multiplicity is inferred for selected nodes and arcs based on the current model estimation. During the M-step, the model parameters are updated based on the estimated multiplicities. This process continues until either the model has converged (the relative convergence tolerance can be set using `-em-conv-eps` (default = 0.001) or until the maximum number of EM iterations is reached (this value can be set using the `-em-max-iter` (default = 25) option. Note that Detox always uses exact inference computations on subgraphs in stage 2.

Stage 2 generates following files:

  * `model.node.st2` and `model.edge.st2` that contain the parameters of the fitted models to the k-mer or q-mer data for the nodes and arcs respectively. These files are prerequisites for stage 3.
  * `histogram.node.st2.gnuplot`, `histogram.node.st2.dat`, `histogram.edge.st2.gnuplot` and `histogram.edge.st2.dat`. These are auxiliary files that can be used to plot the model fit to the histogram data. It can be used to manually check if the fit is good.  You will need the `gnuplot` for that. Simply run `gnuplot histogram.node.st2.gnuplot` and `gnuplot histogram.edge.st2.gnuplot` to generate the files `histogram.node.st2spectrum.pdf` and `histogram.edge.st2spectrum.pdf`, respectively.

The following command line options are relevant to set the **initial** parameter estimates to the k-mer model. It is generally not necessary to provide them. Only when the EM algorithm does not converge to the correct solution, manually setting these these initial estimates might help.

* `-mm-coverage` (default = auto) to provide the average coverage of the non-repeated nodes/arcs of the dataset.
* `-mm-err-cov` (default = 1) to provide the average coverage of erroneous nodes or arcs. 
* `-mm-odf` (default = 1.5) to provide the overdispersion factor of the Negative Binomial distribution. This is the ratio of the variance and the mean.

The number of mixture model components is controlled using the `-mm-components` option. By default, we use 6 components: one to model the erroneous nodes (resp. arcs) and 5 components to model multiplicities 1 to 5.

It is important to note that during the E-step, the multiplicities of nodes and arcs are inferred using the CRF methodology with exact inference on subgraphs. Therefore, the `-crf-nb-size` option to specify the neighborhood size is relevant. In general, a value larger than 0 (i.e., using the CRF methodology) reduces the required number of EM iterations and proves more robust against poor choices of initial model parameter values.

One should only re-run stage 2 to re-train the mixture models. It suffices to remove model files: `rm *.st2`.

See the Detox readme for more details on the parameters relevant for the Detox method.

#### 3. Stage 3: Estimating abundancy thresholds and generating contigs

In stage 3 the abundancy thresholds for including edges in the assembly are computed and then the contigs are generated according to these thresholds.

The thresholds for extension and their corresponding misassembly probabilities are stored in file `thresholds.list` and the final contigs in file `output.fa`.

 ### 4. Other options

 * `-num-threads` to specify the number of threads to use, when using exact inference on subgraphs. By default, the number of threads equals the number of CPU cores or twice that number for CPUs that support hyperthreading.
 * `-help` to display the help page.
 * `-version` to display the version number.
