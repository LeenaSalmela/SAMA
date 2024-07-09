This is the original Detox readme. MAGA code is based on the Detox code so we include it here for reference.

# Detox
Accurate determination of node and arc multiplicities in a de Bruijn graph using Conditional Random Fields

## Prerequisites
In order to build Detox you need to the following software or library packages:

  * CMake version 2.6.3 or higher
  * GCC version 4.7 or higher
  * Google's [sparsehash](https://github.com/sparsehash/sparsehash)
  * ZLIB (optional)
  * recent [Boost](https://www.boost.org/) libraries (for approximate inference)
  * GMP library (for approximate inference)

For approximate inference computations Detox relies on an adaptation of the Loopy belief propagation implementation of the [libDAI](https://bitbucket.org/jorism/libdai/src/master/) C++ libary. 
The classes from libDAI that are necessary for Loopy Belief Propagation with our own additions of more recent message passing schemes are provided in this repository.

In order to run Detox you also need to install and compile [BCALM 2](https://github.com/GATB/bcalm).

## Compilation

Clone the Detox Github repository:

```bash
git clone https://github.com/biointec/detox
```

Next, compile the C++ code as follows:

```bash
cd detox
mkdir build
cd build
cmake ..
make
```

After compilation, the `detox` binary can be found in:

```bash
detox/build/src/detox
```

## Usage

We use BCALM 2 to build a compacted de Bruijn graph from sequencing data. Next, Detox is used to infer its the node and arc multiplicities.
### 1. Generating the de Bruijn graph using BCALM 2
The first step is to prepare a manifest file, e.g. named `reads.mf`, in which the input FASTQ files are listed. For example, suppose the reads are stored in two FASTQ files named `reads_1.fastq` and `reads_2.fastq`, then the `reads.mf` file would simply list these input files, one file per line:

```
reads_1.fastq
reads_2.fastq
```
The set of reads must be in FASTQ (`.fq`/`.fastq`) format.

Using BCALM 2, one can then create a de Bruijn graph as follows:

```bash
bcalm -in reads.mf -kmer-size 21 -abundance-min 1
```

Depending on the size of the genome, and the number of CPU cores of your machine, this process may take minutes (for bacterial genomes) or many hours (for human genome scale genomes). The command line option `-abundance-min 1` instructs BCALM 2 to only retain all k-mers that occurin the input reads. Setting this parameter to a larger value will remove a large fraction of the sequencing errors already, but we observe better performance of our own de Bruijn graph cleaning procedure when not preprocessing the data. Remark that BCALM 2's default value for `-abundance-min` is 2.

When BCALM 2 is finished, it produces the file `reads.unitigs.fa`, in which the de Bruijn graph's unitig node sequences as well as their connectivity are encoded in a FASTA-format. This file serves as input for Detox.

### 2. Inferring the node/arc multiplicities using Detox

The pipeline takes as input a de Bruijn graph in the FASTA-format of BCALM 2 (`reads.unitigs.fa`) as well as the original read set (`reads.mf`). Detox can then be run as follows:

```bash
detox -use-qual -abundance-min 1 -crf-nb-size 5 -approx-inf reads.unitigs.fa reads.mf
```

The two mandatory input arguments are `reads.unitigs.fa` and `reads.mf`, in that order. All other flags or command line options should be listed prior to these two input arguments.

The flag `-use-qual` instructs Detox to weigh the k-mer counts using the Phred quality scores. When BCALM 2 was run using the `-abundance-min`  option, it is best to provide this information also to Detox. When Detox is fitting the error model to the data, it then knows that all k-mers that occur fewer times than the value provided by `abundance-min` are missing. It takes this into account when fitting the model to the data. If you use Detox with unweighted k-mer counts, the value of `-abundance-min` can be automatically inferred.

The `-crf-nb-size` option specifies the size of the neighborhood subgraph around the node/arc whose multiplicity is inferred, whenever detox uses exact inference computations. For smaller, bacterial genomes, a value of 5 can be chosen in order to get the highest possible accuracy. For human genome-scale genomes, a value of 1 is recommended to prevent very long runtimes. In most cases, using a small neighborhood size is enough, however, if regions of extremely high or low coverage span large regions of the de Bruijn graph, it might take a large neighborhood size to obtain a CRF that can infer the correct multiplicities.

Finally, the `-approx-inf` option specifies the use of approximate inference for final multiplicity estimate computations. When this option is left out, final inference computations will use subgraphs with the size specified by the `-crf-nb-size` option.

Upon completion, the estimated node and arc multiplicities are stored in files `estmult.node` and `estmult.edge`, respectively. The file `estmult.node` contains the following fields (tab-separated):

* The node identifier (integer value from 1 to #nodes)
* The estimated multiplicity m (non-negative integer value)
* The log-probability _log[ P(mult. is m) ]_
* The average node coverage (weighed by the Phred quality use if the `-use-qual` flag was used)
* The length of the node expressed as the number of k-mers.

The file `estmult.edge` contains the following fields (tab-separated):

* The node identifier of the node from which the arc is originating
* The node identifier of the destination node
* The estimated multiplicity m (non-negative integer value)
* The log-probability _log[ P(mult. is m) ]_
* The coverage of the arc

Note that the node identifiers in the `estmult.*` files are equal to (the BCALM node identifier + 1). In other words, we number nodes 1, 2, 3, ... whereas BCALM start numbering from 0. A negative node identifier refers to the reverse-complement node.

### 3. Advanced user information 
The Detox pipeline consists of 4 stages. Each stage outputs a number of intermediate files. When re-running Detox (e.g. using different parameter settings), the existence of these intermediate files is checked. If they exist, stages 1, 2 or 3 may be skipped. In order to force re-running stages 1, 2 or 3, simply remove the appropriate files: `rm *.st1` for stage 1,  `rm *.st2` for stage 2 or `rm *.st3` for stage 3. We will now explain the stages in more detail below.

#### Stage 1
The input file `reads.unitigs.fa`, produced by BCALM 2 is read. Even though BCALM 2 provides average k-mer counts for each node (unitig), it does not provide (k+1)-mer counts for the arcs. Also, when using the `-use-qual` flag, Detox needs to weigh the k-mer counts using the Phred quality scores. Note that the Phred score ASCII base value can be set using `-phred-base <value>` (default=33).

So in either case, after reading the de Bruijn graph into memory, Detox needs to stream through all reads in order to establish both the appropriate node and arc counts. The results are written to `reads.unitigs.fa.st1`, a binary file containing the de Bruijn graph annotated with q- or k-mer coverage. This file is a prerequisite for stage 2.

If a file `genome.fasta` is present in current working directory, Detox will assume this contains the reference genome and it will also compute the true node and edge multiplicities. These will be stored in files `truemult.node.st1` and `truemult.edge.st1`, respectively.

In principle, one should only re-run stage 1 when switching between k-mer or q-mer counts.

#### Stage 2
In stage 2, a k-mer and (k+1)-mer mixture model is fitted to the data. By default, Detox selects a number of random nodes and random arcs from the de Bruijn graph on which the model is then trained. The number of nodes and arcs can be specified using the `-em-train-size` (default = 10.000) option. Larger values lead to larger runtimes, but model parameters will be estimated more accurately.

Both models are fitted using the expectation-maximization (EM) algorithm. During the E-step, the multiplicity is inferred for selected nodes and arcs based on the current model estimation. During the M-step, the model parameters are updated based on the estimated multiplicities. This process continues until either the model has converged (the relative convergence tolerance can be set using `-em-conv-eps` (default = 0.001) or until the maximum number of EM iterations is reached (this value can be set using the `-em-max-iter` (default = 25) option. Note that detox always uses exact inference computations on subgraphs in stage 2.

Stage 2 generates following files:

  * `model.node.st2` and `model.edge.st2` that contain the parameters of the fitted models to the k-mer or q-mer data for the nodes and arcs respectively. These files are prerequisites for stage 3.
  * `histogram.node.st2.gnuplot`, `histogram.node.st2.dat`, `histogram.edge.st2.gnuplot` and `histogram.edge.st2.dat`. These are auxiliary files that can be used to plot the model fit to the histogram data. It can be used to manually check if the fit is good.  You will need the `gnuplot` for that. Simply run `gnuplot histogram.node.st2.gnuplot` and `gnuplot histogram.edge.st2.gnuplot` to generate the files `histogram.node.st2spectrum.pdf` and `histogram.edge.st2spectrum.pdf`, respectively.

The following command line options are relevant to set the **initial** parameter estimates to the q-mer or k-mer models. It is generally not necessary to provide them. Only when the EM algorithm does not converge to the correct solution, manually setting these these initial estimates might help.

* `-mm-coverage` (default = auto) to provide the average coverage of the non-repeated nodes/arcs of the dataset.
* `-mm-err-cov` (default = 1) to provide the average coverage of erroneous nodes or arcs. 
* `-mm-odf` (default = 1.5) to provide the overdispersion factor of the Negative Binomial distribution. This is the ratio of the variance and the mean.

The number of mixture model components is controlled using the `-mm-components` option. By default, we use 6 components: one to model the erroneous nodes (resp. arcs) and 5 components to model multiplicities 1 to 5.

It is important to note that during the E-step, the multiplicities of nodes and arcs are inferred using the CRF methodology with exact inference on subgraphs. Therefore, the `-crf-nb-size` option to specify the neighborhood size is relevant. In general, a value larger than 0 (i.e., using the CRF methodology) reduces the required number of EM iterations and proves more robust against poor choices of initial model parameter values.

One should only re-run stage 2 to re-train the mixture models. It suffices to remove model files: `rm *.st2`.

#### Stage 3

In stage 3, spurious nodes and arcs are removed from the de Bruijn graph. This is done iteratively. Based on the coverage model computed in stage 2, a maximum coverage value is selected below which a node or arc is considered as a low coverage node/arc `Cl`. Low coverage nodes and arcs are considered for possible removal. 
In a first pass, low coverage nodes/arcs are only removed when the conservation of flow property holds in all nodes in a subgraph around that node/arc (subgraph size set by `-crf-nb-size`). We split this pass up into several iterations: first only nodes/arcs with coverage below `1/4 * Cl` are considered, next nodes/arcs with a coverage below `2/4 Cl` etc. until `Cl`. After each iteration selected nodes/arcs are removed and the de Bruijn graph is contracted further where possible. In a second pass, the CRF of a subgraph is constructed for each (remaining) low coverage node/arc and the low coverage node/arc is only removed when its CRF-based multiplicity estimate is 0. This pass is again split up in sevaral iterations and selected nodes/arcs are removed and the de Bruijn graph contracted after each iteration.

After graph cleaning, the coverage model is re-estimated under a fixed 0-multiplicity distribution with low weight.  

Stage 3 generates following files:

  * `model.node.st3` and `model.edge.st3` that contain the parameters of the newly fitted models to the k-mer or q-mer data for the nodes and arcs respectively. These files are prerequisites for stage 4.
  * `histogram.node.st3.gnuplot`, `histogram.node.st3.dat`, `histogram.edge.st3.gnuplot` and `histogram.edge.st3.dat`.
  * `reads.unitigs.fa.st3` a binary file containing the new cleaned and compacted de Bruijn graph and node/arc coverage.
  * if `genome.fasta` is present `truemult.node` and `truemult.edge` contain the true multiplicities of the nodes/arcs remaining in the de Bruijn graph.

#### Stage 4
By default, in stage 4, the node and arc multiplicities are computed for **all** nodes and arcs in the cleaned de Bruijn graph.

The most accurate multiplicity estimates can be obtained in feasible runtime when approximate inference is used to estimate multiplicities based on one CRF for the complete de Bruijn graph (use option `-approx-inf`).
If you want to use exact inference on subgraphs, the most important command line option is `-crf-nb-size` to specify the neighborhood size. It is possible use a different neighborhood size in stage 2 (to train to the models) and stage 3 (to infer the multiplicities).

Stae 4 generates the following files: 
* `estmult.node` and `estmult.edge` contain the estimated multiplicities.
* `truemult.node` and `truemult.edge` contain the true multplicities of the cleaned and compacted de Bruijn graph, when `genome.fasta` is present.
* `Cytograph.full.nodes` and `Cytograph.full.arcs` can be used to visualise the de Bruijn graph in Cytoscape (see further).

#### Other parameters that are relevant for CRF inference:

  * `-crf-margin` (default = 2) specifies the number of alternative multiplicities (one-sided) used within the probabilistic factors. For example, for the default value of 2, if the estimated multiplicity of a node is 3 based on its coverage, then the probabilities will be computed for multiplicities {1,2,3,4,5}. If the estimated multiplicity of a node is 0, then probabilities will be computed for multiplicities {0,1,2}.
  * `-crf-flow` (default 1e7) defines the probability ratio for the scenario in which the conservation of flow holds and the scenario in which it does not hold. The higher this value, the stronger the CRF model will impose the conservation of flow rule. If you have a dataset with a very high coverage (much larger than _50x_) you might set this parameter to a higher value, as multiplicity distributions that are far apart require a stronger conservation of flow rule.
  * `-crf-max-fact` (default 1e6) defines the maximum factor size that is allowed during the variable elimination algorithm. If this threshold is exceeded, the algorithm with reduce the `-crf-nb-size` (only for the node/arc for which this threshold was exceeded).
 
 #### Changing the message passing scheme used for Loopy Belief Propagation approximate inference
 
We have tested which of the many parameters for Loopy Belief Propagation (LBP) in libDAI achieved the best results in terms of convergence and speed and set default values accordingly. If you do want to change these parameters, this can be done through a file called `libdai.props` (detox auto-detects a file with this name). 
The file has the following tab-separated format (example will just set the default values):
```
maxiter 500
tol 1e-3
logdomain 1
maxtime 1800
updates SEQMAX0L
weightdecay 1
resinit MESSAGE
damping 0.0
```

* `maxiter` (int) sets the maximum number of LBP iterations
* `tol` (double) sets the tolerance on how much the messages can change between iterations to assume convergence
* `logdomain` (0 or 1) perform message passing computations in the log-domain (advised)
* `maxtime` (int) number of seconds that LBP is run before breaking off computations if convergence is not yet achieved
* `updates` (string) the message passing scheme used (`PARALL`: parallel updates (not advised), `SEQFIX`: fixed order sequential (not advised), `SEQRND`: random order sequential (not advised), `SEQMAX`: maximum residual belief updates, `SEQMAX0L`: maximum residual belief updates with lookahead 0, `SPLASH`: splash belief propagation that can be run multi-threaded (not fully supported yet))
* `weightdecay` (0 or 1) use weight decay
* `resinit` (string) residual initialisation method (only in combination with `SEQMAX0L`) `MESSAGE`: initial pass to compute messages, `UNIFORM` residual upper bound of Sutton and McCallum (2012).
* `damping` (double between 0.0 and 1.0) use dampend belief updates (0.0: no damping)

From our tests Maximum Residual Updates with Lookahead Zero and weight decay proved the best message passing scheme in terms of convergence and speed. 

 ### 4. Other options

 * `-num-threads` to specify the number of threads to use, when using exact inference on subgraphs. By default, the number of threads equals the number of CPU cores or twice that number for CPUs that support hyperthreading.
 * `-help` to display the help page.
 * `-version` to display the version number.

 ### 5. Checking the fit of the model produced in stage 2

If you want to check the fit to the histogram of the model that Detox fits in stage 2 and 3, this can be done visually using [Gnuplot](http://www.gnuplot.info/).
Intermediate files `histogram.node.stx.gnuplot`, `histogram.node.stx.dat`, `histogram.edge.stx.gnuplot` and `histogram.edge.stx.dat` can be used to plot the histogram overlaid with the estimated model. Just run
```
gnuplot histogram.node.stx.gnuplot
```
This will produce a pdf output of the plot.
Detox auto-determines a good scope for the plot, but you can manually change this in the `.gnuplot` files if need be.

#### What to do when the fit seems bad

We have determined what we think are the best default settings for the parameters Detox uses, based on tests on several datasets. However, if your visual inspection of the model fit seems bad, you can try to manually alter some of these parameter values:

##### Check if convergence was obtained

By default, the maximum number of EM-iterations is 25. You will be notified by Detox when model fitting finished with this many iterations. If you see that there was still a relatively large change in parameter values between iteration 24 and 25, it can be beneficial to rerun stage 2 with a higher number of EM-iterations. Option `-em-max-iter [iters]` sets the maximum number of EM-iterations.

##### Adjusting the initial estimates for the parameters of the error distribution and of the multiplicity 1 distribution for the fit in stage 2

If you see that peaks of the estimated distributions do not correspond with the correct peaks in the histogram you can initialize the means of the error- and the multiplicity 1 distribution yourself (`-mm-coverage`, resp. `-mm-err-cov`). Another option is to initialize the multiplicity 1 mean to the average read coverage of your dataset (if this is known to you). If the fit of the distributions seems too wide or too tight, try a different initialization of the overdispersion factor (`-mm-odf`).

 ### 6. Visualizing the de Bruijn graph

You can visualize the final de Bruijn graph using [Cytoscape](https://cytoscape.org/)

The graph is provided as two files `Cytograph.full.nodes` and `Cytograph.full.arcs`. Within Cytoscape you first import the `.arcs` file with `File > Import > Network from file...` and then you import the `.nodes` file with `File > Import > Table from file...`

Additionally, you can rerun Detox with a `-vis-subgraph n` option to export a subgraph of the de Bruijn graph around node `n` (size of the subgraph is specified with the `-crf-nb-size` option. For more information about how to obtain these Cytoscape files and how to visualise them, see the [visualisation folder](visualisation/). There you will also find a small example.
