HogWild++ Experiment Code
===============================================================================

HogWild++ aims to improve the scalability of HogWild! stochastic gradient
descent algorithm on multi-socket (NUMA) machines. In HogWild++, to reduce
inter-socket communication and coherence miss, worker threads are grouped into
several "clusters" with configurable size, and in each cluster threads share
the same model vector.  Each cluster only exchanges model information with its
neighbour(s) to reduce communication overhead, and no centralized model vector
is maintained.  For more details about this algorithm please refer to the
following paper:

```
HogWild++: A New Mechanism for Decentralized Asynchronous Stochastic Gradient Descent
Huan Zhang, Cho-Jui Hsieh and Venkatesh Akella
```

In our implementation, clusters are organized in a logical directional ring,
and there is a single token passing along the ring. The cluster holding the
token is able to communicate with the next cluster on the ring and exchange
model information, and other clusters keep updating their own models. If you
are interested in implementing other (potentially more efficient) topologies
please read the `How to modify` section below.

HogWild++ is based on the HogWild! v03a code, available [here](http://i.stanford.edu/hazy/victor/Hogwild/).

How to build
----------------------

HogWild++ requires `libnuma`. On Debian based systems you can install `libnuma`
using the following command:

```
sudo apt-get install libnuma-dev
```

Then run `make` to build, and the following binaries will be built in `bin` folder:

* svm: Basic SVM example in HogWild!, unchanged, and it is used as the
  performance baseline.

* numasvm: Implementataion of SVM using the HogWild++ algorithm.

* convert: Convert TSV files into binary files. Assumes the rows and columns in
  the TSV file are indexed starting at 0.

* convert_matlab: Converts TSV files into binary files. Assumes the rows and
  columns are indexed starting at 1.

* unconvert: Converts a binary file into a TSV file. The TSV file will be
  indexed starting at 0.

Data Preparation
----------------------

Datasets can be downloaded from [LIBSVM website](http://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary.html)
however you need to convert them to the TSV (Tab Separated Value) format that HogWild! uses.
We provide a script `covert2hogwild.py` to convert LIBSVM format to TSV format.
To reduce data loading time, you should also convert TSV to binary format.

The following commands show how to download and prepare the RCV1 dataset.
Note that for RCV1 we swapped the downloaded training and test set because the "test set"
is actually larger.
```
mkdir data && cd data
# prepare the training set
wget https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/rcv1_test.binary.bz2
bunzip2 rcv1_test.binary.bz2
python ../convert2hogwild.py rcv1_test.binary rcv1_train.tsv
../bin/convert rcv1_train.tsv rcv1_train.bin
# prepare the test set
wget https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/rcv1_train.binary.bz2
bunzip2 rcv1_train.binary.bz2
python ../convert2hogwild.py rcv1_train.binary rcv1_test.tsv
../bin/convert rcv1_test.tsv rcv1_test.bin
# Done
cd ..

```


We have prepared binary files used in the experiments of our paper.
You can download these datasets here:

[http://jaina.cs.ucdavis.edu/datasets/classification_compressed/](http://jaina.cs.ucdavis.edu/datasets/classification_compressed/)

You only need to download .bin.xz files. To save downloading time these files
are compressed. Please decompress them using the `xz` utility before use.


How to run
----------------------

HogWild++ adds three new parameters to HogWild! executable. 

* `cluster_size`: the number of threads in each cluster (referred as `c` in our
  paper). Threads within the same cluster share the same model vector. When the
  cluster size is the same as the total number of threads, HogWild++ is
  equvalent to HogWild!. The best cluster size depends on dataset
  characteristics. As a starting point, you can set the cluster size as the
  total number of physical cores on a single socket.

* `update_delay`: after one cluster finishes synchronization with the next
  cluster, it waits for `update_delay` more updates (inner SGD iterations)
  before passing the token to the next cluster. This delay is referred as
  `\tau_0` in our paper.  If you have less than 4 clusters, 64 is the
  recommended value. If you have more clusters, you can try a smaller update
  delay.

* `tolerence`: When one cluster writes its model delta to the next cluster, it
  will skip some model elements if they only changed very little. The
  threshold is given by the `tolerance` parameter. Usually you should set it
  between 1e-2 to 1e-6.

The following command runs the RCV1 dataset prepared above for 150 epochs with
40 threads, with a cluster size of 10, step size of 5e-01, step decay of 0.928
and update delay of 64:

```
bin/numasvm --epoch 150 --binary 1 --stepinitial 5e-01 --step_decay 0.928 --update_delay 64 --cluster_size 10 --split 40 data/rcv1_train.bin data/rcv1_test.bin
```

Don't forget to change the number of threads (the --splits argument) and the
cluster size to reflect your hardware configuration. For example, if you
have a dual-socket 20-core machine, you can change splits to 20 and cluster
size to 10 or 5. 

If you specify more threads than the total number of physical cores available,
hyper-threading cores will be used as well. This will maximize the computation
power of your machine. However please note that the number of clusters will
still be calculated using the total number of physical cores.  For example, if
you have a dual-socket 12-core machine (24 threads with hyperthreading) and you
run HogWild++ with a cluster size of 6, and 24 worker threads, only 2 clusters
will be created instead of 4.

We provide two python scripts, `collect_svm.py` and `collect_numasvm.py` to run
the baseline HogWild! and our improved algorithm, HogWild++, using the same
datasets and parameters in our paper.  Before running these scripts, make sure
you have downloaded all datasets and put them inside the `data` directory.  To
run a quick experiment on your machine, just run `python collect_svm.py` to
collect HogWild! results and `python collect_numasvm.py` to collect HogWild++
results. Results of runs will be saved in `svm_mmdd-hhmmss` and
`numasvm_mmdd-hhmmss` folders where `mmdd-hhmmss` represents date and time.
You can change the scripts to customize all parameters, like number of threads,
cluster size, step size, etc.

How to modify
----------------------

If you are interested in changing HogWild++ and explore more possibilities,
you can start by reading the following source files:

* `src/numasvm_main.cpp`: this file contains the `main` function for `numasvm`.
  Also, function `CreateNumaClusterRoundRobinRingSVMModel` creates the model
  synchronization ring. If you are interested in changing the synchronization
  topology you should change this function.

* `src/numasvm/svm_exec.hxx`: The SVM solver. Function `ModelUpdate` contains
  the SGD update rule and model synchronization procedure.

* `hazytl/include/hazy/thread/thread_pool-inl.h`: An enhanced thread pool
  implementation supporting CPU topology detection and affinity assignment.

Additional Information
----------------------

If your have any questions or comments, please open an issue on Github,
or send an email to ecezhang@ucdavis.edu. We appreciate your feedback.

