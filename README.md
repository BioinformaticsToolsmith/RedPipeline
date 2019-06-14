# RedPipeline
A pipeline for detecting and clustering repeats

Requirement

These four programs must be installed and accessable globally, i.e. any of them can be executed without providing the full path to the program. You can obtain these prorgram by executing the following commands in a UNIX terminal:

git clone https://github.com/TulsaBioinformaticsToolsmith/Red.git
git clone https://github.com/TulsaBioinformaticsToolsmith/LtrDetector.git
git clone https://github.com/TulsaBioinformaticsToolsmith/Look4TRs.git
git clone https://github.com/TulsaBioinformaticsToolsmith/MeShClust.git

Input

A directory containing FASTA files with ".fa" extension and an ouput directory where intermediate files are stored. 

Output

The most important output by this pipeline is "library.fa," which contains a representative sequence for each family and its possible type e.g. DNA transposon or LINE.

License

Academic use: The software is provided as-is under the GNU GPLv3.

Any restrictions to use for-profit or non-academics: License needed.
