# miRmap - Comprehensive prediction of microRNA target repression strength

The *miRmap* library is a Python library predicting the repression strength of microRNA (miRNA) targets. The model combines:

 - **thermodynamic** features: *ΔG duplex*, *ΔG binding*, *ΔG seed duplex*, *ΔG seed binding*, *ΔG open* and *ΔG total*,
 - **evolutionary** features: *BLS* and *PhyloP*,
 - **probabilistic** features: *P.over binomial* and *P.over exact*, and
 - **sequence-based** features: *AU content*, *UTR position* and *3' pairing*.

> **NOTE** This is a re-implementation by the same author of the miRmap library published in 2011 with most of the core algorithm unchanged. Please refer to the [miRmap1](https://git.sr.ht/~vejnar/miRmap1) repository for the old library.

## Online

*miRmap* is available [online](https://mirmap.ezlab.org).

## Download

See [refs](https://git.sr.ht/~vejnar/miRmap/refs) page.

## Citation

If you use miRmap in your research, please cite:

> Charles E. Vejnar and Evgeny M. Zdobnov<br>
> miRmap: Comprehensive prediction of microRNA target repression strength<br>
> Nucleic Acids Research 2012 Dec 1;40(22):11673-83. doi: [10.1093/nar/gks901](https://doi.org/10.1093/nar/gks901)

## Installation

### External dependencies

1. The [Spatt](https://nuel.perso.math.cnrs.fr/spatt) library is necessary for the *P.over exact* feature.

    Download the latest [Spatt](https://nuel.perso.math.cnrs.fr/spatt) tarball (Version 2.1 was successfully tested), then do:

    ```bash
    cd spatt-<version>
    mkdir build
    cd build
    cmake -DWITH_SHARED_LIB=ON ..
    make
    ```

    To create the library at `libspatt2/libspatt2.so`.

2. [PHAST](http://compgen.bscb.cornell.edu/phast) is necessary for the evolutionary features. Compilation instructions of PHAST are available in this [PKGBUILD](https://aur.archlinux.org/cgit/aur.git/tree/PKGBUILD?h=phast).

### Using `pip`

After installing external dependencies, install *miRmap*:
```bash
pip3 install mirmap
```

Python dependencies [ViennaRNA](https://pypi.org/project/ViennaRNA) and [dendropy](https://pypi.org/project/DendroPy) will be installed from [PyPI](https://pypi.org/).

## Example

```python
import mirmap.target

utr_seq = "ATAGACTGTACATTATGAAGAATACCCAGGAAGACTTTGTGACTGTCACTTGCTGCTTTTTCTGCGCTTCAGTAACAAGT"
mirna_seq = "UAGCAGCACGUAAAUAUUGGCG".replace("U", "T")

targets = mirmap.target.find_targets_with_seed(utr_seq, mirna_seq)
print(targets[0].report())
``````
This will return:
```
          36                   57
          |                    |
CAGGAAGACTTTGTGACTGTCACTTGCTGCTTTTTCTGCGCT
                        |||||||
          GCGGTTATAAATGCACGACGAT
```
Then we can calculate the scores of the miRNA target:
```python
import mirmap.if_lib_spatt
import mirmap.scores

if_spatt = mirmap.if_lib_spatt.Spatt("bin/linux_x86_64/libspatt2.so")

scores = mirmap.scores.calc_scores(
    targets[0],
    if_spatt=if_spatt,
)
print(mirmap.scores.report_scores(scores))
```
This will return:
```
 ΔG duplex (kcal/mol)     -13.8
 ΔG binding (kcal/mol)    -11.95
 ΔG open (kcal/mol)       14.03
 ΔG total (kcal/mol)      0.2345
 AU content               0.6574
 UTR position             22.0
 3' pairing               1.0
 TargetScan score         23.66
 Probability (Exact)      0.03813
 Probability (Binomial)   0.006405
 Conservation (BLS)       0.0
 Conservation (PhyloP)    1.0
 miRmap score             -0.3122
```

## License

The *miRmap* library is distributed under the GNU GPL v3 (see /LICENSE).

Copyright © 2011-2024 Charles E. Vejnar
