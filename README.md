# microclustr

R package `microclustr` performs entity resolution for databases with categorical fields using partition-based Bayesian clustering models. It includes a new family of microclustering prior models for random partitions -- exchangeable sequence of clusters (ESC) models, and the traditional Dirichlet and Pitman-Yor process priors.

## Installation

```r
devtools::install_github("cleanzr/microclustr", build_vignettes = TRUE)
```


## Background

Entity resolution (record linkage or de-duplication) is used to join multiple databases to remove duplicate entities. Recent methods tackle the entity resolution problem as a clustering task. While traditional Bayesian random partition models assume that the size of each cluster grows linearly with the number of data points, this assumption is not appropriate for applications such as entity resolution. This problem requires models that yield clusters whose sizes grow sublinearly with the total number of data points -- **the microclustering property**. The `microclustr` package includes random partition models that satisfy the microclustering property and implements entity resolution with categorical data.

## Main functions

The `microclustr` package contains the main function `SampleCluster` that are used to perform entity resolution as a clustering task under various possible random partition priors. The possible choices of priors are:


- "DP": Random partition induced by Dirichlet process (DP), also known as Chinese restaurant process prior.
- "PY": Random partition induced by Pitman-Yor process (PY).
- "ESCNB": Exchangeable sequence of clusters (ESC) model with zero-truncated negative binomial distribution.
- "ESCD": ESC model with Dirichlet distribution.
- "ESCP": ESC model with zero-truncated Poisson distribution.
- "ESCB": ESC model with zero-truncated binomial distribution with fixed maximum cluster size.
- "ESCBshift": ESC model with shifted binomial distribution, non-fixed maximum cluster size. 

Run the following R codes for more details on the main function:
```r
library(microclustr)
?SampleCluster
```

Additionally, the functions `mean_fnr` and `mean_fdr` can be used to evaluate the record linkage performance when ground truth is available.


For more extensive documentation of the use of this package, please refer to the vignette.

```r
?microclustr
browseVignettes("microclustr")
```

## Citation

This package implements the methods introduced in the following papers:

- Betancourt, B., Zanella, G., & Steorts, R. C. (2020). Random partition models for microclustering tasks. Journal of the American Statistical Association, 1-13.

````{verbatim}
@article{betancourt2020random,
  title={Random partition models for microclustering tasks},
  author={Betancourt, Brenda and Zanella, Giacomo and Steorts, Rebecca C},
  journal={Journal of the American Statistical Association},
  pages={1--13},
  year={2020},
  publisher={Taylor \& Francis}
}
````

- Lee, C. J., & Sang, H. (2022). Why the Rich Get Richer? On the Balancedness of Random Partition Models.  Proceedings of the 39th International Conference on Machine Learning (ICML), PMLR 162:12521 - 12541.

````{verbatim}
@InProceedings{lee2022why,
  title={Why the Rich Get Richer? {O}n the Balancedness of Random Partition Models},
  author={Lee, Changwoo J and Sang, Huiyan},
  booktitle={Proceedings of the 39th International Conference on Machine Learning},
  pages={12521--12541},
  year={2022},
  volume={162},
  series={Proceedings of Machine Learning Research},
  publisher={PMLR}
}
````

## Acknowledgements

This work was supported through the National Science Foundation through NSF CAREER Award 1652431 (PI: Steorts).
