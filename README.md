# Identifying Network Hubs with the Partial-Correlation GLASSO

This repository contains code to reproduce the experiments and figures from the paper  
**“Identifying Network Hubs with the Partial-Correlation Graphical LASSO.”**

---

## TODO

* Check whether those all dependente packages are listed to install
* Delete the `raw_expetiments` folder

---

## Installation

```r
# Install dependencies
install.packages(c("remotes","Matrix","glasso","glmnet","igraph",
                   "ggplot2","dplyr","readr","yaml","pbapply","microbenchmark"))

# Install our method
remotes::install_github("PrzeChoj/pcglassoFast")

# Install competitive method
remotes::install_github("JackStorrorCarter/PCGLASSO")
```

---

## Structure

```
├─ experiments/     # Scripts for plots in manuscript
├─ data/            # Raw and processed data
├─ outputs/         # Results, figures, and tables
└─ manuscript/      # Paper source (arXiv.tex) and PDF
```

---

## Citation

If you use this code, please cite our paper and the [**pcglassoFast**](https://github.com/PrzeChoj/pcglassoFast) package.

```
@misc{bogdan2025identifyingnetworkhubspartial,
      title={Identifying Network Hubs with the Partial Correlation Graphical LASSO}, 
      author={Małgorzata Bogdan and Adam Chojecki and Ivan Hejný and Bartosz Kołodziejek and Jonas Wallin},
      year={2025},
      eprint={2508.12258},
      archivePrefix={arXiv},
      primaryClass={math.ST},
      url={https://arxiv.org/abs/2508.12258}, 
}
```
