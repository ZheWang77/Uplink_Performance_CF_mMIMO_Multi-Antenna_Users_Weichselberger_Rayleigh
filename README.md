# Uplink_Performance_CF_mMIMO_Multi-Antenna_Users_Weichselberger_Rayleigh

This is a code package is related to the following scientific article:

Z. Wang, J. Zhang, B. Ai, C. Yuen, and M. Debbah, "Uplink Performance of Cell-Free Massive MIMO With Multi-Antenna Users Over Jointly-Correlated Rayleigh Fading Channels," in IEEE Transactions on Wireless Communications, vol. 21, no. 9, pp. 7391-7406, Sep. 2022.

Available at: https://arxiv.org/abs/2110.04962

The package contains a simulation environment, based on Matlab, that reproduces the numerical results in the article. *We encourage you to also perform reproducible research!*

## Abstract of Article
In this paper, we investigate a cell-free massive MIMO system with both access points (APs) and user equipments (UEs) equipped with multiple antennas over jointly-correlated Rayleigh fading channels. We study four uplink implementations, from fully centralized processing to fully distributed processing, and derive their achievable spectral efficiency (SE) expressions with minimum mean-squared error successive interference cancellation (MMSE-SIC) detectors and arbitrary combining schemes. Furthermore, the global and local MMSE combining schemes are derived based on full and local channel state information (CSI) obtained under pilot contamination, which can maximize the achievable SE for the fully centralized and distributed implementation, respectively. We study a two-layer decoding implementation with an arbitrary combining scheme in the first layer and optimal large-scale fading decoding (LSFD) in the second layer. Besides, we compute novel closed-form SE expressions for the two-layer decoding implementation with maximum ratio (MR) combining. In the numerical results, we compare the SE performance for different implementation levels, combining schemes, and channel models. It is important to note that increasing the number of antennas per UE may degrade the SE performance.

## Content of Code Package

The package includes the source codes applied in this paper.


## License and Referencing

This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article:

```
@ARTICLE{04962,
author={Wang, Zhe and Zhang, Jiayi and Ai, Bo and Yuen, Chau and Debbah, Merouane},
journal={IEEE Trans. Wireless Commun.},
title={Uplink Performance of Cell-Free Massive {MIMO} With Multi-Antenna Users Over Jointly-Correlated {R}ayleigh Fading Channels},
year={2022},
month = {Sep.},
volume={21},
number={9},
pages={7391-7406},
doi={10.1109/TWC.2022.3158353}}
```
