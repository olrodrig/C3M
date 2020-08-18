# Color-Color Curve Method

When a Type II supernova (SN II) is affected by a non-negligible host galaxy reddening Eh(B-V), the quantification of this value is crucial to accurately infer physical parameters of the SN. In [Rodríguez et al. (2014)](https://ui.adsabs.harvard.edu/abs/2014AJ....148..107R/abstract) I shown that for SNe II the B-V versus V-I color-color curve (C3) during the "plateau" phase can be used to estimate Eh(B-V) through the method proposed by [Natali et al. 1994](https://ui.adsabs.harvard.edu/abs/1994A%26A...289..756N/abstract), which was originally developed to estimate interstellar color excesses for open clusters. This method measures Eh(B-V) from the vertical displacement of a linear C3 with respect to an unreddened one (see [Figure 3](https://iopscience.iop.org/article/10.1088/0004-6256/148/6/107#aj501994f3) of [Rodríguez et al. 2014](https://ui.adsabs.harvard.edu/abs/2014AJ....148..107R/abstract)).

In [Rodríguez et al. (2019)](https://ui.adsabs.harvard.edu/abs/2019MNRAS.483.5459R/abstract) I developed a routine to compute Eh(B-V) given a set of BVI light curves evaluated at the same epochs, called **Color-Color Curve Method (C3M)**. This routine assumes a linear unreddened V-I versus B-V C3 with slope of 0.45±0.07 and y-intercept of 0.107±0.053 mag, and adopts the [Fitzpatrick
(1999)](https://ui.adsabs.harvard.edu/abs/1999PASP..111...63F/abstract) extinction curve with Rv=3.1.

For any question, email me at olrodrig@gmail.com

**If you use the C3M code in your work, please cite [Rodríguez et al. (2014)](https://ui.adsabs.harvard.edu/abs/2014AJ....148..107R/abstract) and [Rodríguez et al. (2019)](https://ui.adsabs.harvard.edu/abs/2019MNRAS.483.5459R/abstract).**
