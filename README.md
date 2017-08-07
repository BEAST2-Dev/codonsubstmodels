# Codon Models

Models describing evolution at the codon level allow 
the estimation of measures of the selective forces 
acting on proteins 
([Kosiol et. al. 2000](https://academic.oup.com/mbe/article-abstract/24/7/1464/986344)). 

The current implementation is M0, also called as GY94, first published by 
[Goldman & Yang 1994](https://academic.oup.com/mbe/article-abstract/11/5/725/1008711), 
and late designated to be one of the M-series (M0-M13) proposed by    
[Yang et al. 2000](http://www.genetics.org/content/155/1/431.short)

The relative instantaneous substitution rate _q<sub>ij</sub>_ from codon _i_ to codon _j_ 
can be calculated as:
![Substitution rate from codon i to codon j](./figures/codon-q.png)
