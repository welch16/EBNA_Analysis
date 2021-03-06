



## Association tests among peak sets


Using the annotation into those three categories, we performed
chi-square tests to compare the annotation distribution among the
different sets.




|Set1   |Set2   | statistic| p.value| minus_log10|
|:------|:------|---------:|-------:|-----------:|
|EBNA2  |EBNA3A |  168.9895|  0.0000|     36.6956|
|EBNA2  |EBNA3B |    3.5597|  0.1687|      0.7730|
|EBNA2  |EBNA3C |  229.7907|  0.0000|     49.8984|
|EBNA2  |RBPJ   |   54.0071|  0.0000|     11.7275|
|EBNA3A |EBNA3B |   91.3476|  0.0000|     19.8359|
|EBNA3A |EBNA3C |  448.3130|  0.0000|     97.3499|
|EBNA3A |RBPJ   |  344.8276|  0.0000|     74.8784|
|EBNA3B |EBNA3C |  178.0941|  0.0000|     38.6726|
|EBNA3B |RBPJ   |   44.1427|  0.0000|      9.5855|
|EBNA3C |RBPJ   |  125.8572|  0.0000|     27.3296|

Here we can see that almost all pairs are strongly associated on this
partition (promoters,gene body and intergenic), except EBNA2 and
EBNA3B.



|Set1   |Set2   | statistic| p.value| minus_log10|
|:------|:------|---------:|-------:|-----------:|
|EBNA2  |EBNA3A |  646.6367|  0.0000|    140.4154|
|EBNA2  |EBNA3B |   13.6394|  0.0011|      2.9618|
|EBNA2  |EBNA3C | 1589.5762|  0.0000|         Inf|
|EBNA2  |RBPJ   |  120.6882|  0.0000|     26.2071|
|EBNA3A |EBNA2  |  234.9315|  0.0000|     51.0147|
|EBNA3A |EBNA3B |  185.8449|  0.0000|     40.3557|
|EBNA3A |EBNA3C | 1899.6263|  0.0000|         Inf|
|EBNA3A |RBPJ   |  521.1501|  0.0000|    113.1663|
|EBNA3B |EBNA2  |    4.8207|  0.0898|      1.0468|
|EBNA3B |EBNA3A |  187.1007|  0.0000|     40.6284|
|EBNA3B |EBNA3C |  642.5271|  0.0000|    139.5230|
|EBNA3B |RBPJ   |   65.0352|  0.0000|     14.1222|
|EBNA3C |EBNA2  |  291.4114|  0.0000|     63.2792|
|EBNA3C |EBNA3A |  700.7123|  0.0000|    152.1577|
|EBNA3C |EBNA3B |  281.0352|  0.0000|     61.0260|
|EBNA3C |RBPJ   |  161.2664|  0.0000|     35.0185|
|RBPJ   |EBNA2  |   98.6643|  0.0000|     21.4247|
|RBPJ   |EBNA3A | 1120.3588|  0.0000|    243.2828|
|RBPJ   |EBNA3B |  138.8953|  0.0000|     30.1607|
|RBPJ   |EBNA3C |  689.8437|  0.0000|    149.7977|

This are goodness of fit tests, therefore we are testing if the
probability distributions of the first set could be generated as a
sample from the probability distribution of the second set.

## Conditional on Dnase




|Set1   |Set2   | statistic| p.value| minus_log10|
|:------|:------|---------:|-------:|-----------:|
|EBNA2  |EBNA3A |  346.2600|  0.0000|     75.1894|
|EBNA2  |EBNA3B |    3.1724|  0.2047|      0.6889|
|EBNA2  |EBNA3C |   70.6071|  0.0000|     15.3321|
|EBNA2  |RBPJ   |   37.0462|  0.0000|      8.0445|
|EBNA3A |EBNA3B |  192.7779|  0.0000|     41.8612|
|EBNA3A |EBNA3C |  364.0075|  0.0000|     79.0432|
|EBNA3A |RBPJ   |  549.8193|  0.0000|    119.3917|
|EBNA3B |EBNA3C |   71.2123|  0.0000|     15.4636|
|EBNA3B |RBPJ   |   40.3631|  0.0000|      8.7647|
|EBNA3C |RBPJ   |   35.7305|  0.0000|      7.7588|

Again in thise case, we are seeing that almost all pairs of peak
sets are strongly associated conditional on Dnase binding except for
EBNA2 and EBNA3B again.



|Set1   |Set2   | statistic| p.value| minus_log10|
|:------|:------|---------:|-------:|-----------:|
|EBNA2  |EBNA3A | 1472.2861|  0.0000|    319.7029|
|EBNA2  |EBNA3B |   11.6849|  0.0029|      2.5373|
|EBNA2  |EBNA3C |  626.7066|  0.0000|    136.0876|
|EBNA2  |RBPJ   |   84.6159|  0.0000|     18.3741|
|EBNA3A |EBNA2  |  478.4613|  0.0000|    103.8965|
|EBNA3A |EBNA3B |  381.0112|  0.0000|     82.7355|
|EBNA3A |EBNA3C | 1657.2395|  0.0000|         Inf|
|EBNA3A |RBPJ   |  829.5572|  0.0000|    180.1361|
|EBNA3B |EBNA2  |    4.3568|  0.1132|      0.9461|
|EBNA3B |EBNA3A |  440.6683|  0.0000|     95.6899|
|EBNA3B |EBNA3C |  297.1910|  0.0000|     64.5342|
|EBNA3B |RBPJ   |   59.8863|  0.0000|     13.0041|
|EBNA3C |EBNA2  |   81.6389|  0.0000|     17.7277|
|EBNA3C |EBNA3A |  597.0912|  0.0000|    129.6567|
|EBNA3C |EBNA3B |   99.2293|  0.0000|     21.5474|
|EBNA3C |RBPJ   |   42.5278|  0.0000|      9.2348|
|RBPJ   |EBNA2  |   66.1668|  0.0000|     14.3679|
|RBPJ   |EBNA3A | 1994.3084|  0.0000|         Inf|
|RBPJ   |EBNA3B |  125.2345|  0.0000|     27.1943|
|RBPJ   |EBNA3C |  240.4488|  0.0000|     52.2128|

Similarly this are goodnes of fit tests as above where we are testing
if the conditional distribution on Dnase binding could be generated
from the second set's distribution.


Note: This analysis doesn't consider any spatial information about the
peak sets.

