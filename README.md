# Pyroptosis
Machine learning model based on cell pyroptosis predicts treatment and prognosis of lung adenocarcinoma

#### Chapter 1. Data preparing

- 从[GeneCards网站](https://www.genecards.org/)输入“pyroptosis”，检索并下载802个细胞焦亡相关基因。
- 整理后各数据库中包含内容如下：

| 数据集    | 总样本数 | LUAD | Normal | 病理参数                             | 预后信息 |
| --------- | -------- | ---- | ------ | ------------------------------------ | -------- |
| GSE3141   | 111      | 58   | 0      | -                                    | OS       |
| GSE30219  | 307      | 124  | 14     | Age,Gender,TNM                       | OS/PFS   |
| GSE31210  | 246      | 226  | 20     | Age,Gender,Smoke,Stage,EGFR,KRAS,ALK | OS/PFS   |
| GSE41271  | 275      | 183  | 0      | Age,Gender,Smoke,Stage               | OS/PFS   |
| GSE50081  | 181      | 129  | 0      | Age,Gender,Stage,Smoke,TNM           | OS/PFS   |
| TCGA_LUAD | 594      | 535  | 59     | ALL                                  | OS/PFS   |

#### Chapter 2. DEGs

- DEGs筛选标准为：abs（log2FC）>1 AND P.adj <0.05。
- 采用DESeq2法进行差异基因分析

​	整理后有Normal和LUAD对比的数据如下：

| 数据集                      | LUAD | Normal | Up   | Down |
| --------------------------- | ---- | ------ | ---- | ---- |
| TCGA_LUAD<br />（DESeq2法） | 535  | 59     | 9878 | 2843 |

#### Chapter 3. DEGOIs ( Differential expression genes of interest)

- 取有Normal样本的数据集进行分析，用DEGsup/DEGsdown ∩ GOIs的方法获得DEGOIs。

  整理后有Normal和LUAD对比的数据如下：

| 数据集                 | DEGOIsup | DEGOIsdown |
| ---------------------- | -------- | ---------- |
| TCGA_LUAD （DESeq2法） | 121      | 66         |

PS. 在TCGA_LUAD 的DESeq2文件种顺便绘制了火山图。

#### Chapter 4. ML (Machine_learning)

​	选定TCGA作为训练集，GSE30219、GSE31210、GSE3141、GSE41271和GSE50081作为验证集，并选定RSF+Ridge作为最佳算法。

#### Chapter 5. Kaplan-Meier curve

​	根据Chapter 4中的结果，我们对TCGA（训练集）和GSE30219、GSE31210、GSE3141、GSE41271、GSE50081（验证集）进行Kaplan-Meier curve的绘制。

#### Chapter 6. ROC

​	在Chapter 4中，我们并没有去绘制ROC曲线，在本章中我们进行ROC曲线的绘制。

整理后结果一览：

| 数据集   | AUC-1year | AUC-3year | AUC-5year |
| -------- | --------- | --------- | --------- |
| GSE3141  | 0.687     | 0.725     | 0.783     |
| GSE30219 | 0.83      | 0.727     | 0.757     |
| GSE31210 | 0.815     | 0.667     | 0.685     |
| GSE41271 | 0.707     | 0.693     | 0.699     |
| GSE50081 | 0.635     | 0.669     | 0.69      |
| UCSC_all | 0.73      | 0.699     | 0.698     |



#### Chapter 7. Cox

​	在Chapter 4中，我们绘制了单因素Cox分析的森林图，在本章中，我们对上述6个数据集（其中GSE3141没有临床病理参数，故实际上只有5个数据集）进行多因素Cox回归分析。

整理后结果一览：

| 数据集                           | 单因素Cox       | 多因素Cox       |
| -------------------------------- | --------------- | --------------- |
| GSE3141                          | -               | -               |
| GSE30219                         | RS+T+N+Stage    | RS              |
| GSE31210(去除基因突变比较有意义) | RS+Stage        | RS(0.087)+Stage |
| GSE41271                         | RS+Stage+Gender | RS+Stage        |
| GSE50081                         | RS+T+N+Stage    | Stage           |
| UCSC_all                         | RS+T+N+Stage    | RS+T+N          |



#### Chapter 8. Nomogram-Calibration-DCA

​	只对UCSC数据集中的all样本进行Nomogram、Calibration和DCA曲线绘制；纳入RS+T+N等3个因素。

#### Chapter 9. Immune infiltration

​	 对UCSC数据集进行免疫浸润分析，包括多种算法下免疫细胞浸润水平及ssGSEA算法评估富集的信号通路。

#### Chapter 10. TMB

​	进行TMB分析。

#### Chapter 11. MSI

​	进行MSI分析。

#### Chapter 12. TIDE

​	进行TIDE分析。

#### Chapter 13. Drug

​	进行药物敏感性分析。

#### Chapter 14. Supplymentary

##### Chapter 14.1 Compare of C-index

​	在本节中，我们比较TCGA、GSE30219、GSE31210、GSE41271、GSE50081等数据集（GSE3141中没有其他临床参数，故舍去）中Risk score和其他病理参数的C-index。

##### Chapter 14.2 Drug

​	在Chapter 13中，我们分析了198种化合物的药物敏感性，也输出了很多结果，但是这个图片太多了，需要进行归纳整理，而且Chapter 13中的代码所输出的图片风格和其他章节（如免疫分析章节）不一致，因此在此重新绘图。

##### Chapter 14.3 29_PRGs

​	整理29个差异表达预后相关PRGs的信息和GO/KEGG分析。

##### Chapter 14.4 GSE3141_Bootstrap

​	GSE3141数据集中样本较少，用Bootstrap重抽样法计算C-index为0.643 ± 0.055，与前文计算的C-index无统计学差异。Bootstrap抽样次数为1000次。

##### Chapter 14.5 density

​	在TCGA、GSE3141、GSE30219、GSE31210、GSE41271和GSE50081中绘制Risk score分布的直方图和密度图。

##### Chapter 14.6 importance

​	本研究采用RSF+Ridge作为最优算法，在进行RSF分析时，输出29个预后相关PRGs的贡献度。

##### Chapter 14.7 29_PRGs_unicox

​	绘制29个差异表达预后相关PRGs的单因素Cox分析森林图。
