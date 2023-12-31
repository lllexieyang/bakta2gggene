# bakta2gggene
Extract "cds" annotation of contigs containing specific gene from [bakta](https://github.com/oschwengers/bakta) results (.tsv).  
Generate input file for drawing gene arrow maps with [gggenes](https://cran.r-project.org/web/packages/gggenes/vignettes/introduction-to-gggenes.html).  
  
Version 2.0  
Last change: Dec 5 2023  
Copy right reserved : Lu Yang (yanglu2016@cau.edu.cn)  
  
## 1. Basic usage

**Extracting [gene_name] from all .tsv files in current directory (./)**  

    python3 bakta2gggene.py -g [gene_name]  
**Extracting [gene_name] from input.tsv in current directory (./)**

    python3 bakta2gggene.py -g [gene_name] -i input.tsv
**Extracting [gene_name] from input.tsv in current directory (./) and save results in ./result/**

    python3 bakta2gggene.py -g [gene_name] -i input.tsv -o result
**Extracting [gene_name] from input1.tsv and input2.tsv (could be more) in current directory (./)**

    python3 bakta2gggene.py -g [gene_name] -i input1.tsv input2.tsv
**Extracting [gene_name] from all .tsv files in [directory of bakta tsv results]**

    python3 bakta2gggene.py -g [gene_name] -d [directory of bakta tsv results]
**Extracting [gene_name] from input.tsv files in [directory of bakta tsv results]**

    python3 bakta2gggene.py -g [gene_name] -d [directory of bakta tsv results] -i input.tsv
## 2. Merge results into one and generate gggenes input file  
**If only *-d* is provided, it will process all input data in the specified directory and merge the results.**  
  如果仅提供了 ***-d*** 参数，将处理指定目录中的所有输入数据并合并结果。  

**If *-i* is used to specify specific input data, merging will not occur.**  
  如果使用 ***-i*** 来指定特定的输入数据，则不会执行合并操作。  

**Additionally, if *-m* or *--merge* is explicitly used, it will override the default behavior and trigger the merging process.**  
  此外，如果通过 ***-m*** 或 ***--merge*** 参数指定了合并，则无视默认规则，执行合并。    

**Example：**

     python3 bakta2gggene.py -g [gene_name] -i input1.tsv input2.tsv -m

## 3. Maximum distance from the target gene
**Specify the maximum distance for extracting annotation results relative to the target gene.**  
指定相对于目标基因提取注释结果的最大距离。

**Example：**

    python3 bakta2gggene.py -g [gene_name] -i input1.tsv input2.tsv --merge --max 10000
> Note: default 5kb  

## 4. Exact match of gene name
加-e精准匹配基因名（默认为查找指定基因开头的结果）  

**Extracting gene exactly the same as provided -> mcr-1.1 only**

    python3 bakta2gggene.py -g mcr-1.1 -e 
**or extracting genes start with "mcr-1.1" -> mcr-1.1, mcr-1.10, mcr-1.12 et al.**

    python3 bakta2gggene.py -g mcr-1.1 
> Note: default set is to extract genes start with [gene_name], suitable for gene clusters like xxxA, xxxB, xxxC...

## 5. Extract gene annotations based on a list 
根据列表提取基因的注释信息（比如ARG_list.txt，每行是一个基因的名字）  

**e.g., ARG_list.txt, where each line is a gene name**  

    python3 bakta2gggene.py -i input1.tsv -l [gene_list_name]
> ⚠️ Note: Extract only the annotation information for the genes themselves, without extracting the gene environment.  
>  注意：只提取基因本身的注释信息，不提取基因环境  

### Version Log
Version 1.0  Sep 11 2023  

Version 2.0  Dec 5 2023  
New feature: Extract gene annotation results based on a specified list (genes only).  
2.0版本新增功能：根据指定列表提取基因注释结果（仅基因）
