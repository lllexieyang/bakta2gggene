# bakta2gggene
Extract annotation of contigs containing specific gene from [bakta](https://github.com/oschwengers/bakta) results (.tsv).  
Generate input file for drawing gene arrow maps with [gggenes](https://cran.r-project.org/web/packages/gggenes/vignettes/introduction-to-gggenes.html).  
  
Version 1.0  
Last change: Sep 11 2023  
Copy right reserved : Lu Yang (yanglu2016@cau.edu.cn)  
  
## 1. Basic usage
#### Extracting [gene_name] from all .tsv files in current directory (./)
    python3 bakta2gggene.py -g [gene_name]  
#### Extracting [gene_name] from input.tsv in current directory (./)
    python3 bakta2gggene.py -g [gene_name] -i input.tsv
#### Extracting [gene_name] from input.tsv in current directory (./) and save results in ./result/
    python3 bakta2gggene.py -g [gene_name] -i input.tsv -o result
#### Extracting [gene_name] from input1.tsv and input2.tsv (could be more) in current directory (./)
    python3 bakta2gggene.py -g [gene_name] -i input1.tsv input2.tsv
#### Extracting [gene_name] from all .tsv files in [directory of bakta tsv results]
    python3 bakta2gggene.py -g [gene_name] -d [directory of bakta tsv results]
#### Extracting [gene_name] from input.tsv files in [directory of bakta tsv results]
    python3 bakta2gggene.py -g [gene_name] -d [directory of bakta tsv results] -i input.tsv
## 2. Merge results into one and generate gggenes input file
#### If only *-d* is provided, it will process all input data in the specified directory and merge the results.  
如果仅提供了 ***-d*** 参数，将处理指定目录中的所有输入数据并合并结果。
#### If *-i* is used to specify specific input data, merging will not occur.  
如果使用 ***-i*** 来指定特定的输入数据，则不会执行合并操作。
#### Additionally, if *-m* or *--merge* is explicitly used, it will override the default behavior and trigger the merging process.
此外，如果通过 ***-m*** 或 ***--merge*** 参数指定了合并，则无视默认规则，执行合并。  
#### Example：
     python3 bakta2gggene.py -g [gene_name] -i input1.tsv input2.tsv -m

## 3. Maximum distance from the target gene
#### Specify the maximum distance for extracting annotation results relative to the target gene.
指定相对于目标基因提取注释结果的最大距离。
#### Example：
    python3 bakta2gggene.py -g [gene_name] -i input1.tsv input2.tsv --merge --max 10000
> Note: default 5kb  
