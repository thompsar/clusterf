# ClusterF
*An application for finding and viewing chemical clusters in the Chembridge libraries based on HTS results*

## Getting Started

1. Download a copy of **clusterf** by clicking the **Code** button and selecting "Download ZIP"
    - The easiest way to **install** necessary python dependancies is using [Anaconda](https://www.anaconda.com/products/individual).
 
2. The following **dependancies must be installed** prior to running HTScampaign (do this by launching the `CMD.exe Prompt` from the Anaconda launcher, type the following commands and follow the prompts):
    - For panel: `conda install -c conda-forge panel`

3. The file structure is as follows:  

```
clusterf
|
+--- compound_libraries
|
|
+--- compound_subsets   (csv files containing **at least** a Compound and Category column
|
|
+--- src
|   
|
+--- clusterf.py  (main program to be run)
```

4. Deposit .csv files into the `compound_subsets` folder, ensuring each subset file containts **at least** a  `Compound` and `Category` column, with values for each.

5. Run clusterf in the command line by typing `panel serve clusterf.py --show` to open the application in your default browser.
