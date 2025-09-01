# HistoMapST
HistoMapST is an R pipeline that allows to combine histopathological annotation with spot information from Visium experiment (v1/v2).  

The overall workflow is shown below:  
![HistoMapST flowchart.](/images/HistoMapST_flowchart.png)

The H&E image is used a source for two complimentary pipelines: high resolution/low dimensionality histopathological annotation (i.e. QuPath) and low resolution/high dimensionality spaceranger pipeline. An output from QuPath containing cell centroids and additional metadata is mapped against spot coordinates from spaceranger.  
