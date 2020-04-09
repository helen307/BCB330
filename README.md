# BCB330
[Journal and Notes](https://github.com/helen307/BCB330/wiki)<br/>
# Workflow
[Summarized Report](https://github.com/helen307/BCB330/wiki/final_report)<br/>
# Reports
1. [Data collection](https://htmlpreview.github.io/?https://github.com/helen307/BCB330/blob/master/get_data_mir_expr/normal_tissues/heart/heart.nb.html)<br/>
* This is an example using tissue - heart.

2. [Statistical overview](https://htmlpreview.github.io/?https://github.com/helen307/BCB330/blob/master/mir_expr_table/iid.nb.html)<br/>
* Includes working example of tissue-specific miRNA catalogue.
* Data from the report was all extracted from the analysis here.

3. Shiny app for tissue-specific miRNA-gene interaction
* Run /mir_gene_table/final_shiny_app.Rmd

# Notes
* __IID human gene data__ is too large, please download it and save it under mir_expr_table/data/ to run iid.Rmd
* __mir_expr_table/iid.nb.html__ is a file that combines gtex.Rmd, iid.Rmd, mir_expr_table.Rmd for the sake of a complete report on all the data analysis
  * mir_expr_table/__mir_expr_table.Rmd__: only examines the miRNA expression. It has some helper files: api.R (from mirDIP), mirDIP_task_helper.Rmd, and task_helper.R (to complete 3 functionalities for miRNA expression catalogue)
  * gtex/__gtex.Rmd__: compares miRNA and gtex gene expression
  * mir_expr_table/__iid.Rmd__: compares miRNA and IID gene expression
* __mir_gene_table/__: 
  * __final_combined.Rmd__: combines the miRNA and IID gene expression data into one table. The resulting table is similar to the tissue-specific IID human gene-gene interaction table.
  * __final_shiny_app.Rmd__: contains tissue-specific miRNA-gene interaction table.
* __saved_expr_data/__: contains both novel and known miRNA, along with publication data.
  * All the tables have the same column names: miRNA name, expressivity sum over all the samples, number of samples selected
* __get_data_mir_expr/__: contains all the data collection step. An example file can be found in the heart/ folder, in html.
* __images/__: stores images for the final report.
* __Reference__: Did not keep the references in my journal. All references are stored in iid.bib and heart.bib for the html report.
  
