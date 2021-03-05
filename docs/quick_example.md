# a minimal executive example

Here we will show how to perform a simple but complete analysis after the installation of pipeline.

<!-- 准备数据文件 -->
## data preparation 
<!-- 从网盘上下载数据文件 -->
You can download a minimal dataset for testing from [test_data]()

Inside this dataset, 

## Compose the configuration file

By invoking new_config.py, you will get a copy of the template. then you can modify this config file to guide your analysis.

The format of invoking new_config.py is as following: 

```shell
new_config.py target
```

`target` is the place where your configuration file will be found.

If this `target` is an existing folder, then a default.cfg file will appear in this folder.

## modifiy this configure file

If you have modified the template after installation, then the changes you have made are already in effect on this copy. 

Before you modify the configuration file, it is recommended that you read [instructions on configure file ](./conf.md)


## build the annotation database manually (optional)

During the process, we will create a binary database file for GTF in order to query the genome information faster.
This '.db' database file will share base name with the GTF file and will locate in the same folder.

However, since there are different versions of GTF files for different species.
So creating the database may encounter some unexpected problems.

For this reason, a prudent solution is to create the relevant database manually




## run 

