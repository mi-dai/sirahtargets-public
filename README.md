# Pipeline to select SIRAH target

After cloning this repo,

first create a conda environment called `sirah_target` using the following command:

> $ conda env create -f sirah_target_environment.yml

and install SciServer in the `sirah_target` environment:
> $ conda activate sirah_target      
> (sirah_target)$ git clone  https://github.com/sciserver/SciScript-Python.git  
> (sirah_target)$ cd SciScript-Python/   
> (sirah_target)$ python Install.py  

then add `sirah_target` to your jupyter kernel:

> (sirah_target)$ ipython kernel install --name "sirah_target" --user

define the following enviroment variables in your system or define it later in the notebook:

> export SIRAHPIPE_DIR=<DIR_OF_THIS_REPO>   
> export SFD_DIR=<DIR_OF_SFDDATA>
 
download [sfddata](https://github.com/kbarbary/sfddata) for dust correction:

> wget https://github.com/kbarbary/sfddata/archive/master.tar.gz    
> tar xzf master.tar.gz

A directory "sfddata-master" will be created. Use it as `<DIR_OF_SFDATA>` above.

You should be able to run the `demo_notebook.ipynb` now.