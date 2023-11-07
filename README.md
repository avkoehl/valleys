# getting leafmap to work

Getting leafmap to work in jupyternotebooks has several undocumented
dependencies and issues

need gdal, ipyleaflet, openEXR, localtileserver, jupyter-contrib-nbextensions

brew install gdal  

then based on the version installed 

poetry add GDAL==<version>

rest should be straight forward with pip


