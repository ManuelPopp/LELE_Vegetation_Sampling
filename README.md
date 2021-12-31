# LELE_Vegetation_Sampling
Scripts for vegetation sampling tasks within the LELE Project

# Prerequisites for usage of the iNaturalist API
Install Python3 and pip3
See also https://realpython.com/installing-python/#how-to-install-python-on-linux
$ sudo apt-get update
$ sudo apt-get install python3.8 python3-pip
Set up virtual environment
See also https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/
# install venv (on a Linux system, in this case)
$ sudo apt install python3.8-venv
# create virtual environment
$ python3 -m venv /home/manuel/Nextcloud/LELE_2021/py3/venv
# activate virtual environment
$ source /home/manuel/Nextcloud/LELE_2021/py3/venv/bin/activate
# install pyinaturalist and argparse
$ pip3 install pyinaturalist
$ pip3 install argparse
# deactivate virtual environment
$ deactivate
# verify deactivation of virtual environment
$ which python
