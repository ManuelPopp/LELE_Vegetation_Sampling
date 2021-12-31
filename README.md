# LELE_Vegetation_Sampling
Scripts for vegetation sampling tasks within the LELE Project

# Prerequisites for usage of the iNaturalist API
1) Install Python3 and pip3
See also https://realpython.com/installing-python/#how-to-install-python-on-linux
$ sudo apt-get update
$ sudo apt-get install python3.8 python3-pip
Set up virtual environment
See also https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/
2) install venv (on a Linux system, in this case)
$ sudo apt install python3.8-venv
2a) create virtual environment
$ python3 -m venv /home/manuel/Nextcloud/LELE_2021/py3/venv
2b) activate virtual environment
$ source /home/manuel/Nextcloud/LELE_2021/py3/venv/bin/activate
3) install pyinaturalist and argparse
$ pip3 install pyinaturalist
$ pip3 install argparse
4a) deactivate virtual environment
$ deactivate
4b) verify deactivation of virtual environment
$ which python
