# Prerequisites for using the iNaturalist API
1) Install Python3 and pip3
See also https://realpython.com/installing-python/#how-to-install-python-on-linux
$ sudo apt-get update
$ sudo apt-get install python3.8 python3-pip
Set up virtual environment
See also https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/
2) install venv (on a Linux system, in this case)
$ sudo apt install python3.8-venv
3) create virtual environment
$ python3 -m venv /home/manuel/Nextcloud/LELE_2021/py3/venv
4) activate virtual environment
$ source /home/manuel/Nextcloud/LELE_2021/py3/venv/bin/activate
5) install pyinaturalist and argparse
$ pip3 install pyinaturalist
$ pip3 install argparse
6) deactivate virtual environment
$ deactivate
7) verify deactivation of virtual environment
$ which python
