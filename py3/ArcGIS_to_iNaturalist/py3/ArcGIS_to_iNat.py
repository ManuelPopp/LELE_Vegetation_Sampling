#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 16:34:30 2022

@author: manuel
@version = "1.0.0"
"""

## directory functions
wd = "/home/manuel/Nextcloud/LELE_2021/py3/ArcGIS_to_iNaturalist/"

import os
def dir_app(filename = None):
    if filename == None:
        return os.path.join(wd, "appdata")
    else:
        return os.path.join(wd, "appdata", filename)

import zipfile
import tempfile
import pickle

## define helper functions
def load_pkl(file):
    with open(file, "rb") as f:
            return pickle.load(f)

def save_login(name, usr, pw, key):
    os.makedirs(dir_app(), exist_ok = True)
    zipObj = zipfile.ZipFile(dir_app(name + ".zip"), "w")
    tmp_dir = tempfile.mkdtemp()
    for var, var_name in zip([usr, pw, key], ["usr", "pw", "key"]):
        with open(os.path.join(tmp_dir, var_name + ".pkl"),\
                  "wb") as f:
            pickle.dump(var, f)
        zipObj.write(os.path.join(tmp_dir, var_name + ".pkl"))

def get_login(name):
    zf = zipfile.ZipFile(dir_app(name + ".zip"))
    tempdir = tempfile.mkdtemp()
    out = dict()
    with zipfile.ZipFile(dir_app(name + ".zip")) as zf:
        zf.extractall(tempdir)
        for var_name in ["usr", "pw", "key"]:
            for root, dirs, files in os.walk(tempdir):
                if var_name + ".pkl" in files:
                    file_name = os.path.join(root, var_name + ".pkl")
            out[var_name] = load_pkl(file_name)
    return out["usr"], out["pw"], out["key"]

## AcrGIS Online log-in
usl = None
if os.path.exists(dir_app("arc.zip")):
    while usl not in ["y", "n"]:
        use_saved_login = input("Do you want to use the last saved" +\
                                "ArcGIS Online login credentials (y/n):\n")
        usl = use_saved_login.casefold()
        if usl not in ["y", "n"]:
            print("Invalid input:", usl + "\nUse 'y' (yes) or 'n' (no).")

from cryptography.fernet import Fernet
if usl == "y":
    # read saved login credentials
    arcgis_usr_enc, arcgis_pw_enc, key = get_login("arc")
    fernet = Fernet(key)
    arcgis_usr = fernet.decrypt(arcgis_usr_enc).decode()
    arcgis_pw = fernet.decrypt(arcgis_pw_enc).decode()
else:
    # prompt user input (ArcGIS login credentials)
    arcgis_usr = input("Enter ArcGIS Online user name:\n")
    arcgis_pw = input("Enter ArcGIS Online password:\n")
    save_creds = input("Save login credentials? (y/n):\n")
    sl = save_creds.casefold()
    if sl == "y":
        # save encrypted login credentials
        key = Fernet.generate_key()
        fernet = Fernet(key)
        arcgis_usr_enc = fernet.encrypt(str(arcgis_usr).encode())
        arcgis_pw_enc = fernet.encrypt(str(arcgis_pw).encode())
        save_login("arc", arcgis_usr_enc, arcgis_pw_enc, key)

## ArcGIS log-in
import arcgis, math
from arcgis.gis import GIS
gis = GIS(None, arcgis_usr, arcgis_pw, verify_cert = False)
ArcGIS_content = input("Enter the name of the ArcGIS Online FeatureLayer " +\
                       "containing the observations to upload:\n")

## iNaturalist log-in
usl = None
if os.path.exists(dir_app("nat.zip")):
    while usl not in ["y", "n"]:
        use_saved_login = input("Do you want to use the last saved" +\
                                "iNaturalist login credentials (y/n):\n")
        usl = use_saved_login.casefold()
        if usl not in ["y", "n"]:
            print("Invalid input:", usl + "\nUse 'y' (yes) or 'n' (no).")

from cryptography.fernet import Fernet
if usl == "y":
    # read saved login credentials
    inat_usr_enc, inat_pw_enc, key = get_login("arc")
    fernet = Fernet(key)
    inat_usr = fernet.decrypt(inat_usr_enc).decode()
    inat_pw = fernet.decrypt(inat_pw_enc).decode()
else:
    # prompt user input (iNaturalist login credentials)
    inat_usr = input("Enter iNaturalist user name:\n")
    inat_pw = input("Enter iNaturalist password:\n")
    save_creds = input("Save login credentials? (y/n):\n")
    sl = save_creds.casefold()
    if sl == "y":
        # save login credentials
        key = Fernet.generate_key()
        fernet = Fernet(key)
        inat_usr_enc = fernet.encrypt(str(inat_usr).encode())
        inat_pw_enc = fernet.encrypt(str(inat_pw).encode())
        save_login("nat", inat_usr_enc, inat_pw_enc, key)

## get ArcGIS Online attachments
search_result = gis.content.search(ArcGIS_content, "Feature Layer")
item = search_result[0]
FeatureLayer = item.layers[0]
PointFeatures = FeatureLayer.query(where = "1=1", outFields = "objectid")
oidFieldName = "OBJECTID"
oids = [f.attributes[oidFieldName] for f in PointFeatures.features]

tmp_dir = tempfile.mkdtemp()
for objID in oids:
    attachments_list = FeatureLayer.attachments.get_list(oid = objID)
    if len(attachments_list) >= 1:
        for attmID in range(1, len(attachments_list) + 1):
            obj_dir = os.path.join(tmp_dir, str(objID))
            os.makedirs(obj_dir)
            FeatureLayer.attachments.download(oid = objID,\
                                              attachment_id = attmID,\
                                              save_path = obj_dir)

## iNaturalist log-in
from pyinaturalist import *
token = get_access_token(
    username = inat_usr,
    password = inat_pw
    )
