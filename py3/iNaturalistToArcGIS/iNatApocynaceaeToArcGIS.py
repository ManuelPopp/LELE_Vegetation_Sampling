#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 17:41:57 2022

@author: manuel
"""
import os
import arcpy
from arcgis.gis import GIS
from arcgis import features
#from arcgis.features import FeatureLayer
#from arcgis.features import FeatureLayerCollection
arcpy.env.overwriteOutput = True

#gis = GIS(None, usr, pw, verify_cert = False)
gis = GIS("pro")
ArcGIS_content = "Apocynaceae"
project_ID = "flora-and-fauna-of-lapalala-wilderness"

## get original ArcGIS Online layer
try:
    FeatureLayer = gis.content.get(ArcGIS_content).layers[0]
except:
    print("Searching for ArcGIS content", ArcGIS_content)
    search_result = gis.content.search(ArcGIS_content, "Feature Layer")
    FeatureLayer = search_result[0].layers[0]
else:
    FeatureLayer = gis.content.get(ArcGIS_content).layers[0]

## get iNaturalist data by running script in Python environment with INat API
from pyinaturalist import *
observations = []
page_results = [1]
p = 1
while len(page_results) > 0:
    page = get_observations(project_id = project_ID,\
                            taxon_id = 47362, per_page = 30, page = p)
    page_results = page["results"]
    observations += page_results
    p += 1

iNat_data = list()
for obs in observations:
    if obs["description"] is not None and obs["description"] != "":
        lat = obs["location"][0]
        lon = obs["location"][1]
        id = obs["id"]
        taxon = obs["taxon"]["name"]
        url = obs["uri"]
        quality = obs["quality_grade"]
        iNat_data.append([lat, lon, taxon, url, quality, id])

## create dataframe
import pandas as pd
data = pd.DataFrame(iNat_data,\
                    columns = ["Latitude", "Longitude", "Taxon", "url",\
                               "Quality", "iNaturalist_ID"])

current_features = FeatureLayer.query(out_fields = "iNaturalist_ID",\
                                      as_df = True)
current_features = current_features["iNaturalist_ID"].tolist()
new_observation_IDs = list(set(data["iNaturalist_ID"]) - set(current_features))
new_data = data.loc[data["iNaturalist_ID"].isin(new_observation_IDs)]

## add columns for ArcGIS Online Feature Layer
AGOL_columns = ["OBJECTID", "esrignss_speed", "esrignss_direction",\
                "esrisnsr_azimuth", "esrignss_positionsourcetype",
                "esrignss_receiver", "esrignss_h_rms", "esrignss_v_rms"
                "esrignss_latitude", "esrignss_longitude",
                "esrignss_altitude", "esrignss_pdop", "esrignss_hdop"
                "esrignss_vdop", "esrignss_fixtype", "esrignss_correctionage"]

for column, index in zip(AGOL_columns, range(0, len(AGOL_columns) + 1)):
    new_data.insert(index, column, [None] * len(new_data))

## save as temporary csv
import tempfile
tmp_dir = tempfile.mkdtemp()
tmp_csv = os.path.join(tmp_dir, "tmp_features.csv")

data.to_csv(tmp_csv)

## delete old csv object
tmp_csv_old = gis.content.search("title: tmp_csv_feature_upload")
for item in tmp_csv_old:
    if item:
        item.delete()

## add updates as csv in AGOL
item_props = {"title" : "tmp_csv_feature_upload",
              "description" : "Temporary csv for feature update",
              "tags" : "csv, update"}

csv_item = gis.content.add(item_properties = item_props, data = tmp_csv,\
                           folder = "tmp_files")

## truncate FeatureLayer
#FeatureLayer.manager.truncate()

## append observations to FeatureLayer
param_0 = gis.content.analyze(item = csv_item.id)
FeatureLayer.append(item_id = csv_item.id, upload_format = "csv",\
                    source_info = param_0["publishParameters"])

# remove temp dir
import shutil
shutil.rmtree(tmp_dir)

## delete used csv object
tmp_csv_old = gis.content.search("title: tmp_csv_feature_upload")
for item in tmp_csv_old:
    if item:
        item.delete()
