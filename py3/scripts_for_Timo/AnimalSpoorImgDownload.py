#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 19:28:47 2022

@author: manuel
Note: Script will only work for AnimalSpoor Layer, since Layer IDs are
explicit
"""
## ArcGIS log-in
import os
import arcgis, math
from arcgis.gis import GIS

ArcGIS_content = "AnimalSpoor"
tmp_dir ="C:\\Users\\Manuel\\Desktop\\test\\"

#gis = GIS(None, usr, pw, verify_cert = False)
gis = GIS("pro")

## get ArcGIS Online attachments
search_result = gis.content.search(ArcGIS_content, "Feature Layer")
item = search_result[0]
FeatureLayer = item.layers[0]
FeatureLayer = gis.content.get("0969767b3a7347dcab2cdcea36d33834").layers[0]
SpeciesTable = gis.content.get("55a51dd698224f85a9fdf6daa69b2cbd").tables[0]\
    .query(out_fields = "English_name", as_df = True)
speciesName = SpeciesTable["English_name"].tolist()
import re
speciesName = [re.sub(r"[^\w]", "", n) for n in speciesName]
#FeatureObjectIDs = FeatureLayer.query(where = "1=1", return_ids_only = True)
attachments = FeatureLayer.attachments.search(as_df = True)
FeatureObjectIDs = attachments["PARENTOBJECTID"].drop_duplicates().to_list()

import pickle
def load_pkl(var_name):
    with open(os.path.join(tmp_dir, var_name + ".pkl"), "rb") as f:
            return pickle.load(f)

def save_pkl(var):
    with open(os.path.join(tmp_dir, "progress" + ".pkl"),\
              "wb") as f:
        pickle.dump(var, f)

if os.path.isdir(os.path.join(tmp_dir, "save_progress.pkl")):
    download_completed = load_pkl("progress")
    FeatureObjectIDs = list(set(FeatureObjectIDs) - set(download_completed))
else:
    download_completed = []

without_id = 0
invalid_ID = 0
step = 1
step_max = len(FeatureObjectIDs)
for objID in FeatureObjectIDs:
    print("Checking object", str(step), "of", str(step_max) + "...")
    attachments_list = FeatureLayer.attachments.get_list(oid = objID)
    current_feature = FeatureLayer.query(where = "OBJECTID=" + str(objID))
    if len(attachments_list) > 0:
        df = current_feature.sdf
        comment = df["NotesComments"][0]
        full_comment = comment
        specID = df["Species_ID"][0]
        if specID is not None:
            specID = speciesName[int(specID) - 1]
        if comment is None:
            Warning("Missing ID for object " + str(objID))
            identifier = "No_ID_obj_"# + str(objID)
            without_id += 1
        else:
            comment = re.sub(r"[^\w]", "", df["NotesComments"][0])
            regEx = re.compile("B\dP\d\d\d\d", re.IGNORECASE)
            if regEx.search(comment) is not None:
                identifier = regEx.search(comment)[0]
            else:
                Warning("Invalid object ID: " + str(comment))
                identifier = "Invalid_ID_obj_"# + str(objID)
                invalid_ID += 1
        identifier = str(objID) + "_" + identifier + str(specID)
        obj_dir = os.path.join(tmp_dir, str(identifier))
        os.makedirs(obj_dir, exist_ok = True)
        with open(os.path.join(obj_dir, "comment.txt"), "w") as f:
            f.write(str(full_comment))
        attachment_IDs = [a["id"] for a in attachments_list]
        for attmID in attachment_IDs:
            FeatureLayer.attachments.download(oid = objID,\
                                              attachment_id = attmID,\
                                              save_path = obj_dir)
        download_completed.append(objID)
        save_pkl(download_completed)
        step += 1
print("Downloads completed.")
