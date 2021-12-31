#!/home/manuel/Nextcloud/LELE_2021/py3/venv/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 26 23:37:02 2021

@author: manuel
"""
import argparse
def parseArguments():
    parser = argparse.ArgumentParser()
    # Positional mandatory arguments
    parser.add_argument("iNatID", help = "Unique ID used in the " + \
                        "iNaturalist observation comment.", type = str)
    parser.add_argument("-pr", "--pr",\
                        help = "iNaturalist project ID.", type = str,\
                            default = "lele-herbs-and-grasses")
    args = parser.parse_args()
    return args
if __name__ == "__main__":
    # Parse the arguments
    args = parseArguments()

iNaturalist_ID = args.iNatID
project_ID = args.pr

from pyinaturalist import *
project_observations = get_observations(project_id = project_ID)
observations = project_observations["results"]

target = []
for obs in observations:
    if obs["description"] is not None and obs["description"] != "":
        if obs["description"] in iNaturalist_ID:
            target.append(obs)

if len(target) < 1:
    print("Found no observation containing iNaturalist-ID '" + \
          iNaturalist_ID + "' in description.")
elif len(target) == 1:
    observation = target[0]
    print("Found one observation containing iNaturalist-ID '" + \
          iNaturalist_ID + "' in description:")
    print("Current taxon:", observation["taxon"]["name"])
    print("Observation url:", observation["uri"])
    similar_observations = []
    for sim_obs in observations:
        if sim_obs["taxon"] == observation["taxon"]["name"]:
            similar_observations.append(sim_obs)
    if len(similar_observations) > 1:
        print("Found the following observations currently " + \
              "listed under the same taxon (" + \
                  observation["taxon"]["name"] + "):")
        for sim_observation in similar_observations:
            print("iNaturalist ID:", sim_observation["description"])
            print("url:", sim_observation["uri"])
else:
    print("Found multiple observations containing iNaturalist-ID '" + \
          iNaturalist_ID + "' in description':")
    for observation in target:
        print("Taxon:", observation["taxon"]["name"])
        print("url:", observation["uri"])