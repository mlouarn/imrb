# https://github.com/lowe-lab-ucl/btrack-examples/blob/main/examples/cell_config.json
# https://btrack.readthedocs.io/en/latest/user_guide/configuration.html
# 

import btrack
import pickle
import numpy as np
import matplotlib.pyplot as plt
# from cv2 import imread # doesn't work the same way
# from cellpose.io import imread
from tifffile import imread
from pathlib import Path
import napari
from pathlib import Path
from skimage.measure import regionprops
import pandas as pd

from collections import Counter

pklpath = "./data/stack_pickle/M76_A1_1_0_stack_correction.pkl"
pklpath = "./data/stack_pickle/M20_A1_7__stack_correction.pkl"
pklpath = "./data/stack_pickle/M20_A1_6__stack_correction.pkl" # use M20_A1_6 json

with open(pklpath, "rb") as pklfile:
    masks = pickle.load(pklfile)

## add a unique label for each mask, per frame
detections = []
for t, mask in enumerate(masks):
    props = regionprops(mask)
    for prop in props:
        detections.append({
            "t":t,
            "x":prop.centroid[1],
            "y":prop.centroid[0],
            "maskID":prop.label # do not add the keyword "label" (btrack uses it for something else), renamed to maskID https://github.com/quantumjot/btrack/issues/187
            })
df_detections = pd.DataFrame(detections)
objects_labeled = btrack.io.localizations_to_objects(df_detections)

# get tracks
## generate
objects = btrack.utils.segmentation_to_objects(masks, properties = ("area",))
tracker = btrack.BayesianTracker()
tracker.configure("./btrack_config_M20_A1_6.json")
tracker.append(objects_labeled)
tracker.volume = ((0, 1200), (0, 1600))
tracker.track_interactive(step_size=100)
tracker.optimize()
tracks = tracker.tracks

tracker.export(f"data/btrack_hdf5/{Path(pklpath).stem}.h5", obj_type="obj_type_1")
pytracks = [track.to_dict() for track in tracks]
with open(f"data/btrack_dict/{Path(pklpath).stem}.pkl", "wb") as pklfile:
    pickle.dump(pytracks, pklfile)

btrack.io.export_CSV(f"data/btrack_csv/{Path(pklpath).stem}.csv", tracks)

## or import
with btrack.io.HDF5FileHandler(
  f"./data/btrack_hdf5/{Path(pklpath).stem}.h5", 'r', obj_type='obj_type_1'
) as reader:
  tracks = reader.tracks

def view_napari():
    viewer = napari.Viewer()
    tifpath = f"data/original/{Path(pklpath).stem}.tif"
    img = imread(tifpath)
    data, properties, graph = btrack.utils.tracks_to_napari(tracks) # data, properties, graph = tracker.to_napari()
    viewer.add_image(img)
    viewer.add_labels(cellpose_masks)
    viewer.add_tracks(data, properties=properties, graph=graph) #
view_napari()

## use napari to see the tracks
viewer = napari.Viewer()
tifpath = f"data/original/{Path(pklpath).stem}.tif"
img = imread(tifpath)
data, properties, graph = btrack.utils.tracks_to_napari(tracks) # data, properties, graph = tracker.to_napari()
viewer.add_image(img)
viewer.add_labels(masks)
viewer.add_tracks(data, properties=properties, graph=graph) #

import pandas as pd
df_tracks = pd.DataFrame(columns=["ID","t","x","y","parent","root"])
for track in tracks:
    for mask_id in range(len(track.to_dict().t)):
        t = track.t[mask_id]
        x = track.x[mask_id]
        y = track.y[mask_id]
        df_tracks.loc[len(df_tracks)] = [track.ID, t, x, y, track.parent, track.root]






# for track in tracks:
#     plt.plot(track.x, track.y)
# plt.show()




import pycellin as pc
import pycellin.graph.properties as pc_core
import collections

with open("data/btrack_tracks.pkl", "rb") as pklfile:
    pytracks = pickle.load(pklfile)

list_t = [d["t"] for d in pytracks]
list_maxt = [max(t) for t in list_t]
list_mint = [min(t) for t in list_t]

def btrack_to_pycellin(btracks:list[collections.OrderedDict]):
    lineages = {}
    for track in btracks:
        lin = pc.CellLineage()
        trackid = track["ID"]
        nodes = []
        for i in range(len(track.t)): # (1, {"frame": 0, "cell_ID": 1, "lineage_ID": 1, "cell_x": 10, "cell_y": 0}),
            node = (i+1, {"frame": track.t, "lineage_ID": trackid, "cell_x": track.x[i], "cell_y": track.y})



pytracks[0]

lin1 = pc.CellLineage()  # division
lin1.add_nodes_from(
    [
        (1, {"frame": 0, "cell_ID": 1, "lineage_ID": 1, "cell_x": 10, "cell_y": 0}),
        (2, {"frame": 1, "cell_ID": 2, "lineage_ID": 1, "cell_x": 20, "cell_y": 0}),
        (3, {"frame": 2, "cell_ID": 3, "lineage_ID": 1, "cell_x": 30, "cell_y": 0}),
        (4, {"frame": 2, "cell_ID": 4, "lineage_ID": 1, "cell_x": 40, "cell_y": 0}),
        (5, {"frame": 3, "cell_ID": 5, "lineage_ID": 1, "cell_x": 50, "cell_y": 0}),
        (6, {"frame": 4, "cell_ID": 6, "lineage_ID": 1, "cell_x": 60, "cell_y": 0}),
    ]
)
lin1.add_edges_from([(1, 2), (2, 3), (2, 4), (3, 5), (5, 6)])
lin1.graph["lineage_ID"] = 1

lin2 = pc.CellLineage()  # gap
lin2.add_nodes_from(
    [
        (7, {"frame": 0, "cell_ID": 7, "lineage_ID": 2, "cell_x": 0, "cell_y": 10}),
        (8, {"frame": 1, "cell_ID": 8, "lineage_ID": 2, "cell_x": 0, "cell_y": 20}),
        (9, {"frame": 3, "cell_ID": 9, "lineage_ID": 2, "cell_x": 0, "cell_y": 30}),
    ]
)
lin2.add_edges_from([(7, 8), (8, 9)])
lin2.graph["lineage_ID"] = 2

lin3 = pc.CellLineage()  # single-node lineage
lin3.add_node(10, frame=0, cell_ID=10, lineage_ID=3, cell_x=0, cell_y=0)
lin3.graph["lineage_ID"] = 3


model = pc.Model(
    data=pc.Data({1: lin1, 2: lin2, 3: lin3}),
    props_metadata=pc.PropsMetadata(
        {
            "frame": pc_core.create_frame_property(),
            "cell_x": pc_core.create_cell_coord_property(axis="x", unit="pixel"),
            "cell_y": pc_core.create_cell_coord_property(axis="y", unit="pixel"),
            "cell_ID": pc_core.create_cell_id_property(),
            "lineage_ID": pc_core.create_lineage_id_property(),
        }
    ),
    reference_time_property="frame",
)
print(model)

model.add_cycle_data()

model.model_metadata.get_standard_metadata()
model.model_metadata.name = "test_model_from_scratch"
model.model_metadata.provenance = "LX"
model.model_metadata.time_unit = "frame"
model.model_metadata.pixel_width = 1
model.model_metadata.pixel_height = 1
model.model_metadata.space_unit = "pixel"

model.model_metadata.segmentation = "ilastik pixel classifier"
model.model_metadata.tracking = "manual"
model.model_metadata.objective = "100X oil"

new_lin_id = model.add_lineage()

# model.add_cell(lid=1, cid=1, time_value=0)
# model.add_cell(lid=1, cid=2, time_value=1)
# model.add_cell(lid=1, cid=3, time_value=2)
# model.add_cell(lid=1, cid=4, time_value=2)
# model.add_cell(lid=1, cid=5, time_value=3)

# model.add_link(source_cid=1, target_cid=2, source_lid=1)
# model.add_link(source_cid=2, target_cid=3, source_lid=1)
# model.add_link(source_cid=2, target_cid=4, source_lid=1)
# model.add_link(source_cid=3, target_cid=5, source_lid=1)

model.set_time_step()

new_lin = model.get_cell_lineage_from_ID(1)
new_lin.plot()

pc.export_TrackMate_XML(model, "data/btrack_toy.xml",
                        units={
                            "space": model.get_space_unit(),
                            "time": model.get_time_unit(),
                        },
                        propagate_cycle_props=True,)