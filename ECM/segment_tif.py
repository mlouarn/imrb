from cellpose import models, io, plot
import matplotlib.pyplot as plt
import pickle
from glob import glob
import numpy as np
from pathlib import Path
import tifffile
from PIL import Image, ImageSequence, TiffImagePlugin
from roifile import roiread
import zipfile
import os

model = models.CellposeModel(gpu=True, pretrained_model = "cpsam_v2")
SLICE_ROIS_DIRPATH = "data/slice_rois" # used to be data/rois
STACK_ROIS_DIRPATH = "data/stack_rois" # used to be data/zip

def segment_img(img:np.ndarray, basename:str):
    result = model.eval(img, flow_threshold=0.2, cellprob_threshold=0.0, niter=None)
    with open(f"data/slice_pickle/{basename}.pkl", 'wb') as picklefile: # used to be pickle instead of slice_pickle
        pickle.dump(result, picklefile)
    masks, flows, styles = result
    io.save_rois(masks, f"{SLICE_ROIS_DIRPATH}/{basename}") # saved as basename_rois.zip
    # fig = plt.figure(figsize=(24,10))
    # plot.show_segmentation(fig, img, masks, flows[0])
    # plt.tight_layout()
    # plt.savefig(f"data/panel/{basename}_panel.png")
    # plt.close()
    return masks


for tifpath in glob("data/original/*.tif"):
    pklpath = f"data/stack_pickle/{Path(tifpath).stem}.pkl"
    if not os.path.exists(pklpath): 
        print(tifpath)
        stack = io.imread(tifpath) # shape (image number, x, y) grey no additional channels
        list_masks = [] 
        for n in range(stack.shape[0]):
            slice = stack[n,:,:]
            basename = f"{Path(tifpath).stem}_{n:0>3}"
            print(basename)
            masks = segment_img(slice, basename) # saves slice : pickle and rois.zip
            list_masks.append(masks)
            ## commented cause no need to store the individual images with the overlay
            # roi = roiread(f"data/rois/{basename}_rois.zip")
            # overlays = [r.tobytes() for r in roi]
            # tifffile.imwrite(f"data/slice/{basename}.tif", slice, imagej=True, metadata={"Overlays":overlays})
        stack_masks = np.stack(list_masks)
        with open(pklpath, "wb") as pklfile:
            pickle.dump(stack_masks, pklfile)

        ## aggregate all roi zips of slides into one for a film (to import in ImageJ)
        print("aggregate")
        basename = Path(tifpath).stem
        zipoutput_path = f"{STACK_ROIS_DIRPATH}/{basename}.zip" # used to be data/zip
        size_stack = io.imread(tifpath).shape[0]

        with zipfile.ZipFile(zipoutput_path, "w") as zipout:
            for slice_id in range(size_stack):
                roispath = f"{SLICE_ROIS_DIRPATH}/{basename}_{slice_id:0>3}_rois.zip"
                with zipfile.ZipFile(roispath, "r") as zipinput:
                    for roiname in zipinput.namelist():
                        with zipinput.open(roiname) as roifile:
                            content = roifile.read()
                            roi_id = roiname[:-4]
                            roi_newname= f"{slice_id+1:0>4}-{roi_id:0>4}-{roi_id:0>4}.roi"
                            zipout.writestr(roi_newname, content)


# take the stack_pickle and output the np.array in a tif image in stack_masks
# for tifpath in glob("./data/stack_pickle/M20_A1_6__stack_correction.pkl"):
tifpath = "./data/stack_pickle/M20_A1_6__stack_correction.pkl"
with open(tifpath, "rb") as pklfile:
    masks = pickle.load(pklfile)
print(masks.shape)
plt.imshow(masks[0,:,:])
plt.show()
npypath = f"./data/stack_masks/{Path(tifpath).stem}.npy"
np.save(npypath, masks)
# with tifffile.TiffWriter(tifmaskpath) as tifhandle:
#     for slice_id in range(masks.shape[0]):
#         mask = masks[slice_id,:,:]
#         tifhandle.write(mask)
# tifffile.imwrite(f"./data/stack_masks/{Path(tifpath).stem}.tif", masks)


# aggregate all slice pkl files into one stack of cellpose labels, save to pkl in stack_pickle
# for tifpath in glob("./data/original/M20_A1_6__stack_correction.tif"):
#     base_pklpath = Path(tifpath).stem
#     list_masks = []
#     for pklpath in glob(f"data/slice_pickle/{base_pklpath}_*"):
#         with open(pklpath, "rb") as pklfile:
#             list_masks.append(pickle.load(pklfile)[0])
#     stack_masks = np.stack(list_masks)
#     with open(f"data/stack_pickle/{base_pklpath}.pkl", "wb") as pklfile:
#         pickle.dump(stack_masks, pklfile)  




# save tif stack with metadata for each slice -> not useful we only aggregate the rois
# def save_tif_stack(tifpath:str):
#     im = Image.open(tifpath)
#     frames = []
#     for i, frame in enumerate(ImageSequence.Iterator(im)):
#         basepath = f"{Path(tifpath).stem}_{i:0>3}"
#         img = Image.open(f"data/slice/{basepath}.tif")
#         frame = frame.convert(frame.mode)
#         info = img.tag_v2
#         frame.encoderinfo = {'tiffinfo': info}
#         frames.append(frame)
#     outpath = f"data/stack/{Path(tifpath).stem}_masks.tif"
#     with open(outpath, "w+b") as fp:
#         with TiffImagePlugin.AppendingTiffWriter(fp) as tf:
#             for frame in frames:
#                 frame.encoderconfig = ()
#                 TiffImagePlugin._save(frame, tf, outpath)
#                 tf.newFrame()

# for tifpath in glob("data/M*tif"):
#     save_tif_stack(tifpath)

