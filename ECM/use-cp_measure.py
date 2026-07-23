import numpy as np
from cp_measure.featurizer import featurize
from tifffile import imread
import matplotlib.pyplot as plt 
import pandas as pd

# img = np.random.default_rng(42).random((1, 240, 240))
# mask = np.zeros((1, 240, 240), dtype=np.int32)
# mask[0, 50:100, 50:100] = 1
# mask[0, 150:200, 150:200] = 2
# data, columns, rows = featurize(img, mask)
# df1 = pd.DataFrame(data=data, columns=columns, index=rows)

col_stats = ['Area', 'BoundingBoxArea', 'ConvexArea', 'EquivalentDiameter', 'Perimeter', 'MajorAxisLength', 
             'MinorAxisLength', 'Eccentricity', 'Orientation', 'Center_X', 'Center_Y', 'BoundingBoxMinimum_X', 
             'BoundingBoxMaximum_X', 'BoundingBoxMinimum_Y', 'BoundingBoxMaximum_Y', 'FormFactor', 'Extent', 'Solidity', 
             'Compactness', 'EulerNumber', 'MaximumRadius', 'MeanRadius', 'MedianRadius', 'FilledArea', 'SpatialMoment_0_0', 
             'SpatialMoment_0_1', 'SpatialMoment_0_2', 'SpatialMoment_0_3', 'SpatialMoment_1_0', 'SpatialMoment_1_1', 
             'SpatialMoment_1_2', 'SpatialMoment_1_3', 'SpatialMoment_2_0', 'SpatialMoment_2_1', 'SpatialMoment_2_2', 
             'SpatialMoment_2_3', 'CentralMoment_0_0', 'CentralMoment_0_1', 'CentralMoment_0_2', 'CentralMoment_0_3', 
             'CentralMoment_1_0', 'CentralMoment_1_1', 'CentralMoment_1_2', 'CentralMoment_1_3', 'CentralMoment_2_0', 
             'CentralMoment_2_1', 'CentralMoment_2_2', 'CentralMoment_2_3', 'NormalizedMoment_0_0', 'NormalizedMoment_0_1', 
             'NormalizedMoment_0_2', 'NormalizedMoment_0_3', 'NormalizedMoment_1_0', 'NormalizedMoment_1_1', 'NormalizedMoment_1_2', 
             'NormalizedMoment_1_3', 'NormalizedMoment_2_0', 'NormalizedMoment_2_1', 'NormalizedMoment_2_2', 'NormalizedMoment_2_3', 
             'NormalizedMoment_3_0', 'NormalizedMoment_3_1', 'NormalizedMoment_3_2', 'NormalizedMoment_3_3', 'HuMoment_0', 
             'HuMoment_1', 'HuMoment_2', 'HuMoment_3', 'HuMoment_4', 'HuMoment_5', 'HuMoment_6', 'InertiaTensor_0_0', 
             'InertiaTensor_0_1', 'InertiaTensor_1_0', 'InertiaTensor_1_1', 'InertiaTensorEigenvalues_0', 'InertiaTensorEigenvalues_1', 
             'PerimeterCrofton', 'Zernike_0_0', 'Zernike_1_1', 'Zernike_2_0', 'Zernike_2_2', 'Zernike_3_1', 'Zernike_3_3', 
             'Zernike_4_0', 'Zernike_4_2', 'Zernike_4_4', 'Zernike_5_1', 'Zernike_5_3', 'Zernike_5_5', 'Zernike_6_0', 
             'Zernike_6_2', 'Zernike_6_4', 'Zernike_6_6', 'Zernike_7_1', 'Zernike_7_3', 'Zernike_7_5', 'Zernike_7_7', 
             'Zernike_8_0', 'Zernike_8_2', 'Zernike_8_4', 'Zernike_8_6', 'Zernike_8_8', 'Zernike_9_1', 'Zernike_9_3', 
             'Zernike_9_5', 'Zernike_9_7', 'Zernike_9_9', 'MinFeretDiameter', 'MaxFeretDiameter', 'Intensity_IntegratedIntensity__ch0', 
             'Intensity_MeanIntensity__ch0', 'Intensity_StdIntensity__ch0', 'Intensity_MinIntensity__ch0', 'Intensity_MaxIntensity__ch0', 
             'Intensity_MassDisplacement__ch0', 'Intensity_LowerQuartileIntensity__ch0', 'Intensity_MedianIntensity__ch0', 
             'Intensity_MADIntensity__ch0', 'Intensity_UpperQuartileIntensity__ch0', 'Location_CenterMassIntensity_X__ch0', 
             'Location_CenterMassIntensity_Y__ch0', 'Location_CenterMassIntensity_Z__ch0', 'Location_MaxIntensity_X__ch0', 
             'Location_MaxIntensity_Y__ch0', 'Location_MaxIntensity_Z__ch0', 'Intensity_IntegratedIntensityEdge__ch0', 
             'Intensity_MeanIntensityEdge__ch0', 'Intensity_StdIntensityEdge__ch0', 'Intensity_MinIntensityEdge__ch0', 
             'Intensity_MaxIntensityEdge__ch0', 'AngularSecondMoment_3_00_256__ch0', 'Contrast_3_00_256__ch0', 'Correlation_3_00_256__ch0', 
             'Variance_3_00_256__ch0', 'InverseDifferenceMoment_3_00_256__ch0', 'SumAverage_3_00_256__ch0', 'SumVariance_3_00_256__ch0', 
             'SumEntropy_3_00_256__ch0', 'Entropy_3_00_256__ch0', 'DifferenceVariance_3_00_256__ch0', 'DifferenceEntropy_3_00_256__ch0', 
             'InfoMeas1_3_00_256__ch0', 'InfoMeas2_3_00_256__ch0', 'AngularSecondMoment_3_01_256__ch0', 'Contrast_3_01_256__ch0', 
             'Correlation_3_01_256__ch0', 'Variance_3_01_256__ch0', 'InverseDifferenceMoment_3_01_256__ch0', 'SumAverage_3_01_256__ch0', 
             'SumVariance_3_01_256__ch0', 'SumEntropy_3_01_256__ch0', 'Entropy_3_01_256__ch0', 'DifferenceVariance_3_01_256__ch0', 
             'DifferenceEntropy_3_01_256__ch0', 'InfoMeas1_3_01_256__ch0', 'InfoMeas2_3_01_256__ch0', 'AngularSecondMoment_3_02_256__ch0',
               'Contrast_3_02_256__ch0', 'Correlation_3_02_256__ch0', 'Variance_3_02_256__ch0', 'InverseDifferenceMoment_3_02_256__ch0',
                 'SumAverage_3_02_256__ch0', 'SumVariance_3_02_256__ch0', 'SumEntropy_3_02_256__ch0', 'Entropy_3_02_256__ch0', 
                 'DifferenceVariance_3_02_256__ch0', 'DifferenceEntropy_3_02_256__ch0', 'InfoMeas1_3_02_256__ch0', 'InfoMeas2_3_02_256__ch0', 
                 'AngularSecondMoment_3_03_256__ch0', 'Contrast_3_03_256__ch0', 'Correlation_3_03_256__ch0', 'Variance_3_03_256__ch0', 
                 'InverseDifferenceMoment_3_03_256__ch0', 'SumAverage_3_03_256__ch0', 'SumVariance_3_03_256__ch0', 'SumEntropy_3_03_256__ch0', 
                 'Entropy_3_03_256__ch0', 'DifferenceVariance_3_03_256__ch0', 'DifferenceEntropy_3_03_256__ch0', 'InfoMeas1_3_03_256__ch0', 
                 'InfoMeas2_3_03_256__ch0', 'Granularity_1__ch0', 'Granularity_2__ch0', 'Granularity_3__ch0', 'Granularity_4__ch0', 'Granularity_5__ch0', 
                 'Granularity_6__ch0', 'Granularity_7__ch0', 'Granularity_8__ch0', 'Granularity_9__ch0', 'Granularity_10__ch0', 'Granularity_11__ch0', 
                 'Granularity_12__ch0', 'Granularity_13__ch0', 'Granularity_14__ch0', 'Granularity_15__ch0', 'Granularity_16__ch0', 
                 'RadialDistribution_FracAtD_1of4__ch0', 'RadialDistribution_MeanFrac_1of4__ch0', 'RadialDistribution_RadialCV_1of4__ch0', 
                 'RadialDistribution_FracAtD_2of4__ch0', 'RadialDistribution_MeanFrac_2of4__ch0', 'RadialDistribution_RadialCV_2of4__ch0', 
                 'RadialDistribution_FracAtD_3of4__ch0', 'RadialDistribution_MeanFrac_3of4__ch0', 'RadialDistribution_RadialCV_3of4__ch0', 
                 'RadialDistribution_FracAtD_4of4__ch0', 'RadialDistribution_MeanFrac_4of4__ch0', 'RadialDistribution_RadialCV_4of4__ch0', 
                 'RadialDistribution_ZernikeMagnitude_0_0__ch0', 'RadialDistribution_ZernikePhase_0_0__ch0', 'RadialDistribution_ZernikeMagnitude_1_1__ch0',
                   'RadialDistribution_ZernikePhase_1_1__ch0', 'RadialDistribution_ZernikeMagnitude_2_0__ch0', 'RadialDistribution_ZernikePhase_2_0__ch0', 
                   'RadialDistribution_ZernikeMagnitude_2_2__ch0', 'RadialDistribution_ZernikePhase_2_2__ch0', 'RadialDistribution_ZernikeMagnitude_3_1__ch0',
                     'RadialDistribution_ZernikePhase_3_1__ch0', 'RadialDistribution_ZernikeMagnitude_3_3__ch0', 'RadialDistribution_ZernikePhase_3_3__ch0', 
                     'RadialDistribution_ZernikeMagnitude_4_0__ch0', 'RadialDistribution_ZernikePhase_4_0__ch0', 'RadialDistribution_ZernikeMagnitude_4_2__ch0', 
                     'RadialDistribution_ZernikePhase_4_2__ch0', 'RadialDistribution_ZernikeMagnitude_4_4__ch0', 'RadialDistribution_ZernikePhase_4_4__ch0', 
                     'RadialDistribution_ZernikeMagnitude_5_1__ch0', 'RadialDistribution_ZernikePhase_5_1__ch0', 'RadialDistribution_ZernikeMagnitude_5_3__ch0',
                       'RadialDistribution_ZernikePhase_5_3__ch0', 'RadialDistribution_ZernikeMagnitude_5_5__ch0', 'RadialDistribution_ZernikePhase_5_5__ch0',
                         'RadialDistribution_ZernikeMagnitude_6_0__ch0', 'RadialDistribution_ZernikePhase_6_0__ch0', 'RadialDistribution_ZernikeMagnitude_6_2__ch0',
                           'RadialDistribution_ZernikePhase_6_2__ch0', 'RadialDistribution_ZernikeMagnitude_6_4__ch0', 'RadialDistribution_ZernikePhase_6_4__ch0', 
                           'RadialDistribution_ZernikeMagnitude_6_6__ch0', 'RadialDistribution_ZernikePhase_6_6__ch0', 'RadialDistribution_ZernikeMagnitude_7_1__ch0', 
                           'RadialDistribution_ZernikePhase_7_1__ch0', 'RadialDistribution_ZernikeMagnitude_7_3__ch0', 'RadialDistribution_ZernikePhase_7_3__ch0', 
                           'RadialDistribution_ZernikeMagnitude_7_5__ch0', 'RadialDistribution_ZernikePhase_7_5__ch0', 'RadialDistribution_ZernikeMagnitude_7_7__ch0',
                             'RadialDistribution_ZernikePhase_7_7__ch0', 'RadialDistribution_ZernikeMagnitude_8_0__ch0', 
                             'RadialDistribution_ZernikePhase_8_0__ch0', 'RadialDistribution_ZernikeMagnitude_8_2__ch0', 
                             'RadialDistribution_ZernikePhase_8_2__ch0', 'RadialDistribution_ZernikeMagnitude_8_4__ch0', 
                             'RadialDistribution_ZernikePhase_8_4__ch0', 'RadialDistribution_ZernikeMagnitude_8_6__ch0', 
                             'RadialDistribution_ZernikePhase_8_6__ch0', 'RadialDistribution_ZernikeMagnitude_8_8__ch0', 
                             'RadialDistribution_ZernikePhase_8_8__ch0', 'RadialDistribution_ZernikeMagnitude_9_1__ch0', 
                             'RadialDistribution_ZernikePhase_9_1__ch0', 'RadialDistribution_ZernikeMagnitude_9_3__ch0', 
                             'RadialDistribution_ZernikePhase_9_3__ch0', 'RadialDistribution_ZernikeMagnitude_9_5__ch0', 
                             'RadialDistribution_ZernikePhase_9_5__ch0', 'RadialDistribution_ZernikeMagnitude_9_7__ch0', 
                             'RadialDistribution_ZernikePhase_9_7__ch0', 'RadialDistribution_ZernikeMagnitude_9_9__ch0', 
                             'RadialDistribution_ZernikePhase_9_9__ch0']
imgs = imread("./data/original/M20_A1_6__stack_correction.tif")
masks = np.load("./data/stack_masks/M20_A1_6__stack_correction.npy")

col_stats.append("Frame")
df = pd.DataFrame(columns=col_stats)
for slice_id in range(masks.shape[0]):
    print(f"{slice_id}/{masks.shape[0]}")
    img = imgs[[0],:,:]
    mask = masks[[0],:,:]
    data, columns, rows = featurize(img, mask)
    df_slice = pd.DataFrame(data=data, columns=columns, index=[row[2] for row in rows])
    df_slice["Frame"] = slice_id
    df = pd.concat([df, df_slice])

df.to_csv('./data/cp_measure_csv/M20_A1_6__stack_correction.csv')

df.loc[(df["Center_X"]>=675.490291) & (df["Center_X"]<=675.490292), "Center_X"]
df.loc[(df["Center_X"]==675.490291), "Center_X"]
df.loc[44, "Center_X"]

import btrack
with btrack.io.HDF5FileHandler(
  f"./data/btrack_hdf5/M20_A1_6__stack_correction.h5", 'r', obj_type='obj_type_1') as reader:
  tracks = reader.tracks

from skimage.measure import regionprops
prop = regionprops(masks[0])[0]
imgs = imread("./data/original/M20_A1_6__stack_correction.tif")
box_coord = prop.bbox
box = imgs[0,:,:]
data, columns, rows = featurize(img, mask)

from cp_measure.core.measureobjectsizeshape import get_feret
mask = np.zeros((50, 50), dtype=np.int32)
mask[5:-6, 5:-6] = 1
get_feret(mask, None)


coordinates = prop.coords.tolist()
coordinates.append(prop.coords[0])
plt.plot(coordinates)
plt.show()