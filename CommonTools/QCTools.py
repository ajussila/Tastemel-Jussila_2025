'''
This package is focused on taking the tools which we have developed in the Ren lab for QC on imaging data,
specicically chromatin tracing and applying it to the outputs of an experiment. IT draws heavily from the scripts we
adopted from Bogdan Bintu's work in the 2018 Science paper.
'''

# external packages
import sys, os, glob
import pickle
import numpy as np
import matplotlib.pyplot as plt
from tqdm.notebook import tqdm
#Internal packages
#add path
workbookDir = os.getcwd()
sys.path.append(r'C:\Users\Administrator\Analysis\ChromatinTracingPipeline_Underdev\CommonTools')
import IOTools as io
import FittingTools as ft
import AlignmentTools as at

'''
This function takes an image from the grab_block function containing only one peak, and computes a local signal
to noise ratio. Here, I am defining the signal to be the median value of the top 10 pixels and the noise to be the median pixel value in the
image given image. The ratio is of the mean signal minus the mean noise, divided by the noise standard deviation.
'''
def single_spot_local_snr(data, max_spot_area=500):

    # flatten the data
    data_2d = list(map(np.ravel, data))
    data_1d = np.ravel(data_2d)
    signal_absolute = single_spot_brightness(data) # get the peak height from top x pixels

    # compute the noise floor as the median pizel value.
    noise_inds = np.argsort(data_1d)[:-max_spot_area] # assuming the spot area won't significantly exceed max_spot_area pixels.
    noise_brightness = list([data_1d[i] for i in noise_inds])
    noise_floor = np.nanmedian(noise_brightness)

    #subtract the noise floor from the signal for a more accurate height of the peak.
    signal_mean = signal_absolute-noise_floor

    # compute the stdev in the noise (all data minus the top 80 px as default).
    noise_stdev = np.nanstd(noise_brightness)

    # return the ratio of signal mean over deviation in the noise.
    return round(signal_mean/noise_stdev, 3)

'''
This function takes an image from the grab_block function containing only one peak and identifies the maximum height of the peak
by ranking the pixel values and taking the median of the top 10 brightest spots.
'''
def single_spot_brightness(data, pixels=50):

    # flatten the data and then take the top 10 brightest pixels. Use their median value as the signal value.
    data_2d = list(map(np.ravel, data))
    data_1d = np.ravel(data_2d)

    # take the median of the top 30 pixels. Median will be robust to artifacts muddying the top values.
    top50_brightest = list(np.argsort(data_1d)[-pixels:])
    signal_absolute = np.nanmedian( list([data_1d[i] for i in top50_brightest] ))

    return signal_absolute

'''
This function takes an image from the grab_block function containing only one peak andd fits a gaussian to compute the area of the
peak in the center of this image. This will be done by computing walking outward from the center of the gaussian until we hit
the half max value in 3 directions, then taking the average to compute an average radius.
'''
def single_spot_area(data, pixels=50, max_spot_area=500):

    #get the half maximum value
    brightness_abs = single_spot_brightness(data, pixels)

    # flatten the data and then take the top 10 brightest pixels. Use their median value as the signal value.
    data_2d = list(map(np.ravel, data))
    data_1d = np.ravel(data_2d)

    # compute the noise floor as the median pizel value.
    noise_brightness = np.argsort(data_1d)[:-max_spot_area] # assuming the spot area won't significantly exceed max_spot_area pixels.
    noise_floor = np.nanmedian(noise_brightness)

    #subtract the noise floor from the signal for a more accurate height of the peak.
    signal_mean = brightness_abs-noise_floor

    # half max is the half height of the peak above the noise floor of the gaussian.
    half_max = signal_mean/2+noise_floor

    area = len(np.where(data_1d > half_max)[0])

    return area

'''
This function is a supporting function for QC stuff which takes some simple experiment information and returns a list containing FOV assignments
of each of the candidates for later use when we analyze which images each candidate is in.
'''
def get_fov_ids(master_folder, analysis_folder, files, fovs):

    # parse out all spots and add them to a list
    dic_spots = []
    for cand_fl in files:
        dic_spot,dic_drift = pickle.load(open(cand_fl,'rb'), encoding='latin1')
        dic_spots.append(dic_spot)

    # append fov number to the list a number of times equal to the number of chromosomes in that FOV.
    fov_ids = list()
    for idx,fov in enumerate(dic_spots):
        for i in fov:
            fov_ids.append(idx)

    fov_ids = np.array(fov_ids)

    return fov_ids

'''
This function initialized important variables for the QC process. It takes just the folder names and finds the other information on its own.
Primarily, it grabs filenames for all dax images, fovs, and folders containing valid data.
'''
def grab_filenames(master_folders, analysis_folder):

    # grab all folder names with valid form (no dapi/h0, etc)
    folders = [folder for folder in glob.glob(master_folders[0]+os.sep+'*')
                                if os.path.isdir(folder)
                                    and os.path.basename(folder)[0]=='H' # if it starts with H
                                    and np.any([pattern in os.path.basename(folder) for pattern in ['R']]) # if it contains any values in ['R']
                                    and np.all([pattern not in os.path.basename(folder)for pattern in ['H0','dapi']])] # if it doesnt contain any elements in [ 'H0', 'dapi']

    # grab all dax images in each folder
    num_daxs = np.array([len(glob.glob(folder+os.sep+'*.dax')) for folder in folders])
    folders = np.array(folders)[num_daxs==np.max(num_daxs)] # only keep the folders that contain the correct number of dax images.

    #sort them by the order they were hybridized in
    folders = np.array(folders)[np.argsort(list(map(io.hybe_number,folders)))]

    # get the list of fov filenames (shared across all imaging rounds) e.g. Conv_zscan_00.dax
    fovs = np.sort(list(map(os.path.basename,glob.glob(folders[0]+os.sep+'*.dax'))))

    folders_keep = list(folders[:])

    # get all filenames for fovs which have identified candidates already
    files = glob.glob(analysis_folder+os.sep+'*__current_cand.pkl')
    fovs = [[os.path.basename(fl).split('__')[0]+'.dax'for fl in files]]


    return fovs, files, folders_keep


def run_QC_single_spot(center, dax_im, roi_size, spot_size):

    if not np.any(np.isnan(center)):
        #grab the data from the center.
        data = ft.grab_block(dax_im, center, [roi_size]*3)

        # compute the values for the given spot
        snr = single_spot_local_snr(data, max_spot_area=spot_size)
        bright = single_spot_brightness(data, 3*roi_size)
        area = single_spot_area(data, pixels=3*roi_size, max_spot_area=spot_size)

    else:
        snr = np.nan
        bright = np.nan
        area = np.nan

    return snr, bright, area

def run_standard_QC(master_folders, analysis_folder, zxys_f, files, fovs, folders_keep, keep, zxy_pix_size=[200,109,109], roi_size=15, suffix='.npy', spot_size=500):

    # got fov_id tags for all candidates.
    fov_ids = get_fov_ids(master_folders[0], analysis_folder, files, fovs)[keep]

    # initialize lists.
    snrs, brights, areas, color, round_ids = [], [], [], [], []

    zxys_f = np.array(zxys_f) # pull from the filtered data where we remove bad spots.

    for fov_number in tqdm(np.unique(fov_ids)):

        # Get selected chromosomal positions corresponding to the fov_number passed into the script.
        fov_id = fov_number
        fov_filename=fovs[0][fov_id]

        # gather dax objects for all imaging rounds (split by color)
        signal_dax_list, signal_color_tags, bead_dax_list, bead_color_tags = io.get_ims_fov_daxmap(folders_keep, fov_filename, col_tags = None, pad = 0)

        # get all spots for a given FOV
        fov_chrs = zxys_f[(np.where(fov_ids == fov_id)[0])]
        # go through all images and find all candidates in that image
        for idx,dax_im in tqdm(list(enumerate(signal_dax_list[:42]))):

            # go through all chromosomes and grab the block for this dax image in each chromosome
            for chrom in (fov_chrs):
                # carry an array that had round ids for each spot.
                round_ids.append(idx)

                # save color for each spot/round.
                color.append(idx%3)
                chr_round = np.array(chrom[idx,:])
                center = np.round(np.divide(chr_round,zxy_pix_size),0)

                # run single spot QCs and then append.
                snr, bright, area = run_QC_single_spot(center, dax_im, roi_size, spot_size)
                snrs.append(snr)
                brights.append(bright)
                areas.append(area)

    # make the folder for QC plot to be saved into
    QC_folder = analysis_folder+os.sep+'QC_folder'+os.sep
    if not os.path.exists(QC_folder):
        os.makedirs(QC_folder)

    # save out all data into the QC_folder along with coordinates used to create these files using the suffix given.
    save_loc = QC_folder+os.sep+r'snrs'+suffix
    np.save(save_loc, snrs)
    save_loc = QC_folder+os.sep+r'brights'+suffix
    np.save(save_loc, brights)
    save_loc = QC_folder+os.sep+r'areas'+suffix
    np.save(save_loc, areas)
    save_loc = QC_folder+os.sep+r'colors'+suffix
    np.save(save_loc, color)
    save_loc = QC_folder+os.sep+r'coords'+suffix
    np.save(save_loc, zxys_f)

    # save out parameter file
    import datetime

    param_fl = open(QC_folder+os.sep+r'params'+suffix[:-4]+r'.txt', 'w')
    param_fl.write("Parameters for "+suffix+ " data run\n")
    datetime_object = datetime.datetime.now()
    param_fl.write("Completed on: " + str(datetime_object)+"\n")
    param_fl.write("Max spot area: "+str(spot_size)+"\n")
    param_fl.write("ROI Size: " + str(roi_size) + " pixels\n")
    param_fl.close()

    print("Parameter file saved at: " + QC_folder+os.sep+r'params'+suffix[:-4]+r'.txt' )

    return snrs, brights, areas, color

# loads in a previous QC run given a suffix or loads in default if no suffix is given. Includes coordinates used to generate this QC result.
def load_qc_data(QC_folder, suffix='.npy'):

    save_loc = QC_folder+os.sep+r'snrs'+suffix
    snrs = np.load(save_loc)
    save_loc = QC_folder+os.sep+r'brights'+suffix
    brights = np.load(save_loc)
    save_loc = QC_folder+os.sep+r'areas'+suffix
    areas = np.load(save_loc)
    save_loc = QC_folder+os.sep+r'colors'+suffix
    colors = np.load(save_loc)
    save_loc = QC_folder+os.sep+r'coords'+suffix
    zxys_f = np.load(save_loc)

    return snrs, brights, areas, colors, zxys_f


def get_spot_vecs(master_folders, analysis_folder, zxys_f, files, fovs, folders_keep, keep, zxy_pix_size=[200,109,109], roi_size=15, suffix='.npy', spot_size=500):

    # got fov_id tags for all candidates.
    fov_ids = get_fov_ids(master_folders[0], analysis_folder, files, fovs)[keep]

    # initialize lists.
    vectors = []

    zxys_f = np.array(zxys_f) # pull from the filtered data where we remove bad spots.

    for fov_number in (np.unique(fov_ids)):

        # Get selected chromosomal positions corresponding to the fov_number passed into the script.
        fov_id = fov_number
        fov_filename=fovs[0][fov_id]

        # gather dax objects for all imaging rounds (split by color)
        signal_dax_list, signal_color_tags, bead_dax_list, bead_color_tags = io.get_ims_fov_daxmap(folders_keep, fov_filename, col_tags = None, pad = 0)

        # get all spots for a given FOV
        fov_chrs = zxys_f[(np.where(fov_ids == fov_id)[0])]
        # go through all images and find all candidates in that image
        for idx,dax_im in (list(enumerate(signal_dax_list[:42]))):

            # go through all chromosomes and grab the block for this dax image in each chromosome
            for chrom in (fov_chrs):

                chr_round = np.array(chrom[idx,:])
                center = np.round(np.divide(chr_round,zxy_pix_size),0)
                if np.isnan(center[0]):
                    continue
                # run single spot QCs and then append.
                data = ft.grab_block(dax_im, center, [roi_size]*3)
                data_2d = list(map(np.ravel, data))
                data_1d = np.ravel(data_2d)
                vectors.append(data_1d)

    vectors = np.array(vectors)
    saveloc = analysis_folder+os.sep+r'spot_vecs.npy'
    np.save(saveloc, vectors)

    return None

'''
This function takes two sets of chromosomes and finds the sets of common
chromosomes (within 3 pixels) and returns the set of unique ones for each data
set and then the shared ones.

Arguments:
data1 - chr positions for one set of chromosomes
data1 - chr positions to compare to data_1
threshold - pixel distance to consider two chrs as overlapping.

Returns:
overlap - set of chromosomes which are considered close enough for overlap
unique_data1, unique_data2 - set of chromsomes unique to one set or the other.
'''
def find_overlapping_chrs(data1, data2, data1_cids, data2_cids, fovs, threshold=10):

    import SimulateData as sd
    import copy

    # find fov assignments for all chrs in each dataset
    data1_fovs = fovs_from_cids(data1_cids, fovs)
    data2_fovs = fovs_from_cids(data2_cids, fovs)
    fov_assignment = []
    idx_pairs = []

    grouping, centers, min_dists = [],[],[]

    # go through each FOV and compare only the points within a given FOV
    for fov in data1_fovs.keys():

        # find the centers of the points in a given FOV
        data1_centers = sd.find_chr_centers(data1[data1_fovs[fov]])
        data2_centers = sd.find_chr_centers(data2[data2_fovs[fov]])

        data2_centers_placeholder = copy.deepcopy(data2_centers)

        # iterate through all points in data1 and search for a point within 3 pixels.
        for idx, point in enumerate(data1_centers):

            # create a list form of fov assignments for the end.
            fov_assignment.append(fov)

            dists = [np.linalg.norm(point-x) for x in data2_centers ]

            lowest = np.argmin(dists)
            min_dists.append(dists[lowest])

            # if the lowest distance is below threshold, add these to the final subset.
            # then remove that point it matched with from data2.
            if dists[lowest] <= threshold:

                # since data2_centers is shrinking we need to index based on a state_dict
                # version of the list to save out the index where we got it from.
                dists_no_remove = [np.linalg.norm(point-x) for x in data2_centers_placeholder ]
                fov_idx = np.where(dists_no_remove == dists[lowest])[0][0]
                idx1 = data1_fovs[fov][idx]
                idx2 = data2_fovs[fov][fov_idx]
                idx_pairs.append((idx1, idx2))

                centers.append(point)
                data2_centers = np.delete(data2_centers, lowest, axis=0)
                grouping.append(0)

            #otherwise add the point to the data1 list of uniques
            else:
                centers.append(point)
                grouping.append(1)

        # any remaining points are unmatched and unique to data2.
        # Extend FOV assignment and grouping accordingly.
        centers.extend(data2_centers)
        grouping.extend([2]*len(data2_centers))
        fov_assignment.extend([fov]*len(data2_centers))

    return grouping, centers, min_dists, fov_assignment, idx_pairs

'''
This function takes two lists of chromosomes along with the indices of where paired
chromosomes are located in the lists in a tuple format, i.e. (idx for data1,
idx for data2).

The function then goes through each point in the chromosome and classifies them
as 0-4, corresponding to:
0 - same spot in both (separated by < threshold nm)
1 - different spot (separated by > threshold nm)
2 - NaN changed to a spot (data1->data2)
3 - A spot changed to NaN (data1->data2)
4 - both data sets found a NaN

Returns:
groupings - list of all chromosomes in order of processing containing lists of
values for each spot from 0 to 4. (2d numpy array)
'''
def find_overlapping_spots(data1, data2, chr_pairs, threshold):

    groupings = []
    dists = []

    # go through all paired chromosomes
    for i, (idx1, idx2) in enumerate(chr_pairs):

        groupings.append([])

        # check each spot in that given chromosome
        for j in range(len(data1[idx1])):

            # if the first point is name and the second is not, add a 2.
            if np.isnan(data1[idx1][j][0]) and not np.isnan(data2[idx2][j][0]):
                groupings[i].append(2)
                continue
            # if the first point is not a NaN and the second is, add a 3.
            elif not np.isnan(data1[idx1][j][0]) and  np.isnan(data2[idx2][j][0]):
                groupings[i].append(3)
                continue
            # if the first point is name and the second is not, add a 2.
            elif np.isnan(data1[idx1][j][0]) and np.isnan(data2[idx2][j][0]):
                groupings[i].append(4)
                continue

            # if none of the values are nan, compute distance and check if it is below threshold nm
            dist = np.linalg.norm(data1[idx1][j]-data2[idx2][j])
            dists.append(dist)

            if dist < threshold:
                groupings[i].append(0)
            else:
                groupings[i].append(1)

    return groupings, dists

'''
Takes the cell ids and makes a dictionary storing the indices of the chromosomes
contained within each FOV with the FOV number as the key.
'''
def fovs_from_cids(cids, num_fovs):

    fovs = {}

    # initialize fov dictionary for each dataset
    for i in range(num_fovs):
        fovs[i] = []

    for idx, cell_id in enumerate(cids):
        fovs[cell_id//1000].append(idx)

    return fovs
