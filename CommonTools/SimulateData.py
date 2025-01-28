'''
This package will have scripts associated with generating simulated data points for
chromosome traces and comparing the 'ground truth' data with algorithmically identified
points to assess improvements to the scripts.

Everything in here assumes coordinates are ordered [z, x, y] and that an image is
70 x 2048 x 2048 pixels.
'''

import numpy as np
from scipy.stats import truncnorm
import copy
from scipy.spatial.distance import cdist

'''
Converts a pixel location to a coordinate
'''
def pixel_to_coord(coord, pix_sz=[200,109,109]):
    return np.multiply(coord, pix_sz)

'''
Converts a coordinate location to a pixel index (subpixel)
'''
def coord_to_pixel(coord, pix_sz=[200,109,109]):
    return np.divide(coord,pix_sz)

'''
This is a weird algorithm that I found which allows me to partition the sq_dist
in an unbiased manner into "partitions" components (z,x,y) of random sizes which always
sum up to the correct distance.

We then return the square root of each component (pixel displacements)
'''
def parition_dist(dist, paritions=3):

    # the sq_dist is made up of a sum of x y and z squared so we partition this into the respective components
    sq_dist= dist**2

    nums = np.random.uniform(0,sq_dist,paritions-1) # generate paritions-1 uniform vals.

    # add 0 and sq_dist to the beginning and end of the sorted list.
    nums = np.concatenate([np.array([0]), sorted(nums), np.array([sq_dist]) ])

    # take pairwise differences to generate a list of paritions values summing to sq_dist
    partition = [nums[i+1]-nums[i] for i in range(len(nums)-1)]

    # return coordinate displacements for z,x,y from the sqrt of the partitions of sq_dist
    return coord_to_pixel(np.sqrt(partition))

'''
Takes a chr_center location (pixel location), number of points in a given chromosome,
the average displacement (in nm), and the stdev (in nm) of points from their nearet neighbors.

Returns a generated chromosome with num_pts and the z,x,y coordinates of each of
those points stored as pixel locations. (subpixel)
'''
def generate_chr(chr_center, num_pts, mean, stdev, im_size=[70,2048,2048]):

    # initialize array of points with z,x,y in each entry
    chr_pts = []

    for pt in chr_pts:
        pt = np.zeros(3)

    # begin iteration as a deviation from chr center location.
    prev_location = chr_center

    # for each point, displace it randomly in z,x,y
    for pt in range(num_pts):

        # generate a displacement distance (in nm)
        dist = truncnorm.rvs(0,2000,mean, stdev)
        # find z,x,y displacement in pixels
        dz,dx,dy = parition_dist(dist)
        disp = np.array([dz,dx,dy])

        # compute the signs (+/-) of each displacement
        signs = np.random.choice([-1,1],3)
        disp = np.multiply(disp,signs)

        current_coord = prev_location+disp

        # check if the current coordinate exceeds limits of image. If it does,
        # add the negative of the value just added and use that instead.
        for i in range(len(current_coord)):
            if current_coord[i]<=0. or current_coord[i]>im_size[i]-1.:
                current_coord[i] -= 2*disp[i] # correct in the opposite direction for that coordinate.

        # append the new coordinate
        chr_pts.append(current_coord)

        prev_location = current_coord

    return np.array(chr_pts)

'''
Given a set of center locations, generate chromosomes centered at those points
with n segments.
'''
def generate_FOV(center_locs, num_segments, mean, stdev, im_sz=[70,2048,2048]):

    fov_chrs = np.array([generate_chr(center_locs[i], num_segments, mean, stdev, im_sz) for i in range(len(center_locs))])

    return fov_chrs

'''
Given an integer, num_centers, and an image size, generates uniformly distr.
points across that image in z,x,y space.

Returns the array of centers.
'''
def generate_centers(num_centers, im_sz=[70,2048,2048]):

    centers = []

    for c in range(num_centers):
        z = np.random.uniform(0,im_sz[0])
        x = np.random.uniform(0,im_sz[1])
        y = np.random.uniform(0,im_sz[2])

        centers.append(np.array([z,x,y]))

    return np.array(centers)

'''
This function takes a set of coordinates and enforces a DE to it, removing a random
portion of the coordinates, replacing them with nan values. Requires a set of coordinates
and a detection efficiency (eff) from 0 to 1.
'''
def apply_detection_eff(coords, eff):

    new_coords = copy.deepcopy(coords)

    for i,chr in enumerate(coords):
        for j,coord in enumerate(chr):
            p = np.random.uniform(0,1)
            if p > eff:
                new_coords[i,j] = [np.nan]*3

    return new_coords


'''
'''
def estimate_errors(ground_truth, alg_results, thresh = 3):

    pairs = find_paired_chrs(ground_truth, alg_results)

    p = [ pair[2] for pair in pairs ]
    print(p)

    for gt, alg, dist in pairs:

        if dist < thresh:
            continue

    return None


'''
'''
def find_paired_chrs(ground_truth, alg_results):

    # find chr centers for each chromosome
    true_centers = find_chr_centers(ground_truth)
    results_centers = find_chr_centers(alg_results)

    print(ground_truth.shape)
    print(results_centers.shape)

    pairs = []

    for i in range(len(results_centers)):
        # go through each center identified in results and find the closest
        # ground truth point in the list.
        dists = cdist([results_centers[i]], true_centers)[0]

        loc = np.argmin(dists)

        print(ground_truth[loc].shape)
        print(results_centers[i].shape)

        # take the pair identified and save as a tuple
        dist = np.linalg.norm(ground_truth[loc]-results_centers[i], axis=0)
        print(dist)
        break
        pairs.append([ground_truth[loc], results_centers[i], dist])

    return pairs

'''
This function just finds the chromosome center for the coordinates and
returns a list of their positions.
'''
def find_chr_centers(coordinates):

    chr_centers = [ np.nanmean(chr,axis=0) for chr in coordinates ]

    return np.array(chr_centers)
