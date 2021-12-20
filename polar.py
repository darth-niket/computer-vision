#!/usr/local/bin/python3
#
# Authors: Niket Malihalli:nmalihal
#          Nitanshu Joshi: nitjoshi
#          Sravani Wayangankar: swayang
#
# Ice layer finder
# Based on skeleton code by D. Crandall, November 2021
#
from PIL import Image
from numpy import *
from scipy.ndimage import filters
import sys
import imageio
import numpy as np
import random
import copy
# calculate "Edge strength map" of an image                                                                                                                                      
def edge_strength(input_image):
    grayscale = array(input_image.convert('L'))
    filtered_y = zeros(grayscale.shape)
    filters.sobel(grayscale,0,filtered_y)
    return sqrt(filtered_y**2)

# draw a "line" on an image (actually just plot the given y-coordinates
#  for each x-coordinate)
# - image is the image to draw on
# - y_coordinates is a list, containing the y-coordinates and length equal to the x dimension size
#   of the image
# - color is a (red, green, blue) color triple (e.g. (255, 0, 0) would be pure red
# - thickness is thickness of line in pixels
#

##Emission Probability
def emissionProbTab(edge_strength):
    N = edge_strength.shape[0]
    T = edge_strength.shape[1]
    Probability = np.zeros((N, T))
    for i in range(N):
        for j in range(T):
            Probability[i][j]=edge_strength[i][j]/sum(edge_strength[:,j])
    return Probability

#Transition Probability
def transitionProbTab(edge_strength):
    N = edge_strength.shape[0]
    T = edge_strength.shape[1]
    Probability = np.zeros((N, T))
    Probability = [[(N - abs(i - j)) / math.sqrt(i + j) for j in range(1, N + 1)] for i in range(1, N + 1)]
    Probability = array(Probability) / array(Probability).sum(axis=0)[:, newaxis]
    return Probability

#Viterbi algorithm
def viterbi(edge_strength, transitionProb, cood_x, cood_y):
    N = edge_strength.shape[0]
    T = edge_strength.shape[1]
    emissionProb = emissionProbTab(edge_strength)
    print(transitionProb)
    vit = np.zeros((N, T))
    bparray = np.zeros((N, T))
    for s in range(N):
        vit[s][0] = emissionProb[s][0]
        bparray[s][0] = 0
    if (cood_x <= -1) and (cood_y <= -1):
        for t in range(1, T):
            for s in range(N):
                maxarray = []
                for i in range(N):
                    maxarray.append(vit[i][t-1] * (transitionProb[i][s]*30000)* emissionProb[s][t])
                vit[s][t] = max(maxarray)
                bparray[s][t] = argmax(maxarray)
    else:
        for t in range(1, T):
            for s in range(N):
                maxarray = []
                for i in range(N):
                    if(cood_x==i) and (cood_y==s):
                        transitionProb[i][s]=1
                        maxarray.append(vit[i][t-1] * transitionProb[i][s]* emissionProb[s][t])
                    else:
                        maxarray.append(vit[i][t-1] * (transitionProb[i][s]*30000)* emissionProb[s][t])
                vit[s][t] = max(maxarray)
                bparray[s][t] = argmax(maxarray)
    temparray = []
    for s in range(N):
        temparray.append(vit[s][T-1])
    bestpathprob = max(temparray)
    bestpathpointer = argmax(temparray)
    bestpath = [0]*T
    for i in range(T-2,-1,-1):
        bestpath[i] = int(bestpathpointer)
        bestpathpointer = bparray[int(bestpathpointer)][i]
    return bestpath

############

def draw_boundary(image, y_coordinates, color, thickness):
    for (x, y) in enumerate(y_coordinates):
        for t in range( int(max(y-int(thickness/2), 0)), int(min(y+int(thickness/2), image.size[1]-1 )) ):
            image.putpixel((x, t), color)
    return image
def draw_asterisk(image, pt, color, thickness):
    for (x, y) in [ (pt[0]+dx, pt[1]+dy) for dx in range(-3, 4) for dy in range(-2, 3) if dx == 0 or dy == 0 or abs(dx) == abs(dy) ]:
        if 0 <= x < image.size[0] and 0 <= y < image.size[1]:
            image.putpixel((x, y), color)
    return image

# Save an image that superimposes three lines (simple, hmm, feedback) in three different colors
# (yellow, blue, red) to the filename
def write_output_image(filename, image, simple, hmm, feedback, feedback_pt):
    new_image = draw_boundary(image, simple, (255, 255, 0), 2)
    new_image = draw_boundary(new_image, hmm, (0, 0, 255), 2)
    new_image = draw_boundary(new_image, feedback, (255, 0, 0), 2)
    new_image = draw_asterisk(new_image, feedback_pt, (255, 0, 0), 2)
    imageio.imwrite(filename, new_image)


# main program
#
if __name__ == "__main__":
    if len(sys.argv) != 6:
        raise Exception("Program needs 5 parameters: input_file airice_row_coord airice_col_coord icerock_row_coord icerock_col_coord")
    input_filename = sys.argv[1]
     # load in image
    input_image = Image.open(input_filename).convert('RGB')
    image_array = array(input_image.convert('L'))
    gt_airice = [ int(i) for i in sys.argv[2:4] ]
    gt_icerock = [ int(i) for i in sys.argv[4:6] ]
    # load in image
    input_image = Image.open(input_filename).convert('RGB')
    image_array = array(input_image.convert('L'))
    # compute edge strength mask -- in case it's helpful. Feel free to use this.
    edge_strength = edge_strength(input_image)
    imageio.imwrite('edges.png', uint8(255 * edge_strength / (amax(edge_strength))))
    # You'll need to add code here to figure out the results! For now,
    # just create some random lines.
    
    #SIMPLE BOUNDARIES
    n = edge_strength.shape[1]
    eprob = emissionProbTab(edge_strength)
    airice_simple = []
    icerock_simple = []
    for i in range(n):
        listnp = eprob[:,i].tolist()
        firstmaxval = max(listnp)
        firstmax = listnp.index(firstmaxval)
        listnp[firstmax] = 0
        while(1):
            secondmaxval = max(listnp)
            secondmax = listnp.index(secondmaxval)
            if abs(firstmax - secondmax) >= 15:
                break
            listnp[secondmax] = 0
        valuemin=min(firstmax, secondmax)
        valuemax=max(firstmax, secondmax)
        airice_simple.append(valuemin)
        icerock_simple.append(valuemax)
    N = edge_strength.shape[0]
    transitionProb = transitionProbTab(edge_strength)


    #HMM CASES
    airice_hmm = viterbi(edge_strength,transitionProb, -1, -1)
    edge_strength_icerock = copy.deepcopy(edge_strength)
    for i in range(len(airice_hmm)):
        for j in range(1,23 ):
            edge_strength_icerock[airice_hmm[i] + j] = [0]*len(edge_strength_icerock[airice_hmm[i] + j])
            edge_strength_icerock[airice_hmm[i] - j] = [0]*len(edge_strength_icerock[airice_hmm[i] + j])  
    transitionProb1=transitionProbTab(edge_strength_icerock)
    for i in range(len(airice_hmm)):
        counter=0
        for j in range(1,23):
            transitionProb1[airice_hmm[i] + j] = [0]*len(transitionProb1[airice_hmm[i] + j])
        while((airice_hmm[i] - counter) >=0 ):
            transitionProb1[airice_hmm[i] - counter] = [0]*len(transitionProb1[airice_hmm[i] -counter])
            counter +=1
            
    icerock_hmm = viterbi(edge_strength_icerock, transitionProb1, -1, -1)

    #FEEDBACK CASES
    airice_feedback=viterbi(edge_strength, transitionProb, gt_airice[0], gt_airice[1])
    icerock_feedback=viterbi(edge_strength_icerock, transitionProb1, gt_icerock[0], gt_icerock[1])
    

    # Now write out the results as images and a text file
    write_output_image("air_ice_output.png", input_image, airice_simple, airice_hmm, airice_feedback, gt_airice)
    write_output_image("ice_rock_output.png", input_image, icerock_simple, icerock_hmm, icerock_feedback, gt_icerock)
    with open("layers_output.txt", "w") as fp:
        for i in (airice_simple, airice_hmm, airice_feedback, icerock_simple, icerock_hmm, icerock_feedback):
            fp.write(str(i) + "\n")