#PROBLEM 2 : Ice tracking

The solution for this problem helps solve global warming as best as it could by analysing the echograms and tracing the air-ice and ice-rock boundaries.

Initial Data Available:
1. m × n echogram
2. echo gram converted to edge strength by the code
3. Test images

Part1 : Using simple bayes net

Approach:
Given a image it is converted to edge-strength by the edge_strength function which returns a mxn array. We can treat these values as emission probabilities and proceed.
We mark the air-ice boundary and ice-rock boundary by visiting each coloumn one at a time.
For the first coloumn we find the maximum value of edge strength and note down its index. We then set the value at this index to zero and find the second maximum value and its index. If the differences between the indices are 15 pixels apart we set the minimum of the two index value as the air-ice boundary and the maximum of the two index value to ice-rock boundary. If the indices are not 15 pixels apart we mark the value at the second index to zero and find the next maximum value and its index. The below snippet of code demostrates how this is handled

****
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
****

Sample output is named as simpleresult.png and is placed in the repository.

Problems faced:
-->
Major problem faced during this implementation is that we assumed that the air-ice boundary will always have a higher pixel strength than the ice-rock boundary. This resulted in the line tracing the ice-rock boundary in some cases while finding the air-ice boundary. 
To solve this the approach of minimum and maximum index as mentioned in the Approach part was implemented. (Given that the air-ice boundary is always above the ice-rock boundary)

Part2:
-->

This part proved to be the most trickiest part as it involved calculation of transition probabilites and implementing the viterbi algorithm


Approach:
-->
To implement the viterbi algorithm we have used the psuedo code which can be found in the repository named viterbialgo.png

The first step was to calculate the emission probabilities. To create an emission probability table of size mxn where each element takes the value of the edge strength at that point divided by the sum of edge strengths of the coloumn.
    e[i][j] = edge_strength[i][j]/sum(edge_strength[:,j])

Then we calculate the transition probability which is given by (total rows -abs(diff between two rows))/sqrt(sum of two rows). If two sates are close to each other then the state is assigned with higher probability else with lower probability.

We then take this values and pass it to the viterbi algorith that returns the best path. 
Inside the viterbi:
--> 
We create two tables, the viterbi table and the backpointer table.
The viterbi table stores the maximized values calculated by multiplying transition probabilty, emission probabilty and the value in the viterbi stored at the previous state.

We then backtrack from the last coloumn to the first to get the maximum value which corresponds to the best path obtained

For tracing the ice rock boundary we make the elements that are above the air-ice boundary and 23 pixels below the air-ice boundary to zero in the edge_strength table, emission table and transition probabilty table.
The code snipet below handles this

****
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

****
A sample of the result is in the repository which is named as hmmresult.png


Problems faced:
--> 
During multiplication of transition probability, the value of trathensition probabilty starts approaching zero as we move along the coloumns and sometimes equates to zero which results in the boundary being traced at the roof of the image. To avoid this a bias term of 30000 is multiplied with the transition probabilty to make it not approach zero. 


Part 3:

In part 3 we incoporate the human feed back where the user gives a correction point which is a (x,y) coordinate on the image.
A simple way to solve this is to make the (x,y) in the transition table as 1 or maximize it so that the previous state is forced to approach this state in the next iteration.
Another way of solving this problem is to set all the values above the x coordinate in the yth row to zero and all the values below the x coordinate in the yth row to zero.

A result sample named feedbackresult.png is available in the repository.


#References

1. https://medium.com/analytics-vidhya/hidden-markov-model-part-1-of-the-hmm-series-3f7fea28a08

2. Discussed the problem with Deepan Elango, Akshay Venkatesh Murthy.

3. Slides from the lectures.






