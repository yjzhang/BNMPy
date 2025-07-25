import os
import numpy as np
import pandas as pd


class ProbabilisticBN(object):

    def __init__(self, numberOfNodes, linkages, numberOfFunctions, functions, probabilities, initialNodeValues,
                 outputFilePath='', nodeDict=None):
        
        self.N = numberOfNodes
        self.varF = np.array(linkages)
        self.nf = np.array( numberOfFunctions )
        self.cij = np.array( probabilities )
        self.F = np.array(functions, dtype=np.int8)
        self.nodes = np.array(initialNodeValues, dtype=np.int8)

        self.networkHistory = np.array([[node for node in self.nodes]])
        self.nf_max = np.max( self.nf )
        self.K = []
        
        # helper array to obtained correc indexes for functions
        temp = np.cumsum( self.nf ) 
        self.cumsum = np.insert(temp,0,0)
        
        self.Nf = np.sum( self.nf ) 
        self.isConstanNode = np.full( self.Nf , False)    

        # Obtain the number of inputs per function
        for i in range( self.Nf ):
            self.K.append(0)
            for linkage in self.varF[ i ]:
                if (linkage == -1):
                    break
                self.K[i] += 1
        self.K = np.array(self.K)
        
        # determine if a node is constant
        for i in range( self.Nf ) :
            if ( self.K[i] == 0 )  :
                self.isConstanNode[i] = True
        
        self.outputFilePath = outputFilePath

        if ( outputFilePath != '' ) :
            self.initializeOutput()
        
        # Add nodeDict attribute for gene name to index mapping
        if nodeDict is None:
            nodeDict = {i: i for i in range(self.N)}
        self.nodeDict = nodeDict        

        # old connectivity matrices - for mutations
        self.old_varF = None
        self.old_F = None
        self.old_nf = None
        self.old_cij = None
        self.old_cumsum = None

        # determine if a node is constant (for efficient lookup)
        self.isConstantNode = np.zeros(self.N, dtype=bool)
        for i in range(self.N):
            if self.nf[i] == 1:
                func_idx = self.cumsum[i]
                if np.all(self.varF[func_idx] == -1):
                    self.isConstantNode[i] = True
                    
        # Initialize cumulative probabilities
        self.update_cumulative_probabilities()
        
        # Track knockdown nodes 
        self.knockdown_nodes = {}  # {node_idx: (target_value, efficacy)}

    def update_cumulative_probabilities(self):
        """
        Update cumulative probabilities for each node's functions.
        This is used for efficient function selection during simulation.
        """
        self.cij_cumsum = []
        
        for i in range(self.N):
            # Get probabilities for this node's functions
            probs = self.cij[i, :self.nf[i]]
            
            # Validate probabilities
            valid_probs = probs[probs >= 0]

            if len(valid_probs) == 0 or np.sum(valid_probs) <= 1e-5:
                # Use uniform distribution
                valid_probs = np.ones(self.nf[i]) / self.nf[i]
            else:
                # normalize probabilities
                valid_probs = valid_probs / np.sum(valid_probs)
            
            # Update Cij
            self.cij[i, :self.nf[i]] = valid_probs

            # Compute cumulative probabilities
            cumsum = np.cumsum(valid_probs)
            
            # Ensure last value is exactly 1.0
            cumsum[-1] = 1.0
            self.cij_cumsum.append(cumsum)
            
    def buildK(self):
        "This rebuilds the K array and related attributes."
        self.K = []
        
        # Rebuild K array
        for i in range(self.Nf):
            self.K.append(0)
            for linkage in self.varF[i]:
                if (linkage == -1):
                    break
                self.K[i] += 1
        self.K = np.array(self.K)
        
        # Rebuild isConstanNode array
        self.isConstanNode = np.full(self.Nf, False)
        for i in range(self.Nf):
            if (self.K[i] == 0):
                self.isConstanNode[i] = True
        
        # Update isConstantNode
        self.isConstantNode = np.zeros(self.N, dtype=bool)
        for i in range(self.N):
            if self.nf[i] == 1:
                func_idx = self.cumsum[i]
                if np.all(self.varF[func_idx] == -1):
                    self.isConstantNode[i] = True

    def setInitialValues(self, initialNodeValues):
        "Sets the initial values of the probabilistic boolean network."
        self.nodes = np.array(initialNodeValues, dtype=np.int8)

    def setInitialValue(self, key, value):
        "Sets a particular node to a given initial value, where the key is indexed in nodeDict."
        ind = self.nodeDict[key]
        self.nodes[ind] = value

    def knockout(self, key, value):
        "Sets a specific node to be permanently fixed to a given value."
        if self.old_varF is None:
            self.old_varF = self.varF.copy()
            self.old_F = self.F.copy()
            self.old_nf = self.nf.copy()
            self.old_cij = self.cij.copy()
            self.old_cumsum = self.cumsum.copy()
            
        # Set the initial value of the node
        self.setInitialValue(key, value)
        
        # Get the node index
        node_idx = self.nodeDict[key]
        
        # Set the number of functions for this node to 1
        self.nf[node_idx] = 1
        
        # Update cumsum array
        temp = np.cumsum(self.nf)
        self.cumsum = np.insert(temp, 0, 0)
        
        # Update Nf
        self.Nf = np.sum(self.nf)
        
        # Set the probability for this node's single function to 1.0
        for i in range(len(self.cij[node_idx])):
            if i == 0:
                self.cij[node_idx, i] = 1.0
            else:
                self.cij[node_idx, i] = -1.0
                
        # Find the function index for this node
        func_idx = self.cumsum[node_idx]
        
        # Set the connectivity of this node to -1 (indicating constant)
        self.varF[func_idx, :] = -1
        
        # Set the truth table entry to the constant value
        self.F[func_idx, 0] = value
        for i in range(1, len(self.F[func_idx])):
            self.F[func_idx, i] = -1
            
        # Rebuild K array and related attributes
        self.buildK()

        # Set this node as constant
        self.isConstantNode[node_idx] = True

    def knockdown(self, key, value, efficacy=1.0):
        """
        Sets a specific node with partial efficacy for experimental perturbation.
        E.g., when using a drug to inhibit a gene, the efficacy is the probability of the drug being effective.
        
        Parameters:
        -----------
        key : str
            Node name (must be in nodeDict)
        value : int
            Target value (0 or 1)
        efficacy : float, optional (default=1.0)
            Efficacy of the perturbation (0-1). 
            - For inhibitors (value=0): probability of achieving 0 = efficacy
            - For stimuli (value=1): probability of achieving 1 = efficacy
        """
        if efficacy >= 1.0:
            # Full efficacy - use regular knockout
            self.knockout(key, value)
        else:
            # Partial efficacy - track as knockdown
            node_idx = self.nodeDict[key]
            
            # Set initial value to target
            self.setInitialValue(key, value)
            
            # Store knockdown information
            self.knockdown_nodes[node_idx] = (value, efficacy)
            
            # Mark node as not constant (it has probabilistic behavior)
            self.isConstantNode[node_idx] = False

    def undoKnockouts(self):
        "Undoes all knockouts and knockdowns. Does not change initial values, however."
        if self.old_varF is not None:
            self.varF = self.old_varF
            self.F = self.old_F
            self.nf = self.old_nf
            self.cij = self.old_cij
            self.cumsum = self.old_cumsum
            self.Nf = np.sum(self.nf)
            
            # Reset stored originals
            self.old_varF = None
            self.old_F = None
            self.old_nf = None
            self.old_cij = None
            self.old_cumsum = None
            
            # Rebuild K array and related attributes
            self.buildK()
            
            # Update cumulative probabilities
            self.update_cumulative_probabilities()
        
        # Clear knockdown nodes
        self.knockdown_nodes = {}

    def initializeOutput(self):
        file = open(self.outputFilePath, 'w')
        stringToWrite = ''
        for i in range(self.N):
            stringToWrite += 'Node_{},'.format(i + 1)
        # replace last comma with newline
        stringToWrite = stringToWrite[:-1] + '\n'
        file.write(stringToWrite)
        file.write(self.stateToWrite())
        file.close()

    def stateToWrite(self):
        stringToWrite = ''
        for node in self.nodes:
            stringToWrite += (str(node) + ',')
        stringToWrite = stringToWrite[:-1] + '\n'

        return stringToWrite

    def writeNetworkHistory(self):
        string = ''
        for timestep in self.networkHistory:
            for node in timestep:
                string += (str(node) + ',')
            string = string[:-1] + '\n'
        file = open(self.outputFilePath, 'w')
        file.write(string)
        file.close()

    def _select_function(self, node_idx: int) -> int:
        """
        Select a function for a node based on its probability distribution.
        Uses pre-computed cumulative probabilities for efficiency.
        
        Parameters:
        -----------
        node_idx : int
            Index of the node
            
        Returns:
        --------
        int
            Selected function index relative to node's function block
        """
        if self.nf[node_idx] == 1:
            return 0
            
        r = np.random.random()
        cumsum = self.cij_cumsum[node_idx]
        return np.searchsorted(cumsum, r)

    def update(self, iterations=1):
        """
        Update the network directly.
        Select a function for each node based on its probability distribution.
        Then, update the network state based on the selected function.

        Parameters:
        -----------
        iterations : int, optional (default=1)
            Number of iterations to simulate
        """
        
        y = np.zeros( (iterations + 1, self.N ) , dtype=np.int8  )
        temp = np.array(  [2**i for i in range(self.N-1, -1, -1) ] ) 
        
        y[0] = self.nodes
        
        for itr in range(iterations):
            for i in range(self.N):
                if self.isConstantNode[i]:
                    y[itr+1][i] = y[itr][i]
                    continue
                
                # Check if this is a knockdown node
                if i in self.knockdown_nodes:
                    target_value, efficacy = self.knockdown_nodes[i]
                    # Apply probabilistic knockdown
                    if np.random.random() < efficacy:
                        y[itr+1][i] = target_value
                    else:
                        y[itr+1][i] = 1 - target_value
                else:
                    # Normal function-based update
                    func_offset = self._select_function(i)
                    idx = self.cumsum[i] + func_offset
                    
                    fInput = 0
                    for j in range(self.K[idx]):    
                        fInput += (y[itr][ self.varF[idx,j]]) * temp[ j - self.K[idx]  ]
                        
                    y[itr+1][i] = self.F[idx,fInput]

        self.nodes = y[-1] # newNodes

        return y
    
    # update the Boolean network with noise, derived from update functions (Boris)
    def update_noise(self, p, iterations=1):
        """
        Update the network with noise parameter p over a given number of iterations.
        
        This method simulates the network's evolution over time with added noise.
        Noise randomly flips node states with probability p after normal network update,
        but is never applied to constant nodes (knockouts).
        
        Parameters:
        -----------
        p : float between 0 and 1
            Probability of applying noise to each non-constant node at each iteration
        iterations : int, optional (default=1)
            Number of iterations to simulate
            
        Returns:
        --------
        y : numpy.ndarray
            Matrix of node states for each iteration (shape: iterations+1 × N)
        """
        y = np.zeros((iterations + 1, self.N), dtype=np.int8)
        temp = np.array([2**i for i in range(self.N-1, -1, -1)])
        
        y[0] = self.nodes

        for itr in range(iterations):
            # First, perform normal network update
            for i in range(self.N):
                if self.isConstantNode[i]:
                    # For constant nodes (knockouts), just copy their value
                    y[itr+1][i] = y[itr][i]
                elif i in self.knockdown_nodes:
                    # For knockdown nodes, apply probabilistic behavior
                    target_value, efficacy = self.knockdown_nodes[i]
                    if np.random.random() < efficacy:
                        y[itr+1][i] = target_value
                    else:
                        y[itr+1][i] = 1 - target_value
                else:
                    # For normal nodes, select a function based on probabilities
                    func_offset = self._select_function(i)
                    idx = self.cumsum[i] + func_offset

                    # Apply the selected function to update the node
                    fInput = 0
                    for j in range(self.K[idx]):
                        fInput += (y[itr][self.varF[idx,j]]) * temp[j - self.K[idx]]
                    
                    y[itr+1][i] = self.F[idx,fInput]
            
            # Then, optionally apply noise (only to non-constant nodes)
            if p > 0:
                for i in range(self.N):
                    if not self.isConstantNode[i] and np.random.rand() < p:
                        # Flip the bit for this node
                        y[itr+1][i] = 1 - y[itr+1][i]
        
        self.nodes = y[-1]  # Update the network state to the final iteration
        return y

    def getRealization(self):
        return (self.F, self.varF)

    def getMeanConnectivity(self):
        return sum(self.K) / len(self.K)

    def getMaxConnectivity(self):
        return max(self.K)

    def getBias(self):
        num0 = 0
        num1 = 0
        for function in self.F:
            for value in function:
                if value == 1:
                    num1 += 1
                elif value == 0:
                    num0 += 1
        return num1 / (num0 + num1)
   
    def getTrajectory( self ) : 
        return  self.networkHistory

    def copy(self) -> 'ProbabilisticBN':
        """
        Create a deep copy of the PBN.
        Required for optimization to avoid modifying original network.
        
        Returns:
        --------
        ProbabilisticBN
            A deep copy of this network
        """
        new_pbn = ProbabilisticBN(
            self.N,
            self.varF.copy(),
            self.nf.copy(),
            self.F.copy(),
            self.cij.copy(),
            self.nodes.copy(),
            nodeDict=self.nodeDict.copy()
        )
        
        # Copy additional attributes
        new_pbn.isConstantNode = self.isConstantNode.copy()
        new_pbn.K = self.K.copy()
        new_pbn.cumsum = self.cumsum.copy()
        
        return new_pbn
