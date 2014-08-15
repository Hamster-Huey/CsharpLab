using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra.Single;
using MathNet.Numerics.Distributions;
using Matrix = MathNet.Numerics.LinearAlgebra.Matrix<float>;
using Evd = MathNet.Numerics.LinearAlgebra.Factorization.Evd<float>;
namespace CsharpLab
{
    public class Functions
    {

        /* Desired Functions: 
         * Connect Components
         * Extract Hand
         * Get Shape Descriptors
         * Get Covariance Matrix (in Math.Net Possibly?)
         * reshape
         * getTfeatures
         */

        /// <summary>
        /// Returns a 2d int array with number of connect components and indices of connected components
        /// </summary>
        /// <param name="array"></param>
        /// <param name="width"></param>
        /// <param name="height"></param>
        /// <returns></returns>
        public static List<BWConnReturn> BWConnComp(DenseMatrix array, int width, int height)
        {
            List<ConnNode> nodes = new List<ConnNode>();
            List<ConnNode> roots = new List<ConnNode>();

            List<BWConnReturn> returnArray = new List<BWConnReturn>();
            //Make each pixel a separate node
            for (int i = 0; i < array.RowCount; i++)
            {
                for (int j = 0; j < array.ColumnCount; j++)
                {

                    ConnNode pixNode = new ConnNode(i * array.ColumnCount + j);
                    if (array[i, j] == 1)
                    {
                        pixNode.setWhite();
                    }
                    nodes.Add(pixNode);
                }
            }

            //iterate through the nodes
            foreach (ConnNode node in nodes)
            {

                //check if node's white
                if (node.white == 1)
                {
                    //first check if root
                    if (node.parent == null)
                        roots.Add(node);

                    //now check up,down,left,right
                    //
                    //up
                    //first check if in top row
                    if (!(node.id < width))
                    {
                        ConnNode newNode = (ConnNode)nodes[node.id - width];
                        //inspected already
                        if (!(newNode.inspected))
                        {
                            //check if white
                            if (newNode.white == 1)
                            {
                                node.AddChild(newNode);
                                newNode.inspected = true;
                            }
                        }
                    }

                    //down
                    //first check if in bottom row
                    if (!(node.id >= (width * (height - 1))))
                    {
                        ConnNode newNode = (ConnNode)nodes[node.id + width];
                        //check if parent or inspected already
                        if (!(newNode.inspected))
                        {
                            //check if white
                            if (newNode.white == 1)
                            {
                                node.AddChild(newNode);
                                newNode.inspected = true;
                            }
                        }
                    }

                    //left
                    //first check if in left column
                    if (node.id % width != 0)
                    {
                        ConnNode newNode = (ConnNode)nodes[node.id - 1];
                        //check if parent or inspected already
                        if (!(newNode.inspected))
                        {
                            //check if white
                            if (newNode.white == 1)
                            {
                                node.AddChild(newNode);
                                newNode.inspected = true;
                            }
                        }
                    }

                    //right
                    //first check if in right column
                    if ((node.id + 1) % width != 0)
                    {
                        ConnNode newNode = (ConnNode)nodes[node.id + 1];
                        //check if parent or inspected already
                        if (!(newNode.inspected))
                        {
                            //check if white
                            if (newNode.white == 1)
                            {
                                node.AddChild(newNode);
                                newNode.inspected = true;
                            }
                        }
                    }
                }
                node.inspected = true;
            }
            //Go through each root and find size and indices
            foreach (ConnNode newnode in roots)
            {
                //make new BWConnReturn
                BWConnReturn values = new BWConnReturn();

                //Get size of tree
                values.size = newnode.GetSize();

                //now get indices
                values.indices = newnode.GetIndices();

                //now add returntype to array
                returnArray.Add(values);

            }
            return returnArray;
        } //checked
        public static List<DenseMatrix> ExtractHand(float[,] recording)
        {
            int x = recording.GetLength(0);
            int y = recording.GetLength(1);


            List<DenseMatrix> depthSequence = new List<DenseMatrix>();
            List<DenseMatrix> depthSequenceOriginal = new List<DenseMatrix>();
            //[frame num, frame]
            List<DenseMatrix> threshFrames = new List<DenseMatrix>();
            //transposed frames of recording
            for (int i = 0; i < x; i++)
            {
                float[] recordFrame = new float[217088];
                for (int j = 0; j < 217088; j++)
                    recordFrame[j] = recording[i, j];
                DenseMatrix frame = new DenseMatrix(512, 424, recordFrame);
                threshFrames.Add((DenseMatrix)frame.Transpose());
            }
            int numFrames = x;

            //Getting the Background Frame
            DenseMatrix backFrame = (DenseMatrix)threshFrames[0];

            //Variables for thresholding

            int numPlast = 0;
            int startThresh = 400;
            int startUpThresh = 200;
            int upThresh = 0;
            int downThresh = 0;
            int maxThresh = 0;
            int minThresh = 0;
            BWConnReturn maxClusterFinal = new BWConnReturn();
            //Thresholding each frame one by one 
            for (int frameNum = 1; frameNum < numFrames; frameNum++)
            {
                bool bad = false;

                //Background subtraction
                DenseMatrix depthFrame = (DenseMatrix)threshFrames[frameNum];
                depthSequenceOriginal.Add(depthFrame);
                depthFrame = backFrame - depthFrame; //Can't find Absolute value function
                DenseMatrix depthFrameTest = Threshold(depthFrame, 3000, 100);
                List<BWConnReturn> CC = BWConnComp(depthFrameTest, 512, 424);

                //ArrayList is a list of objects that each contain clusters with the indeces and number of pixels
                //Goal 1: Find the cluster with the most pixels and store that object
                BWConnReturn maxCluster = new BWConnReturn();
                for (int i = 0; i < CC.Count; i++)
                {
                    if (((BWConnReturn)CC[i]).size > maxCluster.size)
                    {
                        maxCluster = (BWConnReturn)CC[i];
                    }
                }

                //Goal 2: set the indices to 1 and everything else to 0
                depthFrameTest.Clear();
                for (int i = 0; i < maxCluster.indices.Count; i++)
                {
                    //get row and column from index
                    int row = (int)((int)maxCluster.indices[i] / depthFrame.ColumnCount);
                    int column = (int)((int)maxCluster.indices[i] % depthFrame.ColumnCount);

                    //set to 1
                    depthFrameTest[row, column] = 1;
                }

                //Now layer the original depthFrame over the new silhouette 
                depthFrameTest.PointwiseMultiply((DenseMatrix)depthSequenceOriginal[frameNum-1], depthFrame);
                
                //Feedback loop that applies thresholds until satisfied with the result
                bool goodClump = false;
                int iterations = 0;
                while (goodClump == false)
                {
                    int threshold = 0;

                    //Get Hand from First Frame
                    if (frameNum == 1)
                    {
                        List<int> goodIndices = findIndices(depthFrame, startThresh);
                        threshold = goodIndices.OfType<int>().Min();
                        downThresh = threshold - 200;
                        upThresh = threshold + startUpThresh;
                    }
                    else if (bad) //Re-Threshold Frames if first try doesn't work
                    {
                        if (startThresh > 1500)
                        {
                            List<int> goodIndices = findIndices(depthFrame, (int)(maxThresh + minThresh) / 2);
                            threshold = goodIndices.OfType<int>().Min();
                        }
                        else
                        {
                            List<int> goodIndices = findIndices(depthFrame, startThresh);
                            threshold = goodIndices.OfType<int>().Min();
                            startUpThresh = 200;
                            goodClump = true;
                        }
                        downThresh = threshold - 200;
                    }

                    //Threshold
                    depthFrameTest = Threshold(depthFrame, upThresh, downThresh);

                    //Connect Components
                    CC = BWConnComp(depthFrameTest, 512, 424);
                    //Find Max+IDs
                    maxClusterFinal = new BWConnReturn();
                    for (int i = 0; i < CC.Count; i++)
                    {
                        if (((BWConnReturn)CC[i]).size > maxClusterFinal.size)
                        {
                            maxClusterFinal = (BWConnReturn)CC[i];
                        }
                    }

                    //get Difference
                    int numP = maxClusterFinal.size;
                    int difference = numPlast - numP;

                    //Feedback steps for the first frame
                    if (frameNum == 1)
                    {
                        if ((numP > 4500) && (numP < 5500))
                        {
                            goodClump = true; //Decent sized clump
                            numPlast = numP;
                        }
                        else if (numP > 5500)
                            startUpThresh -= 10;
                        else if (numP < 4500)
                            startThresh += 20;
                    }
                    //Feedback for other frames based on difference in size of frame before it, also not too small 
                    else if ((Math.Abs(difference) > 600) || (numP < 1600))
                    {
                        bad = true;
                        //Check to see if anything significant was caught
                        if ((difference > 0) || (numP < 1600))
                        {
                            maxThresh += 2;
                            upThresh = maxThresh;
                        }
                        else
                        {
                            maxThresh -= 2;
                            upThresh = maxThresh;
                        }
                        //Gotta stop sometime
                        iterations++;
                        if (iterations == 3000)
                            goodClump = true;
                    }
                    else
                    {
                        goodClump = true;
                        numPlast = numP;
                    }

                }

                //Set indices to 1
                depthFrame.Clear();
                for (int i = 0; i < maxClusterFinal.indices.Count; i++)
                {
                    //get row and column from index
                    int row = (int)((int)maxClusterFinal.indices[i] / depthFrame.ColumnCount);
                    int column = (int)((int)maxClusterFinal.indices[i] % depthFrame.ColumnCount);

                    //set to 1
                    depthFrame[row, column] = 1;
                }

                //Todo: Fill the Holes
                //
                //

                //Layer Depths onto silhouette to extract the max and min depths
                DenseMatrix depthFrameLayered = new DenseMatrix(424, 512);
                ((DenseMatrix)depthSequenceOriginal[frameNum]).PointwiseMultiply(depthFrame, depthFrameLayered);
                List<int> greaterThanZeros = findIndices(depthFrameLayered, 0);
                maxThresh = greaterThanZeros.OfType<int>().Max();
                minThresh = greaterThanZeros.OfType<int>().Min();

                startThresh = (int)((maxThresh + minThresh) / 2);
                upThresh = maxThresh;
                downThresh = minThresh - 300;

                depthSequence.Add(depthFrame);
            }


            return depthSequence;
        }
        public static DenseMatrix Threshold(DenseMatrix inMatrix, int upThreshold, int downThreshold)
        {
            //Within Threshold => 1, Outside => 0
            for (int i = 0; i < (inMatrix.ColumnCount); i++)
            {
                for (int j = 0; j < inMatrix.RowCount; j++)
                {
                    if ((downThreshold <= Math.Abs(inMatrix[j, i])) && (Math.Abs(inMatrix[j, i]) <= upThreshold))
                    {
                        inMatrix[j, i] = 1;
                    }
                    else { inMatrix[j, i] = 0; }
                }
            }

            return inMatrix;
        }
        public static List<int> findIndices(DenseMatrix inMatrix, int lowerBound) //checked
        {
            List<int> indices = new List<int>();
            for (int i = 0; i < inMatrix.RowCount; i++)
            {
                for (int j = 0; j < inMatrix.ColumnCount; j++)
                {
                    if (inMatrix[i, j] > lowerBound)
                        indices.Add(i * inMatrix.ColumnCount + j);
                }
            }


            return indices;
        }
        /// <summary>
        /// Returns N, S, E, W frames for a given BW frame
        /// </summary>
        /// <param name="inFrame"></param>
        /// <returns></returns>
        public static List<DenseMatrix> getNSEWfeatures(DenseMatrix inFrame)
        {
            List<DenseMatrix> returnList = new List<DenseMatrix>();
            //Assume frame is already logical/binary
            int x = inFrame.ColumnCount;
            int y = inFrame.RowCount;
           
            //Initialize Matrices
            DenseMatrix dN = new DenseMatrix(y, x);
            DenseMatrix dS = new DenseMatrix(y, x);
            DenseMatrix dE = new DenseMatrix(y, x);
            DenseMatrix dW = new DenseMatrix(y, x);

            //get indices of shape of hand
            List<int> indices = findIndices(inFrame, 0);

            //sort indices into x and y lists
            List<int> xs = new List<int>();
            List<int> ys = new List<int>();

            for (int i = 0; i < indices.Count; i++)
            {

                ys.Add((int)((int)indices[i] / inFrame.ColumnCount));
                xs.Add((int)((int)indices[i] % inFrame.ColumnCount));
            }

            //get x and y min/max
            int ymin = ys.OfType<int>().Min();
            int ymax = ys.OfType<int>().Max();
            int xmin = xs.OfType<int>().Min();
            int xmax = xs.OfType<int>().Max();

            //loops

            //DW
            for (int i = ymin; i <= ymax; i++)
            {
                int count = 0;
                for (int j = xmin; j <= xmax; j++)
                {
                    if (inFrame[i, j] == 1)
                    {
                        count = count + 1;
                        dW[i, j] = count;
                    }
                    else
                        count = -1;
                }
            }

            //DE
            for (int i = ymin; i <= ymax; i++)
            {
                int count = 0;
                for (int j = xmax; j >= xmin; j--)
                {
                    if (inFrame[i, j] == 1)
                    {
                        count = count + 1;
                        dE[i, j] = count;
                    }
                    else
                        count = -1;
                }
            }

            //DN
            for (int i = xmin; i <= xmax; i++)
            {
                int count = 0;
                for (int j = ymin; j <= ymax; j++)
                {
                    if (inFrame[j, i] == 1)
                    {
                        count = count + 1;
                        dN[j, i] = count;
                    }
                    else
                        count = 0;
                }
            }

            //DS
            for (int i = xmin; i <= xmax; i++)
            {
                int count = 0;
                for (int j = ymax; j >= ymin; j--)
                {
                    if (inFrame[j, i] == 1)
                    {
                        count = count + 1;
                        dS[j, i] = count;
                    }
                    else
                        count = 0;
                }
            }

            returnList.Add(dN);
            returnList.Add(dS);
            returnList.Add(dE);
            returnList.Add(dW);

            return returnList;
        }//checked

        public static DenseMatrix getdiagboundary(DenseMatrix inFrame)
        {
            DenseMatrix output = new DenseMatrix(inFrame.RowCount, inFrame.ColumnCount);

            for (int diagIdx = -1 * inFrame.RowCount + 1; diagIdx <= inFrame.ColumnCount - 1; diagIdx++)
            {
                MathNet.Numerics.LinearAlgebra.Vector<float> tempDiag = diag(inFrame, diagIdx);
                DenseMatrix tempDiagSum = new DenseMatrix(1, tempDiag.Count);
                int count = 0;

                for (int ii = 0; ii < tempDiag.Count; ii++)
                {
                    if (tempDiag[ii] == 1)
                    {
                        count = count + 1;
                        tempDiagSum[0, ii] = count;
                    }
                    else
                        count = 0;
                }

                if (diagIdx <= 0)
                {
                    int offset = 0;
                    //iterate through diagonal
                    for (int i = -diagIdx; (i < inFrame.RowCount) && (offset < inFrame.ColumnCount); i++)
                    {
                        output[i, offset] = tempDiagSum[0, offset];
                        offset++;
                    }
                }
                else
                {
                    int offset = 0;
                    for (int i = diagIdx; (i < inFrame.ColumnCount) && (offset < inFrame.RowCount); i++)
                    {
                        output[offset, i] = tempDiagSum[0, offset];
                        offset++;
                    }
                }
            }


            return output;
        }//checked

        public static DenseMatrix fliplr(DenseMatrix inFrame)
        {
            for (int i = 0; i < inFrame.RowCount; i++)
            {
                for (int j = 0; j < (int)(inFrame.ColumnCount / 2); j++)
                {
                    inFrame[i, j] += inFrame[i, inFrame.ColumnCount - 1 - j];
                    inFrame[i, inFrame.ColumnCount - 1 - j] = inFrame[i, j] - inFrame[i, inFrame.ColumnCount - 1 - j];
                    inFrame[i, j] -= inFrame[i, inFrame.ColumnCount - 1 - j];
                }
            }
            return inFrame;
        } //checked

        public static DenseMatrix flipud(DenseMatrix inFrame)
        {
            for (int i = 0; i < inFrame.ColumnCount; i++)
            {
                for (int j = 0; j < (int)(inFrame.RowCount / 2); j++)
                {
                    inFrame[j, i] += inFrame[inFrame.RowCount - 1 - j, i];
                    inFrame[inFrame.RowCount - j - 1, i] = inFrame[j, i] - inFrame[inFrame.RowCount - 1 - j, i];
                    inFrame[j, i] -= inFrame[inFrame.RowCount - 1 - j, i];
                }
            }
            return inFrame;
        }//checked

        /// <summary>
        /// return an array of diagonal values of input matrix given an offset
        /// </summary>
        /// <param name="inFrame"></param>
        /// <param name="diag_idx"></param>
        /// <returns></returns>
        public static MathNet.Numerics.LinearAlgebra.Vector<float> diag(DenseMatrix inFrame, int diag_idx)
        {
            ///Idea: Given offset, create new matrix from specified offset and get diagonal
            //Trying to mimic Matlab's function

            DenseMatrix newMatrix = null;
            if (diag_idx <= 0)
            {
                //Along length
                newMatrix = (DenseMatrix)inFrame.SubMatrix(-diag_idx, inFrame.RowCount + diag_idx, 0, inFrame.ColumnCount);
            }
            else if (diag_idx > 0)
            {
                newMatrix = (DenseMatrix)inFrame.SubMatrix(0, inFrame.RowCount, diag_idx, inFrame.ColumnCount - diag_idx);
            }
            MathNet.Numerics.LinearAlgebra.Vector<float> diagonals = newMatrix.Diagonal();

            return diagonals;
        } //checked
        /// <summary>
        /// Take in a CxT matrix (horizontal slice of recording)
        /// </summary>
        /// <param name="inFrame"></param>
        /// <returns></returns>
        public static int[][] getTfeatures(List<DenseMatrix> inFrames)
        {
            //Assume frame is already logical/binary
            int x = ((DenseMatrix)inFrames[0]).ColumnCount;
            int y = ((DenseMatrix)inFrames[0]).RowCount;
            int bignum = x * y;

            int[][] returnList = new int[2][];
            //Initialize Matrices
            int[] dT_Plus = new int[bignum * inFrames.Count];
            int[] dT_Minus = new int[bignum * inFrames.Count];

            int[] flatRecording = new int[bignum * inFrames.Count];

            int count = 0;
            //flatten out recording
            for (int i = 0; i < inFrames.Count; i++)
            {
                for (int j = 0; j < y; j++)
                {
                    for (int k = 0; k < x; k++)
                    {
                        flatRecording[count] = (int)((DenseMatrix)inFrames[i])[j, k];
                        count++;
                    }
                }
            }

            //Now index and stuff
            //Plus
            for (int i = 0; i < bignum; i++)
            {
                int c = 0;
                int jj = i;
                while (jj < bignum * inFrames.Count)
                {
                    if (flatRecording[jj] == 1)
                    {
                        c++;
                        dT_Plus[jj] = c;
                    }
                    else
                        c = 0;
                    jj += bignum;
                }
            }

            //Minus
            for (int i = bignum * inFrames.Count - 1; i >= bignum * inFrames.Count - bignum; i--)
            {
                int c = 0;
                int jj = i;
                while (jj >= 0)
                {
                    if (flatRecording[jj] == 1)
                    {
                        c++;
                        dT_Minus[jj] = c;
                    }
                    else
                        c = 0;
                    jj -= bignum;
                }
            }
            returnList[0] = dT_Plus;
            returnList[1] = dT_Minus;

            return returnList;
        } //checked

        public static float[,] getCovFeaturesFromSequencePath(ushort[,] recording, int dataSize)
        {
            //Convert into float array
            float[,] fRecording = new float[dataSize, recording.GetLength(1)];
            for (int i = 0; i < dataSize; i++)
            {
                for (int j = 0; j < recording.GetLength(1); j++)
                {
                    fRecording[i, j] = (float)recording[i, j];
                }
            }

            //Get List of Matrices for extracted hand sequence
            List<DenseMatrix> cleanSequence = ExtractHand(fRecording);

            //compute time features
            List<DenseMatrix> dE_all = new List<DenseMatrix>();
            List<DenseMatrix> dW_all = new List<DenseMatrix>();
            List<DenseMatrix> dS_all = new List<DenseMatrix>();
            List<DenseMatrix> dN_all = new List<DenseMatrix>();
            List<DenseMatrix> dNE_all = new List<DenseMatrix>();
            List<DenseMatrix> dNW_all = new List<DenseMatrix>();
            List<DenseMatrix> dSE_all = new List<DenseMatrix>();
            List<DenseMatrix> dSW_all = new List<DenseMatrix>();

            for (int frameNum = 0; frameNum < cleanSequence.Count; frameNum++)
            {
                DenseMatrix BW = (DenseMatrix)cleanSequence[frameNum];
                List<DenseMatrix> NSEW = getNSEWfeatures(BW);

                DenseMatrix dNW = null;
                DenseMatrix dNE = null;
                DenseMatrix dSW = null;
                DenseMatrix dSE = null;
                BW.PointwiseMultiply(getdiagboundary(BW), dNW);
                BW.PointwiseMultiply(fliplr(getdiagboundary(fliplr(BW))), dNE);
                BW.PointwiseMultiply(flipud(getdiagboundary(flipud(BW))), dSW);
                dSE = getdiagboundary(flipud(fliplr(BW)));
                dSE = flipud(fliplr(dSE));
                dSE.PointwiseMultiply(BW, dSE);

                dE_all.Add(NSEW[2]);
                dN_all.Add(NSEW[0]);
                dW_all.Add(NSEW[3]);
                dS_all.Add(NSEW[1]);
                dNE_all.Add(dNE);
                dNW_all.Add(dNW);
                dSE_all.Add(dSE);
                dSW_all.Add(dSW);
            }
            int[][] dtFeatures = getTfeatures(cleanSequence);
            int[] dT_plus = dtFeatures[0];
            int[] dT_minus = dtFeatures[1];
            int[] X = new int[((DenseMatrix)cleanSequence[0]).RowCount * ((DenseMatrix)cleanSequence[0]).ColumnCount * cleanSequence.Count];
            int[] Y = new int[((DenseMatrix)cleanSequence[0]).RowCount * ((DenseMatrix)cleanSequence[0]).ColumnCount * cleanSequence.Count];
            int[] T = new int[((DenseMatrix)cleanSequence[0]).RowCount * ((DenseMatrix)cleanSequence[0]).ColumnCount * cleanSequence.Count];

            //Get X
            for (int i = 0; i < ((DenseMatrix)cleanSequence[0]).RowCount * cleanSequence.Count; i++)
            {
                for (int j = 0; j < ((DenseMatrix)cleanSequence[0]).ColumnCount; j++)
                {
                    X[i * ((DenseMatrix)cleanSequence[0]).ColumnCount - 1 + j] = j + 1;
                }
            }

            //Get Y
            for (int u = 0; u < cleanSequence.Count; u++)
            {
                for (int i = 0; i < ((DenseMatrix)cleanSequence[0]).RowCount; i++)
                {
                    for (int j = 0; j < ((DenseMatrix)cleanSequence[0]).ColumnCount; j++)
                    {
                        Y[(u * ((DenseMatrix)cleanSequence[0]).RowCount * ((DenseMatrix)cleanSequence[0]).ColumnCount - 1 + i * ((DenseMatrix)cleanSequence[0]).ColumnCount - 1 + j)] = i + 1;
                    }
                }
            }

            //Get T
            for (int u = 0; u < cleanSequence.Count; u++)
            {
                for (int i = 0; i < ((DenseMatrix)cleanSequence[0]).RowCount; i++)
                {
                    for (int j = 0; j < ((DenseMatrix)cleanSequence[0]).ColumnCount; j++)
                    {
                        T[(u * ((DenseMatrix)cleanSequence[0]).RowCount * ((DenseMatrix)cleanSequence[0]).ColumnCount - 1 + i * ((DenseMatrix)cleanSequence[0]).ColumnCount - 1 + j)] = u;
                    }
                }
            }
            ///LOL I hope this works
            //SO I need to reshape every matrix/array into a vertical line...but only the pixels that are logical



            DenseMatrix F = null;
            for (int frame_num = 0; frame_num < recording.GetLength(0); frame_num++)
            {

                //Find every hand pixel index for each frame
                DenseMatrix curFrame = (DenseMatrix)cleanSequence[frame_num];
                List<int> n = findIndices(curFrame, 0);
                //sort indices into x and y lists
                List<int> xs = new List<int>();
                List<int> ys = new List<int>();

                for (int i = 0; i < n.Count; i++)
                {

                    ys.Add((int)((int)n[i] / curFrame.ColumnCount));
                    xs.Add((int)((int)n[i] % curFrame.ColumnCount));
                }
                //create new 13xindices DenseMatrix to represent the frame
                DenseMatrix features = new DenseMatrix(13, n.Count);

                //now add all of the features
                for (int i = 0; i < n.Count; i++)
                {
                    //X Y and T are all Linear thingies
                    features[0, i] = X[curFrame.RowCount * curFrame.ColumnCount * frame_num + (int)n[i]];
                    features[1, i] = Y[curFrame.RowCount * curFrame.ColumnCount * frame_num + (int)n[i]];
                    features[2, i] = T[curFrame.RowCount * curFrame.ColumnCount * frame_num + (int)n[i]];
                    features[3, i] = dN_all[frame_num][(int)ys[i], (int)xs[i]];
                    features[4, i] = dS_all[frame_num][(int)ys[i], (int)xs[i]];
                    features[5, i] = dE_all[frame_num][(int)ys[i], (int)xs[i]];
                    features[6, i] = dW_all[frame_num][(int)ys[i], (int)xs[i]];
                    features[7, i] = dNE_all[frame_num][(int)ys[i], (int)xs[i]];
                    features[8, i] = dNW_all[frame_num][(int)ys[i], (int)xs[i]];
                    features[9, i] = dSE_all[frame_num][(int)ys[i], (int)xs[i]];
                    features[10, i] = dSW_all[frame_num][(int)ys[i], (int)xs[i]];
                    features[11, i] = dT_plus[curFrame.RowCount * curFrame.ColumnCount * frame_num + (int)n[i]];
                    features[12, i] = dT_minus[curFrame.RowCount * curFrame.ColumnCount * frame_num + (int)n[i]];
                }
                if (F == null)
                {
                    F = features;
                }
                else
                {
                    F = (DenseMatrix)F.Append(features);
                }
            }

            //Now we should have F
            F = (DenseMatrix)F.Transpose();
            float[,] Farray = F.ToArray();
            float[,] Covariance = vanillaLogCov.Covariance(Farray);
            return Covariance;
        }
    }
    class vanillaLogCov
    {

        //put variables here
        DenseMatrix cov1, cov2;
        Matrix diag, logcov1, logcov2;
        Evd evd; //factorization storage
        float distance;

        public vanillaLogCov(float[,] matrix1, float[,] matrix2)
        {
            //deep copy the covariances
            //densematrix allocates new memory for the matrix

            cov1 = DenseMatrix.OfArray(Covariance(matrix1));
            cov2 = DenseMatrix.OfArray(Covariance(matrix2));
            evd = cov1.Evd();
            diag = evd.D;
            //could make this a function but thats for later cleanup.... < . <
            //should be square matrix so matrix lengths dont matter
            for (int i = 0; i < diag.RowCount; i++)
            {
                diag[i, i] = (float)Math.Log(Math.Abs((double)diag[i, i]));
            }
            logcov1 = evd.EigenVectors * diag * evd.EigenVectors.Transpose();

            evd = cov2.Evd();
            diag = evd.D;
            //could make this a function but thats for later cleanup.... < . <
            //should be square matrix so matrix lengths dont matter
            for (int i = 0; i < diag.RowCount; i++)
            {
                diag[i, i] = (float)Math.Log(Math.Abs((double)diag[i, i]));
            }
            logcov2 = evd.EigenVectors * diag * evd.EigenVectors.Transpose();

            //compute distance
            distance = 0; //initialize to empty
            for (int i = 0; i < logcov2.RowCount; i++)
            {
                for (int j = 0; j < logcov2.RowCount; j++)
                {
                    distance += (float)Math.Pow(logcov1[i, j] - logcov2[i, j], 2);
                }
            }
            distance = (float)Math.Sqrt(distance);
        }

        public float getDistance()
        {
            return distance;
        }
        private static float[,] SubtractMeans(float[,] matrix)
        {
            var x = matrix.GetLength(0);
            var y = matrix.GetLength(1);
            var result = new float[x, y];
            
            for (int i = 0; i < y; i++)
            {
                var tmp = 0d;
                for (int j = 0; j < x; j++)
                    tmp += matrix[j,i];
                var mean = tmp / x;
                for (int j = 0; j < x; j++)
                    result[j,i] = matrix[j,i] - (float)mean;
            }
            return result;
        } //Uhh yea
        public static float[,] Transpose(float[,] matrix)
        {
            var x = matrix.GetLength(0);
            var y = matrix.GetLength(1);
            var result = new float[y, x];
            for (int i = 0; i < x; i++)
                for (int j = 0; j < y; j++)
                    result[j, i] = matrix[i, j];
            return result;
        }
        //public unsafe static float[,] Multiply(float[,] matrix1, float[,] matrix2)
        //{
        //    var x = matrix1.GetLength(0);
        //    var y = matrix2.GetLength(1);
        //    var z = matrix1.GetLength(1);
        //    if (matrix1.GetLength(1) != matrix2.GetLength(0))
        //        throw new InvalidOperationException("Can't multiply");
        //    var result = new float[x, y];
        //    fixed (float* p1 = matrix1)
        //    fixed (float* p2 = matrix2)
        //    fixed (float* p3 = result)
        //    {
        //        for (int i = 0; i < x; i++)
        //            for (int j = 0; j < y; j++)
        //            {
        //                float tmp = p3[i * y + j];
        //                for (int k = 0; k < x; k++)
        //                    tmp += p1[i * z + k] * p2[k * y + j];
        //                p3[i * y + j] = tmp;
        //            }
        //    }
        //    return result;
        //}
        public static float[,] Covariance(float[,] matrix)
        {
            var n = SubtractMeans(matrix);
            DenseMatrix matN = DenseMatrix.OfArray(n);
            var t = matN.Transpose();
            var matM = t.Multiply(matN);
            var m = matM.ToArray();
            var x = m.GetLength(0);
            var y = m.GetLength(1);
            
            for (int i = 0; i < x; i++)
                for (int j = 0; j < y; j++)
                {
                    var tmp = m[i, j];
                    tmp /= y - 1;
                    m[i, j] = tmp;
                }
            
            //Need to flatten
            return m;
        }
    }
    public class ConnNode
    {
        public ConnNode parent;

        public List<ConnNode> children;
        public int id;
        public int white;
        public bool inspected;

        public ConnNode(int num)
        {
            this.id = num;
            this.parent = null;
            this.children = new List<ConnNode>();
            this.white = 0;
            inspected = false;
        }

        public void AddChild(ConnNode child)
        {
            children.Add(child);
            child.AddParent(this);
        }

        public void setWhite()
        {
            white = 1;
        }

        public void AddParent(ConnNode tree)
        {
            this.parent = tree;
        }

        public int GetParent()
        {
            return parent.id;
        }

        public List<ConnNode> GetChildren()
        {

            return children;
        }

        public int GetSize()
        {
            int length = 1;
            // Check LRM so left then right then middle
            foreach (ConnNode child in children)
            {
                length += child.GetSize();
            }
            return length;
        }

        public List<int> GetIndices()
        {
            List<int> indices = new List<int>();

            foreach (ConnNode child in children)
            {
                indices.AddRange(child.GetIndices());
            }
            indices.Add(this.id);
            return indices;

        }

    }

    public class BWConnReturn
    {
        public int size;
        public List<int> indices;

        public BWConnReturn()
        {
            size = 0;
            indices = null;
        }
    }

}
