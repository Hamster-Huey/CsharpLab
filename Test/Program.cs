using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;
using System.Threading.Tasks;
using CsharpLab;
using MathNet.Numerics.LinearAlgebra.Single;
using Matrix = MathNet.Numerics.LinearAlgebra.Matrix<float>;
using Evd = MathNet.Numerics.LinearAlgebra.Factorization.Evd<float>;
namespace Test
{
    class Program
    {
        static unsafe void Main(string[] args)
        {
            //ConnNode tree1 = new ConnNode(1);
            //ConnNode tree2 = new ConnNode(2);
            //ConnNode tree3 = new ConnNode(3);
            //ConnNode tree4 = new ConnNode(4);
            //ConnNode tree5 = new ConnNode(5);
            //ConnNode tree6 = new ConnNode(6);
            //ConnNode tree7 = new ConnNode(7);
            //ConnNode tree8 = new ConnNode(8);
            //ConnNode tree9 = new ConnNode(9);
            //ConnNode tree10= new ConnNode(10);
            //ConnNode tree11= new ConnNode(11);
            //ConnNode tree12= new ConnNode(12);
            //ArrayList haha = new ArrayList();
            //Console.Write("Hello \n");
            //tree1.AddChild(tree2);
            //tree1.AddChild(tree3);
            //tree1.AddChild(tree4);
            //tree2.AddChild(tree5);
            //tree2.AddChild(tree6);
            //tree5.AddChild(tree7);
            //tree7.AddChild(tree8);
            //tree7.AddChild(tree9);
            //tree4.AddChild(tree10);
            //tree4.AddChild(tree11);
            //tree11.AddChild(tree12);
            //Console.WriteLine(tree2.GetSize());
            //foreach (int index in tree1.GetIndices())
            //{
            //    Console.WriteLine(index);
            //}

            //int[] frame = new int[49]{0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,};
            //ArrayList connects = Functions.BWConnComp(frame, 7, 7);
            //foreach(BWConnReturn thingy in connects)
            //{
            //    Console.WriteLine(thingy.size);
            //    foreach(int index in thingy.indices)
            //    {
            //        Console.WriteLine(index);
            //    }
            //}


            float[,] nana = new float[,] { 
                                           { 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0}, 
                                           { 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0}, 
                                           { 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1}, 
                                           { 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1}, 
                                           { 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0}, 
                                           { 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0} 
            };
            float[,] nan0 = new float[,] { { 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0 } };
            DenseMatrix m = DenseMatrix.OfArray(nana);
            var T = CsharpLab.Functions.BWConnComp(m,11,6);
#if DEBUG
                Console.WriteLine("Press enter to close...");
            Console.ReadLine();
#endif
        }
    }
}
