using System;
using System.Collections;
using System.Collections.Generic;
namespace xorc
{
	class MainClass
	{
		public int inputN = 2;
		public int outputN = 1;
		public int populationSize = 100;
		public List<Network> networks;
		public Random rand = new Random();
		public static void Main (string[] args)
		{
			//testNetwork ();
		}
		public static void populationSetup()
		{
			for(int i=0;i<populationSize;i++)
			{
				Network net = new Network();

				//add inputs
				for(int j=0;j<inputN;j++)
				{
					Neuron n = new Neuron();
					net.addNeuron(n);
					net.addInput(j);
				}

				//add outputs
				for(int j=0;j<outputN;j++)
				{
					Neuron n = new Neuron();
					net.addNeuron(n);
					net.addOutput(inputN + j);
				}

				//add edges to outputs
				for(int j=0;j<inputN;j++)
				{
					for(int k=0;k<outputN;k++)
					{
						net.addEdge(j, k+inputN, rand.NextDouble() * 4.0 - 2.0);
					}
				}
			}


		}


		public static void testNetwork()
		{
			Network a = new Network ();
			for (int i=0; i<7; i++) {
				Neuron n = new Neuron ();
				a.addNeuron (n);
			}

			a.addInput (0);
			a.addInput (1);
			a.addInput (2);
			a.addOutput (3);
			bool x = a.addEdge (1, 4, 0.1);
			x = a.addEdge (2, 4, 0.2);
			x = a.addEdge (1, 5, 0.3);
			x = a.addEdge (0, 5, 0.4);
			x = a.addEdge (4, 5, 0.5);
			x = a.addEdge (5, 6, 0.6);
			x = a.addEdge (6, 3, 0.7);

			List<double> inputs = new List<double> ();
			inputs.Add (-1);
			inputs.Add (0.5);
			inputs.Add (1.0);
			List<double> outputs = a.calculateOutput (inputs);
			for (int i=0; i<outputs.Count; i++) {
				Console.WriteLine ("i:" + outputs [i]);
			}
			Console.WriteLine ("done");
		}
		
		
		private double xorfitness(double inputA, double inputB, double output)
		{
			double answer = inputA + inputB;
			if (answer > 1.01) {
				answer = 0.0;
			}
			return 1.0 / (1.0 + Math.Abs (answer - output));
		}
	}
}
